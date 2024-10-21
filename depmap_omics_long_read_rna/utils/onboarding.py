import json
import logging
import os
from collections import OrderedDict

import pandas as pd
from nebelung.terra_workspace import TerraWorkspace
from nebelung.utils import expand_dict_columns, type_data_frame
from pandera.typing import DataFrame as TypedDataFrame

from depmap_omics_long_read_rna.types import (
    DeliveryBams,
    GumboClient,
    ModelsAndChildren,
    OnboardingSamples,
    SamplesForGumbo,
    SamplesMaybeInGumbo,
    SamplesWithCDSIDs,
    SamplesWithMetadata,
    SamplesWithShortReadMetadata,
    SeqTable,
    VersionedSamples,
)
from depmap_omics_long_read_rna.utils.gcp import (
    copy_to_cclebams,
    get_objects_metadata,
    update_sample_file_urls,
)
from depmap_omics_long_read_rna.utils.utils import (
    df_to_model,
    model_to_df,
    send_slack_message,
)
from gumbo_gql_client import omics_sequencing_insert_input


def do_onboard_samples(
    gcp_project_id: str,
    unaligned_gcs_destination_bucket: str,
    unaligned_gcs_destination_prefix: str,
    aligned_gcs_destination_bucket: str,
    aligned_gcs_destination_prefix: str,
    terra_workspace: TerraWorkspace,
    gumbo_client: GumboClient,
    dry_run: bool = False,
) -> None:
    # keep track of filtering and quality control checks performed
    stats = OrderedDict()
    report = OrderedDict()

    # get the sequencing and profile tables from Gumbo
    models = model_to_df(gumbo_client.get_models_and_children(), ModelsAndChildren)
    seq_table = explode_and_expand_models(models)
    seq_table_lr = seq_table.loc[seq_table["datatype"].eq("long_read_rna")]

    delivery_bams = terra_workspace.get_entities("delivery_bam", DeliveryBams)

    samples = (
        delivery_bams[
            [
                "delivery_bam_id",
                "model_id",
                "aligned_bam",
                "aligned_bai",
                "delivery_bam",
                "delivery_bam_size",
                "delivery_bam_crc32c",
                "delivery_bam_updated_at",
            ]
        ]
        .dropna()  # remove samples that haven't been aligned yet
        .rename(columns={"delivery_bam_id": "sequencing_id"})
        .reset_index(drop=True)
    )

    # start tracking issues to store in Gumbo seq table
    samples["issue"] = pd.Series([set()] * len(samples))
    samples["blacklist"] = False

    # compare file sizes to filter out samples that are in Gumbo already
    samples = check_already_in_gumbo(
        samples, seq_table_lr, size_col_name="unaligned_bam_size"
    )
    stats["n not yet in Gumbo"] = (~samples["already_in_gumbo"]).sum()
    samples = samples.loc[~samples["already_in_gumbo"]].drop(columns="already_in_gumbo")
    report["not yet in Gumbo"] = samples

    # if len(samples) == 0:
    #     send_slack_message(
    #         os.getenv("SLACK_WEBHOOK_URL_ERRORS"),
    #         os.getenv("SLACK_WEBHOOK_URL_STATS"),
    #         stats,
    #         report,
    #         terra_workspace,
    #         dry_run,
    #     )
    #     return

    # join metadata to current samples
    samples = join_metadata(samples, seq_table_lr)
    samples = join_short_read_metadata(samples, seq_table)

    # check that BAM file sizes are above minimum threshold
    # samples, blacklisted = check_file_sizes(samples)
    # stats["n with BAM file too small"] = blacklisted.sum()
    # report["BAM file too small"] = samples.loc[blacklisted]

    # copy files to our own bucket
    unaligned_sample_files = copy_to_cclebams(
        samples,
        bam_bai_colnames=["delivery_bam"],
        gcp_project_id=gcp_project_id,
        gcs_destination_bucket=unaligned_gcs_destination_bucket,
        gcs_destination_prefix=unaligned_gcs_destination_prefix,
        dry_run=dry_run,
    )

    aligned_sample_files = copy_to_cclebams(
        samples,
        bam_bai_colnames=["aligned_bam", "aligned_bai"],
        gcp_project_id=gcp_project_id,
        gcs_destination_bucket=aligned_gcs_destination_bucket,
        gcs_destination_prefix=aligned_gcs_destination_prefix,
        dry_run=dry_run,
    )

    # replace URLs with ones for our bucket whenever the copy operation succeeded
    samples = update_sample_file_urls(
        samples, unaligned_sample_files, bam_bai_colnames=["delivery_bam"]
    )
    samples = update_sample_file_urls(
        samples,
        aligned_sample_files,
        bam_bai_colnames=["aligned_bam", "aligned_bai"],
    )

    missing_files = pd.Series(
        samples[["delivery_bam", "aligned_bam", "aligned_bai"]].isna().any(axis=1)
    )
    samples.loc[missing_files, "issue"] = samples.loc[missing_files, "issue"].apply(
        lambda x: x.union({"couldn't copy BAM/BAI"})
    )
    samples.loc[missing_files, "blacklist"] = True

    stats["n with failed file copies"] = missing_files.sum()
    report["failed copies"] = samples.loc[missing_files]

    # get object metadata for the aligned BAMs we just copied
    copied_aligned_bams = get_objects_metadata(samples["aligned_bam"]).rename(
        columns={
            "url": "aligned_bam",
            "crc32c": "bam_crc32c_hash",
            "size": "bam_size",
            "gcs_obj_updated_at": "update_time",
        }
    )

    samples_complete = samples.merge(copied_aligned_bams, how="left", on="aligned_bam")

    # rename and set some columns for Gumbo
    gumbo_samples = apply_col_map(samples_complete)

    # increment version numbers for samples with profile IDs already in seq table
    gumbo_samples = increment_sample_versions(gumbo_samples, seq_table)

    # upload the samples to the Gumbo sequencing table
    gumbo_samples = upload_to_gumbo(gumbo_client, gumbo_samples, dry_run)
    stats["n successfully uploaded"] = len(gumbo_samples)
    report["successfully uploaded"] = gumbo_samples

    # upload the samples to Terra
    upsert_terra_samples(terra_workspace, samples_complete, dry_run)

    # send_slack_message(
    #     os.getenv("SLACK_WEBHOOK_URL_ERRORS"),
    #     os.getenv("SLACK_WEBHOOK_URL_STATS"),
    #     stats,
    #     report,
    #     terra_workspace,
    #     dry_run,
    # )


def explode_and_expand_models(
    models: TypedDataFrame[ModelsAndChildren],
) -> TypedDataFrame[SeqTable]:
    seq_table = models.copy()

    for c in ["model_conditions", "omics_profiles", "omics_sequencings"]:
        seq_table = seq_table.explode(c).reset_index(drop=True)
        seq_table = expand_dict_columns(seq_table, name_columns_with_parent=False)

    seq_table = seq_table.dropna(subset=["model_condition_id", "profile_id"])

    seq_table["is_main_sequencing_id"] = seq_table["main_sequencing_id"].eq(
        seq_table["sequencing_id"]
    )

    seq_table[["blacklist_omics", "blacklist"]] = (
        seq_table[["blacklist_omics", "blacklist"]].astype("boolean").fillna(False)
    )

    return type_data_frame(seq_table.drop(columns="main_sequencing_id"), SeqTable)


def join_metadata(
    samples: TypedDataFrame[SamplesMaybeInGumbo],
    seq_table: TypedDataFrame[SeqTable],
) -> TypedDataFrame[SamplesWithMetadata]:
    metadata = seq_table.loc[
        seq_table["datatype"].eq("long_read_rna")
        & ~seq_table["blacklist"]
        & ~seq_table["blacklist_omics"],
        [
            "model_id",
            "model_condition_id",
            "profile_id",
            "cell_line_name",
            "stripped_cell_line_name",
        ],
    ]

    # need to arbitrarily de-dup since there isn't a way of arbitrarily picking a
    # profile when more than one might belong to an model condition
    metadata = metadata.drop_duplicates(subset="model_id")

    samples_annot = samples.merge(metadata, how="left", on="model_id")

    return type_data_frame(samples_annot, SamplesWithMetadata)


def join_short_read_metadata(
    samples: TypedDataFrame[SamplesWithMetadata], seq_table: TypedDataFrame[SeqTable]
) -> TypedDataFrame[SamplesWithShortReadMetadata]:
    # reproduce logic of `makeDefaultModelTable` in depmap_omics_upload
    source_priority = [
        "BROAD",
        "DEPMAP",
        "IBM",
        "CCLE2",
        "SANGER",
        "PRISM",
        "CCLF",
        "CHORDOMA",
        "",
    ]

    sp_df = pd.DataFrame(
        {
            "source": source_priority,
            "source_priority": list(range(len(source_priority))),
        }
    )

    sr = seq_table.loc[
        seq_table["model_id"].isin(samples["model_id"])
        & seq_table["datatype"].eq("rna")
        & ~seq_table["blacklist"]
        & ~seq_table["blacklist_omics"]
    ].copy()

    sr["source"] = sr["source"].fillna("")
    sr = sr.merge(sp_df, how="left", on="source")

    assert sr["source_priority"].notna().all()

    sr_rna = sr.sort_values("source_priority").groupby("model_id").nth(0)

    sr_rna = sr_rna[["model_id", "profile_id"]].rename(
        columns={"profile_id": "sr_profile_id"}
    )

    samples_annot = samples.merge(sr_rna, how="left", on="model_id")

    return type_data_frame(samples_annot, SamplesWithShortReadMetadata)


def check_already_in_gumbo(
    samples: TypedDataFrame[OnboardingSamples],
    seq_table: TypedDataFrame[SeqTable],
    size_col_name: str,
) -> TypedDataFrame[SamplesMaybeInGumbo]:
    """
    Mark samples that are already in Gumbo by comparing stored file sizes.

    :param samples: the data frame of samples
    :param seq_table: the data frame of Gumbo Omics sequencing data
    :param size_col_name: the Gumbo column containing file sizes
    :return: the samples data frame with `already_in_gumbo` column
    """

    logging.info("Checking Gumbo for existing samples...")

    samples["already_in_gumbo"] = (
        samples["delivery_bam_size"]
        .isin(seq_table[size_col_name].dropna().astype("int64"))
        .astype("bool")
    )

    return type_data_frame(samples, SamplesMaybeInGumbo)


def apply_col_map(
    samples: TypedDataFrame[SamplesWithShortReadMetadata],
) -> TypedDataFrame[SamplesForGumbo]:
    """
    Rename the columns in the samples data frame to their corresponding names in Gumbo.

    :param samples: the data frame of samples
    :return: the data frame with just the relevant renamed columns
    """

    logging.info("Renaming columns for Gumbo...")

    gumbo_samples = samples.copy()
    gumbo_samples["source"] = "DEPMAP"
    gumbo_samples["expected_type"] = "long_read_rna"
    gumbo_samples["sequencing_date"] = gumbo_samples["delivery_bam_updated_at"]
    gumbo_samples["update_time"] = gumbo_samples["delivery_bam_updated_at"]
    gumbo_samples["stranded"] = True

    # select columns for Gumbo
    gumbo_samples = gumbo_samples.rename(
        columns={
            "delivery_bam": "unaligned_bam_filepath",
            "delivery_bam_size": "unaligned_bam_size",
            "delivery_bam_crc32c": "unaligned_bam_crc32c_hash",
            "aligned_bam": "bam_filepath",
            "aligned_bai": "bai_filepath",
            "aligned_bam_size": "bam_size",
            "aligned_bam_crc32c": "bam_crc32c_hash",
        }
    ).loc[
        :,
        [
            "sequencing_id",
            "bam_filepath",
            "bai_filepath",
            "bam_crc32c_hash",
            "bam_size",
            "unaligned_bam_filepath",
            "unaligned_bam_crc32c_hash",
            "unaligned_bam_size",
            "profile_id",
            "update_time",
            "sequencing_date",
            "source",
            "expected_type",
            "issue",
            "blacklist",
        ],
    ]

    # concat issues for each sample into single string
    gumbo_samples["issue"] = gumbo_samples["issue"].apply(
        lambda x: "; ".join(x) if len(x) > 0 else pd.NA
    )

    return type_data_frame(gumbo_samples, SamplesForGumbo)


def increment_sample_versions(
    samples: TypedDataFrame[SamplesForGumbo], seq_table: TypedDataFrame[SeqTable]
) -> TypedDataFrame[VersionedSamples]:
    """
    For sample profile IDs that already exist in the Gumbo profiles table, increment
    the version numbers beyond the default value of 1.

    :param samples: the data frame of samples
    :param seq_table: the data frame of Gumbo sequence data
    :return: the samples data frame with optionally incremented version numbers
    """

    logging.info("Incrementing sample version numbers...")

    # count distinct profile IDs in existing seq table from among those observed in the
    # current batch of samples
    existing_profile_id_counts = (
        seq_table.loc[
            seq_table["expected_type"].eq("long_read_rna")
            & seq_table["profile_id"].isin(samples["profile_id"]),
            "profile_id",
        ]
        .value_counts()
        .to_frame(name="version_n")
        .reset_index(names=["profile_id"])
    )

    # increment version numbers in samples if applicable
    versioned_samples = samples.merge(
        existing_profile_id_counts, how="left", on="profile_id"
    )
    versioned_samples["version_n"] = versioned_samples["version_n"].fillna(0)
    versioned_samples["version"] = 1 + versioned_samples["version_n"]
    versioned_samples = versioned_samples.drop(columns="version_n")

    return type_data_frame(versioned_samples, VersionedSamples)


def upload_to_gumbo(
    gumbo_client: GumboClient, samples: TypedDataFrame[VersionedSamples], dry_run: bool
) -> TypedDataFrame[VersionedSamples]:
    """
    Upload the samples data frame to the Gumbo sequencing table.

    :param gumbo_client: an instance of the Gumbo GraphQL client
    :param samples: the data frame of samples
    :param dry_run: whether to skip updates to external data stores
    :return: the samples data frame that was uploaded to Gumbo
    """

    # match historical behavior to exclude bad records except for too-small
    # BAM files (so that other types of bad records like missing GCS
    # objects that might get resolved later still get processed the next time
    # this automation is run)
    samples_to_upload = samples.loc[
        ~samples["blacklist"] | samples["issue"].eq("BAM file too small")
    ]

    if dry_run:
        logging.info(f"(skipping) Inserting {len(samples_to_upload)} samples")
        return samples_to_upload

    logging.info(f"Inserting {len(samples_to_upload)} samples to omics_sequencing")

    if len(samples) > 0:
        objects = df_to_model(samples_to_upload, omics_sequencing_insert_input)
        res = gumbo_client.insert_omics_sequencings(
            username="depmap-omics-long-read-rna", objects=objects
        )
        affected_rows = res.insert_omics_sequencing.affected_rows  # type: ignore
        logging.info(f"Inserted {affected_rows} samples to Gumbo")

    return samples_to_upload


def upsert_terra_samples(
    terra_workspace: TerraWorkspace,
    samples: TypedDataFrame[SamplesWithShortReadMetadata],
    dry_run: bool,
) -> None:
    """
    Upsert sample and participant data to Terra data tables.

    :param terra_workspace: a TerraWorkspace instance
    :param samples: the data frame of samples
    :param dry_run: whether to skip updates to external data stores
    """

    logging.info(f"Upserting data to {terra_workspace.workspace_name} data tables")

    renames = {
        "sequencing_id": "entity:sample_id",
        "model_id": "participant_id",
        "delivery_bam": "ubam",
        "aligned_bam": "aligned_bam",
        "aligned_bai": "aligned_bai",
        "cell_line_name": "cell_line_name",
        "stripped_cell_line_name": "cell_line_name_stripped",
        "model_condition_id": "model_condition_id",
        "profile_id": "profile_id",
        "sr_profile_id": "short_read_profile_id",
    }

    terra_samples = samples.rename(columns=renames).loc[:, renames.values()]

    terra_samples["participant"] = terra_samples["participant_id"].apply(
        lambda x: json.dumps({"entityType": "participant", "entityName": x})
    )

    # upsert participants
    participants = (
        terra_samples[["participant_id"]]
        .rename(columns={"participant": "entity:participant_id"})
        .drop_duplicates()
    )

    if dry_run:
        logging.info(f"(skipping) Upserting {len(participants)} participants")
    else:
        terra_workspace.upload_entities(participants)

    # upsert samples
    if dry_run:
        logging.info(f"(skipping) Upserting {len(terra_samples)} samples")
    else:
        terra_workspace.upload_entities(terra_samples.drop(columns="participant_id"))

    # upsert the join table between participants and samples
    participant_samples = (
        terra_samples.groupby("participant_id")["entity:sample_id"]
        .agg(list)
        .apply(
            lambda x: json.dumps([{"entityType": "sample", "entityName": y} for y in x])
        )
        .reset_index()
        .rename(
            columns={
                "participant": "entity:participant_id",
                "entity:sample_id": "samples",
            }
        )
    )

    if dry_run:
        logging.info(
            f"(skipping) Upserting {len(participant_samples)} participant-samples"
        )
    else:
        terra_workspace.upload_entities(participant_samples)
