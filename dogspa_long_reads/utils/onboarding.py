import json
import logging
import os
from collections import OrderedDict

import pandas as pd
from nebelung.terra_workspace import TerraWorkspace
from nebelung.utils import expand_dict_columns
from pandera.typing import DataFrame as TypedDataFrame

from dogspa_long_reads.types import (
    IdentifiedSrcBams,
    ModelsAndChildren,
    SamplesForGumbo,
    SamplesMaybeInGumbo,
    SamplesWithCDSIDs,
    SamplesWithMetadata,
    SamplesWithShortReadMetadata,
    SeqTable,
    VersionedSamples,
)
from dogspa_long_reads.utils.gcp import (
    check_file_sizes,
    copy_to_cclebams,
    update_sample_file_urls,
)
from dogspa_long_reads.utils.utils import (
    assign_hashed_uuids,
    df_to_model,
    model_to_df,
    send_slack_message,
    uuid_to_base62,
)
from gumbo_gql_client import GumboClient, omics_sequencing_insert_input


def do_onboard_samples(
    gcp_project_id: str,
    gcs_destination_bucket: str,
    gcs_destination_prefix: str,
    uuid_namespace: str,
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

    samples = pd.DataFrame()

    # start tracking issues to store in Gumbo seq table
    samples["issue"] = pd.Series([set()] * len(samples))
    samples["blacklist"] = False

    # compare file sizes to filter out samples that are in Gumbo already
    samples = check_already_in_gumbo(samples, seq_table, size_col_name="bam_size")
    stats["n not yet in Gumbo"] = (~samples["already_in_gumbo"]).sum()
    report["not yet in Gumbo"] = samples.loc[~samples["already_in_gumbo"]]

    if len(samples) == 0:
        send_slack_message(
            os.getenv("SLACK_WEBHOOK_URL_ERRORS"),
            os.getenv("SLACK_WEBHOOK_URL_STATS"),
            stats,
            report,
            terra_workspace,
            dry_run,
        )
        return

    # join metadata to current samples
    samples = join_metadata(samples, seq_table)
    samples = join_short_read_metadata(samples, seq_table)

    # check that BAM file sizes are above minimum threshold
    samples, blacklisted = check_file_sizes(samples)
    stats["n with BAM file too small"] = blacklisted.sum()
    report["BAM file too small"] = samples.loc[blacklisted]

    # assign sequencing IDs
    samples = assign_cds_ids(samples, uuid_namespace)

    # copy files to our own bucket
    sample_files = copy_to_cclebams(
        samples, gcp_project_id, gcs_destination_bucket, gcs_destination_prefix, dry_run
    )

    # replace URLs with ones for our bucket whenever the copy operation succeeded
    samples, blacklisted = update_sample_file_urls(samples, sample_files)
    stats["n with failed file copies"] = blacklisted.sum()
    report["failed copies"] = samples.loc[blacklisted]

    # rename and set some columns for Gumbo
    gumbo_samples = apply_col_map(samples.loc[~samples["already_in_gumbo"]])

    # increment version numbers for samples with profile IDs already in seq table
    gumbo_samples = increment_sample_versions(gumbo_samples, seq_table)

    # upload the samples to the Gumbo sequencing table
    gumbo_samples = upload_to_gumbo(gumbo_client, gumbo_samples, dry_run)
    stats["n successfully uploaded"] = len(gumbo_samples)
    report["successfully uploaded"] = gumbo_samples

    # upload the samples to Terra
    upsert_terra_samples(terra_workspace, samples, dry_run)

    send_slack_message(
        os.getenv("SLACK_WEBHOOK_URL_ERRORS"),
        os.getenv("SLACK_WEBHOOK_URL_STATS"),
        stats,
        report,
        terra_workspace,
        dry_run,
    )


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

    seq_table[["blacklist_omics", "blacklist"]] = seq_table[
        ["blacklist_omics", "blacklist"]
    ].fillna(False)

    return TypedDataFrame[SeqTable](seq_table.drop(columns="main_sequencing_id"))


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

    samples_annot = samples.merge(metadata, how="left", on="model_id")

    return TypedDataFrame[SamplesWithMetadata](samples_annot)


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

    main_seq_ids = (
        sr.loc[
            sr["is_main_sequencing_id"],
            ["model_id", "sequencing_id", "source_priority"],
        ]
        .sort_values("source_priority")
        .groupby("model_id")
        .nth(0)
        .rename(columns={"sequencing_id": "main_sequencing_id"})
    )

    samples_annot = samples.merge(main_seq_ids, how="left", on="model_id")

    sr_rna = sr.sort_values("source_priority").groupby("model_id").nth(0)

    sr_rna = sr_rna[["model_id", "profile_id", "bai_filepath", "bam_filepath"]].rename(
        columns={
            "profile_id": "sr_profile_id",
            "bai_filepath": "sr_bai_filepath",
            "bam_filepath": "sr_bam_filepath",
        }
    )

    samples_annot = samples_annot.merge(sr_rna, how="left", on="model_id")

    return TypedDataFrame[SamplesWithShortReadMetadata](samples_annot)


def assign_cds_ids(
    samples: TypedDataFrame[SamplesWithShortReadMetadata], uuid_namespace: str
) -> TypedDataFrame[SamplesWithCDSIDs]:
    """
    Assign a "CDS-" ID to each sample by hashing relevant columns.

    :param samples: the data frame of samples
    :param uuid_namespace: a namespace for UUIDv3 IDs
    :return: the data frame with CDS IDs
    """

    logging.info("Assigning CDS IDs...")

    samples_w_ids = assign_hashed_uuids(
        samples,
        uuid_namespace=uuid_namespace,
        uuid_col_name="sequencing_id",
        subset=[
            "model_id",
            "crc32c",
            "size",
            "gcs_obj_updated_at",
        ],
    )

    # convert UUIDs to base62 and truncate to match existing ID format
    samples_w_ids["sequencing_id"] = "CDS-" + samples_w_ids["sequencing_id"].apply(
        uuid_to_base62
    ).str.slice(-6)

    return TypedDataFrame[SamplesWithCDSIDs](samples_w_ids)


def check_already_in_gumbo(
    samples: TypedDataFrame[IdentifiedSrcBams],
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
        samples["size"]
        .isin(seq_table[size_col_name].dropna().astype("int64"))
        .astype("bool")
    )

    return TypedDataFrame[SamplesMaybeInGumbo](samples)


def apply_col_map(
    samples: TypedDataFrame[SamplesWithCDSIDs],
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
    gumbo_samples["sequencing_date"] = gumbo_samples["gcs_obj_updated_at"]
    gumbo_samples["stranded"] = True

    # rename columns for Gumbo
    gumbo_samples = gumbo_samples.rename(
        columns={
            "bam_url": "bam_filepath",
            "size": "bam_size",
            "gcs_obj_updated_at": "update_time",
            "crc32c": "bam_crc32c_hash",
        }
    ).loc[
        :,
        [
            "sequencing_id",
            "bam_filepath",
            "profile_id",
            "bam_size",
            "update_time",
            "sequencing_date",
            "source",
            "bam_crc32c_hash",
            "expected_type",
            "issue",
            "blacklist",
        ],
    ]

    # concat issues for each sample into single string
    gumbo_samples["issue"] = gumbo_samples["issue"].apply(
        lambda x: "; ".join(x) if len(x) > 0 else pd.NA
    )

    return TypedDataFrame[SamplesForGumbo](gumbo_samples)


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

    return TypedDataFrame[VersionedSamples](versioned_samples)


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
        res = gumbo_client.insert_omics_sequencings(username="dogspa", objects=objects)
        affected_rows = res.insert_omics_sequencing.affected_rows  # type: ignore
        logging.info(f"Inserted {affected_rows} samples to Gumbo")

    return samples_to_upload


def upsert_terra_samples(
    tw: TerraWorkspace, samples: TypedDataFrame[SamplesWithCDSIDs], dry_run: bool
) -> None:
    """
    Upsert sample and participant data to Terra data tables.

    :param tw: a TerraWorkspace instance
    :param samples: the data frame of samples
    :param dry_run: whether to skip updates to external data stores
    """

    logging.info(f"Upserting data to {tw.workspace_name} data tables")

    renames = {
        "sequencing_id": "entity:sample_id",
        "model_id": "participant_id",
        "bam_url": "LR_bam_filepath",
        "cell_line_name": "CellLineName",
        "stripped_cell_line_name": "StrippedCellLineName",
        "model_condition_id": "ModelCondition",
        "profile_id": "LongReadProfileID",
        "sr_profile_id": "ShortReadProfileID",
        "sr_bai_filepath": "SR_bai_filepath",
        "sr_bam_filepath": "SR_bam_filepath",
        "main_sequencing_id": "MainSequencingID",
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
        tw.upload_entities(participants)

    # upsert samples
    if dry_run:
        logging.info(f"(skipping) Upserting {len(terra_samples)} samples")
    else:
        tw.upload_entities(terra_samples.drop(columns="participant_id"))

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
        tw.upload_entities(participant_samples)
