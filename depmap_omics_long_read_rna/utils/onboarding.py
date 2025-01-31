import json
import logging
from urllib.parse import quote_plus

import pandas as pd
import requests
from nebelung.terra_workspace import TerraWorkspace
from nebelung.utils import type_data_frame
from pandera.typing import DataFrame as TypedDataFrame
from pd_flatten import pd_flatten

from depmap_omics_long_read_rna.types import (
    DeliveryBams,
    GumboClient,
    ModelsAndChildren,
    OnboardingSamples,
    SamplesForGumbo,
    SamplesMaybeInGumbo,
    SamplesWithMetadata,
    SeqTable,
    ShortReadMetadata,
    VersionedSamples,
)
from depmap_omics_long_read_rna.utils.gcp import (
    copy_to_cclebams,
    get_objects_metadata,
    update_sample_file_urls,
)
from depmap_omics_long_read_rna.utils.utils import get_secret_from_sm, model_to_df
from gumbo_gql_client import (
    omics_sequencing_insert_input,
    omics_sequencing_obj_rel_insert_input,
    onboarding_job_insert_input,
    onboarding_sample_arr_rel_insert_input,
    onboarding_sample_insert_input,
    sequencing_alignment_insert_input,
    sequencing_alignment_obj_rel_insert_input,
)


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
    """
    Onboard the latest uBAMs and aligned BAMs to Gumbo and Terra.

    :param gcp_project_id: a GCP project ID to use for billing
    :param unaligned_gcs_destination_bucket: the GCS bucket name for uBAMs
    :param unaligned_gcs_destination_prefix: the GCS bucket prefix for uBAMs
    :param aligned_gcs_destination_bucket: the GCS bucket name for aligned BAMs
    :param aligned_gcs_destination_prefix: the GCS bucket prefix for aligned BAMs
    :param terra_workspace: a TerraWorkspace instance
    :param gumbo_client: an instance of the Gumbo GraphQL client
    :param dry_run: whether to skip updates to external data stores
    """

    # get onboarding workspace metadata from Gumbo (since long reads onboarding is
    # "special", most of the information we'd get there is simply hardcoded in this
    # repo)
    onboarding_workspace = gumbo_client.get_onboarding_workspace()
    excluded_terra_sample_ids = onboarding_workspace.records[
        0
    ].excluded_terra_sample_ids
    min_file_size = onboarding_workspace.records[0].min_file_size

    onboarding_job = onboarding_job_insert_input(
        onboarding_workspace_id=onboarding_workspace.records[0].id,
        succeeded=False,
    )

    try:
        onboarding_job = prep_onboarding_job(
            onboarding_job,
            gcp_project_id,
            unaligned_gcs_destination_bucket,
            unaligned_gcs_destination_prefix,
            aligned_gcs_destination_bucket,
            aligned_gcs_destination_prefix,
            terra_workspace,
            excluded_terra_sample_ids,
            min_file_size,
            gumbo_client,
            dry_run,
        )
        onboarding_job.succeeded = True
    except BrokenPipeError as e:
        logging.error(e)
    finally:
        onboarding_job = persist_onboarding_job(gumbo_client, onboarding_job, dry_run)
        send_slack_message(onboarding_job, dry_run)


def prep_onboarding_job(
    onboarding_job: onboarding_job_insert_input,
    gcp_project_id: str,
    unaligned_gcs_destination_bucket: str,
    unaligned_gcs_destination_prefix: str,
    aligned_gcs_destination_bucket: str,
    aligned_gcs_destination_prefix: str,
    terra_workspace: TerraWorkspace,
    excluded_terra_sample_ids: list[str],
    min_file_size: int,
    gumbo_client: GumboClient,
    dry_run: bool = False,
) -> onboarding_job_insert_input:
    """
    Onboard the latest uBAMs and aligned BAMs to Gumbo and Terra.

    :param onboarding_job: the onboarding job
    :param gcp_project_id: a GCP project ID to use for billing
    :param unaligned_gcs_destination_bucket: the GCS bucket name for uBAMs
    :param unaligned_gcs_destination_prefix: the GCS bucket prefix for uBAMs
    :param aligned_gcs_destination_bucket: the GCS bucket name for aligned BAMs
    :param aligned_gcs_destination_prefix: the GCS bucket prefix for aligned BAMs
    :param terra_workspace: a TerraWorkspace instance
    :param excluded_terra_sample_ids: a list of delivery BAM/sample IDs to exclude
    :param min_file_size: the minimum delivery BAM/CRAM file size
    :param gumbo_client: an instance of the Gumbo GraphQL client
    :param dry_run: whether to skip updates to external data stores
    """

    samples = get_delivery_bams(terra_workspace)
    onboarding_job.n_samples = len(samples)
    onboarding_job.n_samples_succeeded = 0
    onboarding_job.n_samples_failed = 0

    # check if any Terra sample IDs are in the exclusion list and drop them now
    is_excluded = samples["sequencing_id"].isin(excluded_terra_sample_ids)
    onboarding_job.n_samples_excluded = is_excluded.sum()
    samples = samples.loc[~is_excluded]

    # get the sequencing and profile tables from Gumbo
    models = model_to_df(gumbo_client.get_models_and_children(), ModelsAndChildren)
    seq_table = explode_and_expand_models(models)

    # filter so that sequencing alignment-related columns are just the ones for the
    # GP-delivered uBAM
    seq_table_lr = seq_table.loc[
        seq_table["datatype"].eq("long_read_rna")
        & seq_table["sequencing_alignment_source"].eq("GP")
    ]

    # compare file sizes to filter out samples that are in Gumbo already
    samples = check_already_in_gumbo(
        samples, seq_table_lr, size_col_name="cram_bam_size"
    )
    onboarding_job.n_samples_new = (~samples["already_in_gumbo"]).sum()

    if onboarding_job.n_samples_new == 0:
        logging.info("No new samples")
        return onboarding_job

    samples = samples.loc[~samples["already_in_gumbo"]].drop(columns="already_in_gumbo")

    # join metadata to current samples
    samples = join_metadata(samples, seq_table_lr)

    # check that BAM file sizes are above minimum threshold
    samples = check_file_sizes(samples, min_file_size)

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
        lambda x: x.union({"couldn't copy CRAM/CRAI/BAM/BAI"})
    )

    # get object metadata for the aligned BAMs we just copied
    copied_aligned_bams = get_objects_metadata(samples["aligned_bam"]).rename(
        columns={
            "url": "aligned_bam",
            "crc32c": "bam_crc32c_hash",
            "size": "bam_size",
            "gcs_obj_updated_at": "update_time",
        }
    )

    samples_complete = type_data_frame(
        samples.merge(copied_aligned_bams, how="left", on="aligned_bam"),
        SamplesWithMetadata,
    )

    # rename and set some columns for Gumbo
    gumbo_samples = apply_col_map(samples_complete)

    # increment version numbers for samples with profile IDs already in seq table
    gumbo_samples = increment_sample_versions(gumbo_samples, seq_table)

    # prepare the various Gumbo records for sequencings and alignments
    onboarding_job = prep_onboarding_samples(onboarding_job, samples=gumbo_samples)

    # upload the samples to the Terra `sample` data table
    upsert_terra_samples(terra_workspace, samples_complete, dry_run)

    return onboarding_job


def get_delivery_bams(
    terra_workspace: TerraWorkspace,
) -> TypedDataFrame[OnboardingSamples]:
    """
    Get a data frame of samples (e.g. delivery BAMs) from the Terra workspace.

    :param terra_workspace: a TerraWorkspace instance
    :return: a data frame of samples to onboard
    """

    delivery_bams = terra_workspace.get_entities("delivery_bam", DeliveryBams)
    samples = (
        delivery_bams.loc[
            :,
            [
                "delivery_bam_id",
                "model_id",
                "aligned_bam",
                "aligned_bai",
                "delivery_bam",
                "delivery_bam_size",
                "delivery_bam_crc32c",
                "delivery_bam_updated_at",
            ],
        ]
        .dropna()  # remove samples that haven't been aligned yet
        .rename(columns={"delivery_bam_id": "sequencing_id"})
        .reset_index(drop=True)
    )

    # start tracking issues to store in Gumbo seq table
    samples["issue"] = pd.Series([set()] * len(samples))

    return type_data_frame(samples, OnboardingSamples)


def explode_and_expand_models(
    models: TypedDataFrame[ModelsAndChildren],
) -> TypedDataFrame[SeqTable]:
    """
    Unnest columns from the GraphQL call for Gumbo models and their childen.

    :param models: a data frame of models with nested profile/sequencing/etc. data
    :return: a wide version of the data frame without nesting
    """

    seq_table = pd_flatten(models, name_columns_with_parent=False)

    seq_table[["blacklist_omics", "blacklist"]] = (
        seq_table[["blacklist_omics", "blacklist"]].astype("boolean").fillna(False)
    )

    seq_table = seq_table.rename(
        columns={
            "id": "sequencing_alignment_id",
            "url": "cram_bam_url",
            "index_url": "crai_bai_url",
            "size": "cram_bam_size",
        }
    )

    return type_data_frame(seq_table, SeqTable, remove_unknown_cols=True)


def join_metadata(
    samples: TypedDataFrame[SamplesMaybeInGumbo],
    seq_table: TypedDataFrame[SeqTable],
) -> TypedDataFrame[SamplesWithMetadata]:
    """
    Join metadata for long read samples from Gumbo.

    :param samples: the data frame of samples
    :param seq_table: the data frame of Gumbo Omics sequencing data
    :return: the samples with additional metadata from Gumbo
    """

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


def check_file_sizes(
    samples: TypedDataFrame[SamplesWithMetadata],
    min_file_size: int,
) -> TypedDataFrame[SamplesWithMetadata]:
    """
    Check whether BAM file sizes are above configured minimum threshold.

    :param samples: the data frame of samples
    :param min_file_size: the minimum file size
    :return: the samples data frame with issue column filled out for too-small BAM/CRAM
    files
    """

    logging.info("Checking BAM file sizes...")

    bam_too_small = samples["delivery_bam_size"] < min_file_size
    samples.loc[bam_too_small, "issue"] = samples.loc[bam_too_small, "issue"].apply(
        lambda x: x.union({"BAM/CRAM file too small"})
    )

    return type_data_frame(samples, SamplesWithMetadata)


def check_already_in_gumbo(
    samples: TypedDataFrame[OnboardingSamples],
    seq_table: TypedDataFrame[SeqTable],
    size_col_name: str,
) -> TypedDataFrame[SamplesMaybeInGumbo]:
    """
    Mark samples that are already in Gumbo by comparing stored file sizes.
    Get non-blacklisted short read RNA samples to map to long read ones.

    :param samples: the data frame of samples
    :param seq_table: the data frame of Gumbo Omics sequencing data
    :param size_col_name: the Gumbo column containing file sizesa data frame mapping model IDs to short read profile and sequencing IDs
    :return: the samples data frame with `already_in_gumbo` column
    """

    logging.info("Checking Gumbo for existing samples...")

    seq_table_ubam = seq_table.loc[seq_table["reference_genome"].isna()]

    samples["already_in_gumbo"] = (
        samples["delivery_bam_size"]
        .isin(seq_table_ubam[size_col_name].dropna().astype("int64"))
        .astype("bool")
    )

    return type_data_frame(samples, SamplesMaybeInGumbo)


def apply_col_map(
    samples: TypedDataFrame[SamplesWithMetadata],
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
            "delivery_bam": "unaligned_bam_url",
            "delivery_bam_size": "unaligned_bam_size",
            "delivery_bam_crc32c": "unaligned_bam_crc32c_hash",
            "aligned_bam": "bam_url",
            "aligned_bai": "bai_url",
            "aligned_bam_size": "bam_size",
            "aligned_bam_crc32c": "bam_crc32c_hash",
        }
    ).loc[
        :,
        [
            "sequencing_id",
            "bam_url",
            "bam_crc32c_hash",
            "bam_size",
            "bai_url",
            "unaligned_bam_url",
            "unaligned_bam_crc32c_hash",
            "unaligned_bam_size",
            "profile_id",
            "update_time",
            "sequencing_date",
            "source",
            "expected_type",
            "issue",
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
    :param seq_table: the data frame of Gumbo Omics sequencing data
    :return: the samples data frame with optionally incremented version numbers
    """

    logging.info("Incrementing sample version numbers...")

    # count distinct profile IDs in existing seq table from among those observed in the
    # current batch of samples
    existing_profile_id_counts = (
        seq_table.drop_duplicates(subset="sequencing_id")
        .loc[
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


def upsert_terra_samples(
    terra_workspace: TerraWorkspace,
    samples: TypedDataFrame[SamplesWithMetadata],
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


def do_join_short_read_data(
    terra_workspace: TerraWorkspace,
    short_read_terra_workspace: TerraWorkspace,
    gumbo_client: GumboClient,
    dry_run: bool,
) -> None:
    """
    Join columns from the short read RNA workspace to the long read one by mapping on
    common model condition IDs.

    :param terra_workspace: the long read TerraWorkspace instance
    :param short_read_terra_workspace: the short read TerraWorkspace instance
    :param gumbo_client: an instance of the Gumbo GraphQL client
    :param dry_run: whether to skip updates to external data stores
    """

    # get the current long read sample data from Terra
    lr_samples = terra_workspace.get_entities("sample")
    lr_samples["model_id"] = lr_samples["participant"].apply(lambda x: x["entityName"])

    # get the current short read sample data from Terra
    sr_samples = short_read_terra_workspace.get_entities("sample")
    sr_samples = (
        sr_samples.loc[
            :,
            [
                "sample_id",
                "star_junctions",
                "fusion_predictions",
                "fusion_predictions_abridged",
                "rsem_genes",
                "rsem_genes_stranded",
                "rsem_isoforms",
                "rsem_isoforms_stranded",
            ],
        ]
        .rename(
            columns={
                "star_junctions": "sr_star_junctions",
                "fusion_predictions": "sr_fusion_predictions",
                "fusion_predictions_abridged": "sr_fusion_predictions_abridged",
                "rsem_genes": "sr_rsem_genes",
                "rsem_genes_stranded": "sr_rsem_genes_stranded",
                "rsem_isoforms": "sr_rsem_isoforms",
                "rsem_isoforms_stranded": "sr_rsem_isoforms_stranded",
            }
        )
        .astype("string")
    )

    # get the short read sequencing records from Gumbo
    models = model_to_df(gumbo_client.get_models_and_children(), ModelsAndChildren)
    sr_sequencings = explode_and_expand_models(models)
    sr_sequencings = sr_sequencings.loc[
        sr_sequencings["expected_type"].eq("rna"),
        [
            "model_condition_id",
            "profile_id",
            "sequencing_id",
            "datatype",
            "blacklist",
            "blacklist_omics",
            "source",
            "version",
        ],
    ].drop_duplicates(subset="sequencing_id")

    # match a short read sample to each long read sample
    sr_metadata = choose_matched_short_read_sample(
        sr_samples, lr_samples, sr_sequencings
    )

    sr_metadata = sr_metadata.merge(
        sr_samples,
        how="inner",
        left_on="sr_sample_id",
        right_on="sample_id",
    ).drop(columns="sample_id")

    assert bool(~sr_metadata["model_condition_id"].duplicated().any())

    # update long read sample data table in Terra with short read metadata and workflow
    # outputs
    samples_upsert = (
        lr_samples[["sample_id", "model_condition_id"]]
        .merge(sr_metadata, how="inner", on="model_condition_id")
        .drop(columns="model_condition_id")
    )

    if dry_run:
        logging.info(f"(skipping) Upserting {len(samples_upsert)} samples")
    else:
        terra_workspace.upload_entities(samples_upsert)


def choose_matched_short_read_sample(
    sr_samples: pd.DataFrame,
    lr_samples: pd.DataFrame,
    sr_sequencings: TypedDataFrame[SeqTable],
) -> TypedDataFrame[ShortReadMetadata]:
    """
    Choose short read RNA samples to map to long read ones.

    :param sr_samples: the short read workspace sample data table
    :param lr_samples: the long read workspace sample data table
    :param sr_sequencings: the data frame of Gumbo Omics sequencing data
    :return: a data frame mapping common model condition IDs to short read profile and
    sequencing IDs
    """

    # prioritize certain sample source
    source_priority = [
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

    # collect all valid short read samples in Gumbo that are present in both short and
    # long read workspaces
    sr = sr_sequencings.loc[
        sr_sequencings["sequencing_id"].isin(sr_samples["sample_id"])
        & sr_sequencings["model_condition_id"].isin(lr_samples["model_condition_id"])
        & ~sr_sequencings["blacklist"]
        & ~sr_sequencings["blacklist_omics"]
    ].copy()

    # populate the `source_priority` column
    sr["source"] = sr["source"].fillna("")
    sr = sr.merge(sp_df, how="left", on="source")
    assert sr["source_priority"].notna().all()

    # pick a short read sample for each model condition, using the more recent `version`
    # to break ties
    sr_choices = (
        sr.sort_values(
            ["model_condition_id", "source_priority", "version"],
            ascending=[True, True, False],
        )
        .groupby("model_condition_id")
        .nth(0)
    )

    sr_choices = sr_choices[
        ["model_condition_id", "profile_id", "sequencing_id"]
    ].rename(columns={"profile_id": "sr_profile_id", "sequencing_id": "sr_sample_id"})

    return type_data_frame(sr_choices, ShortReadMetadata)


def prep_onboarding_samples(
    onboarding_job: onboarding_job_insert_input,
    samples: TypedDataFrame[VersionedSamples],
) -> onboarding_job_insert_input:
    """
    Create related records for the onboarding job that will be uploaded, including ones
    for new omics_sequencings to be created.

    :param onboarding_job: the onboarding job Pydantic model instances for the workspace
    :param samples: the data frame of samples
    :return: the updated onboarding job
    """

    onboarding_samples = []

    for x in samples.to_dict(orient="records"):
        # make the onboarding_sample record (whether we're creating an omics_sequencing
        # record or not)
        onboarding_sample = onboarding_sample_insert_input.model_validate(
            {
                "terra_sample_id": x["sequencing_id"],
                "omics_profile_id": x["profile_id"],
                "issue": x["issue"],
            }
        )

        if x["issue"] is None:
            # new sample is valid, so create sequencing_alignment and omics_sequencing
            # records (model_validate will extract the fields from x it needs for each)
            onboarding_sample.sequencing_alignment = (
                sequencing_alignment_obj_rel_insert_input(
                    data=sequencing_alignment_insert_input.model_validate(x)
                )
            )

            onboarding_sample.sequencing_alignment.data.omics_sequencing = (
                omics_sequencing_obj_rel_insert_input(
                    data=omics_sequencing_insert_input.model_validate(x)
                )
            )

            onboarding_job.n_samples_succeeded += 1  # pyright: ignore

        else:
            onboarding_job.n_samples_failed += 1  # pyright: ignore

        onboarding_samples.append(onboarding_sample)

    onboarding_job.onboarding_samples = onboarding_sample_arr_rel_insert_input(
        data=onboarding_samples
    )

    return onboarding_job


def persist_onboarding_job(
    gumbo_client: GumboClient,
    onboarding_job: onboarding_job_insert_input,
    dry_run: bool,
) -> onboarding_job_insert_input:
    """
    Upload the onboarding job and related records to Gumbo.

    :param gumbo_client: an instance of the Gumbo GraphQL client
    :param onboarding_job: the onboarding job Pydantic model instances for the workspace
    :param dry_run: whether to skip updates to external data stores
    :return: the updated onboarding job
    """

    logging.info(
        "Onboarding job to persist: "
        + str(
            onboarding_job.model_dump(
                include={
                    "succeeded",
                    "n_samples",
                    "n_samples_new",
                    "n_samples_succeeded",
                    "n_samples_failed",
                }
            )
        )
    )

    if dry_run:
        logging.info(f"(skipping) Inserting onboarding job")
        return onboarding_job

    res = gumbo_client.insert_onboarding_job(
        gumbo_client.username, object=onboarding_job
    )

    job_id = res.insert_onboarding_job_one.id  # pyright: ignore
    logging.info(f"Inserted onboarding_job {job_id}")
    onboarding_job.id = job_id

    # update the corresponding `omics_profile` records' status field.
    update_profile_status(gumbo_client, onboarding_job)

    return onboarding_job


def update_profile_status(
    gumbo_client: GumboClient, onboarding_job: onboarding_job_insert_input
) -> None:
    """
    Update `omics_profile` records for samples that have been successfully created (or
    not).

    :param gumbo_client: an instance of the Gumbo GraphQL client
    :param onboarding_job: the onboarding job Pydantic model instances for the workspace
    """

    if onboarding_job.n_samples_succeeded > 0:  # pyright: ignore
        # set profiles' status to "Done"
        profile_ids = [
            str(x.omics_profile_id)
            for x in onboarding_job.onboarding_samples.data  # pyright: ignore
            if x.issue is None
        ]

        gumbo_client.set_status(gumbo_client.username, profile_ids, status="Done")

    if onboarding_job.n_samples_failed > 0:  # pyright: ignore
        # set profiles' status to a relevant failure term
        profile_ids_fail_small = [
            str(x.omics_profile_id)
            for x in onboarding_job.onboarding_samples.data  # pyright: ignore
            if x.issue is not None and x.issue == "BAM/CRAM file too small"
        ]

        if len(profile_ids_fail_small) > 0:
            gumbo_client.set_status(
                gumbo_client.username,
                profile_ids_fail_small,
                status="Fail - CDS file size",
            )

        profile_ids_fail_sm_id = [
            str(x.omics_profile_id)
            for x in onboarding_job.onboarding_samples.data  # pyright: ignore
            if x.issue is not None and x.issue == "unmappable SM ID"
        ]

        if len(profile_ids_fail_sm_id) > 0:
            gumbo_client.set_status(
                gumbo_client.username,
                profile_ids_fail_sm_id,
                status="Fail - unmappable SM ID",
            )

        profile_ids_fail_other = [
            str(x.omics_profile_id)
            for x in onboarding_job.onboarding_samples.data  # pyright: ignore
            if x.issue is not None
            and x.issue != "BAM/CRAM file too small"
            and x.issue != "unmappable SM ID"
        ]

        if len(profile_ids_fail_other) > 0:
            gumbo_client.set_status(
                gumbo_client.username,
                profile_ids_fail_other,
                status="Fail - CDS QC Other",
            )


def send_slack_message(
    onboarding_job: onboarding_job_insert_input, dry_run: bool
) -> None:
    """
    Send a message to a Slack channel with information about the onboarding job.

    :param onboarding_job: the onboarding job
    :param dry_run: whether to skip sending the Slack message
    """

    if dry_run:
        logging.info("(skipping) Sending results to Slack channel...")
        return

    logging.info("Sending results to Slack channel...")

    # get the Slack webhook URL for the results channel
    webhook_url = get_secret_from_sm(
        "projects/201811582504/secrets/onboarding-webhook-url/versions/latest"
    )

    # construct URLs for pre-filtered onboarding_sample records
    base_url = (
        "https://app.forestadmin.com"
        "/GumboUI/Production/DepMap/data/onboarding_sample/index"
    )

    blocks = []

    filter_params = {
        "type": "and",
        "conditions": [
            {
                "operator": "is",
                "value": onboarding_job.id,
                "fieldName": "onboarding_job",
                "subFieldName": "id",
                "embeddedFieldName": None,
            }
        ],
    }

    filter_encoded = quote_plus(json.dumps(filter_params))
    url = f"{base_url}?filter={filter_encoded}"

    blocks.append(
        {
            "type": "rich_text",
            "elements": [
                {
                    "type": "rich_text_section",
                    "elements": [
                        {
                            "type": "emoji",
                            "name": "white_check_mark"
                            if onboarding_job.n_samples_failed == 0
                            else "exclamation",
                        },
                        {
                            "type": "link",
                            "text": onboarding_job.onboarding_workspace_id,
                            "url": url,
                        },
                    ],
                }
            ],
        }
    )

    fields = [
        {"type": "mrkdwn", "text": "*New*"},
        {"type": "plain_text", "text": str(onboarding_job.n_samples_new)},
    ]

    if onboarding_job.n_samples_succeeded > 0:  # pyright: ignore
        fields.extend(
            [
                {"type": "mrkdwn", "text": "*Succeeded*"},
                {"type": "plain_text", "text": str(onboarding_job.n_samples_succeeded)},
            ]
        )

    if onboarding_job.n_samples_failed > 0:  # pyright: ignore
        fields.extend(
            [
                {"type": "mrkdwn", "text": "*Failed*"},
                {"type": "plain_text", "text": str(onboarding_job.n_samples_failed)},
            ]
        )

    blocks.append({"type": "section", "fields": fields})

    r = requests.post(webhook_url, data=json.dumps({"blocks": blocks}))
    r.raise_for_status()
