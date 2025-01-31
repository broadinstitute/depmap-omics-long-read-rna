import json
import logging

import pandas as pd
from nebelung.terra_workspace import TerraWorkspace
from nebelung.utils import type_data_frame
from pandera.typing import DataFrame as TypedDataFrame

from depmap_omics_long_read_rna.types import (
    GumboClient,
    ModelsAndChildren,
    SamplesWithMetadata,
    VersionedSamples,
)
from depmap_omics_long_read_rna.utils.gcp import (
    copy_to_cclebams,
    get_objects_metadata,
    update_sample_file_urls,
)
from depmap_omics_long_read_rna.utils.metadata import (
    choose_matched_short_read_sample,
    explode_and_expand_models,
)
from depmap_omics_long_read_rna.utils.utils import model_to_df
from gumbo_gql_client import (
    omics_sequencing_insert_input,
    omics_sequencing_obj_rel_insert_input,
    sequencing_alignment_insert_input,
    sequencing_alignment_obj_rel_insert_input,
)


def prep_onboarding_job(
    gcp_project_id: str,
    aligned_gcs_destination_bucket: str,
    aligned_gcs_destination_prefix: str,
    terra_workspace: TerraWorkspace,
    gumbo_client: GumboClient,
    dry_run: bool = False,
) -> None:
    """
    Onboard the latest uBAMs and aligned BAMs to Gumbo and Terra.

    :param gcp_project_id: a GCP project ID to use for billing
    :param aligned_gcs_destination_bucket: the GCS bucket name for aligned BAMs
    :param aligned_gcs_destination_prefix: the GCS bucket prefix for aligned BAMs
    :param terra_workspace: a TerraWorkspace instance
    :param gumbo_client: an instance of the Gumbo GraphQL client
    :param dry_run: whether to skip updates to external data stores
    """

    # get the sequencing and profile tables from Gumbo
    models = model_to_df(gumbo_client.get_models_and_children(), ModelsAndChildren)
    samples = explode_and_expand_models(models)
    samples_lr = samples.loc[samples["datatype"].eq("long_read_rna")]

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
        samples,
        aligned_sample_files,
        bam_bai_colnames=["aligned_bam", "aligned_bai"],
    )

    missing_files = pd.Series(
        samples[["aligned_bam", "aligned_bai"]].isna().any(axis=1)
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

    unaligned_samples = (
        samples.loc[
            :,
            [
                "sequencing_id",
                "profile_id",
                "issue",
                "unaligned_bam_url",
                "unaligned_bam_crc32c_hash",
                "unaligned_bam_size",
            ],
        ]
        .rename(
            columns={
                "sequencing_id": "terra_sample_id",
                "unaligned_bam_url": "url",
                "unaligned_bam_crc32c_hash": "crc32c_hash",
                "unaligned_bam_size": "size",
            }
        )
        .assign(sequencing_alignment_source="GP", reference_genome=None)
    )

    aligned_samples = (
        samples.loc[
            :,
            [
                "sequencing_id",
                "profile_id",
                "issue",
                "bam_url",
                "bam_crc32c_hash",
                "bam_size",
            ],
        ]
        .rename(
            columns={
                "sequencing_id": "terra_sample_id",
                "bam_url": "url",
                "bam_crc32c_hash": "crc32c_hash",
                "bam_size": "size",
            }
        )
        .assign(sequencing_alignment_source="CDS", reference_genome="hg38")
    )

    samples_long = pd.concat([unaligned_samples, aligned_samples])

    onboarding_samples = []

    for x in samples_long.to_dict(orient="records"):
        # make the onboarding_sample record (whether we're creating an omics_sequencing
        # record or not)
        onboarding_sample = onboarding_sample_insert_input.model_validate(
            {
                "terra_sample_id": x["terra_sample_id"],
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

            if x["reference_genome"] is None:
                # don't double count because we're technically onboarding both the
                # unaligned, GP-delivered file and the aligned CDS one
                onboarding_job.n_samples_succeeded += 1  # pyright: ignore

        elif x["reference_genome"] is None:
            onboarding_job.n_samples_failed += 1  # pyright: ignore

        onboarding_samples.append(onboarding_sample)

    onboarding_job.onboarding_samples = onboarding_sample_arr_rel_insert_input(
        data=onboarding_samples
    )

    return onboarding_job
