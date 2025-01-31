import logging
import os
from pathlib import Path

import pandas as pd
from nebelung.terra_workflow import TerraWorkflow
from nebelung.terra_workspace import TerraWorkspace
from nebelung.utils import type_data_frame
from pandera.typing import DataFrame as TypedDataFrame

from depmap_omics_long_read_rna.types import (
    DeliveryBams,
    ObjectMetadata,
    SamplesWithCDSIDs,
)
from depmap_omics_long_read_rna.utils.gcp import list_blobs
from depmap_omics_long_read_rna.utils.utils import assign_hashed_uuids, uuid_to_base62


def do_upsert_delivery_bams(
    gcs_source_bucket: str,
    gcs_source_glob: str,
    uuid_namespace: str,
    terra_workspace: TerraWorkspace,
    dry_run: bool,
) -> None:
    """
    Search for uBAMs delivered to a GCS bucket by GP and populate a `sample` data table
    on Terra.

    :param gcs_source_bucket: the GCS bucket where GP delivers uBAMs
    :param gcs_source_glob: a glob expression to search for uBAMs in the bucket
    :param uuid_namespace: a namespace for generated UUIDv3s
    :param terra_workspace: a TerraWorkspace instance
    :param dry_run: whether to skip updates to external data stores
    """

    # get delivered (u)BAM file metadata
    src_bams = list_blobs(gcs_source_bucket, glob=gcs_source_glob)
    bams = make_delivery_bam_df(src_bams)

    # generate sample/sequencing/CDS IDs
    bams = assign_cds_ids(bams, uuid_namespace)

    # upsert to Terra data table
    bam_ids = bams.pop("cds_id")
    bams.insert(0, "entity:sample_id", bam_ids)

    if dry_run:
        logging.info(f"(skipping) Upserting {len(bams)} delivery_bams")
        return

    terra_workspace.upload_entities(bams)


def make_delivery_bam_df(
    bams: TypedDataFrame[ObjectMetadata],
) -> TypedDataFrame[DeliveryBams]:
    """
    Prepare data frame for the Terra delivery_bam data table using inventory of GCS
    blobs.

    :param bams: a data frame of delivered uBAMs
    :return: a data frame of delivered uBAM metadata
    """

    bams_w_ids = bams.copy()

    # extract the Gumbo model ID from the BAM filenames
    bams_w_ids["model_id"] = (
        bams_w_ids["url"]
        .apply(lambda x: Path(x).name)
        .str.extract(r"^(ACH-[A-Z 0-9]+)")
    )

    if bool(bams_w_ids["model_id"].isna().any()):
        raise ValueError("There are BAM files not named with model IDs (ACH-*)")

    bams_w_ids = bams_w_ids.rename(columns={"url": "bam_url"})

    bams_w_ids["bam_id"] = bams_w_ids["bam_url"].apply(os.path.basename)

    bams_w_ids = bams_w_ids.rename(
        columns={
            "bam_url": "delivery_bam",
            "crc32c": "delivery_bam_crc32c",
            "size": "delivery_bam_size",
            "gcs_obj_updated_at": "delivery_bam_updated_at",
        }
    )

    return type_data_frame(bams_w_ids, DeliveryBams)


def assign_cds_ids(
    samples: TypedDataFrame[DeliveryBams], uuid_namespace: str
) -> TypedDataFrame[SamplesWithCDSIDs]:
    """
    Assign a "CDS-" ID to each sample by hashing relevant columns.

    :param samples: the data frame of samples
    :param uuid_namespace: a namespace for generated UUIDv3s
    :return: the data frame with CDS IDs
    """

    logging.info("Assigning CDS IDs...")

    samples_w_ids = assign_hashed_uuids(
        samples,
        uuid_namespace=uuid_namespace,
        uuid_col_name="cds_id",
        subset=["delivery_bam_crc32c", "delivery_bam_size", "delivery_bam_updated_at"],
    )

    # convert UUIDs to base62 and truncate to match existing ID format
    samples_w_ids["cds_id"] = "CDS-" + samples_w_ids["cds_id"].apply(
        uuid_to_base62
    ).str.slice(-6)

    return type_data_frame(samples_w_ids, SamplesWithCDSIDs)


def do_delta_align_delivery_bams(
    terra_workspace: TerraWorkspace, terra_workflow: TerraWorkflow, dry_run: bool
) -> None:
    """
    Identify delivered uBAMs that haven't been aligned yet and submit a job to do that.

    :param terra_workspace: a TerraWorkspace instance
    :param terra_workflow: a TerraWorkflow instance for the alignment method
    :param dry_run: whether to skip updates to external data stores
    """

    bams = terra_workspace.get_entities("sample", DeliveryBams)

    if "aligned_bam" not in bams.columns:
        bams["aligned_bam"] = pd.NA

    bams_to_align = bams.loc[bams["aligned_bam"].isna()]

    if len(bams_to_align) == 0:
        logging.info("No new BAMs to align")
        return

    # get statuses of submitted entity workflow statuses
    submittable_entities = terra_workspace.check_submittable_entities(
        entity_type="sample",
        entity_ids=bams_to_align["sample_id"],
        terra_workflow=terra_workflow,
        resubmit_n_times=1,
        force_retry=False,
    )

    logging.info(f"Submittable entities: {submittable_entities}")

    if len(submittable_entities["failed"]) > 0:
        raise RuntimeError("Some entities have failed too many times")

    # don't submit jobs for entities that are currently running, completed, or failed
    # too many times
    bams_to_align = bams_to_align.loc[
        bams_to_align["sample_id"].isin(
            submittable_entities["unsubmitted"].union(submittable_entities["retryable"])
        )
    ]

    if dry_run:
        logging.info("(skipping) Submitting align_long_reads job")
        return

    bam_set_id = terra_workspace.create_entity_set(
        entity_type="sample",
        entity_ids=bams_to_align["sample_id"],
        suffix="align",
    )

    terra_workspace.submit_workflow_run(
        terra_workflow=terra_workflow,
        entity=bam_set_id,
        etype="sample_set",
        expression="this.samples",
        use_callcache=True,
        use_reference_disks=False,
        memory_retry_multiplier=1.5,
    )
