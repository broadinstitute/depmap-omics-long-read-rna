import logging
import os
from pathlib import Path

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
    src_bams = list_blobs(bucket_name=gcs_source_bucket, glob=gcs_source_glob)
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

    # extract the profile ID from the BAM filenames
    bams_w_ids["omics_profile_id"] = (
        bams_w_ids["url"]
        .apply(lambda x: Path(x).name)
        .str.extract(r"^(PR-[A-Z a-z 0-9]+)")
    )

    if bool(bams_w_ids["omics_profile_id"].isna().any()):
        raise ValueError("There are BAM files not named with profile IDs (PR-*.bam)")

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
