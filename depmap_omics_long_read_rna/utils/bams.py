import logging
import os
from pathlib import Path

import pandas as pd
from firecloud import api as firecloud_api
from nebelung.terra_workspace import TerraWorkspace
from nebelung.utils import call_firecloud_api, type_data_frame
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
    uuid_namespace: str,
    gcs_source_glob: str,
    terra_workspace: TerraWorkspace,
) -> None:
    # get delivered BAM file metadata
    src_bams = list_blobs(gcs_source_bucket, glob=gcs_source_glob)
    bams = id_bams(src_bams)

    # generate sample/sequencing/CDS IDs
    bams = assign_cds_ids(bams, uuid_namespace)

    # upsert to Terra data table
    bam_ids = bams.pop("cds_id")
    bams.insert(0, "entity:delivery_bam_id", bam_ids)
    terra_workspace.upload_entities(bams)


def id_bams(bams: TypedDataFrame[ObjectMetadata]) -> TypedDataFrame[DeliveryBams]:
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
    :param uuid_namespace: a namespace for UUIDv3 IDs
    :return: the data frame with CDS IDs
    """

    logging.info("Assigning CDS IDs...")

    samples_w_ids = assign_hashed_uuids(
        samples,
        uuid_namespace=uuid_namespace,
        uuid_col_name="cds_id",
        subset=[
            "delivery_bam_crc32c",
            "delivery_bam_size",
            "delivery_bam_updated_at",
        ],
    )

    # convert UUIDs to base62 and truncate to match existing ID format
    samples_w_ids["cds_id"] = "CDS-" + samples_w_ids["cds_id"].apply(
        uuid_to_base62
    ).str.slice(-6)

    return type_data_frame(samples_w_ids, SamplesWithCDSIDs)


def do_delta_align_delivery_bams(terra_workspace: TerraWorkspace) -> None:
    bams = terra_workspace.get_entities("delivery_bam", DeliveryBams)

    if "aligned_bam" not in bams.columns:
        bams["aligned_bam"] = pd.NA

    bams_to_align = bams.loc[bams["aligned_bam"].isna()]

    if len(bams_to_align) == 0:
        logging.info("No new BAMs to align")
        return

    bam_set_id = terra_workspace.create_entity_set(
        entity_type="delivery_bam",
        entity_ids=bams_to_align["delivery_bam_id"],
        suffix="align",
    )

    call_firecloud_api(
        firecloud_api.create_submission,
        wnamespace=terra_workspace.workspace_namespace,
        workspace=terra_workspace.workspace_name,
        cnamespace=terra_workspace.workspace_namespace,  # happens to be the same
        config="align_long_reads",
        entity=bam_set_id,
        etype="delivery_bam_set",
        expression="this.delivery_bams",
        use_callcache=True,
        use_reference_disks=False,
        memory_retry_multiplier=1.5,
    )
