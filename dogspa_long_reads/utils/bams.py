import logging
import os
from pathlib import Path

import pandas as pd
from nebelung.terra_workflow import TerraWorkflow
from nebelung.terra_workspace import TerraWorkspace
from nebelung.utils import type_data_frame
from pandera.typing import DataFrame as TypedDataFrame

from dogspa_long_reads.types import ObjectMetadata, TerraBams
from dogspa_long_reads.utils.gcp import list_blobs


def do_upsert_delivery_bams(
    gcs_source_bucket: str, gcs_source_glob: str, terra_workspace: TerraWorkspace
) -> None:
    legacy_bams_csv = Path(
        os.path.join(
            Path(os.path.dirname(__file__)).parent.parent, "data", "legacy_bams.csv"
        )
    ).resolve()

    # get delivered BAM file metadata
    src_bams = list_blobs(gcs_source_bucket, glob=gcs_source_glob)
    bams = id_bams(src_bams, legacy_bams_csv)

    upsert_bams = bams.rename(columns={"bam_id": "entity:bam_id"})[
        ["entity:bam_id", "bam", "model_id", "size", "crc32c", "gcs_obj_updated_at"]
    ]

    terra_workspace.upload_entities(upsert_bams)


def id_bams(
    bams: TypedDataFrame[ObjectMetadata], legacy_bams_f: Path
) -> TypedDataFrame[TerraBams]:
    bams_w_ids = bams.copy()

    bams_w_ids["model_id"] = (
        bams_w_ids["url"]
        .apply(lambda x: Path(x).name)
        .str.extract(r"^(ACH-[A-Z 0-9]+)")
    )

    if bool(bams_w_ids["model_id"].isna().any()):
        raise ValueError("There are BAM files not named with model IDs (ACH-*)")

    bams_w_ids = bams_w_ids.rename(columns={"url": "bam_url"})

    legacy_bams = pd.read_csv(legacy_bams_f)
    bams_w_ids = pd.concat([bams_w_ids, legacy_bams])

    bams_w_ids["bam_id"] = bams_w_ids["bam_url"].apply(os.path.basename)
    assert ~bams_w_ids["bam_id"].duplicated().any()

    bams_w_ids = bams_w_ids.rename(columns={"bam_url": "bam"})

    return type_data_frame(bams_w_ids, TerraBams)


def do_delta_index_delivery_bams(
    terra_workspace: TerraWorkspace, terra_workflow: TerraWorkflow
) -> None:
    bams = terra_workspace.get_entities("bam", TerraBams)

    if "bai" not in bams.columns:
        bams["bai"] = pd.NA

    bams_to_index = bams.loc[bams["bai"].isna()]

    if len(bams_to_index) == 0:
        logging.info("No new BAMs to index")
        return

    bam_set_id = terra_workspace.create_entity_set(
        entity_type="bam", entity_ids=bams_to_index["bam_id"], suffix="index_bam"
    )

    terra_workspace.submit_workflow_run(
        terra_workflow,
        entity=bam_set_id,
        etype="bam_set",
        expression="this.bams",
        use_callcache=True,
        use_reference_disks=False,
        memory_retry_multiplier=1.2,
    )
