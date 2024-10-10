import os

import pandas as pd
from dotenv import load_dotenv
from nebelung.terra_workspace import TerraWorkspace

from dogspa_long_reads.types import GumboClient, ModelsAndChildren
from dogspa_long_reads.utils.onboarding import explode_and_expand_models
from dogspa_long_reads.utils.utils import model_to_df

load_dotenv()

terra_workspace = TerraWorkspace(
    workspace_namespace="broad-firecloud-ccle", workspace_name="Long_Read_Omics"
)

terra_samples = terra_workspace.get_entities("sample")

pilot = pd.read_csv("./data/pilot.csv")
pb_pre = pd.read_csv("./data/pb_pre.csv")
tag = pd.read_csv("./data/tag.csv")

gumbo_client = GumboClient(
    url=os.environ["HASURA_URL"],
    username="snp_str_qc",
    headers={"X-Hasura-Admin-Secret": os.environ["HASURA_ADMIN_SECRET"]},
)

models = model_to_df(gumbo_client.get_models_and_children(), ModelsAndChildren)
models = explode_and_expand_models(models)
models = models.loc[models["datatype"].eq("long_read_rna")]
models["smid_ordered"] = models["smid_ordered"].str.strip()

delivery_bams_pilot = pilot.merge(
    models.loc[models["smid_ordered"].isna()], how="left", on="model_id"
)

assert ~delivery_bams_pilot["model_id"].duplicated().any()
assert len(delivery_bams_pilot) == len(pilot)

delivery_bams_pb_pre_wdups = pb_pre.merge(
    models.rename(columns={"smid_ordered": "sm_id"}),
    how="left",
    on=["model_id"],
)

delivery_bams_pb_pre = delivery_bams_pb_pre_wdups.loc[
    delivery_bams_pb_pre_wdups["sm_id_y"].isna()
    | delivery_bams_pb_pre_wdups["sm_id_x"].eq(delivery_bams_pb_pre_wdups["sm_id_y"])
]

delivery_bams_pb_pre_bad_smid = delivery_bams_pb_pre_wdups.loc[
    ~delivery_bams_pb_pre_wdups["bam_url"].isin(delivery_bams_pb_pre["bam_url"])
]

delivery_bams_pb_pre = pd.concat([delivery_bams_pb_pre, delivery_bams_pb_pre_bad_smid])

assert ~delivery_bams_pb_pre["model_id"].duplicated().any()
assert len(delivery_bams_pb_pre) == len(pb_pre)

delivery_bams_tag = tag.merge(
    models.dropna(subset="smid_ordered"), how="left", on="model_id"
)

assert ~delivery_bams_tag["model_id"].duplicated().any()
assert len(delivery_bams_tag) == len(tag)

delivery_bams = pd.concat(
    [
        delivery_bams_pilot[
            [
                "bam_url",
                "crc32c",
                "size",
                "gcs_obj_updated_at",
                "model_id",
                "profile_id",
            ]
        ],
        delivery_bams_pb_pre[
            [
                "bam_url",
                "crc32c",
                "size",
                "gcs_obj_updated_at",
                "model_id",
                "profile_id",
            ]
        ],
        delivery_bams_tag[
            [
                "bam_url",
                "crc32c",
                "size",
                "gcs_obj_updated_at",
                "model_id",
                "profile_id",
            ]
        ],
    ]
).rename(
    columns={
        "bam_url": "delivery_bam",
        "crc32c": "delivery_bam_crc32c",
        "size": "delivery_bam_size",
        "gcs_obj_updated_at": "delivery_bam_updated_at",
    }
)

assert ~delivery_bams["profile_id"].duplicated().any()

delivery_bams = delivery_bams.merge(
    terra_samples[["sample_id", "LongReadProfileID"]].rename(
        columns={"LongReadProfileID": "profile_id"}
    ),
    how="left",
    on="profile_id",
)

assert delivery_bams["sample_id"].notna().all()

aligned_bams = (
    pd.read_csv("./data/legacy_aligned_bams.tsv", sep="\t")
    .dropna()
    .rename(
        columns={"minimap2_bam": "aligned_bam", "minimap2_bam_index": "aligned_bai"}
    )
)

delivery_bams = delivery_bams.merge(aligned_bams, how="left", on="model_id")

assert ~delivery_bams["model_id"].duplicated().any()

delivery_bams["delivery_bam_id"] = delivery_bams["delivery_bam"].apply(os.path.basename)
bam_ids = delivery_bams.pop("delivery_bam_id")
delivery_bams.insert(0, "entity:delivery_bam_id", bam_ids)

terra_workspace.upload_entities(delivery_bams)
