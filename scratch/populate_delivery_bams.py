import os

import pandas as pd
from dotenv import load_dotenv
from nebelung.terra_workspace import TerraWorkspace

from dogspa_long_reads.types import GumboClient, ModelsAndChildren
from dogspa_long_reads.utils.bams import assign_cds_ids
from dogspa_long_reads.utils.onboarding import explode_and_expand_models
from dogspa_long_reads.utils.utils import model_to_df

load_dotenv()

terra_workspace = TerraWorkspace(
    workspace_namespace="broad-firecloud-ccle", workspace_name="Long_Read_Omics"
)

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

delivery_bams["delivery_bam_id"] = delivery_bams["delivery_bam"].apply(os.path.basename)
delivery_bams = assign_cds_ids(delivery_bams, "00000000-0000-0000-0000-000000000000")

aligned_bams = (
    pd.read_csv("./data/legacy_aligned_bams.tsv", sep="\t")
    .dropna()
    .rename(
        columns={"minimap2_bam": "aligned_bam", "minimap2_bam_index": "aligned_bai"}
    )
)

last_batch = {
    "ACH-000040",
    "ACH-000212",
    "ACH-000246",
    "ACH-000261",
    "ACH-000280",
    "ACH-000288",
    "ACH-000455",
    "ACH-000622",
    "ACH-000634",
    "ACH-000783",
    "ACH-000906",
    "ACH-000913",
    "ACH-000946",
    "ACH-001067",
    "ACH-001329",
    "ACH-001386",
    "ACH-001415",
    "ACH-001538",
    "ACH-001606",
    "ACH-001607",
    "ACH-001608",
    "ACH-001609",
    "ACH-001611",
    "ACH-001619",
    "ACH-001622",
    "ACH-001623",
    "ACH-001652",
    "ACH-001673",
    "ACH-001687",
    "ACH-001843",
    "ACH-001844",
    "ACH-001850",
    "ACH-001853",
    "ACH-001960",
    "ACH-001961",
    "ACH-002070",
    "ACH-002533",
    "ACH-002680",
    "ACH-002681",
    "ACH-003130",
    "ACH-003193",
    "ACH-003233",
    "ACH-003235",
    "ACH-003237",
    "ACH-003337",
}

aligned_bams = aligned_bams.loc[~aligned_bams["model_id"].isin(last_batch)]

delivery_bams = delivery_bams.merge(aligned_bams, how="left", on="model_id")

assert ~delivery_bams["model_id"].duplicated().any()

bam_ids = delivery_bams.pop("delivery_bam_id")
delivery_bams.insert(0, "entity:delivery_bam_id", bam_ids)

terra_workspace.upload_entities(delivery_bams)
