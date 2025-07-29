import tomllib

import pandas as pd
from nebelung.terra_workspace import TerraWorkspace
from nebelung.utils import type_data_frame

from depmap_omics_long_read_rna.types import DeliveryBams
from depmap_omics_long_read_rna.utils.delivery_bams import assign_cds_ids
from depmap_omics_long_read_rna.utils.gcp import get_objects_metadata

pd.set_option("display.max_columns", 30)
pd.set_option("display.max_colwidth", 50)
pd.set_option("display.max_info_columns", 30)
pd.set_option("display.max_info_rows", 20)
pd.set_option("display.max_rows", 20)
pd.set_option("display.max_seq_items", None)
pd.set_option("display.width", 200)
pd.set_option("expand_frame_repr", True)
pd.set_option("mode.chained_assignment", "warn")

# use same config loading as when calling the module CLI
with open("config.toml", "rb") as f:
    config = tomllib.load(f)

terra_delivery_workspace = TerraWorkspace(
    workspace_namespace=config["terra"]["delivery_workspace_namespace"],
    workspace_name=config["terra"]["delivery_workspace_name"],
)

"""
SELECT
    omics_profile.id,
    omics_profile.smid_ordered,
    omics_profile.model_condition_id,
    omics_profile.sm_id_matched,
    model_condition.model_id,
    model.cell_line_name,
    model.stripped_cell_line_name,
    model.growth_pattern,
    model.ccle_name,
    model.cell_line_aliases,
    depmap_model_type.subtype
FROM
    omics_profile
        INNER JOIN model_condition ON omics_profile.model_condition_id = model_condition.id
        INNER JOIN model ON model_condition.model_id = model.id
        INNER JOIN depmap_model_type ON model.depmap_model_type_id = depmap_model_type.id
WHERE
    omics_profile.id IN ('PR-qx8jgr', 'PR-dg3WGe', 'PR-TDyUCV', 'PR-uh2aqm', 'PR-QNteGy', 'PR-x2MpVH', 'PR-rad8fr', 'PR-J7KMX7', 'PR-HhsB29', 'PR-P6ePwA', 'PR-ANwkaE', 'PR-TJRxuv', 'PR-xJDUEA', 'PR-WP8rrP', 'PR-es46WJ', 'PR-uBh2jB', 'PR-3m8DpR', 'PR-G6NVZP', 'PR-EWGEfV', 'PR-gnKZ6h', 'PR-s3PXZd', 'PR-xJDGKV', 'PR-QNdycn', 'PR-uBmmCn', 'PR-e7WZgh', 'PR-Zt6vPz', 'PR-BM8Xfs', 'PR-c7BZjk', 'PR-gncaFN', 'PR-bDZ8f9', 'PR-cCEb2X', 'PR-8defUr', 'PR-BN39kt', 'PR-zDVcme', 'PR-TPeXbb', 'PR-DyEzaK', 'PR-8ufpPy', 'PR-kqnUhV', 'PR-NFvH9r', 'PR-g6r4pp', 'PR-4hHvvZ', 'PR-8yxHqx', 'PR-pGHcyY', 'PR-VJQbkf', 'PR-NH4CRN', 'PR-JJHZWm', 'PR-PVpzMb', 'PR-ZthTvx', 'PR-8fBxMm', 'PR-Ueatdv', 'PR-DsKEB4', 'PR-hBHRvg', 'PR-Jnayyu', 'PR-BzKhfE', 'PR-gAbz3w', 'PR-CTmE4r', 'PR-wJpnmw', 'PR-Md4MhQ', 'PR-8PhQXY', 'PR-wrnwtD', 'PR-xZDbTj', 'PR-Mdmep4', 'PR-nNRvzg', 'PR-Mrv7KF', 'PR-paJDNr', 'PR-4TX4HP', 'PR-NCsNCz', 'PR-BAVuf3', 'PR-tcyuUs', 'PR-bE6WTC', 'PR-Ya3pwR', 'PR-8TWXuR', 'PR-zJVmvd', 'PR-7wPPaV', 'PR-txzd88', 'PR-xf48KE', 'PR-uApUtm', 'PR-TaEGx2', 'PR-b6jE42', 'PR-WKF46v', 'PR-peyZeM', 'PR-uxnfg3', 'PR-J3NnJA', 'PR-aAUwXm', 'PR-G73Je6', 'PR-bRzbgW', 'PR-UVPGmp', 'PR-pUVZN3', 'PR-Y3d7p6', 'PR-pkX2pT', 'PR-ahbNrp', 'PR-RG722P', 'PR-eNXKFg', 'PR-xJ4Bqe', 'PR-d2nasX', 'PR-qBxsCC', 'PR-ZGkPBA', 'PR-WZQymY', 'PR-jGvwFu', 'PR-BUFjQa');
"""

profiles = pd.read_csv("./data/cclf100profiles.tsv", sep="\t")
inventory = pd.read_csv("./data/CCLF_SU2C_First100_GEA_LongReadRNAseq.tsv", sep="\t")
mapping = pd.read_csv("./data/first100_SU2C_longreadrnaseq_with_names.tsv", sep="\t")

df = (
    profiles.merge(
        inventory,
        how="inner",
        left_on="sm_id_matched",
        right_on="Cell Line Name",
    )
    .merge(mapping, how="inner", left_on="sm_id_matched", right_on="cell_line_name")
    .loc[
        :,
        [
            "model_id",
            "cell_line_name_x",
            "stripped_cell_line_name",
            "model_condition_id",
            "id",
            "smid_ordered",
            "sm_id_matched",
            "growth_pattern",
            "cell_line_aliases",
            "subtype",
            "Culture system",
            "Flask Conditions",
            "Culture media",
            "bam",
        ],
    ]
    .rename(
        columns={
            "cell_line_name_x": "cell_line_name",
            "id": "omics_profile_id",
            "Culture system": "culture_system",
            "Flask Conditions": "flask_conditions",
            "Culture media": "culture_media",
        }
    )
)

df["file_id"] = df["bam"].str.extract(r"\/([^/]+).merged")

blobs = get_objects_metadata(df["bam"])
print(blobs["size"].duplicated().any())

bams = type_data_frame(
    df.merge(blobs, how="inner", left_on="bam", right_on="url")
    .loc[:, ["omics_profile_id", "url", "crc32c", "size", "gcs_obj_updated_at"]]
    .rename(
        columns={
            "url": "delivery_bam",
            "crc32c": "delivery_bam_crc32c",
            "size": "delivery_bam_size",
            "gcs_obj_updated_at": "delivery_bam_updated_at",
        }
    ),
    DeliveryBams,
)

bams = assign_cds_ids(bams, config["uuid_namespace"])

sample_ids = bams.pop("cds_id")
bams.insert(0, "entity:sample_id", sample_ids)

existing_samples = terra_delivery_workspace.get_entities("sample")
bams = bams.loc[
    ~bams["delivery_bam_crc32c"].isin(existing_samples["delivery_bam_crc32c"])
]

if len(bams) > 0:
    terra_delivery_workspace.upload_entities(bams)
