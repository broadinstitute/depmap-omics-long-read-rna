import re
from io import StringIO

import pandas as pd
from google.cloud import storage

from depmap_omics_long_read_rna.utils.gcp import list_blobs

pd.set_option("display.max_columns", 30)
pd.set_option("display.max_colwidth", 100)
pd.set_option("display.max_info_columns", 30)
pd.set_option("display.max_info_rows", 20)
pd.set_option("display.max_rows", 20)
pd.set_option("display.max_seq_items", None)
pd.set_option("display.width", 200)
pd.set_option("expand_frame_repr", True)
pd.set_option("mode.chained_assignment", "warn")


storage_client = storage.Client()


bams = list_blobs(
    "fc-aaf4de93-c104-45c4-a01a-a036869119c6",
    glob="PB_preprocess/SM-*/merge/*.sorted.bam",
)

csvs = list_blobs(
    "fc-aaf4de93-c104-45c4-a01a-a036869119c6", glob="PB_preprocess/SM-*/*.csv"
)


def read_gcs_csv(url: str) -> pd.DataFrame:
    blob = storage.Blob.from_string(url, client=storage_client)
    csv_str = blob.download_as_string().decode()
    return pd.read_csv(StringIO(csv_str))


metadata = pd.concat([read_gcs_csv(x) for x in csvs["url"]])
metadata.columns = [
    "sample_id",
    "kinnex_adapter",
    "isoseq_primer",
    "ach_id",
    "sm_id",
    "rin_score",
    "rin",
]
metadata = (
    metadata.drop(columns=["kinnex_adapter", "isoseq_primer", "rin_score", "rin"])
    .drop_duplicates()
    .rename(columns={"ach_id": "model_id", "batch_id": "sm_id"})
)

bams["sm_id"] = bams["url"].str.extract(r"/(SM-[^/]+)/")
bams["sample_id"] = bams["url"].str.extract(r".+/([^./]+)")


def canonize(x: str) -> str:
    y = x.lower()
    y = re.sub(r"\.[0-9]+$", "", y)
    y = re.sub(r"[-_ /]", "", y)
    return y


bams["sample_id_canon"] = bams["sample_id"].apply(canonize)
metadata["sample_id_canon"] = metadata["sample_id"].apply(canonize)

bams_dedup = (
    bams.sort_values("gcs_obj_updated_at", ascending=False)
    .groupby("sample_id_canon")
    .nth(0)
)

metadata_dedup = metadata.drop_duplicates(subset="sample_id_canon")

df = bams_dedup.merge(
    metadata_dedup, on=["sm_id", "sample_id_canon"], how="outer", indicator=True
)

unmapped_sample_id_canon = df.loc[
    df["sample_id_canon"].duplicated(keep=False), "sample_id_canon"
]

df = df.loc[df["_merge"].eq("both")]

df2 = bams.merge(metadata, on=["sm_id", "sample_id_canon"], how="outer", indicator=True)
df2 = df2.loc[df2["sample_id_canon"].isin(unmapped_sample_id_canon)]
df2_dedup = (
    df2.sort_values("gcs_obj_updated_at", ascending=False)
    .groupby("sample_id_canon")
    .nth(0)
)

pb_pre = pd.concat([df, df2_dedup]).rename(columns={"url": "bam_url"})
pb_pre["size"] = pb_pre["size"].astype("int64")
pb_pre = pb_pre[
    ["model_id", "sm_id", "bam_url", "crc32c", "size", "gcs_obj_updated_at"]
]
pb_pre.to_csv("./data/pb_pre.csv", index=False)
