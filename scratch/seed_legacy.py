import re
from io import StringIO

import pandas as pd
from google.cloud import storage

pd.set_option("display.max_columns", 30)
pd.set_option("display.max_colwidth", 100)
pd.set_option("display.max_info_columns", 30)
pd.set_option("display.max_info_rows", 20)
pd.set_option("display.max_rows", 20)
pd.set_option("display.max_seq_items", None)
pd.set_option("display.width", 200)
pd.set_option("expand_frame_repr", True)
pd.set_option("mode.chained_assignment", "warn")


def list_blobs(
    bucket_name: str, prefix: str | None = None, glob: str | None = None
) -> pd.DataFrame:
    """
    Get the names and sizes of existing blobs in a GCS bucket.

    :param bucket_name: the name of the GCS bucket
    :param prefix: an optional prefix for listing
    :param glob: an optional glob for listing
    :return: a data frame of object names and sizes
    """

    storage_client = storage.Client()

    if prefix is not None and glob is not None:
        raise ValueError("At most one of `prefix` and `glob` can be specified")
    elif prefix is not None:
        pages = storage_client.list_blobs(
            bucket_or_name=bucket_name,
            prefix=prefix,
            delimiter="/",
            fields="items(name,crc32c,size,updated),nextPageToken",
        ).pages
    elif glob is not None:
        pages = storage_client.list_blobs(
            bucket_or_name=bucket_name,
            match_glob=glob,
            fields="items(name,crc32c,size,updated),nextPageToken",
        ).pages
    else:
        pages = storage_client.list_blobs(
            bucket_or_name=bucket_name,
            fields="items(name,crc32c,size,updated),nextPageToken",
        ).pages

    blobs = []

    for page in pages:
        blobs.extend(
            [
                {
                    "url": "gs://" + bucket_name + "/" + x.name,
                    "crc32c": x.crc32c,
                    "size": x.size,
                    "gcs_obj_updated_at": x.updated,
                }
                for x in page
            ]
        )

    df = pd.DataFrame(blobs)

    df["gcs_obj_updated_at"] = df["gcs_obj_updated_at"].astype("datetime64[ns, UTC]")

    return df


bams = list_blobs(
    "fc-aaf4de93-c104-45c4-a01a-a036869119c6",
    glob="PB_preprocess/SM-*/merge/*.sorted.bam",
)

bams = pd.DataFrame(bams)

csvs = list_blobs(
    "fc-aaf4de93-c104-45c4-a01a-a036869119c6", glob="PB_preprocess/SM-*/*.csv"
)

storage_client = storage.Client()


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
    "batch_id",
    "rin_score",
    "rin",
]
metadata = (
    metadata.drop(columns=["kinnex_adapter", "isoseq_primer", "rin_score", "rin"])
    .drop_duplicates()
    .rename(columns={"ach_id": "model_id"})
)

bams["batch_id"] = bams["url"].str.extract(r"/(SM-[^/]+)/")
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
    metadata_dedup, on=["batch_id", "sample_id_canon"], how="outer", indicator=True
)

unmapped_sample_id_canon = df.loc[
    df["sample_id_canon"].duplicated(keep=False), "sample_id_canon"
]

df = df.loc[df["_merge"].eq("both")]

df2 = bams.merge(
    metadata, on=["batch_id", "sample_id_canon"], how="outer", indicator=True
)
df2 = df2.loc[df2["sample_id_canon"].isin(unmapped_sample_id_canon)]
df2_dedup = (
    df2.sort_values("gcs_obj_updated_at", ascending=False)
    .groupby("sample_id_canon")
    .nth(0)
)

legacy = pd.concat([df, df2_dedup]).rename(columns={"url": "bam_url"})
legacy["size"] = legacy["size"].astype("int64")
legacy = legacy[["model_id", "bam_url", "crc32c", "size", "gcs_obj_updated_at"]]

pilot_bams = list_blobs(
    "fc-aaf4de93-c104-45c4-a01a-a036869119c6", glob="pilot/merge/*.sorted.bam"
)

pilot_bams = pilot_bams.rename(columns={"url": "bam_url"})
pilot_bams["model_id"] = pilot_bams["bam_url"].str.extract(r"(ACH-.+).sorted.bam$")

legacy = pd.concat((legacy, pilot_bams))
legacy["gcs_obj_updated_at"] = legacy["gcs_obj_updated_at"].dt.date.astype("str")

legacy.to_csv("./data/legacy_bams.csv", index=False)
