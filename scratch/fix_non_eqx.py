import tomllib
from typing import Iterable

import pandas as pd
from google.cloud import storage
from nebelung.terra_workspace import TerraWorkspace

from depmap_omics_long_read_rna.utils.gcp import copy_to_cclebams

pd.set_option("display.max_columns", 30)
pd.set_option("display.max_colwidth", 50)
pd.set_option("display.max_info_columns", 30)
pd.set_option("display.max_info_rows", 20)
pd.set_option("display.max_rows", 20)
pd.set_option("display.max_seq_items", None)
pd.set_option("display.width", 200)
pd.set_option("expand_frame_repr", True)
pd.set_option("mode.chained_assignment", "warn")


def get_objects_metadata(urls: Iterable[str]) -> pd.DataFrame:
    storage_client = storage.Client(project="depmap-omics")

    blobs = {}

    with storage_client.batch(raise_exception=False):
        for url in urls:
            blob = storage.Blob.from_string(url, client=storage_client)
            bucket = storage_client.bucket(blob.bucket.name)
            blob = bucket.get_blob(blob.name)
            blobs[url] = blob

    metadata = [
        {"url": k, "crc32c": v.crc32c, "size": v.size, "gcs_obj_updated_at": v.updated}
        for k, v in blobs.items()
    ]

    return pd.DataFrame(metadata)


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
        pages = storage_client.list_blobs(  # pyright: ignore
            bucket_or_name=bucket_name,
            prefix=prefix,
            delimiter="/",
            fields="items(name,crc32c,size,updated),nextPageToken",
        ).pages
    elif glob is not None:
        pages = storage_client.list_blobs(  # pyright: ignore
            bucket_or_name=bucket_name,
            match_glob=glob,
            fields="items(name,crc32c,size,updated),nextPageToken",
        ).pages
    else:
        pages = storage_client.list_blobs(  # pyright: ignore
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

    return pd.DataFrame(blobs)


# use same config loading as when calling the module CLI
with open("config.toml", "rb") as f:
    config = tomllib.load(f)

workspace = TerraWorkspace(
    workspace_namespace=config["terra"]["workspace_namespace"],
    workspace_name=config["terra"]["workspace_name"],
)

samples = workspace.get_entities("sample")

blobs = get_objects_metadata(urls=samples["aligned_cram_bam"].dropna())
missing_blobs = blobs.loc[blobs["size"].isna()]

all_blobs = list_blobs(
    bucket_name="fc-secure-d068ecea-bdcf-4e38-a449-edcd8933a71c",
    glob="submissions/final-outputs/**/CDS-*.bam",
)

missing_blobs["cds_id"] = missing_blobs["url"].str.extract(r"(CDS-.{6})")
all_blobs["cds_id"] = all_blobs["url"].str.extract(r"(CDS-.{6})")

missing_blobs = missing_blobs.merge(all_blobs, how="left", on="cds_id")

fixed_samples = (
    missing_blobs.copy()
    .rename(columns={"cds_id": "sample_id", "url_y": "aligned_cram_bam"})
    .loc[:, ["sample_id", "aligned_cram_bam"]]
)

fixed_samples["aligned_crai_bai"] = fixed_samples["aligned_cram_bam"].str.replace(
    ".bam", ".bai", regex=False
)

sample_files = copy_to_cclebams(
    fixed_samples,
    gcp_project_id="depmap-omics",
    gcs_destination_bucket="cclebams",
    gcs_destination_prefix="long_read_rna_hg38",
    dry_run=False,
)

blobs = get_objects_metadata(urls=samples["aligned_cram_bam"].dropna())

[
    print(x)
    for x in "update sequencing_alignment set crc32c_hash='"
    + blobs["crc32c"]
    + "', size="
    + blobs["size"].astype("Int64").astype("string")
    + " where url='"
    + blobs["url"]
    + "';"
]
