import logging
import pathlib
from urllib.parse import urlunsplit

import pandas as pd
import tomllib
from google.cloud import storage
from nebelung.terra_workspace import TerraWorkspace
from tqdm import tqdm

pd.set_option("display.max_columns", 30)
pd.set_option("display.max_colwidth", 50)
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

    return pd.DataFrame(blobs)


def copy_to_cclebams(
    samples: pd.DataFrame,
    gcp_project_id: str,
    gcs_destination_bucket: str,
    gcs_destination_prefix: str,
    dry_run: bool,
) -> pd.DataFrame:
    """
    Copy all BAM files in the samples data frame to our bucket

    :param samples: the data frame of samples
    :param gcp_project_id: a GCP project ID to use for billing
    :param gcs_destination_bucket: the name of the destination bucket
    :param gcs_destination_prefix: an object prefix for the copied BAMs
    :param dry_run: whether to skip updates to external data stores
    :return: a data frame of files we attempted to copy
    """

    logging.info("Copying files to our bucket...")

    storage_client = storage.Client(project=gcp_project_id)

    # collect all CRAM/BAM/CRAI/BAI files to copy
    sample_files = samples.melt(
        id_vars="sequencing_id",
        value_vars=["crai_or_bai_url", "cram_or_bam_url"],
        var_name="url_kind",
        value_name="url",
    ).dropna()

    # all copied files will have same destination bucket and prefix
    dest_bucket = storage_client.bucket(
        gcs_destination_bucket, user_project=gcp_project_id
    )
    prefix = gcs_destination_prefix.strip("/")

    # keep track of copy attempts
    copy_results = []

    logging.info(f"Checking {len(sample_files)} files to copy...")

    # can't use rewrite in a batch context, so do plain iteration
    for r in tqdm(sample_files.itertuples(index=False), total=len(sample_files)):
        url = r[sample_files.columns.get_loc("url")]

        try:
            # construct the source blob
            src_blob = storage.Blob.from_string(url, client=storage_client)

            # construct the destination blob (named by sample ID)
            dest_file_ext = pathlib.Path(str(src_blob.name)).suffix
            dest_obj_key = (
                "/".join(
                    [
                        prefix,
                        str(r[sample_files.columns.get_loc("sequencing_id")]),
                    ]
                )
                + dest_file_ext
            )
            dest_blob = storage.Blob(dest_obj_key, bucket=dest_bucket)

            # GCS rewrite operation is instantaneous if location and storage class match
            dest_blob.storage_class = src_blob.bucket.get_blob(
                src_blob.name
            ).storage_class

            new_url = urlunsplit(("gs", dest_bucket.name, dest_obj_key, "", ""))

            if dest_blob.exists():
                logging.info(f"{dest_blob.name} already exists in {dest_bucket.name}")
                copy_results.append({"url": url, "new_url": new_url, "copied": True})

            else:
                if dry_run:
                    logging.info(
                        f"(skipping) Copying {dest_blob.name} to {dest_bucket.name}"
                    )
                else:
                    logging.info(f"Copying {dest_blob.name} to {dest_bucket.name}")
                    rewrite_blob(src_blob, dest_blob)

                copy_results.append({"url": url, "new_url": new_url, "copied": True})

        except Exception as e:
            logging.error(f"Error copying {url} to {dest_bucket}: {e}")
            copy_results.append({"url": url, "new_url": None, "copied": False})

    sample_files = sample_files.merge(pd.DataFrame(copy_results), how="inner", on="url")

    return sample_files


def rewrite_blob(src_blob: storage.Blob, dest_blob: storage.Blob) -> None:
    """
    Use Blob.rewrite to copy blob to another location

    :param src_blob: a storage Blob to copy from
    :param dest_blob: a storage Blob to copy to
    """

    logging.info(f"Copying {dest_blob.name}...")

    # assume the rewrite will take multiple requests and keep checking as long as
    # the request returns a rewrite token
    token = None

    while True:
        token, bytes_rewritten, total_bytes = dest_blob.rewrite(
            source=src_blob, token=token, timeout=60
        )

        logging.info(
            "Copied {bytes_rewritten} / {total_bytes} bytes ({prop}%)".format(
                bytes_rewritten=bytes_rewritten,
                total_bytes=total_bytes,
                prop=round(bytes_rewritten / total_bytes * 100),
            )
        )

        if token is None:
            break


logger = logging.getLogger()
logger.addHandler(logging.StreamHandler())
logger.setLevel(logging.INFO)

config = {}

with open("./config.toml", "rb") as f:
    config.update(tomllib.load(f))

delivery_workspace = TerraWorkspace(
    workspace_namespace=config["terra"]["delivery_workspace_namespace"],
    workspace_name=config["terra"]["delivery_workspace_name"],
    owners=[],
)

samples = delivery_workspace.get_entities("sample")

alignments = samples[["sample_id", "aligned_bam", "aligned_bai"]].dropna()

alignments["aligned_bam_dest"] = (
    "gs://cclebams/long_read_rna_hg38/" + alignments["sample_id"] + ".bam"
)
alignments["aligned_bai_dest"] = (
    "gs://cclebams/long_read_rna_hg38/" + alignments["sample_id"] + ".bai"
)

existing_blobs = list_blobs(bucket_name="cclebams", prefix="long_read_rna_hg38/")

to_copy = alignments.loc[
    ~(
        alignments["aligned_bam_dest"].isin(existing_blobs["url"])
        & alignments["aligned_bai_dest"].isin(existing_blobs["url"])
    )
].rename(
    columns={
        "sample_id": "sequencing_id",
        "aligned_bam": "cram_or_bam_url",
        "aligned_bai": "crai_or_bai_url",
    }
)

copy_to_cclebams(
    samples=to_copy,
    gcp_project_id="depmap-omics",
    gcs_destination_bucket="cclebams",
    gcs_destination_prefix="long_read_rna_hg38",
    dry_run=False,
)

existing_blobs = list_blobs(bucket_name="cclebams", prefix="long_read_rna_hg38/")

to_copy_meta = to_copy[["sequencing_id", "aligned_bam_dest", "aligned_bai_dest"]].merge(
    existing_blobs, how="inner", left_on="aligned_bam_dest", right_on="url"
)

to_copy_meta["reference_genome"] = "hg38"
to_copy_meta["sequencing_alignment_source"] = "CDS"

to_copy_meta = to_copy_meta.rename(
    columns={
        "sequencing_id": "omics_sequencing_id",
        "aligned_bam_dest": "url",
        "aligned_bai_dest": "index_url",
        "crc32c": "crc32c_hash",
    }
)

to_copy_meta.to_csv("/tmp/to_copy_meta.csv", index=False)
