import logging
import pathlib
from typing import Optional, Tuple
from urllib.parse import urlunsplit

import pandas as pd
from google.cloud import storage
from pandera.typing import DataFrame as TypedDataFrame
from tqdm import tqdm

from dogspa_long_reads.types import (
    CopiedSampleFiles,
    ObjectMetadata,
    SamplesWithCDSIDs,
    SamplesWithShortReadMetadata,
)


def list_blobs(
    bucket_name: str, prefix: str | None = None, glob: str | None = None
) -> TypedDataFrame[ObjectMetadata]:
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

    if len(blobs) == 0:
        return TypedDataFrame[ObjectMetadata](
            pd.DataFrame(ObjectMetadata.example(size=0))
        )

    df = pd.DataFrame(blobs)

    # drop time component
    df["gcs_obj_updated_at"] = (
        df["gcs_obj_updated_at"].astype("datetime64[ns, UTC]").dt.date.astype("str")
    )

    return TypedDataFrame[ObjectMetadata](df)


def check_file_sizes(
    samples: TypedDataFrame[SamplesWithShortReadMetadata],
) -> Tuple[TypedDataFrame[SamplesWithShortReadMetadata], pd.Series]:
    """
    Check whether BAM file sizes are above configured minimum threshold.

    :param samples: the data frame of samples
    :return: the samples data frame with issue and blacklist columns filled out for
    too-small BAM files and a series indicating blacklisted rows
    """

    logging.info("Checking BAM file sizes...")

    bam_too_small = samples["size"] < 0  # TODO: use 2e9 again

    samples.loc[bam_too_small, "issue"] = samples.loc[bam_too_small, "issue"].apply(
        lambda x: x.union({"BAM file too small"})
    )
    samples.loc[bam_too_small, "blacklist"] = True

    return TypedDataFrame[SamplesWithShortReadMetadata](samples), bam_too_small


def copy_to_cclebams(
    samples: TypedDataFrame[SamplesWithCDSIDs],
    gcp_project_id: str,
    gcs_destination_bucket: str,
    gcs_destination_prefix: str,
    dry_run: bool,
) -> TypedDataFrame[CopiedSampleFiles]:
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
        value_vars=["bai", "bam"],
        var_name="url_kind",
        value_name="url",
    )

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

    return TypedDataFrame[CopiedSampleFiles](sample_files)


def rewrite_blob(src_blob: storage.Blob, dest_blob: storage.Blob) -> None:
    """
    Use Blob.rewrite to copy blob to another location

    :param src_blob: a storage Blob to copy from
    :param dest_blob: a storage Blob to copy to
    """

    logging.info(f"Copying {dest_blob.name}...")

    # assume the rewrite will take multiple requests and keep checking as long as
    # the request returns a rewrite token
    token: Optional[str] = None

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


def update_sample_file_urls(
    samples: TypedDataFrame[SamplesWithCDSIDs],
    sample_files: TypedDataFrame[CopiedSampleFiles],
) -> Tuple[TypedDataFrame[SamplesWithCDSIDs], pd.Series]:
    """
    Replace BAM URLs with new ones used in `copy_to_depmap_omics_bucket`

    :param samples: the data frame of samples
    :param sample_files: a data frame of files we attempted to copy
    :return: the samples data frame with issue and blacklist columns filled out for
    rows with files we couldn't copy and a series indicating blacklisted rows
    """

    logging.info("Updating GCS file URLs...")

    samples_updated = samples

    sample_file_urls = sample_files.loc[sample_files["copied"], ["bam_url", "new_url"]]
    samples_updated = samples_updated.merge(sample_file_urls, how="left", on="bam_url")
    samples_updated["bam_url"] = samples_updated["new_url"]
    samples_updated = samples_updated.drop(columns=["new_url"])

    missing_files = pd.Series(samples_updated[["bam_url"]].isna().any(axis=1))
    samples_updated.loc[missing_files, "issue"] = samples_updated.loc[
        missing_files, "issue"
    ].apply(lambda x: x.union({"couldn't copy BAM"}))
    samples_updated.loc[missing_files, "blacklist"] = True

    return TypedDataFrame[SamplesWithCDSIDs](samples_updated), missing_files
