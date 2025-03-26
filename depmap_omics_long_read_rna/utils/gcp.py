import pandas as pd
from google.cloud import storage
from nebelung.utils import type_data_frame
from pandera.typing import DataFrame as TypedDataFrame

from depmap_omics_long_read_rna.types import ObjectMetadata


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

    if len(blobs) == 0:
        return TypedDataFrame[ObjectMetadata](
            pd.DataFrame(ObjectMetadata.example(size=0))
        )

    df = pd.DataFrame(blobs)

    # drop time component
    df["gcs_obj_updated_at"] = (
        df["gcs_obj_updated_at"].astype("datetime64[ns, UTC]").dt.date.astype("str")
    )

    df = df.loc[~df["url"].str.contains("2240/")]

    return type_data_frame(df, ObjectMetadata)
