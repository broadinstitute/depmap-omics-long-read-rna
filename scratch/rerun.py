import contextlib
import os
import tempfile
from pathlib import Path
from pprint import pp

import google.auth.transport.requests
import google.oauth2.id_token
import pandas as pd
import pysam
import tqdm
from nebelung.terra_workspace import TerraWorkspace
from pd_flatten import pd_flatten

from depmap_omics_long_read_rna.utils.gcp import get_objects_metadata, list_blobs

pd.set_option("display.max_columns", 30)
pd.set_option("display.max_colwidth", 50)
pd.set_option("display.max_info_columns", 30)
pd.set_option("display.max_info_rows", 20)
pd.set_option("display.max_rows", 20)
pd.set_option("display.max_seq_items", None)
pd.set_option("display.width", 200)
pd.set_option("expand_frame_repr", True)
pd.set_option("mode.chained_assignment", "warn")


def make_delivery_bam_df(bams: pd.DataFrame) -> pd.DataFrame:
    """
    Prepare data frame for the Terra delivery_bam data table using inventory of GCS
    blobs.

    :param bams: a data frame of delivered uBAMs
    :return: a data frame of delivered uBAM metadata
    """

    bams_w_ids = bams.copy()

    # extract the profile ID from the BAM filenames
    bams_w_ids["omics_profile_id"] = (
        bams_w_ids["url"]
        .apply(lambda x: Path(x).name)
        .str.extract(r"^(PR-[A-Z a-z 0-9]+)")
    )

    if bool(bams_w_ids["omics_profile_id"].isna().any()):
        raise ValueError("There are BAM files not named with profile IDs (PR-*.bam)")

    bams_w_ids = bams_w_ids.rename(columns={"url": "bam_url"})

    bams_w_ids = bams_w_ids.rename(
        columns={
            "bam_url": "delivery_bam",
            "crc32c": "delivery_bam_crc32c",
            "size": "delivery_bam_size",
            "gcs_obj_updated_at": "delivery_bam_updated_at",
        }
    )

    return bams_w_ids


def get_gcp_access_token() -> str:
    """
    Get a GCP access token for current credentials. (Some tools like pysam that don't go
    through gcloud native Python APIs need this.)

    :return: the auth/bearer token
    """

    auth_req = google.auth.transport.requests.Request()  # type: ignore

    creds, _ = google.auth.default(
        scopes=["https://www.googleapis.com/auth/cloud-platform"]
    )
    creds.refresh(auth_req)  # type: ignore
    token = creds.token  # type: ignore

    if token is None:
        raise ValueError("GCP auth token cannot be None")

    return token


@contextlib.contextmanager
def tmp_cwd():
    """
    Change the working directory of the context to a temporary directory
    """

    orig_cwd = os.getcwd()

    try:
        with tempfile.TemporaryDirectory() as tmpdirname:
            os.chdir(tmpdirname)
            yield

    finally:
        os.chdir(orig_cwd)


def get_bam_headers(urls):
    # pysam needs GCS auth set as an ENV var
    os.environ["GCS_OAUTH_TOKEN"] = get_gcp_access_token()

    headers = []

    # no way to tell pysam not to download BAI files to the current working directory,
    # so temporarily change the working directory for this process and its children
    with tmp_cwd():
        for url in tqdm.tqdm(urls):
            afile = pysam.AlignmentFile(url, "rb", check_sq=False)
            headers.append({"url": url, "headers": afile.header.to_dict()})

    return pd_flatten(pd.DataFrame(headers), sep="_")


mapping = pd.read_csv(
    "./data/updated_with_merged_unaligned_07022025.tsv",
    sep="\t",
    usecols=[
        "entity:sample_id",
        "omics_profile_id",
        "delivery_bam",
        "updated_bam_path_0702",
    ],
)

pp(
    mapping.loc[mapping["omics_profile_id"].str.contains("PR-LaBpwf")].to_dict(
        orient="records"
    )
)

mapping_old_blobs = get_objects_metadata(mapping["delivery_bam"])
mapping_new_blobs = get_objects_metadata(mapping["updated_bam_path_0702"])

mapping_blobs = mapping.merge(
    mapping_old_blobs, how="left", left_on="delivery_bam", right_on="url"
).merge(
    mapping_new_blobs,
    how="left",
    left_on="updated_bam_path_0702",
    right_on="url",
    suffixes=("_old", "_new"),
)

mapping_blobs["size_rat"] = mapping_blobs["size_new"] / mapping_blobs["size_old"]
mapping_blobs = mapping_blobs.sort_values("size_rat")

mapping_old_headers = get_bam_headers(mapping["delivery_bam"])
mapping_new_headers = get_bam_headers(mapping["updated_bam_path_0702"])

all_headers = pd.concat([mapping_old_headers, mapping_new_headers])
all_headers = all_headers.loc[
    all_headers["headers_PG_ID"].eq("samtools"),
    ["url", "headers_PG_ID", "headers_PG_CL"],
].drop_duplicates()

pp(
    all_headers.loc[all_headers["headers_PG_CL"].str.contains("bcM0001.bc06")].to_dict(
        orient="records"
    )
)

merge_headers = mapping.merge(
    mapping_old_headers.loc[
        mapping_old_headers["headers_PG_ID"].eq("samtools"),
        ["url", "headers_PG_ID", "headers_PG_CL"],
    ].drop_duplicates(),
    how="left",
    left_on="delivery_bam",
    right_on="url",
).merge(
    mapping_new_headers.loc[
        mapping_new_headers["headers_PG_ID"].eq("samtools"),
        ["url", "headers_PG_ID", "headers_PG_CL"],
    ].drop_duplicates(),
    how="left",
    left_on="updated_bam_path_0702",
    right_on="url",
    suffixes=("_old", "_new"),
)

merge_headers["cmd_old"] = merge_headers["headers_PG_CL_old"].str.extract(
    r"\.bam\s+(.*)"
)
merge_headers["cmd_new"] = merge_headers["headers_PG_CL_new"].str.extract(
    r"\.bam\s+(.*)"
)

pp(
    merge_headers.loc[merge_headers["cmd_old"].ne(merge_headers["cmd_new"]),].to_dict(
        orient="records"
    )
)

pp(
    merge_headers.loc[
        merge_headers["omics_profile_id"].isin(["PR-BmQo3y", "PR-LaBpwf", "PR-2P0zxk"])
    ].to_dict(orient="records")
)


pp(
    mapping_old_headers.loc[
        mapping_old_headers["headers_PG_CL"].str.contains("bcM0001.bc06")
    ].to_dict(orient="records")
)

mapping_blobs.to_csv("./data/mapping.csv", index=False)

gcs_source_bucket = "fc-aaf4de93-c104-45c4-a01a-a036869119c6"
gcs_source_glob = "*rerun*/**/*.*.bam"
src_bams = list_blobs(bucket_name=gcs_source_bucket, glob=gcs_source_glob)
bams = make_delivery_bam_df(src_bams)

bams.loc[bams["omics_profile_id"].eq("PR-LaBpwf")].values
bams.loc[~bams["delivery_bam"].isin(mapping["updated_bam_path_0702"])]

bams = bams.loc[
    bams["delivery_bam"].ne(
        "gs://fc-aaf4de93-c104-45c4-a01a-a036869119c6/PDO-34724_rerun_July2025/PR-LaBpwf.merged.unaligned.bam"
    )
    & bams["omics_profile_id"].ne("PR-GvR49r")
]

terra_workspace = TerraWorkspace("broad-firecloud-ccle", "Long_Read_Omics")
existing_samples = terra_workspace.get_entities("sample")

m = bams.merge(
    existing_samples, on="omics_profile_id", how="left", suffixes=("", "_old")
)

m["ratio"] = m["delivery_bam_size"] / m["delivery_bam_size_old"]
m = m.sort_values("ratio")

pp(m.loc[m["omics_profile_id"].duplicated(keep=False)])

m = m.loc[
    :,
    [
        "sample_id",
        "delivery_bam",
        "delivery_bam_updated_at",
        "delivery_bam_crc32c",
        "delivery_bam_size",
    ],
]

terra_workspace.create_entity_set("sample", m["sample_id"], "polya_rerun")

terra_workspace.upload_entities(m)
