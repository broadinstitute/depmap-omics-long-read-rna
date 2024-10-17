import json
import logging
import os
import pathlib
from urllib.parse import urlunsplit

import pandas as pd
from dotenv import load_dotenv
from google.cloud import storage
from nebelung.terra_workspace import TerraWorkspace
from tqdm import tqdm

from depmap_omics_long_read_rna.types import GumboClient, ModelsAndChildren
from depmap_omics_long_read_rna.utils.gcp import get_objects_metadata, rewrite_blob
from depmap_omics_long_read_rna.utils.onboarding import (
    explode_and_expand_models,
    join_metadata,
)
from depmap_omics_long_read_rna.utils.utils import df_to_model, model_to_df
from gumbo_gql_client import omics_sequencing_insert_input

load_dotenv()


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
        value_vars=["bai", "bam"],
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


terra_workspace = TerraWorkspace(
    workspace_namespace="broad-firecloud-ccle", workspace_name="Long_Read_Omics"
)

delivery_bams = terra_workspace.get_entities("delivery_bam")
delivery_bams["delivery_bam_size"] = delivery_bams["delivery_bam_size"].astype("int64")

unaligned_bams = delivery_bams[["sample_id", "delivery_bam"]].dropna()
unaligned_bams["delivery_bai"] = pd.NA
unaligned_bams = unaligned_bams.rename(
    columns={"sample_id": "sequencing_id", "delivery_bam": "bam", "delivery_bai": "bai"}
)

unaligned_bams_res = copy_to_cclebams(
    unaligned_bams,
    gcp_project_id="depmap-omics",
    gcs_destination_bucket="cclebams",
    gcs_destination_prefix="long_read_rna/",
    dry_run=False,
)

assert unaligned_bams_res["copied"].all()

aligned_bams = delivery_bams[["sample_id", "aligned_bam", "aligned_bai"]].dropna()
aligned_bams = aligned_bams.rename(
    columns={"sample_id": "sequencing_id", "aligned_bam": "bam", "aligned_bai": "bai"}
)

aligned_bams_res = copy_to_cclebams(
    aligned_bams,
    gcp_project_id="depmap-omics",
    gcs_destination_bucket="cclebams",
    gcs_destination_prefix="long_read_rna_giab/",
    dry_run=False,
)

assert aligned_bams_res["copied"].all()

gumbo_unaligned = (
    unaligned_bams_res[["sequencing_id", "url", "new_url"]]
    .rename(columns={"url": "delivery_bam", "new_url": "unaligned_bam"})
    .merge(
        delivery_bams[["delivery_bam", "delivery_bam_crc32c", "delivery_bam_size"]],
        how="left",
        on="delivery_bam",
    )
    .rename(
        columns={
            "delivery_bam_crc32c": "unaligned_bam_crc32c_hash",
            "delivery_bam_size": "unaligned_bam_size",
        }
    )
    .drop(columns="delivery_bam")
)

gumbo_aligned = (
    aligned_bams_res[["sequencing_id", "new_url", "url_kind"]]
    .pivot(columns="url_kind", index="sequencing_id", values="new_url")
    .reset_index()
)

gumbo_aligned_blobs = get_objects_metadata(gumbo_aligned["bam"])

gumbo_aligned_blobs = gumbo_aligned_blobs.rename(
    columns={
        "url": "bam",
        "crc32c": "bam_crc32c_hash",
        "size": "bam_size",
        "gcs_obj_updated_at": "sequencing_date",
    }
)

gumbo_aligned_blobs["update_time"] = gumbo_aligned_blobs["sequencing_date"]

gumbo_aligned = gumbo_aligned.merge(gumbo_aligned_blobs, how="inner", on="bam")

gumbo_samples = gumbo_unaligned.merge(gumbo_aligned, how="outer", on="sequencing_id")

gumbo_samples["blacklist"] = False
gumbo_samples["version"] = 1
gumbo_samples["expected_type"] = "long_read_rna"
gumbo_samples["source"] = "DEPMAP"
gumbo_samples["stranded"] = True

gumbo_client = GumboClient(
    url=os.environ["HASURA_URL"],
    username="snp_str_qc",
    headers={"X-Hasura-Admin-Secret": os.environ["HASURA_ADMIN_SECRET"]},
)

gumbo_samples = gumbo_samples.merge(
    delivery_bams[["sample_id", "profile_id"]].rename(
        columns={"sample_id": "sequencing_id"}
    ),
    how="left",
    on="sequencing_id",
)

gumbo_samples = gumbo_samples.rename(
    columns={
        "bam": "bam_filepath",
        "bai": "bai_filepath",
        "unaligned_bam": "unaligned_bam_filepath",
    }
)

objects = df_to_model(gumbo_samples, omics_sequencing_insert_input)
objects[0].model_dump()
res = gumbo_client.insert_omics_sequencings(
    username="depmap-omics-long-read-rna", objects=objects
)

models = model_to_df(gumbo_client.get_models_and_children(), ModelsAndChildren)
seq_table = explode_and_expand_models(models)
# seq_table = seq_table.loc[seq_table["datatype"].eq("long_read_rna")]

terra_samples = gumbo_samples[
    [
        "sequencing_id",
        "unaligned_bam_filepath",
        "bam_filepath",
        "bai_filepath",
        "profile_id",
    ]
].rename(columns={"sequencing_id": "sample_id"})


def join_metadata(samples: pd.DataFrame, seq_table: pd.DataFrame) -> pd.DataFrame:
    metadata = seq_table.loc[
        seq_table["datatype"].eq("long_read_rna")
        & ~seq_table["blacklist"]
        & ~seq_table["blacklist_omics"],
        [
            "model_id",
            "model_condition_id",
            "profile_id",
            "cell_line_name",
            "stripped_cell_line_name",
        ],
    ]

    samples_annot = samples.merge(metadata, how="left", on="profile_id")

    return samples_annot


def join_short_read_metadata(
    samples: pd.DataFrame, seq_table: pd.DataFrame
) -> pd.DataFrame:
    # reproduce logic of `makeDefaultModelTable` in depmap_omics_upload
    source_priority = [
        "BROAD",
        "DEPMAP",
        "IBM",
        "CCLE2",
        "SANGER",
        "PRISM",
        "CCLF",
        "CHORDOMA",
        "",
    ]

    sp_df = pd.DataFrame(
        {
            "source": source_priority,
            "source_priority": list(range(len(source_priority))),
        }
    )

    sr = seq_table.loc[
        seq_table["model_id"].isin(samples["model_id"])
        & seq_table["datatype"].eq("rna")
        & ~seq_table["blacklist"]
        & ~seq_table["blacklist_omics"]
    ].copy()

    sr["source"] = sr["source"].fillna("")
    sr = sr.merge(sp_df, how="left", on="source")

    assert sr["source_priority"].notna().all()

    main_seq_ids = (
        sr.loc[
            sr["is_main_sequencing_id"],
            ["model_id", "sequencing_id", "source_priority"],
        ]
        .sort_values("source_priority")
        .groupby("model_id")
        .nth(0)
        .rename(columns={"sequencing_id": "main_sequencing_id"})
    )

    samples_annot = samples.merge(main_seq_ids, how="left", on="model_id")

    sr_rna = sr.sort_values("source_priority").groupby("model_id").nth(0)

    sr_rna = sr_rna[["model_id", "profile_id", "bai_filepath", "bam_filepath"]].rename(
        columns={
            "profile_id": "sr_profile_id",
            "bai_filepath": "sr_bai_filepath",
            "bam_filepath": "sr_bam_filepath",
        }
    )

    samples_annot = samples_annot.merge(sr_rna, how="left", on="model_id")

    return samples_annot


terra_samples = join_metadata(terra_samples, seq_table)
terra_samples = join_short_read_metadata(terra_samples, seq_table)

terra_samples = terra_samples.rename(columns={"sample_id": "sequencing_id"})


def upsert_terra_samples(
    terra_workspace: TerraWorkspace, terra_samples: pd.DataFrame, dry_run: bool
) -> None:
    """
    Upsert sample and participant data to Terra data tables.

    :param terra_workspace: a TerraWorkspace instance
    :param samples: the data frame of samples
    :param dry_run: whether to skip updates to external data stores
    """

    logging.info(f"Upserting data to {terra_workspace.workspace_name} data tables")

    renames = {
        "sequencing_id": "entity:sample_id",
        "model_id": "participant_id",
        "unaligned_bam_filepath": "LR_unaligned_bam_filepath",
        "bam_filepath": "LR_bam_filepath",
        "bai_filepath": "LR_bai_filepath",
        "cell_line_name": "CellLineName",
        "stripped_cell_line_name": "StrippedCellLineName",
        "model_condition_id": "ModelCondition",
        "profile_id": "LongReadProfileID",
        "sr_profile_id": "ShortReadProfileID",
        "sr_bai_filepath": "SR_bai_filepath",
        "sr_bam_filepath": "SR_bam_filepath",
        "main_sequencing_id": "MainSequencingID",
    }

    terra_samples = terra_samples.rename(columns=renames).loc[:, renames.values()]

    terra_samples["participant"] = terra_samples["participant_id"].apply(
        lambda x: json.dumps({"entityType": "participant", "entityName": x})
    )

    # upsert participants
    participants = (
        terra_samples[["participant_id"]]
        .rename(columns={"participant": "entity:participant_id"})
        .drop_duplicates()
    )

    if dry_run:
        logging.info(f"(skipping) Upserting {len(participants)} participants")
    else:
        terra_workspace.upload_entities(participants)

    # upsert samples
    if dry_run:
        logging.info(f"(skipping) Upserting {len(terra_samples)} samples")
    else:
        terra_workspace.upload_entities(terra_samples.drop(columns="participant_id"))

    # upsert the join table between participants and samples
    participant_samples = (
        terra_samples.groupby("participant_id")["entity:sample_id"]
        .agg(list)
        .apply(
            lambda x: json.dumps([{"entityType": "sample", "entityName": y} for y in x])
        )
        .reset_index()
        .rename(
            columns={
                "participant": "entity:participant_id",
                "entity:sample_id": "samples",
            }
        )
    )

    if dry_run:
        logging.info(
            f"(skipping) Upserting {len(participant_samples)} participant-samples"
        )
    else:
        terra_workspace.upload_entities(participant_samples)


upsert_terra_samples(terra_workspace, terra_samples, dry_run=True)
