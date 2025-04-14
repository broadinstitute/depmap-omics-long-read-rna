import logging
from functools import partial

from nebelung.terra_workspace import TerraWorkspace
from nebelung.utils import type_data_frame
from pandera.typing import DataFrame as TypedDataFrame
from pd_flatten import pd_flatten

from depmap_omics_long_read_rna.types import (
    AlignedSamplesWithObjectMetadata,
    CopiedSampleFiles,
    ExistingAlignments,
    GumboClient,
    NewSequencingAlignments,
)
from depmap_omics_long_read_rna.utils.gcp import copy_to_cclebams, get_objects_metadata
from depmap_omics_long_read_rna.utils.utils import model_to_df
from gumbo_gql_client import sequencing_alignment_insert_input


def onboard_aligned_bams(
    terra_workspace: TerraWorkspace, gumbo_client: GumboClient, dry_run: bool
) -> None:
    """
    TODO

    :param terra_workspace: a TerraWorkspace instance
    :param gumbo_client: an instance of the Gumbo GraphQL client
    :param dry_run: whether to skip updates to external data stores
    """

    samples = terra_workspace.get_entities("sample")
    samples = samples.loc[:, ["sample_id", "aligned_bam", "aligned_bai"]].dropna()

    existing_alignments = model_to_df(
        gumbo_client.long_read_sequencing_alignments(),
        ExistingAlignments,
        mutator=partial(pd_flatten, name_columns_with_parent=False),
    )

    seq_ids_no_cds = set(
        existing_alignments.loc[
            existing_alignments["sequencing_alignment_source"].eq("GP"),
            "omics_sequencing_id",
        ]
    ).difference(
        existing_alignments.loc[
            existing_alignments["sequencing_alignment_source"].eq("CDS"),
            "omics_sequencing_id",
        ]
    )

    samples = samples.loc[samples["sample_id"].isin(list(seq_ids_no_cds))]

    objects_metadata = get_objects_metadata(samples["aligned_bam"])

    samples = type_data_frame(
        samples.merge(
            objects_metadata, how="left", left_on="aligned_bam", right_on="url"
        ).drop(columns=["url", "gcs_obj_updated_at"]),
        AlignedSamplesWithObjectMetadata,
    )

    assert ~bool(samples["size"].isin(existing_alignments["size"]).any())

    sample_files = copy_to_cclebams(
        samples,
        gcp_project_id="depmap-omics",
        gcs_destination_bucket="cclebams",
        gcs_destination_prefix="long_read_rna_hg38",
        dry_run=dry_run,
    )

    samples = update_sample_file_urls(samples, sample_files)

    persist_sequencing_alignments(gumbo_client, samples, dry_run)


def update_sample_file_urls(
    samples: TypedDataFrame[AlignedSamplesWithObjectMetadata],
    sample_files: TypedDataFrame[CopiedSampleFiles],
) -> TypedDataFrame[AlignedSamplesWithObjectMetadata]:
    """
    Replace BAM URLs with new ones used in `copy_to_depmap_omics_bucket`.

    :param samples: the data frame of samples
    :param sample_files: a data frame of files we attempted to copy
    :return: the samples data frame with issue column filled out for rows with files we
    couldn't copy
    """

    logging.info("Updating GCS file URLs...")

    samples_updated = samples.copy()

    for c in ["aligned_bai", "aligned_bam"]:
        sample_file_urls = sample_files.loc[sample_files["copied"], ["url", "new_url"]]
        samples_updated = samples_updated.merge(
            sample_file_urls, how="left", left_on=c, right_on="url"
        )
        samples_updated[c] = samples_updated["new_url"]
        samples_updated = samples_updated.drop(columns=["url", "new_url"])

    return type_data_frame(samples_updated, AlignedSamplesWithObjectMetadata)


def persist_sequencing_alignments(
    gumbo_client: GumboClient,
    samples: TypedDataFrame[AlignedSamplesWithObjectMetadata],
    dry_run: bool,
) -> None:
    """
    Insert many sequencing alignments to Gumbo.

    :param gumbo_client: an instance of the Gumbo GraphQL client
    :param samples: a data frame of prepared sequencing alignments to insert
    :param dry_run: whether to skip updates to external data stores
    """

    new_sequencing_alignments = samples.rename(
        columns={
            "sample_id": "omics_sequencing_id",
            "crc32c": "crc32c_hash",
            "aligned_bam": "url",
            "aligned_bai": "index_url",
        }
    )

    new_sequencing_alignments["reference_genome"] = "hg38"
    new_sequencing_alignments["sequencing_alignment_source"] = "CDS"

    new_sequencing_alignments = type_data_frame(
        new_sequencing_alignments, NewSequencingAlignments
    )

    sequencing_alignment_inserts = [
        sequencing_alignment_insert_input.model_validate(x)
        for x in new_sequencing_alignments.to_dict(orient="records")
    ]

    if dry_run:
        logging.info(
            f"(skipping) Inserting {len(sequencing_alignment_inserts)} "
            "sequencing alignments"
        )
        return

    res = gumbo_client.insert_sequencing_alignments(
        gumbo_client.username, objects=sequencing_alignment_inserts
    )

    affected_rows = res.insert_sequencing_alignment.affected_rows  # pyright: ignore
    logging.info(f"Inserted {affected_rows} sequencing alignments")
