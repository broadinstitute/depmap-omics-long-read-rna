import pandas as pd
from nebelung.terra_workspace import TerraWorkspace
from nebelung.utils import type_data_frame
from pandera.typing import DataFrame as TypedDataFrame
from pd_flatten import pd_flatten

from depmap_omics_long_read_rna.types import (
    AlignmentMetadataLong,
    GumboClient,
    LongReadAlignmentMetadata,
    LongReadTerraSamples,
    ModelsAndChildren,
    ShortReadMetadata,
    ShortReadTerraSamples,
)
from depmap_omics_long_read_rna.utils.utils import model_to_df


def do_refresh_terra_samples(
    terra_workspace: TerraWorkspace,
    short_read_terra_workspace: TerraWorkspace,
    gumbo_client: GumboClient,
) -> None:
    """
    Upsert the sample data table in Terra with the latest long read sample metadata from
    Gumbo and the short read Terra workspace.

    :param terra_workspace: the Terra workspace for long read RNA samples
    :param short_read_terra_workspace: the Terra workspace for short read RNA samples
    :param gumbo_client: an instance of the Gumbo GraphQL client
    """

    models = model_to_df(
        gumbo_client.get_models_and_children(timeout=30.0), ModelsAndChildren
    )

    # make wide data frame of Gumbo sequencing_alignment and related records
    alignments = explode_and_expand_models(models)

    # start constructing a sample data frame
    samples = collect_lr_alignments(alignments)

    # join short read data
    samples = join_sr_metadata(samples, alignments, short_read_terra_workspace)

    terra_workspace.upload_entities(df=samples)


def explode_and_expand_models(
    models: TypedDataFrame[ModelsAndChildren],
) -> TypedDataFrame[AlignmentMetadataLong]:
    """
    Unnest columns from the GraphQL call for Gumbo models and their childen.

    :param models: a data frame of models with nested profile/sequencing/etc. data
    :return: a wide version of the data frame without nesting
    """

    alignments = pd_flatten(models, name_columns_with_parent=False)

    alignments = alignments.dropna(
        subset=[
            "model_id",
            "model_condition_id",
            "omics_profile_id",
            "omics_sequencing_id",
            "sequencing_alignment_id",
        ]
    )

    return type_data_frame(alignments, AlignmentMetadataLong, remove_unknown_cols=True)


def collect_lr_alignments(alignments: TypedDataFrame[AlignmentMetadataLong]):
    """
    Construct a wide data frame of long read RNA sample metadata using data from Gumbo.

    :param alignments: a data frame of sequencing_alignment metadata from Gumbo
    :return: a data frame of long read samples
    """

    # collect sequencing_alignment records for GP-delivered uBAMs and our analysis-ready
    # hg38-aligned BAMs
    lr_alignments = alignments.loc[alignments["datatype"].eq("long_read_rna")]
    lr_alignments_gp = lr_alignments.loc[
        lr_alignments["sequencing_alignment_source"].eq("GP")
    ]
    lr_alignments_cds = lr_alignments.loc[
        lr_alignments["sequencing_alignment_source"].eq("CDS")
    ]

    samples = (
        lr_alignments_gp[
            [
                "omics_sequencing_id",
                "omics_profile_id",
                "model_condition_id",
                "model_id",
                "cell_line_name",
                "stripped_cell_line_name",
                "sequencing_alignment_id",
                "cram_bam_url",
            ]
        ]
        .rename(
            columns={
                "omics_sequencing_id": "sample_id",
                "sequencing_alignment_id": "delivery_sequencing_alignment_id",
                "cram_bam_url": "delivery_cram_bam",
            }
        )
        .merge(
            lr_alignments_cds[
                [
                    "omics_sequencing_id",
                    "sequencing_alignment_id",
                    "cram_bam_url",
                    "crai_bai_url",
                    "reference_genome",
                ]
            ].rename(
                columns={
                    "omics_sequencing_id": "sample_id",
                    "sequencing_alignment_id": "aligned_sequencing_alignment_id",
                    "cram_bam_url": "aligned_cram_bam",
                    "crai_bai_url": "aligned_crai_bai",
                }
            ),
            how="left",
            on="sample_id",
        )
    )

    return type_data_frame(samples, LongReadAlignmentMetadata)


def join_sr_metadata(
    samples: TypedDataFrame[LongReadAlignmentMetadata],
    alignments: TypedDataFrame[AlignmentMetadataLong],
    short_read_terra_workspace: TerraWorkspace,
) -> TypedDataFrame[LongReadTerraSamples]:
    """
    Collect metadata from Gumbo and workflow outputs from Terra relating to short read
    RNA samples and join them to the long read samples.

    :param samples: a data frame of long read samples
    :param alignments: a data frame of sequencing_alignment metadata from Gumbo
    :param short_read_terra_workspace: the Terra workspace for short read RNA samples
    :return: the long read samples with short read sample metaadata joined to it
    """

    # get the current short read sample data from Terra
    sr_terra_samples = get_sr_terra_samples(short_read_terra_workspace)

    # collect short read metadata from Gumbo
    sr_sequencings = alignments.loc[
        alignments["datatype"].eq("rna"),
        [
            "model_condition_id",
            "omics_profile_id",
            "omics_sequencing_id",
            "source",
            "version",
        ],
    ].drop_duplicates(subset="omics_sequencing_id")

    # match a short read sample to each long read sample
    sr_metadata = choose_matched_short_read_sample(samples, sr_sequencings)

    sr_metadata = sr_terra_samples.merge(
        sr_metadata,
        how="inner",
        left_on="sample_id",
        right_on="sr_omics_sequencing_id",
    ).drop(columns="sample_id")

    samples_w_metadata = samples.merge(sr_metadata, how="left", on="model_condition_id")

    return type_data_frame(samples_w_metadata, LongReadTerraSamples)


def get_sr_terra_samples(
    short_read_terra_workspace: TerraWorkspace,
) -> TypedDataFrame[ShortReadTerraSamples]:
    """
    Get a subset of short read RNA sample workflow output URLs from Terra.

    :param short_read_terra_workspace: the Terra workspace for short read RNA samples
    :return: a data frame of shot read output file URLs
    """

    sr_terra_samples = short_read_terra_workspace.get_entities("sample")

    sr_terra_samples = (
        sr_terra_samples.loc[
            :,
            [
                "sample_id",
                "star_junctions",
                "fusion_predictions",
                "fusion_predictions_abridged",
                "rsem_genes",
                "rsem_genes_stranded",
                "rsem_isoforms",
                "rsem_isoforms_stranded",
            ],
        ]
        .rename(
            columns={
                "star_junctions": "sr_star_junctions",
                "fusion_predictions": "sr_fusion_predictions",
                "fusion_predictions_abridged": "sr_fusion_predictions_abridged",
                "rsem_genes": "sr_rsem_genes",
                "rsem_genes_stranded": "sr_rsem_genes_stranded",
                "rsem_isoforms": "sr_rsem_isoforms",
                "rsem_isoforms_stranded": "sr_rsem_isoforms_stranded",
            }
        )
        .astype("string")
    )

    return type_data_frame(sr_terra_samples, ShortReadTerraSamples)


def choose_matched_short_read_sample(
    samples: pd.DataFrame, sr_sequencings: pd.DataFrame
) -> TypedDataFrame[ShortReadMetadata]:
    """
    Choose short read RNA samples to map to long read ones.

    :param samples: a data frame of long read samples
    :param sr_sequencings: the short read workspace sample data table
    :return: a data frame mapping common model condition IDs to short read profile and
    sequencing IDs
    """

    # prioritize certain sample source
    source_priority = [
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

    # collect all valid short read samples in Gumbo that are present in both short and
    # long read workspaces
    sr = sr_sequencings.loc[
        sr_sequencings["model_condition_id"].isin(samples["model_condition_id"])
    ].copy()

    # populate the `source_priority` column
    sr["source"] = sr["source"].fillna("")
    sr = sr.merge(sp_df, how="left", on="source")
    assert sr["source_priority"].notna().all()

    # pick a short read sample for each model condition, using the more recent `version`
    # to break ties
    sr_choices = (
        sr.sort_values(
            ["model_condition_id", "source_priority", "version"],
            ascending=[True, True, False],
        )
        .groupby("model_condition_id")
        .nth(0)
    )

    sr_choices = sr_choices[
        ["model_condition_id", "omics_profile_id", "omics_sequencing_id"]
    ].rename(
        columns={
            "omics_profile_id": "sr_omics_profile_id",
            "omics_sequencing_id": "sr_omics_sequencing_id",
        }
    )

    return type_data_frame(sr_choices, ShortReadMetadata)
