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
    ShortReadTerraSamples,
)
from depmap_omics_long_read_rna.utils.utils import model_to_df


def refresh_terra_samples(
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
        gumbo_client.mapped_short_long_rna_samples(timeout=30.0), ModelsAndChildren
    )

    # make wide data frame of Gumbo sequencing_alignment and related records
    alignments = explode_and_expand_models(models)

    # start constructing a sample data frame
    samples = collect_lr_alignments(alignments)

    # join short read data
    samples = join_sr_data(samples, alignments, short_read_terra_workspace)

    terra_workspace.upload_entities(df=samples)


def explode_and_expand_models(
    models: TypedDataFrame[ModelsAndChildren],
) -> TypedDataFrame[AlignmentMetadataLong]:
    """
    Unnest columns from the GraphQL call for Gumbo models and their childen.

    :param models: a data frame of models with nested model/sequencing/etc. data
    :return: a wide version of the data frame without nesting
    """

    alignments = pd_flatten(models, name_columns_with_parent=False).dropna(
        subset="sequencing_alignment_id"
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


def join_sr_data(
    samples: TypedDataFrame[LongReadAlignmentMetadata],
    alignments: TypedDataFrame[AlignmentMetadataLong],
    short_read_terra_workspace: TerraWorkspace,
) -> TypedDataFrame[LongReadTerraSamples]:
    """
    Collect data from Gumbo for short read RNA samples and join them to the long read
    samples.

    :param samples: a data frame of long read samples
    :param alignments: a data frame of sequencing_alignment metadata from Gumbo
    :param short_read_terra_workspace: the Terra workspace for short read RNA samples
    :return: the long read samples with short read sample metaadata joined to it
    """

    # collect short read metadata from Gumbo
    sr_sequencings = alignments.loc[alignments["datatype"].eq("rna")]

    # get short read outputs from Terra
    sr_terra_samples = get_sr_terra_samples(short_read_terra_workspace)

    # join the best short read sample available to each long read sample
    samples_w_sr = choose_matched_short_read_sample(
        samples, sr_sequencings, sr_terra_samples
    )

    return type_data_frame(samples_w_sr, LongReadTerraSamples)


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
            ["sample_id", "star_junctions"],
        ]
        .rename(columns={"star_junctions": "sr_star_junctions"})
        .astype("string")
    )

    return type_data_frame(sr_terra_samples, ShortReadTerraSamples)


def choose_matched_short_read_sample(
    samples: pd.DataFrame,
    sr_sequencings: TypedDataFrame[AlignmentMetadataLong],
    sr_terra_samples: TypedDataFrame[ShortReadTerraSamples],
) -> TypedDataFrame[LongReadTerraSamples]:
    """
    Choose short read RNA samples to map to long read ones.

    :param samples: a data frame of long read samples
    :param sr_sequencings: a data frame of short read sequencings from Gumbo
    :param sr_terra_samples: a data frame of short read outputs from Terra
    :return: long read samples with short read data joined when possible
    """

    # collect IDs long read side
    lr = samples.loc[
        :, ["sample_id", "model_id", "model_condition_id"]
    ].drop_duplicates()

    # collect IDs (ones needed for joining or unique identification) on short read side
    sr = sr_sequencings.loc[
        :,
        [
            "model_id",
            "model_condition_id",
            # (we're actually choosing a sequencing alignment record on this side)
            "sequencing_alignment_id",
            "sequencing_alignment_source",
            "priority",
        ],
    ].drop_duplicates()

    # join on model_id (model_condition_id might or might not match)
    pairs = lr.merge(sr, how="inner", on="model_id", suffixes=("_lr", "_sr"))

    # add another priority column for whether the LR-SR model conditions match
    pairs["mc_priority"] = 1
    pairs.loc[
        pairs["model_condition_id_lr"].ne(pairs["model_condition_id_sr"]), "mc_priority"
    ] = 2

    # add another priority column for choosing the GP-delivered CRAM/BAM over our
    # analysis-ready one (since at one point we were dropping unmapped reads)
    pairs["sequencing_alignment_priority"] = pd.NA
    pairs.loc[
        pairs["sequencing_alignment_source"].eq("GP"), "sequencing_alignment_priority"
    ] = 1
    pairs.loc[
        pairs["sequencing_alignment_source"].eq("CDS"), "sequencing_alignment_priority"
    ] = 2
    assert bool(pairs["sequencing_alignment_priority"].notna().all())

    # pick the best option within matching model condition when possible, then use the
    # omics_mapping priority, then a choose GP vs. CDS alignment file
    choices = (
        pairs.sort_values(["mc_priority", "priority", "sequencing_alignment_priority"])
        .groupby("sample_id")
        .nth(0)
    )

    # collect SR columns to include in the final LR sample data table
    sr_choices = (
        choices.merge(sr_sequencings, how="inner", on="sequencing_alignment_id")
        .loc[
            :,
            [
                "sample_id",
                "model_condition_id",
                "omics_profile_id",
                "omics_sequencing_id",
                "sequencing_alignment_id",
                "crai_bai_url",
                "cram_bam_url",
            ],
        ]
        .rename(
            columns={
                "model_condition_id": "sr_model_condition_id",
                "omics_profile_id": "sr_omics_profile_id",
                "omics_sequencing_id": "sr_omics_sequencing_id",
                "sequencing_alignment_id": "sr_sequencing_alignment_id",
                "crai_bai_url": "sr_crai_bai",
                "cram_bam_url": "sr_cram_bam",
            }
        )
    )

    sr_choices["sr_file_format"] = (
        sr_choices["sr_cram_bam"].str.rsplit(".", n=1).str.get(1).str.upper()
    )

    # can now finally join SR data to LR data
    samples = samples.merge(sr_choices, how="left", on="sample_id")

    # join SR output files
    samples = samples.merge(
        sr_terra_samples.rename(
            columns={
                "sample_id": "sr_omics_sequencing_id",
                "star_junctions": "sr_star_junctions",
            }
        ),
        how="left",
        on="sr_omics_sequencing_id",
    )

    return type_data_frame(samples, LongReadTerraSamples, remove_unknown_cols=True)
