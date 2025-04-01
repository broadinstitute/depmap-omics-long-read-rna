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
    ShortReadPriorities,
)
from depmap_omics_long_read_rna.utils.utils import model_to_df


def do_refresh_terra_samples(
    terra_workspace: TerraWorkspace,
    gumbo_client: GumboClient,
) -> None:
    """
    Upsert the sample data table in Terra with the latest long read sample metadata from
    Gumbo and the short read Terra workspace.

    :param terra_workspace: the Terra workspace for long read RNA samples
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
    samples = join_sr_metadata(samples, alignments)

    terra_workspace.upload_entities(df=samples)


def explode_and_expand_models(
    models: TypedDataFrame[ModelsAndChildren],
) -> TypedDataFrame[AlignmentMetadataLong]:
    """
    Unnest columns from the GraphQL call for Gumbo models and their childen.

    :param models: a data frame of models with nested profile/sequencing/etc. data
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


def join_sr_metadata(
    samples: TypedDataFrame[LongReadAlignmentMetadata],
    alignments: TypedDataFrame[AlignmentMetadataLong],
) -> TypedDataFrame[LongReadTerraSamples]:
    """
    Collect metadata from Gumbo for short read RNA samples and join them to the long
    read samples.

    :param samples: a data frame of long read samples
    :param alignments: a data frame of sequencing_alignment metadata from Gumbo
    :return: the long read samples with short read sample metaadata joined to it
    """

    # collect short read metadata from Gumbo
    sr_sequencings = alignments.loc[alignments["datatype"].eq("rna")]

    # join the best short read sample available to each long read sample
    sr_priorities = prioritize_short_read_sample(samples, sr_sequencings)
    samples_w_metadata = choose_matched_short_read_sample(samples, sr_priorities)

    return type_data_frame(samples_w_metadata, LongReadTerraSamples)


def prioritize_short_read_sample(
    samples: pd.DataFrame, sr_sequencings: pd.DataFrame
) -> TypedDataFrame[ShortReadPriorities]:
    """
    Add columns to the data frame of short read samples to be used for sorting and
    grouping, and thus allow choosing the "best" matched short read sample for each
    long read sample.

    :param samples: a data frame of long read samples
    :param sr_sequencings: a data frame of short read data from Gumbo
    :return: `sr_sequencings` with `_priority` columns
    """

    # collect short read samples in Gumbo matching long read samples' model IDs
    sr = (
        sr_sequencings.loc[sr_sequencings["model_id"].isin(samples["model_id"])]
        .copy()
        .rename(
            columns={
                "omics_profile_id": "sr_omics_profile_id",
                "omics_sequencing_id": "sr_omics_sequencing_id",
                "sequencing_alignment_id": "sr_sequencing_alignment_id",
                "crai_bai_url": "sr_crai_bai",
                "cram_bam_url": "sr_cram_bam",
            }
        )
    )

    # clean the `*_priority` columns
    for c in ["drug", "source", "expansion_team", "sequencing_alignment_source"]:
        sr[c] = sr[c].str.lower().fillna("")

    # dict values are in order of preference (higher to lower)
    priorities = {
        "drug": ["", "dmso"],
        "stranded": [True, False],
        "source": ["broad", "sanger", "stjude", "other", "dfci", ""],
        "expansion_team": [
            "dmx",
            "tda",
            "gpp",
            "other",
            "cclf",
            "stjude",
            "sanger",
            "none",
            "unknown",
            "",
        ],
        "sequencing_alignment_source": ["gp", "cds"],
    }

    for c, vals in priorities.items():
        priorities_df = pd.DataFrame({c: vals, f"{c}_priority": list(range(len(vals)))})
        sr = sr.merge(priorities_df, how="left", on=c)

    return type_data_frame(sr, ShortReadPriorities, remove_unknown_cols=True)


def choose_matched_short_read_sample(
    samples: pd.DataFrame, sr_priorities: pd.DataFrame
) -> TypedDataFrame[LongReadTerraSamples]:
    """
    Choose short read RNA samples to map to long read ones.

    :param samples: a data frame of long read samples
    :param sr_priorities: a data frame of short read data from
    `prioritize_short_read_sample`
    :return: long read samples with short read data joined when possible
    """

    # collect base long read metadata for joining
    lr = samples.loc[:, ["model_id", "model_condition_id"]].drop_duplicates()

    # get top short read sample within matching model condition
    lr_sr = (
        lr.merge(sr_priorities, how="inner", on=["model_id", "model_condition_id"])
        .sort_values(
            [
                "model_condition_id",
                "drug_priority",
                "stranded_priority",
                "source_priority",
                "expansion_team_priority",
                "sequencing_alignment_source_priority",
                "version",  # tiebreaker (if multiple sequencings per profile)
            ],
            ascending=[True, True, True, True, True, True, False],
        )
        .groupby("model_condition_id")
        .nth(0)
    )

    # LR-SR model conditions are the same in this joined data frame
    lr_sr["sr_model_condition_id"] = lr_sr["model_condition_id"]

    # if we were unable to find a match based on common model condition, go up to the
    # model level and try again
    lr_todo = lr.loc[~lr["model_id"].isin(lr_sr["model_id"])]

    if len(lr_todo) > 0:
        sr_todo = sr_priorities.loc[sr_priorities["model_id"].isin(lr_todo["model_id"])]

        lr_sr_model = (
            lr_todo.merge(
                sr_todo.rename(columns={"model_condition_id": "sr_model_condition_id"}),
                how="inner",
                on="model_id",
            )
            .sort_values(
                [
                    "model_id",
                    "drug_priority",
                    "stranded_priority",
                    "source_priority",
                    "expansion_team_priority",
                    "sequencing_alignment_source_priority",
                    "version",  # tiebreaker (if multiple sequencings per profile)
                ],
                ascending=[True, True, True, True, True, True, False],
            )
            .groupby("model_id")
            .nth(0)
        )

        lr_sr = pd.concat([lr_sr, lr_sr_model], ignore_index=True)

    # can finally join SR data to LR data
    samples = samples.merge(lr_sr, how="left", on=["model_id", "model_condition_id"])

    return type_data_frame(samples, LongReadTerraSamples, remove_unknown_cols=True)
