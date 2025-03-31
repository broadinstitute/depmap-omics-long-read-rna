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
    ShortReadTerraSamples,
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
    Collect metadata from Gumbo and workflow outputs from Terra relating to short read
    RNA samples and join them to the long read samples.

    :param samples: a data frame of long read samples
    :param alignments: a data frame of sequencing_alignment metadata from Gumbo
    :return: the long read samples with short read sample metaadata joined to it
    """

    # collect short read metadata from Gumbo
    sr_sequencings = alignments.loc[alignments["datatype"].eq("rna")]

    # match a short read sample to each long read sample
    sr_priorities = prioritize_short_read_sample(samples, sr_sequencings)
    samples_w_metadata = choose_matched_short_read_sample(samples, sr_priorities)

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
        sr_terra_samples.loc[:, ["sample_id", "delivery_cram_bam"]]
        .rename(columns={"delivery_cram_bam": "sr_cram_bam"})
        .astype("string")
    )

    return type_data_frame(sr_terra_samples, ShortReadTerraSamples)


def prioritize_short_read_sample(
    samples: pd.DataFrame, sr_sequencings: pd.DataFrame
) -> TypedDataFrame[ShortReadPriorities]:
    """
    TODO

    :param samples: a data frame of long read samples
    :param sr_sequencings: a data frame of short read data from Gumbo
    :return: `sr_sequencings` with `_priority` columns
    """

    drug_priority = ["", "dmso"]

    drug_df = pd.DataFrame(
        {
            "drug": drug_priority,
            "drug_priority": list(range(len(drug_priority))),
        }
    )

    stranded = [True, False]

    stranded_df = pd.DataFrame(
        {"stranded": stranded, "stranded_priority": list(range(len(stranded)))}
    )

    source_priority = ["broad", "sanger", "stjude", "other", "dfci", ""]

    source_df = pd.DataFrame(
        {
            "source": source_priority,
            "source_priority": list(range(len(source_priority))),
        }
    )

    expansion_team_priority = [
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
    ]

    expansion_team_df = pd.DataFrame(
        {
            "expansion_team": expansion_team_priority,
            "expansion_team_priority": list(range(len(expansion_team_priority))),
        }
    )

    sequencing_alignment_source_priority = ["gp", "cds"]

    sequencing_alignment_source_df = pd.DataFrame(
        {
            "sequencing_alignment_source": sequencing_alignment_source_priority,
            "sequencing_alignment_source_priority": list(
                range(len(sequencing_alignment_source_priority))
            ),
        }
    )

    # collect all valid short read samples in Gumbo that are present in both short and
    # long read workspaces
    sr = sr_sequencings.loc[sr_sequencings["model_id"].isin(samples["model_id"])].copy()

    # populate the `*_priority` columns
    for c in ["drug", "source", "expansion_team", "sequencing_alignment_source"]:
        sr[c] = sr[c].str.lower().fillna("")

    sr = (
        sr.merge(drug_df, how="left", on="drug")
        .merge(stranded_df, how="left", on="stranded")
        .merge(source_df, how="left", on="source")
        .merge(expansion_team_df, how="left", on="expansion_team")
        .merge(
            sequencing_alignment_source_df, how="left", on="sequencing_alignment_source"
        )
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

    return type_data_frame(sr, ShortReadPriorities, remove_unknown_cols=True)


def choose_matched_short_read_sample(
    samples: pd.DataFrame, sr_priorities: pd.DataFrame
) -> TypedDataFrame[LongReadTerraSamples]:
    """
    Choose short read RNA samples to map to long read ones.

    :param samples: a data frame of long read samples
    :param sr_priorities: TODO
    :return: TODO
    """

    lr = samples[["model_id", "model_condition_id"]].drop_duplicates()

    lr_sr_mc = (
        lr.merge(sr_priorities, how="inner", on=["model_id", "model_condition_id"])
        .sort_values(
            [
                "model_id",
                "model_condition_id",
                "drug_priority",
                "stranded_priority",
                "source_priority",
                "expansion_team_priority",
                "sequencing_alignment_source_priority",
                "version",
            ],
            ascending=[True, True, True, True, True, True, True, False],
        )
        .groupby("model_condition_id")
        .nth(0)
    )

    lr_sr_mc["sr_model_condition_id"] = lr_sr_mc["model_condition_id"]

    lr_todo = lr.loc[~lr["model_id"].isin(lr_sr_mc["model_id"])]
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
                "version",
            ],
            ascending=[True, True, True, True, True, True, False],
        )
        .groupby("model_id")
        .nth(0)
    )

    lr_sr_all = pd.concat([lr_sr_mc, lr_sr_model], ignore_index=True)

    samples = samples.merge(
        lr_sr_all, how="left", on=["model_id", "model_condition_id"]
    )

    return type_data_frame(samples, LongReadTerraSamples, remove_unknown_cols=True)
