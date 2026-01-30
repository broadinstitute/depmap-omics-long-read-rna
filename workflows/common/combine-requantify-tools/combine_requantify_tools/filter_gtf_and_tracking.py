import logging
from pathlib import Path

import pandas as pd

from combine_requantify_tools.types import (
    Gtf,
    Sqanti3Classification,
    TypedDataFrame,
    UpdatedTracking,
)
from combine_requantify_tools.utils import type_data_frame


def filter_tracking(
    tracking_in: Path, squanti_classification: Path, prefix: str
) -> tuple[TypedDataFrame[UpdatedTracking], TypedDataFrame[Sqanti3Classification]]:
    """
    Load the tracking file (output from process_tracking_file) and the Sqanti3
    classification file (output from combine_gtfs.run_sqanti3 task) TODO

    :param tracking_in: output from process_tracking_file
    :param squanti_classification: Sqanti3 classification file
    :param prefix: prefix previously used for gffcompare
    :return: filtered tracking and classification data frames
    """

    # load updated tracking
    logging.info(f"Loading tracking file from {tracking_in}")
    updated_tracking = type_data_frame(pd.read_parquet(tracking_in), UpdatedTracking)

    # collect observed transcript IDs (i.e. "TCONS_*")
    updated_tracking["transcript_id"] = (
        updated_tracking["transcript_id"].str.split("|").str.get(0)
    )
    transcript_ids = updated_tracking["transcript_id"].drop_duplicates()

    # load Sqanti3 classification file
    logging.info(f"Loading Sqanti3 classification file from {squanti_classification}")
    sq = type_data_frame(
        pd.read_table(
            squanti_classification,
            sep="\t",
            na_values=["NA"],
            dtype={
                "isoform": "string",
                "structural_category": "string",
                "associated_transcript": "string",
                "RTS_stage": "boolean",
                "coding": "string",
                "ORF_seq": "string",
                # other columns aren't needed
            },
            usecols=[
                "isoform",
                "structural_category",
                "associated_transcript",
                "RTS_stage",
                "coding",
                "ORF_seq",
            ],
        ),
        Sqanti3Classification,
    )

    # split off known isoforms
    sq_annotated = sq.loc[sq["isoform"].str.startswith("ENST")]

    # split off TCONS + novel_*_catalog + coding
    sq_filtered = sq.loc[
        sq["isoform"].str.startswith(prefix)
        & sq["structural_category"].isin(["novel_not_in_catalog", "novel_in_catalog"])
        & (sq["coding"] == "coding")
    ]

    # use ORF to join known isoform IDs to the filtered records
    merged_sq = sq_filtered.merge(
        sq_annotated.loc[:, ["ORF_seq", "isoform"]],
        on="ORF_seq",
        how="left",
        suffixes=("_tcons", "_enst"),
        indicator=True,
    )

    # identify ORFs in novel but not in annotated
    merged_sq = merged_sq.loc[merged_sq["_merge"].eq("left_only")]

    # TODO
    sq_filtered = sq_filtered.loc[
        sq_filtered["isoform"].isin(merged_sq["isoform_tcons"])
        & ~sq_filtered["RTS_stage"]
        & sq_filtered["isoform"].isin(transcript_ids)
    ]

    # TODO
    sq_existing = sq.loc[
        sq["structural_category"].isin(["full-splice_match"])
        & ~(
            (sq["associated_transcript"].str.startswith("ENST"))
            | (sq["associated_transcript"] == "novel")
        )
        & (sq["coding"] == "coding")
    ]

    sq_filtered_tracking = pd.concat([sq_filtered, sq_existing], ignore_index=True)

    updated_tracking = updated_tracking.loc[
        updated_tracking["transcript_id"].isin(sq_filtered_tracking["isoform"])
    ]

    return type_data_frame(updated_tracking, UpdatedTracking), type_data_frame(
        sq_filtered, Sqanti3Classification
    )


def filter_gtf(
    gtf_in: Path, sq_filtered: TypedDataFrame[Sqanti3Classification]
) -> TypedDataFrame[Gtf]:
    """
    Filter combined GTF to features observed in filtered Sqanti3 classification file.

    :param gtf_in: combined GTF
    :param sq_filtered: filtered Sqanti3 classification from filter_tracking
    :return: filtered GTF data frame
    """

    logging.info(f"Loading GTF from {gtf_in}")

    gtf = type_data_frame(
        pd.read_csv(
            gtf_in,
            sep="\t",
            comment="#",
            header=None,
            names=[
                "seqname",
                "source",
                "feature",
                "start",
                "end",
                "score",
                "strand",
                "frame",
                "attribute",
            ],
            dtype={
                "seqname": "string",
                "source": "string",
                "feature": "string",
                "start": "int64",
                "end": "int64",
                "score": "string",
                "strand": "string",
                "frame": "string",
                "attribute": "string",
            },
        ),
        Gtf,
    )

    # remove features without a strand value (+/-)
    gtf = gtf.loc[gtf["strand"].ne(".")]

    # extract transcript_id from attributes column
    gtf["transcript_id"] = gtf["attribute"].str.extract(r'transcript_id "([^"]+)"')
    assert bool(gtf["transcript_id"].notna().all())

    # filter GTF
    gtf = gtf.loc[gtf["transcript_id"].isin(sq_filtered["isoform"])]
    gtf = gtf.drop(columns=["transcript_id"])

    return type_data_frame(gtf, Gtf)
