from pathlib import Path

import pandas as pd


def filter_tracking(
    tracking_in: Path, squanti_classification: Path, prefix: str
) -> tuple[pd.DataFrame, pd.DataFrame]:
    # Load updated tracking
    updated_tracking = pd.read_parquet(tracking_in)

    transcript_ids = updated_tracking["transcript_id"].drop_duplicates()
    transcript_ids = transcript_ids.str.split("|").str.get(0)

    # Load classification
    sq = pd.read_table(
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
        },
        usecols=[
            "isoform",
            "structural_category",
            "associated_transcript",
            "RTS_stage",
            "coding",
            "ORF_seq",
        ],
        low_memory=False,
    )

    sq_annotated = sq[sq["isoform"].str.startswith("ENST")]

    # Filter TCONS + novel_*_catalog + coding
    sq_filtered = sq[
        sq["isoform"].str.startswith(prefix)
        & sq["structural_category"].isin(["novel_not_in_catalog", "novel_in_catalog"])
        & (sq["coding"] == "coding")
    ]
    sq_annotated_sm = sq_annotated[["ORF_seq", "isoform"]]
    merged_sq = sq_filtered.merge(
        sq_annotated_sm,
        on="ORF_seq",
        how="left",
        suffixes=("_tcons", "_enst"),
        indicator=True,
    )
    merged_sq = merged_sq[
        merged_sq["_merge"] == "left_only"
    ]  # orf in novel and not in annotated

    sq_filtered = sq_filtered[sq_filtered["isoform"].isin(merged_sq["isoform_tcons"])]
    sq_filtered = sq_filtered[sq_filtered["RTS_stage"] == False]  # Filter for RTS_stage
    sq_filtered = sq_filtered[sq_filtered["isoform"].isin(transcript_ids)]

    sq_filtered_previouslyfound = sq[
        sq["structural_category"].isin(["full-splice_match"])
        & ~(
            (sq["associated_transcript"].str.startswith("ENST"))
            | (sq["associated_transcript"] == "novel")
        )
        & (sq["coding"] == "coding")
    ]

    sq_filtered_tracking = pd.concat([sq_filtered, sq_filtered_previouslyfound])

    updated_tracking["transcript_id"] = (
        updated_tracking["transcript_id"].str.split("|").str.get(0)
    )
    updated_tracking = updated_tracking[
        updated_tracking["transcript_id"].isin(sq_filtered_tracking["isoform"])
    ]

    return updated_tracking, sq_filtered


def filter_gtf(
    annotation_filtered_gtf: Path, sq_filtered: pd.DataFrame
) -> pd.DataFrame:
    gtf = pd.read_csv(
        annotation_filtered_gtf,
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
    )

    # Extract transcript_id from attributes column
    gtf["transcript_id"] = gtf["attribute"].str.extract(r'transcript_id "([^"]+)"')
    assert bool(gtf["transcript_id"].notna().all())

    # Filter GTF
    gtf = gtf.loc[gtf["transcript_id"].isin(sq_filtered["isoform"])]
    gtf = gtf.drop(columns=["transcript_id"])

    return gtf
