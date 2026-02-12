import logging
from pathlib import Path

import pandas as pd

from combine_requantify_tools.types import (
    Gtf,
    TranscriptCounts,
    TypedDataFrame,
    UpdatedTracking,
)
from combine_requantify_tools.utils import read_gtf_from_file, type_data_frame


def do_filter_sample_gtf(
    gtf_in: Path, tracking_in: Path, sample_id: str, transcript_counts_in: Path
) -> TypedDataFrame[Gtf]:
    """
    Filter input GTF based nonzero transcripts in a single sample's transcript counts
    file.

    :param gtf_in: output GTF from `process_gtf`
    :param tracking_in: output tracking from `filter_gtf_and_tracking`
    :param sample_id: sample ID associated with this transcript counts file
    :param transcript_counts_in: path to sample's transcript counts file from initial
    Isoquant run
    :return: filtered GTF file
    """

    gtf = read_gtf_from_file(gtf_in)

    # extract transcript and gene IDs from attributes
    gtf["transcript_id"] = gtf["attribute"].str.extract(r'transcript_id "([^"]+)"')
    gtf["gene_id"] = gtf["attribute"].str.extract(r'gene_id "([^"]+)"')

    # filter tracking file for this sample
    logging.info(f"Reading {tracking_in}")
    tracking = type_data_frame(pd.read_parquet(tracking_in), UpdatedTracking)
    tracking = tracking.loc[tracking["sample"].eq(sample_id)]

    # get transcripts with non-zero counts
    logging.info(f"Reading {transcript_counts_in}")
    tc = type_data_frame(
        pd.read_csv(
            transcript_counts_in,
            sep="\t",
            comment="#",
            header=None,
            names=["id1", "count"],
            compression="gzip",
            dtype={"id1": "string", "count": "int64"},
        ),
        TranscriptCounts,
    )

    # filter out zero counts
    tc = tc.loc[tc["count"].gt(0)]

    # combine transcript IDs from tracking and sample's counts
    transcript_ids = pd.concat(
        [tracking["transcript_id"], tc["id1"]], ignore_index=True
    ).drop_duplicates()

    # get transcript and gene IDs for the transcripts we want to keep
    match_tx_gene_id = gtf.loc[
        gtf["transcript_id"].isin(transcript_ids), ["transcript_id", "gene_id"]
    ].drop_duplicates()

    assert bool(match_tx_gene_id.notna().all(axis=None))

    # filter GTF by observed genes and transcripts
    gtf_updated = gtf.loc[
        (gtf["feature"].eq("gene") & (gtf["gene_id"].isin(match_tx_gene_id["gene_id"])))
        | gtf["transcript_id"].isin(transcript_ids)
    ]

    return type_data_frame(gtf_updated, Gtf, remove_unknown_cols=True)
