import logging
import re
from concurrent.futures import ALL_COMPLETED, ThreadPoolExecutor, as_completed, wait
from pathlib import Path

import pandas as pd
from tqdm.auto import tqdm


def do_process_tracking_file(
    tracking_in: Path,
    sample_ids: list[str],
    discovered_transcript_counts: list[str],
    min_count: int,
) -> pd.DataFrame:
    """
    Update the input tracking file with counts from the samples' transcript counts
    files.

    :param tracking_in: input tracking file path
    :param sample_ids: list of sample IDs
    :param discovered_transcript_counts: list of paths to transcript counts files
    :param min_count: min value to use to filter counts
    :return: updated tracking file
    """

    # read the input tracking file
    fixed_cols = ["transcript_id", "loc", "gene_id", "val"]

    logging.info(f"Reading {tracking_in}")
    tracking = pd.read_csv(
        tracking_in,
        sep="\t",
        header=None,
        names=[*fixed_cols, *sample_ids],
        na_values="-",
        dtype="string",
    )

    tracking["gene_id"] = tracking["gene_id"].fillna(".")

    # make data frame long
    tracking = tracking.melt(
        id_vars=fixed_cols, var_name="sample", value_name="id"
    ).dropna()

    tracking["gene_id"] = tracking["gene_id"].replace({".": pd.NA})
    tracking["id1"] = tracking["id"].str.split("|").str.get(1)

    # make mapping from sample ID to transcript count file
    sample_to_tpm_file = {
        extract_sample_id(path): path for path in discovered_transcript_counts
    }

    assert len(set(sample_to_tpm_file.keys()).symmetric_difference(sample_ids)) == 0, (
        "Provided sample IDs and transcript count files don't match"
    )

    def process_single_sample(sample_id: str, tpm_file: str) -> pd.DataFrame:
        """
        Read a single sample's transcript counts file, filter based on `min_count`, and
        return the matching subset of the input tracking file with the sample's `count`
        values joined.

        :param sample_id: a sample ID
        :param tpm_file: the transcript counts file for this sample
        :return: subset of the input tracking file with the sample's `count` values
        """

        tc = pd.read_csv(
            tpm_file,
            sep="\t",
            comment="#",
            header=None,
            names=["id1", "count"],
            compression="gzip",
            dtype={"id1": "string", "count": "int64"},
        )

        tracking_sample = tracking.loc[tracking["sample"].eq(sample_id)]
        tracking_sample = tracking_sample.merge(
            tc, on="id1", how="left", validate="one_to_one"
        )
        tracking_sample = tracking_sample.loc[tracking_sample["count"].ge(min_count)]

        return tracking_sample

    # iterate over transcript count files
    logging.info(f"Reading transcript counts for {len(sample_ids)} samples")

    with ThreadPoolExecutor() as executor:
        futures = [
            executor.submit(process_single_sample, sample_id, tpm_file, min_count)
            for sample_id, tpm_file in sample_to_tpm_file.items()
        ]

        with tqdm(total=len(futures)) as pbar:
            for _ in as_completed(futures):
                pbar.update(1)

        wait(futures, return_when=ALL_COMPLETED)

    logging.info(f"Combining data frames")
    updated_tracking = pd.concat(
        [x.result() for x in futures], ignore_index=True
    ).sort_values(["sample", "id"])

    return updated_tracking


def extract_sample_id(path: str) -> str:
    """
    Extract the sample/sequencing ID from a transcript counts file's name, e.g.
    "CDS-ABCDEF" from "path/CDS-ABCDEF/CDS-ABCDEF.discovered_transcript_tpm.tsv.gz".

    :param path: path to a transcript counts file
    :return: the CDS-* sample ID
    """

    m = re.findall(r"CDS-[A-Za-z0-9]{6}", path)
    assert len(set(m)) == 1, f"Couldn't find single sample ID in '{path}'"
    return m[0]
