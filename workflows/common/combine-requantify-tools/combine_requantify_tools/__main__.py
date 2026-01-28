import logging
import re
from concurrent.futures import ALL_COMPLETED, ThreadPoolExecutor, as_completed, wait
from pathlib import Path
from typing import Annotated, Any

import pandas as pd
import typer
from tqdm.auto import tqdm

pd.set_option("display.max_columns", 30)
pd.set_option("display.max_colwidth", 50)
pd.set_option("display.max_info_columns", 30)
pd.set_option("display.max_info_rows", 20)
pd.set_option("display.max_rows", 20)
pd.set_option("display.max_seq_items", None)
pd.set_option("display.width", 200)
pd.set_option("expand_frame_repr", True)
pd.set_option("mode.chained_assignment", "warn")

app = typer.Typer(rich_markup_mode=None, pretty_exceptions_enable=False)

config: dict[str, Any] = {}


# noinspection PyUnusedLocal
def done(*args, **kwargs):
    logging.info("Done.")


@app.command()
def process_tracking_file(
    tracking_in: Annotated[
        Path,
        typer.Option(
            "path to input tracking file (from combine_gtfs.run_gffcompare task)"
        ),
    ],
    tracking_out: Annotated[Path, typer.Option("path to write output tracking file")],
    sample_ids_list: Annotated[
        Path, typer.Option("path to text file containing list of sample IDs")
    ],
    transcript_counts_file_list: Annotated[
        Path,
        typer.Option(
            "path to text file containing list paths to samples' transcript count files"
        ),
    ],
    min_count: Annotated[
        int, typer.Option("min value of count to filter transcript counts")
    ],
) -> None:
    import pandas as pd

    # read file containing list of sample IDs
    logging.info(f"Reading {sample_ids_list}")
    with open(sample_ids_list, "r") as f:
        sample_ids = f.read().splitlines()

    sample_ids = [x.strip() for x in sample_ids]

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

    # read file containing list of transcript count files
    logging.info(f"Reading {transcript_counts_file_list}")
    with open(transcript_counts_file_list, "r") as f:
        transcript_counts = f.read().splitlines()

    transcript_counts = [x.strip() for x in transcript_counts]

    def extract_sample_id(path: str) -> str:
        """
        Extract the sample/sequencing ID from a transcript counts file's name, e.g.
        "CDS-ABCDEF" from "path/CDS-ABCDEF.discovered_transcript_tpm.tsv.gz".

        :param path: path to a transcript counts file
        :return: the CDS-* sample ID
        """

        m = re.findall(r"CDS-[A-Z a-z 0-9]{6}", path)
        assert len(m) == 1, f"Couldn't find single sample ID in '{path}'"
        return m[0]

    # make mapping from sample ID to transcript count file
    sample_to_tpm_file = {extract_sample_id(path): path for path in transcript_counts}

    assert len(set(sample_to_tpm_file.keys()).symmetric_difference(sample_ids)) == 0, (
        "Provided sample IDs and transcript count files don't match"
    )

    # iterate over transcript count files
    logging.info(f"Reading transcript counts for {len(sample_ids)} samples")

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

    with ThreadPoolExecutor() as executor:
        futures = [
            executor.submit(process_single_sample, sample_id, tpm_file)
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

    updated_tracking.to_parquet(tracking_out, index=False)


@app.callback(result_callback=done)
def main():
    logger = logging.getLogger()
    logger.addHandler(logging.StreamHandler())
    logger.setLevel(logging.INFO)


if __name__ == "__main__":
    app()
