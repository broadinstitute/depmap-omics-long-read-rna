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
            help="path to input tracking file (from combine_gtfs.run_gffcompare task)"
        ),
    ],
    tracking_out: Annotated[
        Path, typer.Option(help="path to write output tracking file")
    ],
    sample_ids_list: Annotated[
        Path, typer.Option(help="path to text file containing list of sample IDs")
    ],
    discovered_transcript_counts_file_list: Annotated[
        Path,
        typer.Option(
            help="path to text file containing paths to samples' "
            "discovered_transcript_counts files (from quantify_lr_rna workflow)"
        ),
    ],
    min_count: Annotated[int, typer.Option(help="min value to use to filter counts")],
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
    logging.info(f"Reading {discovered_transcript_counts_file_list}")
    with open(discovered_transcript_counts_file_list, "r") as f:
        discovered_transcript_counts = f.read().splitlines()

    discovered_transcript_counts = [x.strip() for x in discovered_transcript_counts]

    def extract_sample_id(path: str) -> str:
        """
        Extract the sample/sequencing ID from a transcript counts file's name, e.g.
        "CDS-ABCDEF" from "path/CDS-ABCDEF/CDS-ABCDEF.discovered_transcript_tpm.tsv.gz".

        :param path: path to a transcript counts file
        :return: the CDS-* sample ID
        """

        m = re.findall(r"CDS-[A-Z a-z 0-9]{6}", path)
        assert len(set(m)) == 1, f"Couldn't find single sample ID in '{path}'"
        return m[0]

    # make mapping from sample ID to transcript count file
    sample_to_tpm_file = {
        extract_sample_id(path): path for path in discovered_transcript_counts
    }

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


@app.command()
def filter_gtf_and_tracking(
    updated_tracking: Annotated[
        Path, typer.Option(help="path to updated tracking file (parquet format)")
    ],
    squanti_classification: Annotated[
        Path, typer.Option(help="path to SQANTI3 classification file")
    ],
    annotation_filtered_gtf: Annotated[
        Path, typer.Option(help="path to input filtered GTF file")
    ],
    prefix: Annotated[str, typer.Option(help="transcript prefix (e.g., TCONS)")],
    updated_tracking_out: Annotated[
        Path, typer.Option(help="path to write filtered tracking file")
    ],
    filtered_gtf_out: Annotated[
        Path, typer.Option(help="path to write filtered GTF file")
    ],
) -> None:
    import pandas as pd

    # Load updated tracking
    updated_tracking_df = pd.read_parquet(updated_tracking)
    updated_tracking_nodups = updated_tracking_df[["transcript_id"]].drop_duplicates()
    updated_tracking_nodups["transcript_id"] = (
        updated_tracking_nodups["transcript_id"].str.split("|").str[0]
    )

    # Load classification
    sq = pd.read_csv(squanti_classification, sep="\t")
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
    sq_filtered = sq_filtered[
        sq_filtered["isoform"].isin(updated_tracking_nodups["transcript_id"])
    ]

    sq_filtered_previouslyfound = sq[
        sq["structural_category"].isin(["full-splice_match"])
        & ~(
            (sq["associated_transcript"].str.startswith("ENST"))
            | (sq["associated_transcript"] == "novel")
        )
        & (sq["coding"] == "coding")
    ]

    sq_filtered_tracking = pd.concat([sq_filtered, sq_filtered_previouslyfound])

    updated_tracking_df["transcript_id"] = (
        updated_tracking_df["transcript_id"].str.split("|").str[0]
    )
    updated_tracking_df = updated_tracking_df[
        updated_tracking_df["transcript_id"].isin(sq_filtered_tracking["isoform"])
    ]
    updated_tracking_df.to_csv(updated_tracking_out, sep="\t", index=False)

    # Load GTF
    gtf = pd.read_csv(annotation_filtered_gtf, sep="\t", comment="#", header=None)

    # Extract transcript_id from attributes column
    gtf["transcript_id"] = gtf[8].str.extract(r'transcript_id "([^"]+)"')

    # Filter GTF
    gtf_filtered = gtf[gtf["transcript_id"].isin(sq_filtered["isoform"])].copy()
    gtf_filtered.drop(columns=["transcript_id"], inplace=True)

    # Write output
    gtf_filtered.to_csv(
        filtered_gtf_out, sep="\t", index=False, header=False, quoting=3
    )

    # Add GTF header
    header = """##gff-version 3
##description: evidence-based annotation of the human genome (GRCh38), version 38 (Ensembl 104)
##provider: GENCODE
##contact: gencode-help@ebi.ac.uk
##format: gtf
##date: 2021-03-12"""

    with open(filtered_gtf_out, "r") as f:
        lines = f.readlines()
    with open(filtered_gtf_out, "w") as f:
        f.write(header + "\n")
        f.writelines(lines)


@app.callback(result_callback=done)
def main():
    logger = logging.getLogger()
    logger.addHandler(logging.StreamHandler())
    logger.setLevel(logging.INFO)


if __name__ == "__main__":
    app()
