import csv
import logging
from pathlib import Path
from typing import Annotated, Any

import pandas as pd
import typer

from combine_requantify_tools.filter_gtf_and_tracking import filter_gtf, filter_tracking
from combine_requantify_tools.map_transcript_ids import do_map_transcript_ids
from combine_requantify_tools.process_gtf import do_process_gtf
from combine_requantify_tools.process_tracking_file import do_process_tracking_file
from combine_requantify_tools.utils import read_sample_ids_from_file

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
def map_transcript_ids(
    tracking_in: Annotated[
        Path,
        typer.Option(help="path to input tracking file (from gffcompare)"),
    ],
    gtf_in: Annotated[
        Path,
        typer.Option(help="path to input GTF file (from gffcompare)"),
    ],
    sample_ids_list: Annotated[
        Path, typer.Option(help="path to text file containing list of sample IDs")
    ],
    tracking_out: Annotated[
        Path, typer.Option(help="path to write output tracking file")
    ],
    gtf_out: Annotated[
        Path,
        typer.Option(help="path to  write output GTF file"),
    ],
) -> None:
    # read file containing list of sample IDs
    sample_ids = read_sample_ids_from_file(sample_ids_list)

    tracking_mapped, gtf_mapped = do_map_transcript_ids(tracking_in, gtf_in, sample_ids)

    logging.info(f"Writing {tracking_out}")
    tracking_mapped.to_parquet(tracking_out, index=False)

    logging.info(f"Writing {gtf_out}")
    gtf_mapped.to_csv(
        gtf_out, sep="\t", index=False, header=False, quoting=csv.QUOTE_NONE
    )


@app.command()
def process_tracking_file(
    tracking_in: Annotated[
        Path,
        typer.Option(help="path to input tracking file (from process_tracking_file)"),
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
    tracking_out: Annotated[
        Path, typer.Option(help="path to write output tracking file")
    ],
) -> None:
    # read file containing list of sample IDs
    sample_ids = read_sample_ids_from_file(sample_ids_list)

    # read file containing list of transcript count files
    logging.info(f"Reading {discovered_transcript_counts_file_list}")
    with open(discovered_transcript_counts_file_list, "r") as f:
        discovered_transcript_counts = f.read().splitlines()

    discovered_transcript_counts = [x.strip() for x in discovered_transcript_counts]

    # process the tracking file
    updated_tracking = do_process_tracking_file(
        tracking_in, sample_ids, discovered_transcript_counts, min_count
    )

    updated_tracking.to_parquet(tracking_out, index=False)


@app.command()
def filter_gtf_and_tracking(
    tracking_in: Annotated[
        Path,
        typer.Option(help="path to updated tracking file (from process_tracking_file)"),
    ],
    squanti_classification: Annotated[
        Path, typer.Option(help="path to SQANTI3 classification file")
    ],
    gtf_in: Annotated[Path, typer.Option(help="path to input filtered GTF file")],
    prefix: Annotated[str, typer.Option(help="transcript prefix (e.g. 'TCONS')")],
    tracking_out: Annotated[
        Path, typer.Option(help="path to write filtered tracking file")
    ],
    gtf_out: Annotated[Path, typer.Option(help="path to write filtered GTF file")],
) -> None:
    updated_tracking, sq_filtered = filter_tracking(
        tracking_in, squanti_classification, prefix
    )

    gtf_filtered = filter_gtf(gtf_in, sq_filtered)

    logging.info(f"Writing {tracking_out}")
    updated_tracking.to_parquet(tracking_out, index=False)

    logging.info(f"Writing {gtf_out}")
    gtf_filtered.to_csv(
        gtf_out, sep="\t", index=False, header=False, quoting=csv.QUOTE_NONE
    )


@app.command()
def process_gtf(
    gtf_in: Annotated[Path, typer.Option(help="path to input GTF file")],
    gtf_out: Annotated[Path, typer.Option(help="path to write GTF file")],
) -> None:
    gtf_processed = do_process_gtf(gtf_in)

    logging.info(f"Writing {gtf_out}")
    gtf_processed.to_csv(
        gtf_out, sep="\t", index=False, header=False, quoting=csv.QUOTE_NONE
    )


@app.callback(result_callback=done)
def main():
    logger = logging.getLogger()
    logger.addHandler(logging.StreamHandler())
    logger.setLevel(logging.INFO)


if __name__ == "__main__":
    app()
