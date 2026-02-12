import logging
from pathlib import Path
from typing import Type

import pandas as pd

from combine_requantify_tools.types import (
    Gtf,
    PanderaBaseSchema,
    Sqanti3Classification,
    Tracking,
    TypedDataFrame,
)


def type_data_frame(
    df: pd.DataFrame,
    pandera_schema: Type[PanderaBaseSchema],
    reorder_cols: bool = True,
    remove_unknown_cols: bool = False,
) -> TypedDataFrame[PanderaBaseSchema]:
    """
    Coerce a data frame into one specified by a Pandera schema and optionally reorder
    columns and remove unknown columns.

    :param df: a data frame
    :param pandera_schema: a Pandera schema
    :param reorder_cols: reorder columns as specified in the schema
    :param remove_unknown_cols: remove columns not specified in the schema
    :return: a data frame validated with the provided Pandera schema
    """

    if len(df) == 0:
        # make an empty data frame that conforms to the Pandera schema
        s = pandera_schema.to_schema()

        # `example` doesn't know how to instantiate columns with structured data
        dict_cols = []
        list_cols = []

        for c in s.columns:
            if s.columns[c].dtype.type is dict:
                dict_cols.append(c)
                s = s.remove_columns([c])
            elif s.columns[c].dtype.type is list:
                list_cols.append(c)
                s = s.remove_columns([c])

        df = pd.DataFrame(s.example(size=1))

        if len(dict_cols) > 0:
            for c in dict_cols:
                df[c] = [{}] * len(df)

        if len(list_cols) > 0:
            for c in list_cols:
                df[c] = [[]] * len(df)

        df = df.iloc[:0]
        return TypedDataFrame[pandera_schema](df)

    if not remove_unknown_cols and not reorder_cols:
        # can type and return
        return TypedDataFrame[pandera_schema](df)

    # we need to collect the current columns and schema columns (in original orders)
    current_cols = list(df.columns)
    schema_cols = list(pandera_schema.to_schema().columns.keys())

    if remove_unknown_cols:
        # drop excess columns (if any)
        excess_cols = list(set(current_cols) - set(schema_cols))

        if len(excess_cols) > 0:
            df = df.drop(columns=excess_cols)
            current_cols = list(df.columns)

    # `df` might contain extra columns, but we can still type it now
    df = TypedDataFrame[pandera_schema](df)

    if reorder_cols:
        # put columns in schema order, with extra columns in original order at the end
        all_cols = schema_cols.copy()
        all_cols.extend(current_cols)
        all_cols = list(dict.fromkeys(all_cols))
        df = TypedDataFrame[pandera_schema](df.loc[:, all_cols])

    return df


def read_sample_ids_from_file(sample_ids_list: Path) -> list[str]:
    """
    Read an ordered list of sample IDs from a text file.

    :param sample_ids_list: path to the file
    :return: a list of sample IDs
    """

    logging.info(f"Reading {sample_ids_list}")

    with open(sample_ids_list, "r") as f:
        sample_ids = f.read().splitlines()

    sample_ids = [x.strip() for x in sample_ids]

    return sample_ids


def read_tracking_file(
    tracking_in: Path, sample_ids: None | list[str] = None
) -> TypedDataFrame[Tracking]:
    """
    Read a tracking file as a typed data frame.

    :param tracking_in: path to the tracking file
    :param sample_ids: the names of the sample ID columns
    :return: a data frame of tracking data
    """

    logging.info(f"Reading {tracking_in}")

    if tracking_in.suffix == ".parquet":
        return type_data_frame(pd.read_parquet(tracking_in), Tracking)
    else:
        # we need the names of the sample ID columns since there's no header
        assert sample_ids is not None

        return type_data_frame(
            pd.read_table(
                tracking_in,
                sep="\t",
                header=None,
                names=["transcript_id", "loc", "gene_id", "val", *sample_ids],
                na_values="-",
                dtype="string",
            ),
            Tracking,
        )


def read_gtf_from_file(gtf_in: str | Path) -> TypedDataFrame[Gtf]:
    """
    Read a GTF file as a typed data frame.

    :param gtf_in: path to the GTF file
    :return: a data frame of GTF data
    """

    logging.info(f"Reading {gtf_in}")
    return type_data_frame(
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


def read_sqanti3_classification_file(
    squanti_classification: str | Path,
) -> TypedDataFrame[Sqanti3Classification]:
    """
    Read a Sqanti3 classication file as a typed data frame.

    :param squanti_classification: path to the tracking file
    :return: a data frame of classication data
    """

    logging.info(f"Loading Sqanti3 classification file from {squanti_classification}")
    return type_data_frame(
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
