from typing import TypeVar

import pandas as pd
import pandera.pandas as pa
from pandera.api.dataframe.model_config import BaseConfig
from pandera.typing import DataFrame, Series


class CoercedDataFrame(pa.DataFrameModel):
    class Config(BaseConfig):  # pyright: ignore
        coerce = True


class NonemptyDataFrame(CoercedDataFrame):
    @classmethod
    @pa.dataframe_check
    def not_empty(cls, df: pd.DataFrame) -> bool:
        return len(df) > 0


class TrackingToProcess(CoercedDataFrame):
    transcript_id: Series[pd.StringDtype]
    loc: Series[pd.StringDtype]
    gene_id: Series[pd.StringDtype] = pa.Field(nullable=True)
    val: Series[pd.StringDtype]
    # other columns are named for sample IDs


class TrackingLong(CoercedDataFrame):
    transcript_id: Series[pd.StringDtype]
    loc: Series[pd.StringDtype]
    gene_id: Series[pd.StringDtype] = pa.Field(nullable=True)
    val: Series[pd.StringDtype]
    sample: Series[pd.StringDtype]
    id: Series[pd.StringDtype]
    id1: Series[pd.StringDtype]


class TranscriptCounts(CoercedDataFrame):
    id1: Series[pd.StringDtype]
    count: Series[pd.Int64Dtype]


class TrackingProcessed(TrackingLong):
    count: Series[pd.Int64Dtype]


class UpdatedTracking(CoercedDataFrame):
    transcript_id: Series[pd.StringDtype]
    loc: Series[pd.StringDtype]
    gene_id: Series[pd.StringDtype] = pa.Field(nullable=True)
    val: Series[pd.StringDtype]
    sample: Series[pd.StringDtype]
    id: Series[pd.StringDtype]
    id1: Series[pd.StringDtype]
    count: Series[pd.Int64Dtype]


class Sqanti3Classification(CoercedDataFrame):
    isoform: Series[pd.StringDtype]
    structural_category: Series[pd.StringDtype]
    associated_transcript: Series[pd.StringDtype]
    RTS_stage: Series[pd.BooleanDtype]
    coding: Series[pd.StringDtype]
    ORF_seq: Series[pd.StringDtype] = pa.Field(nullable=True)


class Gtf(NonemptyDataFrame):
    seqname: Series[pd.StringDtype]
    source: Series[pd.StringDtype]
    feature: Series[pd.StringDtype] = pa.Field(isin=["exon", "transcript"])
    start: Series[pd.StringDtype]
    end: Series[pd.StringDtype]
    score: Series[pd.StringDtype]
    strand: Series[pd.StringDtype] = pa.Field(isin=["+", "-", "."])
    frame: Series[pd.StringDtype]
    attribute: Series[pd.StringDtype]


PanderaBaseSchema = TypeVar("PanderaBaseSchema", bound=CoercedDataFrame)
TypedDataFrame = DataFrame
