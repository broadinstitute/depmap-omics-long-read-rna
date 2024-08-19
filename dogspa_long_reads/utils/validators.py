from typing import Any, Literal, Optional, TypeVar

import pandas as pd
import pandera as pa
import psutil
from pandera.api.pandas.model_config import BaseConfig as PaBaseConfig
from pandera.typing import Series
from pydantic import BaseModel


class Workspace(BaseModel):
    namespace: str
    name: str


class GcsBucket(BaseModel):
    bucket: str
    prefix: str


class DogspaConfig(BaseModel):
    workspace: Workspace
    gcp_project: str
    gcs_source: GcsBucket
    gcs_destination: GcsBucket
    uuid_namespace: str = "00000000-0000-0000-0000-000000000000"
    ncpus: int = psutil.cpu_count()
    dry_run: bool = True


class BamUpdatedAt(BaseModel):
    url: str
    rg_updated_at: Optional[str] = None
    new_issue: Optional[str] = None


class CoercedDataFrame(pa.DataFrameModel):
    class Config(PaBaseConfig):
        coerce = True  # convert to indicated dtype upon TypedDataFrame init


class ModelsAndChildren(CoercedDataFrame):
    model_id: Series[pd.StringDtype]
    cell_line_name: Series[pd.StringDtype]
    stripped_cell_line_name: Series[pd.StringDtype]
    model_conditions: Series


class SeqTable(CoercedDataFrame):
    model_id: Series[pd.StringDtype]
    cell_line_name: Series[pd.StringDtype]
    stripped_cell_line_name: Series[pd.StringDtype]
    model_condition_id: Series[pd.StringDtype]
    profile_id: Series[pd.StringDtype]
    datatype: Series[pd.StringDtype]
    legacy_size: Series[pd.Int64Dtype] = pa.Field(nullable=True)
    sequencing_id: Series[pd.StringDtype] = pa.Field(nullable=True)
    is_main_sequencing_id: Series[pd.BooleanDtype]
    blacklist_omics: Series[pd.BooleanDtype]
    blacklist: Series[pd.BooleanDtype] = pa.Field(nullable=True)
    expected_type: Series[pd.StringDtype] = pa.Field(nullable=True)
    bam_filepath: Series[pd.StringDtype] = pa.Field(nullable=True)
    bai_filepath: Series[pd.StringDtype] = pa.Field(nullable=True)


class ObjectMetadata(CoercedDataFrame):
    url: Series[pd.StringDtype]
    crc32c: Series[pd.StringDtype]
    size: Series[pd.Int64Dtype]
    gcs_obj_updated_at: Series[pd.StringDtype]


class IdentifiedSrcBam(CoercedDataFrame):
    bam_url: Series[pd.StringDtype]
    crc32c: Series[pd.StringDtype]
    size: Series[pd.Int64Dtype]
    gcs_obj_updated_at: Series[pd.StringDtype]
    model_id: Series[pd.StringDtype]
    issue: Series[pd.StringDtype] = pa.Field(nullable=True)
    blacklist: Series[pd.BooleanDtype]


class SamplesMaybeInGumbo(IdentifiedSrcBam):
    already_in_gumbo: Series[pd.BooleanDtype]


class SamplesWithShortReadMetadata(IdentifiedSrcBam):
    main_sequencing_id: Series[pd.StringDtype]
    sr_profile_id: Series[pd.StringDtype]
    sr_bam_filepath: Series[pd.StringDtype]
    sr_bai_filepath: Series[pd.StringDtype]


class SamplesWithCDSIDs(IdentifiedSrcBam):
    sequencing_id: Series[pd.StringDtype]


class TerraSample(CoercedDataFrame):
    sample_id: Series[pd.StringDtype]


class SamplesWithProfileIds(ObjectMetadata):
    profile_id: Series[pd.StringDtype] = pa.Field(nullable=True)


class SamplesWithRgUpdatedAt(SamplesWithProfileIds):
    rg_updated_at: Series[pd.StringDtype] = pa.Field(nullable=True)


class CopiedSampleFiles(CoercedDataFrame):
    sequencing_id: Series[pd.StringDtype]
    url_kind: Series[pd.StringDtype]
    url: Series[pd.StringDtype]
    copied: Series[pd.BooleanDtype]


class SamplesForGumbo(CoercedDataFrame):
    sequencing_id: Series[pd.StringDtype]
    sm_id: Series[pd.StringDtype] = pa.Field(nullable=True)
    profile_id: Series[pd.StringDtype] = pa.Field(nullable=True)
    legacy_size: Series[pd.Int64Dtype] = pa.Field(nullable=True)
    update_time: Series[pd.StringDtype] = pa.Field(nullable=True)
    sequencing_date: Series[pd.StringDtype] = pa.Field(nullable=True)
    source: Series[pd.StringDtype]
    crc32c_hash: Series[pd.StringDtype] = pa.Field(nullable=True)
    expected_type: Series[pd.StringDtype]
    issue: Series[pd.StringDtype] = pa.Field(nullable=True)
    blacklist: Series[pd.BooleanDtype]


class VersionedSamples(SamplesForGumbo):
    version: Series[int]


PydanticBaseModel = TypeVar("PydanticBaseModel", bound=BaseModel)
PanderaBaseSchema = TypeVar("PanderaBaseSchema", bound=CoercedDataFrame)
TypedDataFrame = pa.typing.DataFrame
