from typing import TypeVar

import pandas as pd
import pandera as pa
import psutil
from nebelung.types import CoercedDataFrame
from pandera.typing import Series  # pyright: ignore
from pydantic import BaseModel


class Workspace(BaseModel):
    namespace: str
    name: str


class GcsSource(BaseModel):
    bucket: str
    glob: str


class GcsDest(BaseModel):
    bucket: str
    prefix: str


class DogspaConfig(BaseModel):
    workspace: Workspace
    gcp_project: str
    gcs_source: GcsSource
    gcs_destination: GcsDest
    uuid_namespace: str = "00000000-0000-0000-0000-000000000000"
    ncpus: int = psutil.cpu_count()
    dry_run: bool = True


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
    blacklist: Series[pd.BooleanDtype]
    expected_type: Series[pd.StringDtype] = pa.Field(nullable=True)
    bam_filepath: Series[pd.StringDtype] = pa.Field(nullable=True)
    bai_filepath: Series[pd.StringDtype] = pa.Field(nullable=True)


class ObjectMetadata(CoercedDataFrame):
    url: Series[pd.StringDtype] = pa.Field(unique=True)
    crc32c: Series[pd.StringDtype] = pa.Field(unique=True)
    size: Series[pd.Int64Dtype] = pa.Field(unique=True)
    gcs_obj_updated_at: Series[pd.StringDtype]


class IdentifiedSrcBam(CoercedDataFrame):
    bam_url: Series[pd.StringDtype] = pa.Field(unique=True)
    crc32c: Series[pd.StringDtype] = pa.Field(unique=True)
    size: Series[pd.Int64Dtype] = pa.Field(unique=True)
    gcs_obj_updated_at: Series[pd.StringDtype]
    model_id: Series[pd.StringDtype] = pa.Field(unique=True)
    issue: Series  # storing sets in this column, so it's a generic Pandas object dtype
    blacklist: Series[pd.BooleanDtype]


class SamplesMaybeInGumbo(IdentifiedSrcBam):
    already_in_gumbo: Series[pd.BooleanDtype]


class SamplesWithMetadata(SamplesMaybeInGumbo):
    model_condition_id: Series[pd.StringDtype]
    profile_id: Series[pd.StringDtype]


class SamplesWithShortReadMetadata(SamplesWithMetadata):
    main_sequencing_id: Series[pd.StringDtype] = pa.Field(nullable=True)
    sr_profile_id: Series[pd.StringDtype] = pa.Field(nullable=True)
    sr_bam_filepath: Series[pd.StringDtype] = pa.Field(nullable=True)
    sr_bai_filepath: Series[pd.StringDtype] = pa.Field(nullable=True)


class SamplesWithCDSIDs(SamplesWithShortReadMetadata):
    sequencing_id: Series[pd.StringDtype]


class CopiedSampleFiles(SamplesWithCDSIDs):
    new_url: Series[pd.StringDtype]
    copied: Series[pd.BooleanDtype]


class SamplesForGumbo(CoercedDataFrame):
    sequencing_id: Series[pd.StringDtype]
    bam_filepath: Series[pd.StringDtype]
    profile_id: Series[pd.StringDtype]
    legacy_size: Series[pd.Int64Dtype]
    update_time: Series[pd.StringDtype]
    sequencing_date: Series[pd.StringDtype]
    source: Series[pd.StringDtype]
    crc32c_hash: Series[pd.StringDtype]
    expected_type: Series[pd.StringDtype]
    issue: Series[pd.StringDtype] = pa.Field(nullable=True)
    blacklist: Series[pd.BooleanDtype]


class VersionedSamples(SamplesForGumbo):
    version: Series[int]


PydanticBaseModel = TypeVar("PydanticBaseModel", bound=BaseModel)
PanderaBaseSchema = TypeVar("PanderaBaseSchema", bound=CoercedDataFrame)
TypedDataFrame = pa.typing.DataFrame  # pyright: ignore
