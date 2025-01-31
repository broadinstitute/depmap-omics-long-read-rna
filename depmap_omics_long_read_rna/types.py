from typing import Optional, TypeVar

import httpx
import pandas as pd
import pandera as pa
from nebelung.types import CoercedDataFrame
from pandera.typing import Series
from pydantic import BaseModel

from gumbo_gql_client.gumbo_client import GumboClient as AriadneGumboClient


class GumboClient(AriadneGumboClient):
    def __init__(
        self,
        url: str,
        username: str,
        headers: dict[str, str],
        http_client: Optional[httpx.Client] = None,
    ):
        super().__init__(url=url, headers=headers, http_client=http_client)
        self.username = username  # store username on this object for use in mutations


class ModelsAndChildren(CoercedDataFrame):
    model_id: Series[pd.StringDtype]
    cell_line_name: Series[pd.StringDtype]
    stripped_cell_line_name: Series[pd.StringDtype] = pa.Field(nullable=True)
    model_conditions: Series


class SampleMetadata(CoercedDataFrame):
    model_id: Series[pd.StringDtype]
    cell_line_name: Series[pd.StringDtype]
    stripped_cell_line_name: Series[pd.StringDtype]
    model_condition_id: Series[pd.StringDtype]
    omics_profile_id: Series[pd.StringDtype]
    datatype: Series[pd.StringDtype]
    omics_sequencing_id: Series[pd.StringDtype] = pa.Field(nullable=True)
    source: Series[pd.StringDtype] = pa.Field(nullable=True)
    version: Series[pd.Int64Dtype] = pa.Field(nullable=True)
    sequencing_alignment_id: Series[pd.Int64Dtype] = pa.Field(nullable=True)
    crai_bai_url: Series[pd.StringDtype] = pa.Field(nullable=True)
    cram_bam_url: Series[pd.StringDtype] = pa.Field(nullable=True)
    reference_genome: Series[pd.StringDtype] = pa.Field(nullable=True)
    sequencing_alignment_source: Series[pd.StringDtype] = pa.Field(nullable=True)


class ObjectMetadata(CoercedDataFrame):
    url: Series[pd.StringDtype] = pa.Field(unique=True)
    crc32c: Series[pd.StringDtype] = pa.Field(unique=True)
    size: Series[pd.Int64Dtype] = pa.Field(unique=True)
    gcs_obj_updated_at: Series[pd.StringDtype]


class DeliveryBams(CoercedDataFrame):
    delivery_bam: Series[pd.StringDtype] = pa.Field(unique=True)
    delivery_bam_crc32c: Series[pd.StringDtype] = pa.Field(unique=True)
    delivery_bam_size: Series[pd.Int64Dtype] = pa.Field(unique=True)
    delivery_bam_updated_at: Series[pd.StringDtype]
    aligned_bam: Optional[Series[pd.StringDtype]] = pa.Field(nullable=True)
    aligned_bai: Optional[Series[pd.StringDtype]] = pa.Field(nullable=True)
    model_id: Series[pd.StringDtype] = pa.Field(unique=True)


class SamplesWithCDSIDs(DeliveryBams):
    cds_id: Series[pd.StringDtype] = pa.Field(unique=True)


class OnboardingSamples(CoercedDataFrame):
    sequencing_id: Series[pd.StringDtype] = pa.Field(unique=True)
    delivery_bam: Series[pd.StringDtype] = pa.Field(unique=True)
    delivery_bam_crc32c: Series[pd.StringDtype] = pa.Field(unique=True)
    delivery_bam_size: Series[pd.Int64Dtype] = pa.Field(unique=True)
    delivery_bam_updated_at: Series[pd.StringDtype]
    aligned_bam: Series[pd.StringDtype] = pa.Field(unique=True)
    aligned_bai: Series[pd.StringDtype] = pa.Field(unique=True)
    model_id: Series[pd.StringDtype] = pa.Field(unique=True)
    issue: Series  # storing sets in this column, so it's a generic Pandas object dtype


class SamplesMaybeInGumbo(OnboardingSamples):
    already_in_gumbo: Series[pd.BooleanDtype]


class SamplesWithMetadata(OnboardingSamples):
    model_condition_id: Series[pd.StringDtype]
    profile_id: Series[pd.StringDtype]


class CopiedSampleFiles(CoercedDataFrame):
    sequencing_id: Series[pd.StringDtype]
    url_kind: Series[pd.StringDtype]
    new_url: Series[pd.StringDtype]
    url: Series[pd.StringDtype]
    copied: Series[pd.BooleanDtype]


class SamplesForGumbo(CoercedDataFrame):
    sequencing_id: Series[pd.StringDtype]
    unaligned_bam_url: Series[pd.StringDtype]
    unaligned_bam_crc32c_hash: Series[pd.StringDtype]
    unaligned_bam_size: Series[pd.StringDtype]
    bam_url: Series[pd.StringDtype]
    bam_crc32c_hash: Series[pd.StringDtype]
    bam_size: Series[pd.StringDtype]
    bai_url: Series[pd.StringDtype]
    profile_id: Series[pd.StringDtype]
    update_time: Series[pd.StringDtype]
    sequencing_date: Series[pd.StringDtype]
    source: Series[pd.StringDtype]
    expected_type: Series[pd.StringDtype]
    issue: Series[pd.StringDtype] = pa.Field(nullable=True)


class VersionedSamples(SamplesForGumbo):
    version: Series[int]


class ShortReadMetadata(CoercedDataFrame):
    sr_omics_sequencing_id: Series[pd.StringDtype]
    sr_omics_profile_id: Series[pd.StringDtype]
    model_condition_id: Series[pd.StringDtype]


PydanticBaseModel = TypeVar("PydanticBaseModel", bound=BaseModel)
