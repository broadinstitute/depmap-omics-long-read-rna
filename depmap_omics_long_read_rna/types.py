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


class AlignmentMetadataLong(CoercedDataFrame):
    model_id: Series[pd.StringDtype]
    cell_line_name: Series[pd.StringDtype]
    stripped_cell_line_name: Series[pd.StringDtype]
    drug: Series[pd.StringDtype] = pa.Field(nullable=True)
    expansion_team: Series[pd.StringDtype] = pa.Field(nullable=True)
    model_condition_id: Series[pd.StringDtype]
    omics_profile_id: Series[pd.StringDtype]
    datatype: Series[pd.StringDtype]
    stranded: Series[pd.BooleanDtype] = pa.Field(nullable=True)
    omics_sequencing_id: Series[pd.StringDtype] = pa.Field(nullable=True)
    source: Series[pd.StringDtype] = pa.Field(nullable=True)
    version: Series[pd.Int64Dtype] = pa.Field(nullable=True)
    sequencing_alignment_id: Series[pd.Int64Dtype] = pa.Field(nullable=True)
    crai_bai_url: Series[pd.StringDtype] = pa.Field(nullable=True)
    cram_bam_url: Series[pd.StringDtype] = pa.Field(nullable=True)
    reference_genome: Series[pd.StringDtype] = pa.Field(nullable=True)
    sequencing_alignment_source: Series[pd.StringDtype] = pa.Field(nullable=True)


class LongReadAlignmentMetadata(CoercedDataFrame):
    sample_id: Series[pd.StringDtype]
    omics_profile_id: Series[pd.StringDtype]
    model_condition_id: Series[pd.StringDtype]
    model_id: Series[pd.StringDtype]
    cell_line_name: Series[pd.StringDtype]
    stripped_cell_line_name: Series[pd.StringDtype]
    delivery_sequencing_alignment_id: Series[pd.Int64Dtype]
    delivery_cram_bam: Series[pd.StringDtype]
    aligned_sequencing_alignment_id: Series[pd.Int64Dtype] = pa.Field(nullable=True)
    aligned_cram_bam: Series[pd.StringDtype] = pa.Field(nullable=True)
    aligned_crai_bai: Series[pd.StringDtype] = pa.Field(nullable=True)
    reference_genome: Series[pd.StringDtype] = pa.Field(nullable=True)


class ShortReadPriorities(CoercedDataFrame):
    model_id: Series[pd.StringDtype]
    model_condition_id: Series[pd.StringDtype]
    sr_omics_profile_id: Series[pd.StringDtype]
    sr_omics_sequencing_id: Series[pd.StringDtype]
    sr_sequencing_alignment_id: Series[pd.StringDtype]
    sr_crai_bai: Series[pd.StringDtype] = pa.Field(nullable=True)
    sr_cram_bam: Series[pd.StringDtype]
    version: Series[pd.Int64Dtype]
    drug_priority: Series[pd.Int64Dtype]
    stranded_priority: Series[pd.Int64Dtype]
    source_priority: Series[pd.Int64Dtype]
    expansion_team_priority: Series[pd.Int64Dtype]
    sequencing_alignment_source_priority: Series[pd.Int64Dtype]


class LongReadTerraSamples(CoercedDataFrame):
    sample_id: Series[pd.StringDtype] = pa.Field(unique=True)
    omics_profile_id: Series[pd.StringDtype] = pa.Field(unique=True)
    model_condition_id: Series[pd.StringDtype]
    model_id: Series[pd.StringDtype]
    cell_line_name: Series[pd.StringDtype]
    stripped_cell_line_name: Series[pd.StringDtype]
    delivery_sequencing_alignment_id: Series[pd.Int64Dtype] = pa.Field(unique=True)
    delivery_cram_bam: Series[pd.StringDtype] = pa.Field(unique=True)
    aligned_sequencing_alignment_id: Series[pd.Int64Dtype] = pa.Field(nullable=True)
    aligned_cram_bam: Series[pd.StringDtype] = pa.Field(nullable=True)
    aligned_crai_bai: Series[pd.StringDtype] = pa.Field(nullable=True)
    reference_genome: Series[pd.StringDtype] = pa.Field(nullable=True)
    sr_model_condition_id: Series[pd.StringDtype] = pa.Field(nullable=True)
    sr_omics_profile_id: Series[pd.StringDtype] = pa.Field(nullable=True)
    sr_omics_sequencing_id: Series[pd.StringDtype] = pa.Field(nullable=True)
    sr_sequencing_alignment_id: Series[pd.StringDtype] = pa.Field(nullable=True)
    sr_crai_bai: Series[pd.StringDtype] = pa.Field(nullable=True)
    sr_cram_bam: Series[pd.StringDtype] = pa.Field(nullable=True)
    sr_file_format: Series[pd.StringDtype] = pa.Field(
        isin={"CRAM", "BAM"}, nullable=True
    )

    @pa.dataframe_check
    def crai_required(cls, df: pd.DataFrame) -> pd.Series:
        return ~(df["sr_cram_bam"].str.endswith(".cram") & df["sr_crai_bai"].isna())


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
    omics_profile_id: Series[pd.StringDtype] = pa.Field(unique=True)


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


class SamplesWithMetadata(OnboardingSamples):
    model_condition_id: Series[pd.StringDtype]
    profile_id: Series[pd.StringDtype]


class CopiedSampleFiles(CoercedDataFrame):
    sample_id: Series[pd.StringDtype]
    url_kind: Series[pd.StringDtype]
    new_url: Series[pd.StringDtype]
    url: Series[pd.StringDtype]
    copied: Series[pd.BooleanDtype]


class ExistingAlignments(CoercedDataFrame):
    omics_sequencing_id: Series[pd.StringDtype]
    sequencing_alignment_source: Series[pd.StringDtype]
    size: Series[pd.Int64Dtype] = pa.Field(unique=True)


class AlignedSamples(CoercedDataFrame):
    sample_id: Series[pd.StringDtype] = pa.Field(unique=True)
    aligned_bam: Series[pd.StringDtype] = pa.Field(unique=True)
    aligned_bai: Series[pd.StringDtype] = pa.Field(unique=True)


class AlignedSamplesWithObjectMetadata(AlignedSamples):
    crc32c: Series[pd.StringDtype] = pa.Field(unique=True)
    size: Series[pd.Int64Dtype] = pa.Field(unique=True)


class NewSequencingAlignments(CoercedDataFrame):
    omics_sequencing_id: Series[pd.StringDtype] = pa.Field(unique=True)
    url: Series[pd.StringDtype] = pa.Field(unique=True)
    index_url: Series[pd.StringDtype] = pa.Field(unique=True)
    sequencing_alignment_source: Series[pd.StringDtype]
    reference_genome: Series[pd.StringDtype]
    crc32c_hash: Series[pd.StringDtype] = pa.Field(unique=True)
    size: Series[pd.Int64Dtype] = pa.Field(unique=True)


PydanticBaseModel = TypeVar("PydanticBaseModel", bound=BaseModel)
