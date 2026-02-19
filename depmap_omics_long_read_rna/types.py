from typing import Literal, Optional

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
    model_condition_id: Series[pd.StringDtype]
    omics_profile_id: Series[pd.StringDtype]
    omics_sequencing_id: Series[pd.StringDtype]
    datatype: Series[pd.StringDtype]
    priority: Series[pd.Int64Dtype]
    model: Series
    omics_sequencing: Series


class AlignmentMetadataLong(CoercedDataFrame):
    model_id: Series[pd.StringDtype]
    model_condition_id: Series[pd.StringDtype]
    omics_profile_id: Series[pd.StringDtype]
    omics_sequencing_id: Series[pd.StringDtype]
    datatype: Series[pd.StringDtype]
    priority: Series[pd.Int64Dtype]
    cell_line_name: Series[pd.StringDtype]
    stripped_cell_line_name: Series[pd.StringDtype]
    sequencing_alignment_id: Series[pd.Int64Dtype]
    crai_bai_url: Series[pd.StringDtype] = pa.Field(nullable=True)
    cram_bam_url: Series[pd.StringDtype]
    reference_genome: Series[pd.StringDtype] = pa.Field(nullable=True)
    sequencing_alignment_source: Series[pd.StringDtype]
    size: Series[pd.Int64Dtype]


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


class ShortReadTerraSamples(CoercedDataFrame):
    sample_id: Series[pd.StringDtype]
    sr_star_junctions: Series[pd.StringDtype] = pa.Field(nullable=True)


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
    sr_reference_genome: Series[pd.StringDtype] = pa.Field(nullable=True)
    sr_ref_fasta: Series[pd.StringDtype] = pa.Field(nullable=True)
    sr_ref_fasta_index: Series[pd.StringDtype] = pa.Field(nullable=True)
    sr_model_condition_id: Series[pd.StringDtype] = pa.Field(nullable=True)
    sr_omics_profile_id: Series[pd.StringDtype] = pa.Field(nullable=True)
    sr_omics_sequencing_id: Series[pd.StringDtype] = pa.Field(nullable=True)
    sr_sequencing_alignment_id: Series[pd.StringDtype] = pa.Field(nullable=True)
    sr_crai_bai: Series[pd.StringDtype] = pa.Field(nullable=True)
    sr_cram_bam: Series[pd.StringDtype] = pa.Field(nullable=True)
    sr_file_format: Series[pd.StringDtype] = pa.Field(
        isin={"CRAM", "BAM"}, nullable=True
    )
    sr_star_junctions: Series[pd.StringDtype] = pa.Field(nullable=True)

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
    omics_profile_id: Series[pd.StringDtype] = pa.Field(unique=True)


class SamplesWithCDSIDs(DeliveryBams):
    cds_id: Series[pd.StringDtype] = pa.Field(unique=True)


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
    aligned_cram_bam: Series[pd.StringDtype] = pa.Field(unique=True)
    aligned_crai_bai: Series[pd.StringDtype] = pa.Field(unique=True)


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


class DeltaJob(BaseModel):
    workflow_name: str
    entity_type: str
    entity_set_type: str
    entity_id_col: str
    expression: str
    input_cols: set[str] | None = None
    output_cols: set[str] | None = None
    resubmit_n_times: int = 0
    force_retry: bool = False
    use_callcache: bool = True
    use_reference_disks: bool = False
    delete_intermediate_output_files: bool = False
    memory_retry_multiplier: float = 1.0
    per_workflow_cost_cap: float | None = None
    workflow_failure_mode: Literal["NoNewCalls", "ContinueWhilePossible"] = "NoNewCalls"
    user_comment: str | None = None
    max_n_entities: int | None = None
    dry_run: bool = False
