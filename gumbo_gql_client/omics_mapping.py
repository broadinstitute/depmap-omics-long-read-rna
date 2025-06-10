from typing import List, Optional

from .base_model import BaseModel


class OmicsMapping(BaseModel):
    records: List["OmicsMappingRecords"]


class OmicsMappingRecords(BaseModel):
    model_id: Optional[str]
    model_condition_id: Optional[str]
    omics_profile_id: Optional[str]
    omics_sequencing_id: Optional[str]
    priority: Optional[int]
    model: Optional["OmicsMappingRecordsModel"]
    omics_sequencing: Optional["OmicsMappingRecordsOmicsSequencing"]


class OmicsMappingRecordsModel(BaseModel):
    cell_line_name: Optional[str]
    stripped_cell_line_name: Optional[str]


class OmicsMappingRecordsOmicsSequencing(BaseModel):
    sequencing_alignments: List[
        "OmicsMappingRecordsOmicsSequencingSequencingAlignments"
    ]


class OmicsMappingRecordsOmicsSequencingSequencingAlignments(BaseModel):
    sequencing_alignment_id: int
    crai_bai_url: Optional[str]
    cram_bam_url: str
    reference_genome: Optional[str]
    sequencing_alignment_source: str
    size: int


OmicsMapping.model_rebuild()
OmicsMappingRecords.model_rebuild()
OmicsMappingRecordsOmicsSequencing.model_rebuild()
