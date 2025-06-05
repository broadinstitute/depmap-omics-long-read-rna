from typing import List, Optional

from .base_model import BaseModel


class MappedShortLongRnaSamples(BaseModel):
    records: List["MappedShortLongRnaSamplesRecords"]


class MappedShortLongRnaSamplesRecords(BaseModel):
    model_id: Optional[str]
    model_condition_id: Optional[str]
    omics_profile_id: Optional[str]
    omics_sequencing_id: Optional[str]
    datatype: Optional[str]
    priority: Optional[int]
    model: Optional["MappedShortLongRnaSamplesRecordsModel"]
    omics_sequencing: Optional["MappedShortLongRnaSamplesRecordsOmicsSequencing"]


class MappedShortLongRnaSamplesRecordsModel(BaseModel):
    cell_line_name: Optional[str]
    stripped_cell_line_name: Optional[str]


class MappedShortLongRnaSamplesRecordsOmicsSequencing(BaseModel):
    sequencing_alignments: List[
        "MappedShortLongRnaSamplesRecordsOmicsSequencingSequencingAlignments"
    ]


class MappedShortLongRnaSamplesRecordsOmicsSequencingSequencingAlignments(BaseModel):
    sequencing_alignment_id: int
    crai_bai_url: Optional[str]
    cram_bam_url: str
    reference_genome: Optional[str]
    sequencing_alignment_source: str
    size: int


MappedShortLongRnaSamples.model_rebuild()
MappedShortLongRnaSamplesRecords.model_rebuild()
MappedShortLongRnaSamplesRecordsOmicsSequencing.model_rebuild()
