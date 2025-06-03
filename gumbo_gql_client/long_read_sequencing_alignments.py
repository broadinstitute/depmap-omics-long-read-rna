from typing import List

from .base_model import BaseModel


class LongReadSequencingAlignments(BaseModel):
    records: List["LongReadSequencingAlignmentsRecords"]


class LongReadSequencingAlignmentsRecords(BaseModel):
    omics_sequencing_id: str
    sequencing_alignments: List[
        "LongReadSequencingAlignmentsRecordsSequencingAlignments"
    ]


class LongReadSequencingAlignmentsRecordsSequencingAlignments(BaseModel):
    sequencing_alignment_source: str
    size: int


LongReadSequencingAlignments.model_rebuild()
LongReadSequencingAlignmentsRecords.model_rebuild()
