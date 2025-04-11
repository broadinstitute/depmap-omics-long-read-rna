from typing import List

from .base_model import BaseModel


class LongReadSequencingAlignments(BaseModel):
    omics_sequencing: List["LongReadSequencingAlignmentsOmicsSequencing"]


class LongReadSequencingAlignmentsOmicsSequencing(BaseModel):
    sequencing_id: str
    sequencing_alignments: List[
        "LongReadSequencingAlignmentsOmicsSequencingSequencingAlignments"
    ]


class LongReadSequencingAlignmentsOmicsSequencingSequencingAlignments(BaseModel):
    sequencing_alignment_source: str


LongReadSequencingAlignments.model_rebuild()
LongReadSequencingAlignmentsOmicsSequencing.model_rebuild()
