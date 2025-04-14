from typing import List, Optional

from .base_model import BaseModel


class InsertSequencingAlignments(BaseModel):
    set_username: List["InsertSequencingAlignmentsSetUsername"]
    insert_sequencing_alignment: Optional[
        "InsertSequencingAlignmentsInsertSequencingAlignment"
    ]


class InsertSequencingAlignmentsSetUsername(BaseModel):
    username: str


class InsertSequencingAlignmentsInsertSequencingAlignment(BaseModel):
    affected_rows: int


InsertSequencingAlignments.model_rebuild()
