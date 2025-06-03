from typing import List, Optional

from .base_model import BaseModel


class InsertSequencingAlignment(BaseModel):
    set_username: List["InsertSequencingAlignmentSetUsername"]
    insert_sequencing_alignment: Optional[
        "InsertSequencingAlignmentInsertSequencingAlignment"
    ]


class InsertSequencingAlignmentSetUsername(BaseModel):
    username: str


class InsertSequencingAlignmentInsertSequencingAlignment(BaseModel):
    affected_rows: int


InsertSequencingAlignment.model_rebuild()
