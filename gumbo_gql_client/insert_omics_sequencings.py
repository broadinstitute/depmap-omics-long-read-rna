from typing import List, Optional

from .base_model import BaseModel


class InsertOmicsSequencings(BaseModel):
    set_username: List["InsertOmicsSequencingsSetUsername"]
    insert_omics_sequencing: Optional["InsertOmicsSequencingsInsertOmicsSequencing"]


class InsertOmicsSequencingsSetUsername(BaseModel):
    username: str


class InsertOmicsSequencingsInsertOmicsSequencing(BaseModel):
    affected_rows: int


InsertOmicsSequencings.model_rebuild()
