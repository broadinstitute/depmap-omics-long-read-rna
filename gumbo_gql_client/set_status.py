from typing import List, Optional

from .base_model import BaseModel


class SetStatus(BaseModel):
    set_username: List["SetStatusSetUsername"]
    update_omics_profile_many: Optional[
        List[Optional["SetStatusUpdateOmicsProfileMany"]]
    ]


class SetStatusSetUsername(BaseModel):
    username: str


class SetStatusUpdateOmicsProfileMany(BaseModel):
    affected_rows: int


SetStatus.model_rebuild()
