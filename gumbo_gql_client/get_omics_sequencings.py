from typing import List, Optional

from .base_model import BaseModel


class GetOmicsSequencings(BaseModel):
    records: List["GetOmicsSequencingsRecords"]


class GetOmicsSequencingsRecords(BaseModel):
    legacy_size: Optional[int]
    expected_type: Optional[str]
    profile_id: Optional[str]


GetOmicsSequencings.model_rebuild()
