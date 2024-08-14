from typing import List, Optional

from .base_model import BaseModel


class GetOmicsProfiles(BaseModel):
    records: List["GetOmicsProfilesRecords"]


class GetOmicsProfilesRecords(BaseModel):
    profile_id: str
    datatype: Optional[str]
    smid_ordered: Optional[str]
    sm_id_matched: Optional[str]
    smid_returned: Optional[str]


GetOmicsProfiles.model_rebuild()
