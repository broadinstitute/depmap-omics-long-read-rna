from typing import List

from .base_model import BaseModel


class GetOnboardingWorkspace(BaseModel):
    records: List["GetOnboardingWorkspaceRecords"]


class GetOnboardingWorkspaceRecords(BaseModel):
    id: str
    excluded_terra_sample_ids: List[str]
    min_file_size: int


GetOnboardingWorkspace.model_rebuild()
