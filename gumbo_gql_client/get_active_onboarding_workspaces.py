from typing import List

from .base_model import BaseModel


class GetActiveOnboardingWorkspaces(BaseModel):
    records: List["GetActiveOnboardingWorkspacesRecords"]


class GetActiveOnboardingWorkspacesRecords(BaseModel):
    id: str
    excluded_terra_sample_ids: List[str]


GetActiveOnboardingWorkspaces.model_rebuild()
