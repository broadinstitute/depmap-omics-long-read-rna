from typing import List, Optional

from .base_model import BaseModel


class InsertOnboardingJob(BaseModel):
    set_username: List["InsertOnboardingJobSetUsername"]
    insert_onboarding_job_one: Optional["InsertOnboardingJobInsertOnboardingJobOne"]


class InsertOnboardingJobSetUsername(BaseModel):
    username: str


class InsertOnboardingJobInsertOnboardingJobOne(BaseModel):
    id: int


InsertOnboardingJob.model_rebuild()
