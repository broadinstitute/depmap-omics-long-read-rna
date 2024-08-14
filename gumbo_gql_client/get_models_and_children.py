from typing import List, Optional

from .base_model import BaseModel


class GetModelsAndChildren(BaseModel):
    records: List["GetModelsAndChildrenRecords"]


class GetModelsAndChildrenRecords(BaseModel):
    model_id: str
    cell_line_name: Optional[str]
    stripped_cell_line_name: Optional[str]
    model_conditions: List["GetModelsAndChildrenRecordsModelConditions"]


class GetModelsAndChildrenRecordsModelConditions(BaseModel):
    model_condition_id: str
    omics_profiles: List["GetModelsAndChildrenRecordsModelConditionsOmicsProfiles"]


class GetModelsAndChildrenRecordsModelConditionsOmicsProfiles(BaseModel):
    profile_id: str
    main_sequencing_id: Optional[str]
    datatype: Optional[str]
    omics_sequencings: List[
        "GetModelsAndChildrenRecordsModelConditionsOmicsProfilesOmicsSequencings"
    ]


class GetModelsAndChildrenRecordsModelConditionsOmicsProfilesOmicsSequencings(
    BaseModel
):
    sequencing_id: str
    legacy_size: Optional[int]
    sequencing_id: str
    expected_type: Optional[str]
    bam_filepath: Optional[str]
    bai_filepath: Optional[str]


GetModelsAndChildren.model_rebuild()
GetModelsAndChildrenRecords.model_rebuild()
GetModelsAndChildrenRecordsModelConditions.model_rebuild()
GetModelsAndChildrenRecordsModelConditionsOmicsProfiles.model_rebuild()
