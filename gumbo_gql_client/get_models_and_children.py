from typing import Any, List, Optional

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
    blacklist_omics: Optional[bool]
    omics_order_date: Optional[Any]
    smid_ordered: Optional[str]
    smid_returned: Optional[str]
    omics_sequencings: List[
        "GetModelsAndChildrenRecordsModelConditionsOmicsProfilesOmicsSequencings"
    ]


class GetModelsAndChildrenRecordsModelConditionsOmicsProfilesOmicsSequencings(
    BaseModel
):
    bai_filepath: Optional[str]
    bam_filepath: Optional[str]
    blacklist: bool
    expected_type: Optional[str]
    sequencing_id: str
    source: Optional[str]
    unaligned_bam_size: Optional[int]


GetModelsAndChildren.model_rebuild()
GetModelsAndChildrenRecords.model_rebuild()
GetModelsAndChildrenRecordsModelConditions.model_rebuild()
GetModelsAndChildrenRecordsModelConditionsOmicsProfiles.model_rebuild()
