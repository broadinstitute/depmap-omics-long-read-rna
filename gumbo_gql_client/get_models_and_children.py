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
    blacklist: bool
    expected_type: Optional[str]
    sequencing_id: str
    source: Optional[str]
    version: Optional[int]
    sequencing_alignments: List[
        "GetModelsAndChildrenRecordsModelConditionsOmicsProfilesOmicsSequencingsSequencingAlignments"
    ]


class GetModelsAndChildrenRecordsModelConditionsOmicsProfilesOmicsSequencingsSequencingAlignments(
    BaseModel
):
    id: int
    url: str
    index_url: Optional[str]
    size: int
    reference_genome: Optional[str]
    sequencing_alignment_source: str


GetModelsAndChildren.model_rebuild()
GetModelsAndChildrenRecords.model_rebuild()
GetModelsAndChildrenRecordsModelConditions.model_rebuild()
GetModelsAndChildrenRecordsModelConditionsOmicsProfiles.model_rebuild()
GetModelsAndChildrenRecordsModelConditionsOmicsProfilesOmicsSequencings.model_rebuild()
