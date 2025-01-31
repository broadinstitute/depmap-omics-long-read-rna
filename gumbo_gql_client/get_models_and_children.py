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
    omics_profile_id: str
    datatype: Optional[str]
    omics_sequencings: List[
        "GetModelsAndChildrenRecordsModelConditionsOmicsProfilesOmicsSequencings"
    ]


class GetModelsAndChildrenRecordsModelConditionsOmicsProfilesOmicsSequencings(
    BaseModel
):
    omics_sequencing_id: str
    blacklist: bool
    expected_type: Optional[str]
    source: Optional[str]
    version: Optional[int]
    sequencing_alignments: List[
        "GetModelsAndChildrenRecordsModelConditionsOmicsProfilesOmicsSequencingsSequencingAlignments"
    ]


class GetModelsAndChildrenRecordsModelConditionsOmicsProfilesOmicsSequencingsSequencingAlignments(
    BaseModel
):
    sequencing_alignment_id: int
    crai_bai_url: Optional[str]
    cram_bam_url: str
    reference_genome: Optional[str]
    sequencing_alignment_source: str


GetModelsAndChildren.model_rebuild()
GetModelsAndChildrenRecords.model_rebuild()
GetModelsAndChildrenRecordsModelConditions.model_rebuild()
GetModelsAndChildrenRecordsModelConditionsOmicsProfiles.model_rebuild()
GetModelsAndChildrenRecordsModelConditionsOmicsProfilesOmicsSequencings.model_rebuild()
