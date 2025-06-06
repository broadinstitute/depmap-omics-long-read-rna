from typing import List, Optional

from .base_model import BaseModel


class ShortReadOmicsMapping(BaseModel):
    records: List["ShortReadOmicsMappingRecords"]


class ShortReadOmicsMappingRecords(BaseModel):
    model_id: Optional[str]
    model_condition_id: Optional[str]
    omics_sequencing_id: Optional[str]
    priority: Optional[int]


ShortReadOmicsMapping.model_rebuild()
