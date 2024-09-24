from typing import Any, List, Optional

from pydantic import Field

from .base_model import BaseModel
from .enums import (
    audit_user_constraint,
    audit_user_update_column,
    cursor_ordering,
    depmap_model_type_constraint,
    depmap_model_type_update_column,
    media_constraint,
    media_update_column,
    model_condition_constraint,
    model_condition_select_column,
    model_condition_select_column_model_condition_aggregate_bool_exp_bool_and_arguments_columns,
    model_condition_select_column_model_condition_aggregate_bool_exp_bool_or_arguments_columns,
    model_condition_update_column,
    model_constraint,
    model_update_column,
    omics_profile_constraint,
    omics_profile_select_column,
    omics_profile_select_column_omics_profile_aggregate_bool_exp_bool_and_arguments_columns,
    omics_profile_select_column_omics_profile_aggregate_bool_exp_bool_or_arguments_columns,
    omics_profile_update_column,
    omics_sequencing_constraint,
    omics_sequencing_select_column,
    omics_sequencing_select_column_omics_sequencing_aggregate_bool_exp_bool_and_arguments_columns,
    omics_sequencing_select_column_omics_sequencing_aggregate_bool_exp_bool_or_arguments_columns,
    omics_sequencing_update_column,
    order_by,
    patient_constraint,
    patient_update_column,
    snp_fingerprint_comparison_constraint,
    snp_fingerprint_comparison_update_column,
    snp_fingerprint_constraint,
    snp_fingerprint_update_column,
    str_profile_constraint,
    str_profile_update_column,
    task_entity_constraint,
    task_entity_update_column,
    task_result_constraint,
    task_result_select_column,
    task_result_update_column,
    terra_sync_constraint,
    terra_sync_update_column,
)


class Boolean_comparison_exp(BaseModel):
    eq: Optional[bool] = Field(alias="_eq", default=None)
    gt: Optional[bool] = Field(alias="_gt", default=None)
    gte: Optional[bool] = Field(alias="_gte", default=None)
    in_: Optional[List[bool]] = Field(alias="_in", default=None)
    is_null: Optional[bool] = Field(alias="_is_null", default=None)
    lt: Optional[bool] = Field(alias="_lt", default=None)
    lte: Optional[bool] = Field(alias="_lte", default=None)
    neq: Optional[bool] = Field(alias="_neq", default=None)
    nin: Optional[List[bool]] = Field(alias="_nin", default=None)


class Int_comparison_exp(BaseModel):
    eq: Optional[int] = Field(alias="_eq", default=None)
    gt: Optional[int] = Field(alias="_gt", default=None)
    gte: Optional[int] = Field(alias="_gte", default=None)
    in_: Optional[List[int]] = Field(alias="_in", default=None)
    is_null: Optional[bool] = Field(alias="_is_null", default=None)
    lt: Optional[int] = Field(alias="_lt", default=None)
    lte: Optional[int] = Field(alias="_lte", default=None)
    neq: Optional[int] = Field(alias="_neq", default=None)
    nin: Optional[List[int]] = Field(alias="_nin", default=None)


class String_comparison_exp(BaseModel):
    eq: Optional[str] = Field(alias="_eq", default=None)
    gt: Optional[str] = Field(alias="_gt", default=None)
    gte: Optional[str] = Field(alias="_gte", default=None)
    ilike: Optional[str] = Field(alias="_ilike", default=None)
    in_: Optional[List[str]] = Field(alias="_in", default=None)
    iregex: Optional[str] = Field(alias="_iregex", default=None)
    is_null: Optional[bool] = Field(alias="_is_null", default=None)
    like: Optional[str] = Field(alias="_like", default=None)
    lt: Optional[str] = Field(alias="_lt", default=None)
    lte: Optional[str] = Field(alias="_lte", default=None)
    neq: Optional[str] = Field(alias="_neq", default=None)
    nilike: Optional[str] = Field(alias="_nilike", default=None)
    nin: Optional[List[str]] = Field(alias="_nin", default=None)
    niregex: Optional[str] = Field(alias="_niregex", default=None)
    nlike: Optional[str] = Field(alias="_nlike", default=None)
    nregex: Optional[str] = Field(alias="_nregex", default=None)
    nsimilar: Optional[str] = Field(alias="_nsimilar", default=None)
    regex: Optional[str] = Field(alias="_regex", default=None)
    similar: Optional[str] = Field(alias="_similar", default=None)


class audit_user_bool_exp(BaseModel):
    and_: Optional[List["audit_user_bool_exp"]] = Field(alias="_and", default=None)
    not_: Optional["audit_user_bool_exp"] = Field(alias="_not", default=None)
    or_: Optional[List["audit_user_bool_exp"]] = Field(alias="_or", default=None)
    username: Optional["String_comparison_exp"] = None


class audit_user_insert_input(BaseModel):
    username: Optional[str] = None


class audit_user_on_conflict(BaseModel):
    constraint: audit_user_constraint
    update_columns: List[audit_user_update_column]
    where: Optional["audit_user_bool_exp"] = None


class audit_user_order_by(BaseModel):
    username: Optional[order_by] = None


class audit_user_pk_columns_input(BaseModel):
    username: str


class audit_user_set_input(BaseModel):
    username: Optional[str] = None


class audit_user_stream_cursor_input(BaseModel):
    initial_value: "audit_user_stream_cursor_value_input"
    ordering: Optional[cursor_ordering] = None


class audit_user_stream_cursor_value_input(BaseModel):
    username: Optional[str] = None


class audit_user_updates(BaseModel):
    set: Optional["audit_user_set_input"] = Field(alias="_set", default=None)
    where: "audit_user_bool_exp"


class bigint_comparison_exp(BaseModel):
    eq: Optional[int] = Field(alias="_eq", default=None)
    gt: Optional[int] = Field(alias="_gt", default=None)
    gte: Optional[int] = Field(alias="_gte", default=None)
    in_: Optional[List[int]] = Field(alias="_in", default=None)
    is_null: Optional[bool] = Field(alias="_is_null", default=None)
    lt: Optional[int] = Field(alias="_lt", default=None)
    lte: Optional[int] = Field(alias="_lte", default=None)
    neq: Optional[int] = Field(alias="_neq", default=None)
    nin: Optional[List[int]] = Field(alias="_nin", default=None)


class date_comparison_exp(BaseModel):
    eq: Optional[Any] = Field(alias="_eq", default=None)
    gt: Optional[Any] = Field(alias="_gt", default=None)
    gte: Optional[Any] = Field(alias="_gte", default=None)
    in_: Optional[List[Any]] = Field(alias="_in", default=None)
    is_null: Optional[bool] = Field(alias="_is_null", default=None)
    lt: Optional[Any] = Field(alias="_lt", default=None)
    lte: Optional[Any] = Field(alias="_lte", default=None)
    neq: Optional[Any] = Field(alias="_neq", default=None)
    nin: Optional[List[Any]] = Field(alias="_nin", default=None)


class depmap_model_type_bool_exp(BaseModel):
    and_: Optional[List["depmap_model_type_bool_exp"]] = Field(
        alias="_and", default=None
    )
    not_: Optional["depmap_model_type_bool_exp"] = Field(alias="_not", default=None)
    or_: Optional[List["depmap_model_type_bool_exp"]] = Field(alias="_or", default=None)
    depmap_code: Optional["String_comparison_exp"] = None
    lineage: Optional["String_comparison_exp"] = None
    oncotree_code: Optional["String_comparison_exp"] = None
    primary_disease: Optional["String_comparison_exp"] = None
    subtype: Optional["String_comparison_exp"] = None


class depmap_model_type_insert_input(BaseModel):
    depmap_code: Optional[str] = None
    lineage: Optional[str] = None
    oncotree_code: Optional[str] = None
    primary_disease: Optional[str] = None
    subtype: Optional[str] = None


class depmap_model_type_obj_rel_insert_input(BaseModel):
    data: "depmap_model_type_insert_input"
    on_conflict: Optional["depmap_model_type_on_conflict"] = None


class depmap_model_type_on_conflict(BaseModel):
    constraint: depmap_model_type_constraint
    update_columns: List[depmap_model_type_update_column]
    where: Optional["depmap_model_type_bool_exp"] = None


class depmap_model_type_order_by(BaseModel):
    depmap_code: Optional[order_by] = None
    lineage: Optional[order_by] = None
    oncotree_code: Optional[order_by] = None
    primary_disease: Optional[order_by] = None
    subtype: Optional[order_by] = None


class depmap_model_type_pk_columns_input(BaseModel):
    depmap_code: str


class depmap_model_type_set_input(BaseModel):
    depmap_code: Optional[str] = None
    lineage: Optional[str] = None
    oncotree_code: Optional[str] = None
    primary_disease: Optional[str] = None
    subtype: Optional[str] = None


class depmap_model_type_stream_cursor_input(BaseModel):
    initial_value: "depmap_model_type_stream_cursor_value_input"
    ordering: Optional[cursor_ordering] = None


class depmap_model_type_stream_cursor_value_input(BaseModel):
    depmap_code: Optional[str] = None
    lineage: Optional[str] = None
    oncotree_code: Optional[str] = None
    primary_disease: Optional[str] = None
    subtype: Optional[str] = None


class depmap_model_type_updates(BaseModel):
    set: Optional["depmap_model_type_set_input"] = Field(alias="_set", default=None)
    where: "depmap_model_type_bool_exp"


class fingerprintable_samples_bool_exp(BaseModel):
    and_: Optional[List["fingerprintable_samples_bool_exp"]] = Field(
        alias="_and", default=None
    )
    not_: Optional["fingerprintable_samples_bool_exp"] = Field(
        alias="_not", default=None
    )
    or_: Optional[List["fingerprintable_samples_bool_exp"]] = Field(
        alias="_or", default=None
    )
    bai_filepath: Optional["String_comparison_exp"] = None
    bam_filepath: Optional["String_comparison_exp"] = None
    cell_line_name: Optional["String_comparison_exp"] = None
    expected_type: Optional["String_comparison_exp"] = None
    hg_19_bai_filepath: Optional["String_comparison_exp"] = Field(
        alias="hg19_bai_filepath", default=None
    )
    hg_19_bam_filepath: Optional["String_comparison_exp"] = Field(
        alias="hg19_bam_filepath", default=None
    )
    hg_38_crai_filepath: Optional["String_comparison_exp"] = Field(
        alias="hg38_crai_filepath", default=None
    )
    hg_38_cram_filepath: Optional["String_comparison_exp"] = Field(
        alias="hg38_cram_filepath", default=None
    )
    id: Optional["bigint_comparison_exp"] = None
    low_quality: Optional["Boolean_comparison_exp"] = None
    model_condition: Optional["String_comparison_exp"] = None
    model_id: Optional["String_comparison_exp"] = None
    omics_sequencing_blacklist: Optional["Boolean_comparison_exp"] = None
    patient_id: Optional["String_comparison_exp"] = None
    profile_blacklist_omics: Optional["Boolean_comparison_exp"] = None
    profile_id: Optional["String_comparison_exp"] = None
    sequencing_id: Optional["String_comparison_exp"] = None


class fingerprintable_samples_order_by(BaseModel):
    bai_filepath: Optional[order_by] = None
    bam_filepath: Optional[order_by] = None
    cell_line_name: Optional[order_by] = None
    expected_type: Optional[order_by] = None
    hg_19_bai_filepath: Optional[order_by] = Field(
        alias="hg19_bai_filepath", default=None
    )
    hg_19_bam_filepath: Optional[order_by] = Field(
        alias="hg19_bam_filepath", default=None
    )
    hg_38_crai_filepath: Optional[order_by] = Field(
        alias="hg38_crai_filepath", default=None
    )
    hg_38_cram_filepath: Optional[order_by] = Field(
        alias="hg38_cram_filepath", default=None
    )
    id: Optional[order_by] = None
    low_quality: Optional[order_by] = None
    model_condition: Optional[order_by] = None
    model_id: Optional[order_by] = None
    omics_sequencing_blacklist: Optional[order_by] = None
    patient_id: Optional[order_by] = None
    profile_blacklist_omics: Optional[order_by] = None
    profile_id: Optional[order_by] = None
    sequencing_id: Optional[order_by] = None


class fingerprintable_samples_stream_cursor_input(BaseModel):
    initial_value: "fingerprintable_samples_stream_cursor_value_input"
    ordering: Optional[cursor_ordering] = None


class fingerprintable_samples_stream_cursor_value_input(BaseModel):
    bai_filepath: Optional[str] = None
    bam_filepath: Optional[str] = None
    cell_line_name: Optional[str] = None
    expected_type: Optional[str] = None
    hg_19_bai_filepath: Optional[str] = Field(alias="hg19_bai_filepath", default=None)
    hg_19_bam_filepath: Optional[str] = Field(alias="hg19_bam_filepath", default=None)
    hg_38_crai_filepath: Optional[str] = Field(alias="hg38_crai_filepath", default=None)
    hg_38_cram_filepath: Optional[str] = Field(alias="hg38_cram_filepath", default=None)
    id: Optional[int] = None
    low_quality: Optional[bool] = None
    model_condition: Optional[str] = None
    model_id: Optional[str] = None
    omics_sequencing_blacklist: Optional[bool] = None
    patient_id: Optional[str] = None
    profile_blacklist_omics: Optional[bool] = None
    profile_id: Optional[str] = None
    sequencing_id: Optional[str] = None


class float8_comparison_exp(BaseModel):
    eq: Optional[Any] = Field(alias="_eq", default=None)
    gt: Optional[Any] = Field(alias="_gt", default=None)
    gte: Optional[Any] = Field(alias="_gte", default=None)
    in_: Optional[List[Any]] = Field(alias="_in", default=None)
    is_null: Optional[bool] = Field(alias="_is_null", default=None)
    lt: Optional[Any] = Field(alias="_lt", default=None)
    lte: Optional[Any] = Field(alias="_lte", default=None)
    neq: Optional[Any] = Field(alias="_neq", default=None)
    nin: Optional[List[Any]] = Field(alias="_nin", default=None)


class jsonb_cast_exp(BaseModel):
    string: Optional["String_comparison_exp"] = Field(alias="String", default=None)


class jsonb_comparison_exp(BaseModel):
    cast: Optional["jsonb_cast_exp"] = Field(alias="_cast", default=None)
    contained_in: Optional[Any] = Field(alias="_contained_in", default=None)
    contains: Optional[Any] = Field(alias="_contains", default=None)
    eq: Optional[Any] = Field(alias="_eq", default=None)
    gt: Optional[Any] = Field(alias="_gt", default=None)
    gte: Optional[Any] = Field(alias="_gte", default=None)
    has_key: Optional[str] = Field(alias="_has_key", default=None)
    has_keys_all: Optional[List[str]] = Field(alias="_has_keys_all", default=None)
    has_keys_any: Optional[List[str]] = Field(alias="_has_keys_any", default=None)
    in_: Optional[List[Any]] = Field(alias="_in", default=None)
    is_null: Optional[bool] = Field(alias="_is_null", default=None)
    lt: Optional[Any] = Field(alias="_lt", default=None)
    lte: Optional[Any] = Field(alias="_lte", default=None)
    neq: Optional[Any] = Field(alias="_neq", default=None)
    nin: Optional[List[Any]] = Field(alias="_nin", default=None)


class media_bool_exp(BaseModel):
    and_: Optional[List["media_bool_exp"]] = Field(alias="_and", default=None)
    not_: Optional["media_bool_exp"] = Field(alias="_not", default=None)
    or_: Optional[List["media_bool_exp"]] = Field(alias="_or", default=None)
    formulation: Optional["String_comparison_exp"] = None
    media_id: Optional["String_comparison_exp"] = None
    serum_free: Optional["Boolean_comparison_exp"] = None


class media_insert_input(BaseModel):
    formulation: Optional[str] = None
    media_id: Optional[str] = None
    serum_free: Optional[bool] = None


class media_obj_rel_insert_input(BaseModel):
    data: "media_insert_input"
    on_conflict: Optional["media_on_conflict"] = None


class media_on_conflict(BaseModel):
    constraint: media_constraint
    update_columns: List[media_update_column]
    where: Optional["media_bool_exp"] = None


class media_order_by(BaseModel):
    formulation: Optional[order_by] = None
    media_id: Optional[order_by] = None
    serum_free: Optional[order_by] = None


class media_pk_columns_input(BaseModel):
    media_id: str


class media_set_input(BaseModel):
    formulation: Optional[str] = None
    media_id: Optional[str] = None
    serum_free: Optional[bool] = None


class media_stream_cursor_input(BaseModel):
    initial_value: "media_stream_cursor_value_input"
    ordering: Optional[cursor_ordering] = None


class media_stream_cursor_value_input(BaseModel):
    formulation: Optional[str] = None
    media_id: Optional[str] = None
    serum_free: Optional[bool] = None


class media_updates(BaseModel):
    set: Optional["media_set_input"] = Field(alias="_set", default=None)
    where: "media_bool_exp"


class model_bool_exp(BaseModel):
    and_: Optional[List["model_bool_exp"]] = Field(alias="_and", default=None)
    not_: Optional["model_bool_exp"] = Field(alias="_not", default=None)
    or_: Optional[List["model_bool_exp"]] = Field(alias="_or", default=None)
    age: Optional["Int_comparison_exp"] = None
    age_category: Optional["String_comparison_exp"] = None
    ancestry: Optional["String_comparison_exp"] = None
    catalog_number: Optional["String_comparison_exp"] = None
    ccle_line: Optional["Boolean_comparison_exp"] = None
    ccle_name: Optional["String_comparison_exp"] = None
    cell_line_aliases: Optional["String_comparison_exp"] = None
    cell_line_in_stock: Optional["Boolean_comparison_exp"] = None
    cell_line_name: Optional["String_comparison_exp"] = None
    cell_line_ordered_date: Optional["date_comparison_exp"] = None
    cell_line_received: Optional["String_comparison_exp"] = None
    comments: Optional["String_comparison_exp"] = None
    consent_2015: Optional["String_comparison_exp"] = None
    converge_id: Optional["String_comparison_exp"] = None
    cosmic_id: Optional["Int_comparison_exp"] = None
    cultured_drug_resistance: Optional["String_comparison_exp"] = None
    date_cell_line_received: Optional["date_comparison_exp"] = None
    date_first_publication: Optional["date_comparison_exp"] = None
    date_model_derived: Optional["date_comparison_exp"] = None
    date_shared_in_dbgap: Optional["date_comparison_exp"] = None
    dbgap: Optional["String_comparison_exp"] = None
    depmap_model_type_id: Optional["String_comparison_exp"] = None
    derived_outside_us: Optional["Boolean_comparison_exp"] = None
    engineered_model: Optional["String_comparison_exp"] = None
    engineered_model_details: Optional["String_comparison_exp"] = None
    first_publication_link: Optional["String_comparison_exp"] = None
    geo_loc: Optional["String_comparison_exp"] = None
    growth_pattern: Optional["String_comparison_exp"] = None
    hcmi_id: Optional["String_comparison_exp"] = None
    inferred_ethnicity: Optional["String_comparison_exp"] = None
    lineage: Optional["String_comparison_exp"] = None
    medium: Optional["media_bool_exp"] = None
    model_conditions: Optional["model_condition_bool_exp"] = None
    model_conditions_aggregate: Optional["model_condition_aggregate_bool_exp"] = None
    model_data_sharing: Optional["String_comparison_exp"] = None
    model_data_sharing_comments: Optional["String_comparison_exp"] = None
    model_derivation_material: Optional["String_comparison_exp"] = None
    model_id: Optional["String_comparison_exp"] = None
    model_subtype_features: Optional["String_comparison_exp"] = None
    model_transfer: Optional["String_comparison_exp"] = None
    model_transfer_comments: Optional["String_comparison_exp"] = None
    model_transferred_to_stjude: Optional["String_comparison_exp"] = None
    model_type: Optional["String_comparison_exp"] = None
    molecular_subtype: Optional["String_comparison_exp"] = None
    ncit_code: Optional["String_comparison_exp"] = None
    ncit_subtype: Optional["String_comparison_exp"] = None
    new_histological_subtype: Optional["String_comparison_exp"] = None
    new_molecular_subtype: Optional["String_comparison_exp"] = None
    onboarded_doubling_time: Optional["String_comparison_exp"] = None
    onboarded_media: Optional["String_comparison_exp"] = None
    orspid: Optional["String_comparison_exp"] = None
    part_of_prism: Optional["String_comparison_exp"] = None
    patient_id: Optional["String_comparison_exp"] = None
    patient_resistance: Optional["String_comparison_exp"] = None
    patient_response_score: Optional["String_comparison_exp"] = None
    patient_response_score_system: Optional["String_comparison_exp"] = None
    patient_subtype_features: Optional["String_comparison_exp"] = None
    patient_treatment_type: Optional["String_comparison_exp"] = None
    patient_tumor_grade: Optional["String_comparison_exp"] = None
    peddep_line: Optional["Boolean_comparison_exp"] = None
    peddep_subgroup: Optional["String_comparison_exp"] = None
    permission_to_release: Optional["Boolean_comparison_exp"] = None
    plate_coating: Optional["String_comparison_exp"] = None
    primary_diagnosis: Optional["String_comparison_exp"] = None
    primary_disease: Optional["String_comparison_exp"] = None
    primary_or_metastasis: Optional["String_comparison_exp"] = None
    proposed_deliverable: Optional["String_comparison_exp"] = None
    proposed_release_date: Optional["date_comparison_exp"] = None
    public_comments: Optional["String_comparison_exp"] = None
    recurrent: Optional["String_comparison_exp"] = None
    registration_complete: Optional["Boolean_comparison_exp"] = None
    rrid: Optional["String_comparison_exp"] = None
    sample_collection_site: Optional["String_comparison_exp"] = None
    sanger_model_id: Optional["String_comparison_exp"] = None
    sex: Optional["String_comparison_exp"] = None
    sj_compbio_id: Optional["String_comparison_exp"] = None
    source_detail: Optional["String_comparison_exp"] = None
    source_type: Optional["String_comparison_exp"] = None
    stage: Optional["String_comparison_exp"] = None
    staging_system: Optional["String_comparison_exp"] = None
    stated_race: Optional["String_comparison_exp"] = None
    stjude_derived: Optional["Boolean_comparison_exp"] = None
    stripped_cell_line_name: Optional["String_comparison_exp"] = None
    sub_subtype: Optional["String_comparison_exp"] = None
    subtype: Optional["String_comparison_exp"] = None
    the_depmap_model_type: Optional["depmap_model_type_bool_exp"] = None
    tissue_origin: Optional["String_comparison_exp"] = None
    transformed_type: Optional["String_comparison_exp"] = None
    treatment_details: Optional["String_comparison_exp"] = None
    treatment_status: Optional["String_comparison_exp"] = None
    tumor_regression_score: Optional["String_comparison_exp"] = None
    vendor_link: Optional["String_comparison_exp"] = None
    wtsi_master_cell_id: Optional["Int_comparison_exp"] = None


class model_condition_aggregate_bool_exp(BaseModel):
    bool_and: Optional["model_condition_aggregate_bool_exp_bool_and"] = None
    bool_or: Optional["model_condition_aggregate_bool_exp_bool_or"] = None
    count: Optional["model_condition_aggregate_bool_exp_count"] = None


class model_condition_aggregate_bool_exp_bool_and(BaseModel):
    arguments: model_condition_select_column_model_condition_aggregate_bool_exp_bool_and_arguments_columns
    distinct: Optional[bool] = None
    filter: Optional["model_condition_bool_exp"] = None
    predicate: "Boolean_comparison_exp"


class model_condition_aggregate_bool_exp_bool_or(BaseModel):
    arguments: model_condition_select_column_model_condition_aggregate_bool_exp_bool_or_arguments_columns
    distinct: Optional[bool] = None
    filter: Optional["model_condition_bool_exp"] = None
    predicate: "Boolean_comparison_exp"


class model_condition_aggregate_bool_exp_count(BaseModel):
    arguments: Optional[List[model_condition_select_column]] = None
    distinct: Optional[bool] = None
    filter: Optional["model_condition_bool_exp"] = None
    predicate: "Int_comparison_exp"


class model_condition_aggregate_order_by(BaseModel):
    avg: Optional["model_condition_avg_order_by"] = None
    count: Optional[order_by] = None
    max: Optional["model_condition_max_order_by"] = None
    min: Optional["model_condition_min_order_by"] = None
    stddev: Optional["model_condition_stddev_order_by"] = None
    stddev_pop: Optional["model_condition_stddev_pop_order_by"] = None
    stddev_samp: Optional["model_condition_stddev_samp_order_by"] = None
    sum: Optional["model_condition_sum_order_by"] = None
    var_pop: Optional["model_condition_var_pop_order_by"] = None
    var_samp: Optional["model_condition_var_samp_order_by"] = None
    variance: Optional["model_condition_variance_order_by"] = None


class model_condition_arr_rel_insert_input(BaseModel):
    data: List["model_condition_insert_input"]
    on_conflict: Optional["model_condition_on_conflict"] = None


class model_condition_avg_order_by(BaseModel):
    batch_doubling_time: Optional[order_by] = None
    number_vials_available: Optional[order_by] = None
    source_doubling_time: Optional[order_by] = None


class model_condition_bool_exp(BaseModel):
    and_: Optional[List["model_condition_bool_exp"]] = Field(alias="_and", default=None)
    not_: Optional["model_condition_bool_exp"] = Field(alias="_not", default=None)
    or_: Optional[List["model_condition_bool_exp"]] = Field(alias="_or", default=None)
    batch_doubling_time: Optional["Int_comparison_exp"] = None
    cell_characteristics: Optional["String_comparison_exp"] = None
    cell_format: Optional["String_comparison_exp"] = None
    cell_grouping: Optional["String_comparison_exp"] = None
    cell_has_debris: Optional["String_comparison_exp"] = None
    cell_morphology: Optional["String_comparison_exp"] = None
    cell_shape: Optional["String_comparison_exp"] = None
    cell_size: Optional["String_comparison_exp"] = None
    comments: Optional["String_comparison_exp"] = None
    condition_only: Optional["String_comparison_exp"] = None
    contaminated: Optional["Boolean_comparison_exp"] = None
    contamination_details: Optional["String_comparison_exp"] = None
    days_with_drug: Optional["String_comparison_exp"] = None
    dmx_priority: Optional["String_comparison_exp"] = None
    drug: Optional["String_comparison_exp"] = None
    drug_concentration: Optional["String_comparison_exp"] = None
    expansion_completed: Optional["date_comparison_exp"] = None
    expansion_completed_date: Optional["date_comparison_exp"] = None
    expansion_issues: Optional["String_comparison_exp"] = None
    expansion_team: Optional["String_comparison_exp"] = None
    freeze_media: Optional["String_comparison_exp"] = None
    freezerpro_uid: Optional["String_comparison_exp"] = None
    growth_media: Optional["String_comparison_exp"] = None
    initials_status_pic: Optional["String_comparison_exp"] = None
    line_received_for_expansion: Optional["date_comparison_exp"] = None
    measured_survival: Optional["String_comparison_exp"] = None
    medium: Optional["media_bool_exp"] = None
    model: Optional["model_bool_exp"] = None
    model_condition_id: Optional["String_comparison_exp"] = None
    model_id: Optional["String_comparison_exp"] = None
    number_vials_available: Optional["Int_comparison_exp"] = None
    omics_profiles: Optional["omics_profile_bool_exp"] = None
    omics_profiles_aggregate: Optional["omics_profile_aggregate_bool_exp"] = None
    onboarding_myco_order: Optional["date_comparison_exp"] = None
    onboarding_str: Optional["String_comparison_exp"] = None
    onboarding_str_order: Optional["date_comparison_exp"] = None
    parent_model_condition_id: Optional["String_comparison_exp"] = None
    passage_number: Optional["String_comparison_exp"] = None
    plate_coating: Optional["String_comparison_exp"] = None
    prism_notes: Optional["String_comparison_exp"] = None
    project: Optional["String_comparison_exp"] = None
    resistance_mechanism: Optional["String_comparison_exp"] = None
    source: Optional["String_comparison_exp"] = None
    source_doubling_time: Optional["Int_comparison_exp"] = None
    source_growth_pattern: Optional["String_comparison_exp"] = None
    split_recommendation: Optional["String_comparison_exp"] = None
    supplements: Optional["String_comparison_exp"] = None
    thaw_date: Optional["date_comparison_exp"] = None
    to_gpp: Optional["Boolean_comparison_exp"] = None


class model_condition_inc_input(BaseModel):
    batch_doubling_time: Optional[int] = None
    number_vials_available: Optional[int] = None
    source_doubling_time: Optional[int] = None


class model_condition_insert_input(BaseModel):
    batch_doubling_time: Optional[int] = None
    cell_characteristics: Optional[str] = None
    cell_format: Optional[str] = None
    cell_grouping: Optional[str] = None
    cell_has_debris: Optional[str] = None
    cell_morphology: Optional[str] = None
    cell_shape: Optional[str] = None
    cell_size: Optional[str] = None
    comments: Optional[str] = None
    condition_only: Optional[str] = None
    contaminated: Optional[bool] = None
    contamination_details: Optional[str] = None
    days_with_drug: Optional[str] = None
    dmx_priority: Optional[str] = None
    drug: Optional[str] = None
    drug_concentration: Optional[str] = None
    expansion_completed: Optional[Any] = None
    expansion_completed_date: Optional[Any] = None
    expansion_issues: Optional[str] = None
    expansion_team: Optional[str] = None
    freeze_media: Optional[str] = None
    freezerpro_uid: Optional[str] = None
    growth_media: Optional[str] = None
    initials_status_pic: Optional[str] = None
    line_received_for_expansion: Optional[Any] = None
    measured_survival: Optional[str] = None
    medium: Optional["media_obj_rel_insert_input"] = None
    model: Optional["model_obj_rel_insert_input"] = None
    model_condition_id: Optional[str] = None
    model_id: Optional[str] = None
    number_vials_available: Optional[int] = None
    omics_profiles: Optional["omics_profile_arr_rel_insert_input"] = None
    onboarding_myco_order: Optional[Any] = None
    onboarding_str: Optional[str] = None
    onboarding_str_order: Optional[Any] = None
    parent_model_condition_id: Optional[str] = None
    passage_number: Optional[str] = None
    plate_coating: Optional[str] = None
    prism_notes: Optional[str] = None
    project: Optional[str] = None
    resistance_mechanism: Optional[str] = None
    source: Optional[str] = None
    source_doubling_time: Optional[int] = None
    source_growth_pattern: Optional[str] = None
    split_recommendation: Optional[str] = None
    supplements: Optional[str] = None
    thaw_date: Optional[Any] = None
    to_gpp: Optional[bool] = None


class model_condition_max_order_by(BaseModel):
    batch_doubling_time: Optional[order_by] = None
    cell_characteristics: Optional[order_by] = None
    cell_format: Optional[order_by] = None
    cell_grouping: Optional[order_by] = None
    cell_has_debris: Optional[order_by] = None
    cell_morphology: Optional[order_by] = None
    cell_shape: Optional[order_by] = None
    cell_size: Optional[order_by] = None
    comments: Optional[order_by] = None
    condition_only: Optional[order_by] = None
    contamination_details: Optional[order_by] = None
    days_with_drug: Optional[order_by] = None
    dmx_priority: Optional[order_by] = None
    drug: Optional[order_by] = None
    drug_concentration: Optional[order_by] = None
    expansion_completed: Optional[order_by] = None
    expansion_completed_date: Optional[order_by] = None
    expansion_issues: Optional[order_by] = None
    expansion_team: Optional[order_by] = None
    freeze_media: Optional[order_by] = None
    freezerpro_uid: Optional[order_by] = None
    growth_media: Optional[order_by] = None
    initials_status_pic: Optional[order_by] = None
    line_received_for_expansion: Optional[order_by] = None
    measured_survival: Optional[order_by] = None
    model_condition_id: Optional[order_by] = None
    model_id: Optional[order_by] = None
    number_vials_available: Optional[order_by] = None
    onboarding_myco_order: Optional[order_by] = None
    onboarding_str: Optional[order_by] = None
    onboarding_str_order: Optional[order_by] = None
    parent_model_condition_id: Optional[order_by] = None
    passage_number: Optional[order_by] = None
    plate_coating: Optional[order_by] = None
    prism_notes: Optional[order_by] = None
    project: Optional[order_by] = None
    resistance_mechanism: Optional[order_by] = None
    source: Optional[order_by] = None
    source_doubling_time: Optional[order_by] = None
    source_growth_pattern: Optional[order_by] = None
    split_recommendation: Optional[order_by] = None
    supplements: Optional[order_by] = None
    thaw_date: Optional[order_by] = None


class model_condition_min_order_by(BaseModel):
    batch_doubling_time: Optional[order_by] = None
    cell_characteristics: Optional[order_by] = None
    cell_format: Optional[order_by] = None
    cell_grouping: Optional[order_by] = None
    cell_has_debris: Optional[order_by] = None
    cell_morphology: Optional[order_by] = None
    cell_shape: Optional[order_by] = None
    cell_size: Optional[order_by] = None
    comments: Optional[order_by] = None
    condition_only: Optional[order_by] = None
    contamination_details: Optional[order_by] = None
    days_with_drug: Optional[order_by] = None
    dmx_priority: Optional[order_by] = None
    drug: Optional[order_by] = None
    drug_concentration: Optional[order_by] = None
    expansion_completed: Optional[order_by] = None
    expansion_completed_date: Optional[order_by] = None
    expansion_issues: Optional[order_by] = None
    expansion_team: Optional[order_by] = None
    freeze_media: Optional[order_by] = None
    freezerpro_uid: Optional[order_by] = None
    growth_media: Optional[order_by] = None
    initials_status_pic: Optional[order_by] = None
    line_received_for_expansion: Optional[order_by] = None
    measured_survival: Optional[order_by] = None
    model_condition_id: Optional[order_by] = None
    model_id: Optional[order_by] = None
    number_vials_available: Optional[order_by] = None
    onboarding_myco_order: Optional[order_by] = None
    onboarding_str: Optional[order_by] = None
    onboarding_str_order: Optional[order_by] = None
    parent_model_condition_id: Optional[order_by] = None
    passage_number: Optional[order_by] = None
    plate_coating: Optional[order_by] = None
    prism_notes: Optional[order_by] = None
    project: Optional[order_by] = None
    resistance_mechanism: Optional[order_by] = None
    source: Optional[order_by] = None
    source_doubling_time: Optional[order_by] = None
    source_growth_pattern: Optional[order_by] = None
    split_recommendation: Optional[order_by] = None
    supplements: Optional[order_by] = None
    thaw_date: Optional[order_by] = None


class model_condition_obj_rel_insert_input(BaseModel):
    data: "model_condition_insert_input"
    on_conflict: Optional["model_condition_on_conflict"] = None


class model_condition_on_conflict(BaseModel):
    constraint: model_condition_constraint
    update_columns: List[model_condition_update_column]
    where: Optional["model_condition_bool_exp"] = None


class model_condition_order_by(BaseModel):
    batch_doubling_time: Optional[order_by] = None
    cell_characteristics: Optional[order_by] = None
    cell_format: Optional[order_by] = None
    cell_grouping: Optional[order_by] = None
    cell_has_debris: Optional[order_by] = None
    cell_morphology: Optional[order_by] = None
    cell_shape: Optional[order_by] = None
    cell_size: Optional[order_by] = None
    comments: Optional[order_by] = None
    condition_only: Optional[order_by] = None
    contaminated: Optional[order_by] = None
    contamination_details: Optional[order_by] = None
    days_with_drug: Optional[order_by] = None
    dmx_priority: Optional[order_by] = None
    drug: Optional[order_by] = None
    drug_concentration: Optional[order_by] = None
    expansion_completed: Optional[order_by] = None
    expansion_completed_date: Optional[order_by] = None
    expansion_issues: Optional[order_by] = None
    expansion_team: Optional[order_by] = None
    freeze_media: Optional[order_by] = None
    freezerpro_uid: Optional[order_by] = None
    growth_media: Optional[order_by] = None
    initials_status_pic: Optional[order_by] = None
    line_received_for_expansion: Optional[order_by] = None
    measured_survival: Optional[order_by] = None
    medium: Optional["media_order_by"] = None
    model: Optional["model_order_by"] = None
    model_condition_id: Optional[order_by] = None
    model_id: Optional[order_by] = None
    number_vials_available: Optional[order_by] = None
    omics_profiles_aggregate: Optional["omics_profile_aggregate_order_by"] = None
    onboarding_myco_order: Optional[order_by] = None
    onboarding_str: Optional[order_by] = None
    onboarding_str_order: Optional[order_by] = None
    parent_model_condition_id: Optional[order_by] = None
    passage_number: Optional[order_by] = None
    plate_coating: Optional[order_by] = None
    prism_notes: Optional[order_by] = None
    project: Optional[order_by] = None
    resistance_mechanism: Optional[order_by] = None
    source: Optional[order_by] = None
    source_doubling_time: Optional[order_by] = None
    source_growth_pattern: Optional[order_by] = None
    split_recommendation: Optional[order_by] = None
    supplements: Optional[order_by] = None
    thaw_date: Optional[order_by] = None
    to_gpp: Optional[order_by] = None


class model_condition_pk_columns_input(BaseModel):
    model_condition_id: str


class model_condition_set_input(BaseModel):
    batch_doubling_time: Optional[int] = None
    cell_characteristics: Optional[str] = None
    cell_format: Optional[str] = None
    cell_grouping: Optional[str] = None
    cell_has_debris: Optional[str] = None
    cell_morphology: Optional[str] = None
    cell_shape: Optional[str] = None
    cell_size: Optional[str] = None
    comments: Optional[str] = None
    condition_only: Optional[str] = None
    contaminated: Optional[bool] = None
    contamination_details: Optional[str] = None
    days_with_drug: Optional[str] = None
    dmx_priority: Optional[str] = None
    drug: Optional[str] = None
    drug_concentration: Optional[str] = None
    expansion_completed: Optional[Any] = None
    expansion_completed_date: Optional[Any] = None
    expansion_issues: Optional[str] = None
    expansion_team: Optional[str] = None
    freeze_media: Optional[str] = None
    freezerpro_uid: Optional[str] = None
    growth_media: Optional[str] = None
    initials_status_pic: Optional[str] = None
    line_received_for_expansion: Optional[Any] = None
    measured_survival: Optional[str] = None
    model_condition_id: Optional[str] = None
    model_id: Optional[str] = None
    number_vials_available: Optional[int] = None
    onboarding_myco_order: Optional[Any] = None
    onboarding_str: Optional[str] = None
    onboarding_str_order: Optional[Any] = None
    parent_model_condition_id: Optional[str] = None
    passage_number: Optional[str] = None
    plate_coating: Optional[str] = None
    prism_notes: Optional[str] = None
    project: Optional[str] = None
    resistance_mechanism: Optional[str] = None
    source: Optional[str] = None
    source_doubling_time: Optional[int] = None
    source_growth_pattern: Optional[str] = None
    split_recommendation: Optional[str] = None
    supplements: Optional[str] = None
    thaw_date: Optional[Any] = None
    to_gpp: Optional[bool] = None


class model_condition_stddev_order_by(BaseModel):
    batch_doubling_time: Optional[order_by] = None
    number_vials_available: Optional[order_by] = None
    source_doubling_time: Optional[order_by] = None


class model_condition_stddev_pop_order_by(BaseModel):
    batch_doubling_time: Optional[order_by] = None
    number_vials_available: Optional[order_by] = None
    source_doubling_time: Optional[order_by] = None


class model_condition_stddev_samp_order_by(BaseModel):
    batch_doubling_time: Optional[order_by] = None
    number_vials_available: Optional[order_by] = None
    source_doubling_time: Optional[order_by] = None


class model_condition_stream_cursor_input(BaseModel):
    initial_value: "model_condition_stream_cursor_value_input"
    ordering: Optional[cursor_ordering] = None


class model_condition_stream_cursor_value_input(BaseModel):
    batch_doubling_time: Optional[int] = None
    cell_characteristics: Optional[str] = None
    cell_format: Optional[str] = None
    cell_grouping: Optional[str] = None
    cell_has_debris: Optional[str] = None
    cell_morphology: Optional[str] = None
    cell_shape: Optional[str] = None
    cell_size: Optional[str] = None
    comments: Optional[str] = None
    condition_only: Optional[str] = None
    contaminated: Optional[bool] = None
    contamination_details: Optional[str] = None
    days_with_drug: Optional[str] = None
    dmx_priority: Optional[str] = None
    drug: Optional[str] = None
    drug_concentration: Optional[str] = None
    expansion_completed: Optional[Any] = None
    expansion_completed_date: Optional[Any] = None
    expansion_issues: Optional[str] = None
    expansion_team: Optional[str] = None
    freeze_media: Optional[str] = None
    freezerpro_uid: Optional[str] = None
    growth_media: Optional[str] = None
    initials_status_pic: Optional[str] = None
    line_received_for_expansion: Optional[Any] = None
    measured_survival: Optional[str] = None
    model_condition_id: Optional[str] = None
    model_id: Optional[str] = None
    number_vials_available: Optional[int] = None
    onboarding_myco_order: Optional[Any] = None
    onboarding_str: Optional[str] = None
    onboarding_str_order: Optional[Any] = None
    parent_model_condition_id: Optional[str] = None
    passage_number: Optional[str] = None
    plate_coating: Optional[str] = None
    prism_notes: Optional[str] = None
    project: Optional[str] = None
    resistance_mechanism: Optional[str] = None
    source: Optional[str] = None
    source_doubling_time: Optional[int] = None
    source_growth_pattern: Optional[str] = None
    split_recommendation: Optional[str] = None
    supplements: Optional[str] = None
    thaw_date: Optional[Any] = None
    to_gpp: Optional[bool] = None


class model_condition_sum_order_by(BaseModel):
    batch_doubling_time: Optional[order_by] = None
    number_vials_available: Optional[order_by] = None
    source_doubling_time: Optional[order_by] = None


class model_condition_updates(BaseModel):
    inc: Optional["model_condition_inc_input"] = Field(alias="_inc", default=None)
    set: Optional["model_condition_set_input"] = Field(alias="_set", default=None)
    where: "model_condition_bool_exp"


class model_condition_var_pop_order_by(BaseModel):
    batch_doubling_time: Optional[order_by] = None
    number_vials_available: Optional[order_by] = None
    source_doubling_time: Optional[order_by] = None


class model_condition_var_samp_order_by(BaseModel):
    batch_doubling_time: Optional[order_by] = None
    number_vials_available: Optional[order_by] = None
    source_doubling_time: Optional[order_by] = None


class model_condition_variance_order_by(BaseModel):
    batch_doubling_time: Optional[order_by] = None
    number_vials_available: Optional[order_by] = None
    source_doubling_time: Optional[order_by] = None


class model_inc_input(BaseModel):
    age: Optional[int] = None
    cosmic_id: Optional[int] = None
    wtsi_master_cell_id: Optional[int] = None


class model_insert_input(BaseModel):
    age: Optional[int] = None
    age_category: Optional[str] = None
    ancestry: Optional[str] = None
    catalog_number: Optional[str] = None
    ccle_line: Optional[bool] = None
    ccle_name: Optional[str] = None
    cell_line_aliases: Optional[str] = None
    cell_line_in_stock: Optional[bool] = None
    cell_line_name: Optional[str] = None
    cell_line_ordered_date: Optional[Any] = None
    cell_line_received: Optional[str] = None
    comments: Optional[str] = None
    consent_2015: Optional[str] = None
    converge_id: Optional[str] = None
    cosmic_id: Optional[int] = None
    cultured_drug_resistance: Optional[str] = None
    date_cell_line_received: Optional[Any] = None
    date_first_publication: Optional[Any] = None
    date_model_derived: Optional[Any] = None
    date_shared_in_dbgap: Optional[Any] = None
    dbgap: Optional[str] = None
    depmap_model_type_id: Optional[str] = None
    derived_outside_us: Optional[bool] = None
    engineered_model: Optional[str] = None
    engineered_model_details: Optional[str] = None
    first_publication_link: Optional[str] = None
    geo_loc: Optional[str] = None
    growth_pattern: Optional[str] = None
    hcmi_id: Optional[str] = None
    inferred_ethnicity: Optional[str] = None
    lineage: Optional[str] = None
    medium: Optional["media_obj_rel_insert_input"] = None
    model_conditions: Optional["model_condition_arr_rel_insert_input"] = None
    model_data_sharing: Optional[str] = None
    model_data_sharing_comments: Optional[str] = None
    model_derivation_material: Optional[str] = None
    model_id: Optional[str] = None
    model_subtype_features: Optional[str] = None
    model_transfer: Optional[str] = None
    model_transfer_comments: Optional[str] = None
    model_transferred_to_stjude: Optional[str] = None
    model_type: Optional[str] = None
    molecular_subtype: Optional[str] = None
    ncit_code: Optional[str] = None
    ncit_subtype: Optional[str] = None
    new_histological_subtype: Optional[str] = None
    new_molecular_subtype: Optional[str] = None
    onboarded_doubling_time: Optional[str] = None
    onboarded_media: Optional[str] = None
    orspid: Optional[str] = None
    part_of_prism: Optional[str] = None
    patient_id: Optional[str] = None
    patient_resistance: Optional[str] = None
    patient_response_score: Optional[str] = None
    patient_response_score_system: Optional[str] = None
    patient_subtype_features: Optional[str] = None
    patient_treatment_type: Optional[str] = None
    patient_tumor_grade: Optional[str] = None
    peddep_line: Optional[bool] = None
    peddep_subgroup: Optional[str] = None
    permission_to_release: Optional[bool] = None
    plate_coating: Optional[str] = None
    primary_diagnosis: Optional[str] = None
    primary_disease: Optional[str] = None
    primary_or_metastasis: Optional[str] = None
    proposed_deliverable: Optional[str] = None
    proposed_release_date: Optional[Any] = None
    public_comments: Optional[str] = None
    recurrent: Optional[str] = None
    registration_complete: Optional[bool] = None
    rrid: Optional[str] = None
    sample_collection_site: Optional[str] = None
    sanger_model_id: Optional[str] = None
    sex: Optional[str] = None
    sj_compbio_id: Optional[str] = None
    source_detail: Optional[str] = None
    source_type: Optional[str] = None
    stage: Optional[str] = None
    staging_system: Optional[str] = None
    stated_race: Optional[str] = None
    stjude_derived: Optional[bool] = None
    stripped_cell_line_name: Optional[str] = None
    sub_subtype: Optional[str] = None
    subtype: Optional[str] = None
    the_depmap_model_type: Optional["depmap_model_type_obj_rel_insert_input"] = None
    tissue_origin: Optional[str] = None
    transformed_type: Optional[str] = None
    treatment_details: Optional[str] = None
    treatment_status: Optional[str] = None
    tumor_regression_score: Optional[str] = None
    vendor_link: Optional[str] = None
    wtsi_master_cell_id: Optional[int] = None


class model_obj_rel_insert_input(BaseModel):
    data: "model_insert_input"
    on_conflict: Optional["model_on_conflict"] = None


class model_on_conflict(BaseModel):
    constraint: model_constraint
    update_columns: List[model_update_column]
    where: Optional["model_bool_exp"] = None


class model_order_by(BaseModel):
    age: Optional[order_by] = None
    age_category: Optional[order_by] = None
    ancestry: Optional[order_by] = None
    catalog_number: Optional[order_by] = None
    ccle_line: Optional[order_by] = None
    ccle_name: Optional[order_by] = None
    cell_line_aliases: Optional[order_by] = None
    cell_line_in_stock: Optional[order_by] = None
    cell_line_name: Optional[order_by] = None
    cell_line_ordered_date: Optional[order_by] = None
    cell_line_received: Optional[order_by] = None
    comments: Optional[order_by] = None
    consent_2015: Optional[order_by] = None
    converge_id: Optional[order_by] = None
    cosmic_id: Optional[order_by] = None
    cultured_drug_resistance: Optional[order_by] = None
    date_cell_line_received: Optional[order_by] = None
    date_first_publication: Optional[order_by] = None
    date_model_derived: Optional[order_by] = None
    date_shared_in_dbgap: Optional[order_by] = None
    dbgap: Optional[order_by] = None
    depmap_model_type_id: Optional[order_by] = None
    derived_outside_us: Optional[order_by] = None
    engineered_model: Optional[order_by] = None
    engineered_model_details: Optional[order_by] = None
    first_publication_link: Optional[order_by] = None
    geo_loc: Optional[order_by] = None
    growth_pattern: Optional[order_by] = None
    hcmi_id: Optional[order_by] = None
    inferred_ethnicity: Optional[order_by] = None
    lineage: Optional[order_by] = None
    medium: Optional["media_order_by"] = None
    model_conditions_aggregate: Optional["model_condition_aggregate_order_by"] = None
    model_data_sharing: Optional[order_by] = None
    model_data_sharing_comments: Optional[order_by] = None
    model_derivation_material: Optional[order_by] = None
    model_id: Optional[order_by] = None
    model_subtype_features: Optional[order_by] = None
    model_transfer: Optional[order_by] = None
    model_transfer_comments: Optional[order_by] = None
    model_transferred_to_stjude: Optional[order_by] = None
    model_type: Optional[order_by] = None
    molecular_subtype: Optional[order_by] = None
    ncit_code: Optional[order_by] = None
    ncit_subtype: Optional[order_by] = None
    new_histological_subtype: Optional[order_by] = None
    new_molecular_subtype: Optional[order_by] = None
    onboarded_doubling_time: Optional[order_by] = None
    onboarded_media: Optional[order_by] = None
    orspid: Optional[order_by] = None
    part_of_prism: Optional[order_by] = None
    patient_id: Optional[order_by] = None
    patient_resistance: Optional[order_by] = None
    patient_response_score: Optional[order_by] = None
    patient_response_score_system: Optional[order_by] = None
    patient_subtype_features: Optional[order_by] = None
    patient_treatment_type: Optional[order_by] = None
    patient_tumor_grade: Optional[order_by] = None
    peddep_line: Optional[order_by] = None
    peddep_subgroup: Optional[order_by] = None
    permission_to_release: Optional[order_by] = None
    plate_coating: Optional[order_by] = None
    primary_diagnosis: Optional[order_by] = None
    primary_disease: Optional[order_by] = None
    primary_or_metastasis: Optional[order_by] = None
    proposed_deliverable: Optional[order_by] = None
    proposed_release_date: Optional[order_by] = None
    public_comments: Optional[order_by] = None
    recurrent: Optional[order_by] = None
    registration_complete: Optional[order_by] = None
    rrid: Optional[order_by] = None
    sample_collection_site: Optional[order_by] = None
    sanger_model_id: Optional[order_by] = None
    sex: Optional[order_by] = None
    sj_compbio_id: Optional[order_by] = None
    source_detail: Optional[order_by] = None
    source_type: Optional[order_by] = None
    stage: Optional[order_by] = None
    staging_system: Optional[order_by] = None
    stated_race: Optional[order_by] = None
    stjude_derived: Optional[order_by] = None
    stripped_cell_line_name: Optional[order_by] = None
    sub_subtype: Optional[order_by] = None
    subtype: Optional[order_by] = None
    the_depmap_model_type: Optional["depmap_model_type_order_by"] = None
    tissue_origin: Optional[order_by] = None
    transformed_type: Optional[order_by] = None
    treatment_details: Optional[order_by] = None
    treatment_status: Optional[order_by] = None
    tumor_regression_score: Optional[order_by] = None
    vendor_link: Optional[order_by] = None
    wtsi_master_cell_id: Optional[order_by] = None


class model_pk_columns_input(BaseModel):
    model_id: str


class model_set_input(BaseModel):
    age: Optional[int] = None
    age_category: Optional[str] = None
    ancestry: Optional[str] = None
    catalog_number: Optional[str] = None
    ccle_line: Optional[bool] = None
    ccle_name: Optional[str] = None
    cell_line_aliases: Optional[str] = None
    cell_line_in_stock: Optional[bool] = None
    cell_line_name: Optional[str] = None
    cell_line_ordered_date: Optional[Any] = None
    cell_line_received: Optional[str] = None
    comments: Optional[str] = None
    consent_2015: Optional[str] = None
    converge_id: Optional[str] = None
    cosmic_id: Optional[int] = None
    cultured_drug_resistance: Optional[str] = None
    date_cell_line_received: Optional[Any] = None
    date_first_publication: Optional[Any] = None
    date_model_derived: Optional[Any] = None
    date_shared_in_dbgap: Optional[Any] = None
    dbgap: Optional[str] = None
    depmap_model_type_id: Optional[str] = None
    derived_outside_us: Optional[bool] = None
    engineered_model: Optional[str] = None
    engineered_model_details: Optional[str] = None
    first_publication_link: Optional[str] = None
    geo_loc: Optional[str] = None
    growth_pattern: Optional[str] = None
    hcmi_id: Optional[str] = None
    inferred_ethnicity: Optional[str] = None
    lineage: Optional[str] = None
    model_data_sharing: Optional[str] = None
    model_data_sharing_comments: Optional[str] = None
    model_derivation_material: Optional[str] = None
    model_id: Optional[str] = None
    model_subtype_features: Optional[str] = None
    model_transfer: Optional[str] = None
    model_transfer_comments: Optional[str] = None
    model_transferred_to_stjude: Optional[str] = None
    model_type: Optional[str] = None
    molecular_subtype: Optional[str] = None
    ncit_code: Optional[str] = None
    ncit_subtype: Optional[str] = None
    new_histological_subtype: Optional[str] = None
    new_molecular_subtype: Optional[str] = None
    onboarded_doubling_time: Optional[str] = None
    onboarded_media: Optional[str] = None
    orspid: Optional[str] = None
    part_of_prism: Optional[str] = None
    patient_id: Optional[str] = None
    patient_resistance: Optional[str] = None
    patient_response_score: Optional[str] = None
    patient_response_score_system: Optional[str] = None
    patient_subtype_features: Optional[str] = None
    patient_treatment_type: Optional[str] = None
    patient_tumor_grade: Optional[str] = None
    peddep_line: Optional[bool] = None
    peddep_subgroup: Optional[str] = None
    permission_to_release: Optional[bool] = None
    plate_coating: Optional[str] = None
    primary_diagnosis: Optional[str] = None
    primary_disease: Optional[str] = None
    primary_or_metastasis: Optional[str] = None
    proposed_deliverable: Optional[str] = None
    proposed_release_date: Optional[Any] = None
    public_comments: Optional[str] = None
    recurrent: Optional[str] = None
    registration_complete: Optional[bool] = None
    rrid: Optional[str] = None
    sample_collection_site: Optional[str] = None
    sanger_model_id: Optional[str] = None
    sex: Optional[str] = None
    sj_compbio_id: Optional[str] = None
    source_detail: Optional[str] = None
    source_type: Optional[str] = None
    stage: Optional[str] = None
    staging_system: Optional[str] = None
    stated_race: Optional[str] = None
    stjude_derived: Optional[bool] = None
    stripped_cell_line_name: Optional[str] = None
    sub_subtype: Optional[str] = None
    subtype: Optional[str] = None
    tissue_origin: Optional[str] = None
    transformed_type: Optional[str] = None
    treatment_details: Optional[str] = None
    treatment_status: Optional[str] = None
    tumor_regression_score: Optional[str] = None
    vendor_link: Optional[str] = None
    wtsi_master_cell_id: Optional[int] = None


class model_stream_cursor_input(BaseModel):
    initial_value: "model_stream_cursor_value_input"
    ordering: Optional[cursor_ordering] = None


class model_stream_cursor_value_input(BaseModel):
    age: Optional[int] = None
    age_category: Optional[str] = None
    ancestry: Optional[str] = None
    catalog_number: Optional[str] = None
    ccle_line: Optional[bool] = None
    ccle_name: Optional[str] = None
    cell_line_aliases: Optional[str] = None
    cell_line_in_stock: Optional[bool] = None
    cell_line_name: Optional[str] = None
    cell_line_ordered_date: Optional[Any] = None
    cell_line_received: Optional[str] = None
    comments: Optional[str] = None
    consent_2015: Optional[str] = None
    converge_id: Optional[str] = None
    cosmic_id: Optional[int] = None
    cultured_drug_resistance: Optional[str] = None
    date_cell_line_received: Optional[Any] = None
    date_first_publication: Optional[Any] = None
    date_model_derived: Optional[Any] = None
    date_shared_in_dbgap: Optional[Any] = None
    dbgap: Optional[str] = None
    depmap_model_type_id: Optional[str] = None
    derived_outside_us: Optional[bool] = None
    engineered_model: Optional[str] = None
    engineered_model_details: Optional[str] = None
    first_publication_link: Optional[str] = None
    geo_loc: Optional[str] = None
    growth_pattern: Optional[str] = None
    hcmi_id: Optional[str] = None
    inferred_ethnicity: Optional[str] = None
    lineage: Optional[str] = None
    model_data_sharing: Optional[str] = None
    model_data_sharing_comments: Optional[str] = None
    model_derivation_material: Optional[str] = None
    model_id: Optional[str] = None
    model_subtype_features: Optional[str] = None
    model_transfer: Optional[str] = None
    model_transfer_comments: Optional[str] = None
    model_transferred_to_stjude: Optional[str] = None
    model_type: Optional[str] = None
    molecular_subtype: Optional[str] = None
    ncit_code: Optional[str] = None
    ncit_subtype: Optional[str] = None
    new_histological_subtype: Optional[str] = None
    new_molecular_subtype: Optional[str] = None
    onboarded_doubling_time: Optional[str] = None
    onboarded_media: Optional[str] = None
    orspid: Optional[str] = None
    part_of_prism: Optional[str] = None
    patient_id: Optional[str] = None
    patient_resistance: Optional[str] = None
    patient_response_score: Optional[str] = None
    patient_response_score_system: Optional[str] = None
    patient_subtype_features: Optional[str] = None
    patient_treatment_type: Optional[str] = None
    patient_tumor_grade: Optional[str] = None
    peddep_line: Optional[bool] = None
    peddep_subgroup: Optional[str] = None
    permission_to_release: Optional[bool] = None
    plate_coating: Optional[str] = None
    primary_diagnosis: Optional[str] = None
    primary_disease: Optional[str] = None
    primary_or_metastasis: Optional[str] = None
    proposed_deliverable: Optional[str] = None
    proposed_release_date: Optional[Any] = None
    public_comments: Optional[str] = None
    recurrent: Optional[str] = None
    registration_complete: Optional[bool] = None
    rrid: Optional[str] = None
    sample_collection_site: Optional[str] = None
    sanger_model_id: Optional[str] = None
    sex: Optional[str] = None
    sj_compbio_id: Optional[str] = None
    source_detail: Optional[str] = None
    source_type: Optional[str] = None
    stage: Optional[str] = None
    staging_system: Optional[str] = None
    stated_race: Optional[str] = None
    stjude_derived: Optional[bool] = None
    stripped_cell_line_name: Optional[str] = None
    sub_subtype: Optional[str] = None
    subtype: Optional[str] = None
    tissue_origin: Optional[str] = None
    transformed_type: Optional[str] = None
    treatment_details: Optional[str] = None
    treatment_status: Optional[str] = None
    tumor_regression_score: Optional[str] = None
    vendor_link: Optional[str] = None
    wtsi_master_cell_id: Optional[int] = None


class model_updates(BaseModel):
    inc: Optional["model_inc_input"] = Field(alias="_inc", default=None)
    set: Optional["model_set_input"] = Field(alias="_set", default=None)
    where: "model_bool_exp"


class omics_profile_aggregate_bool_exp(BaseModel):
    bool_and: Optional["omics_profile_aggregate_bool_exp_bool_and"] = None
    bool_or: Optional["omics_profile_aggregate_bool_exp_bool_or"] = None
    count: Optional["omics_profile_aggregate_bool_exp_count"] = None


class omics_profile_aggregate_bool_exp_bool_and(BaseModel):
    arguments: omics_profile_select_column_omics_profile_aggregate_bool_exp_bool_and_arguments_columns
    distinct: Optional[bool] = None
    filter: Optional["omics_profile_bool_exp"] = None
    predicate: "Boolean_comparison_exp"


class omics_profile_aggregate_bool_exp_bool_or(BaseModel):
    arguments: omics_profile_select_column_omics_profile_aggregate_bool_exp_bool_or_arguments_columns
    distinct: Optional[bool] = None
    filter: Optional["omics_profile_bool_exp"] = None
    predicate: "Boolean_comparison_exp"


class omics_profile_aggregate_bool_exp_count(BaseModel):
    arguments: Optional[List[omics_profile_select_column]] = None
    distinct: Optional[bool] = None
    filter: Optional["omics_profile_bool_exp"] = None
    predicate: "Int_comparison_exp"


class omics_profile_aggregate_order_by(BaseModel):
    count: Optional[order_by] = None
    max: Optional["omics_profile_max_order_by"] = None
    min: Optional["omics_profile_min_order_by"] = None


class omics_profile_arr_rel_insert_input(BaseModel):
    data: List["omics_profile_insert_input"]
    on_conflict: Optional["omics_profile_on_conflict"] = None


class omics_profile_bool_exp(BaseModel):
    and_: Optional[List["omics_profile_bool_exp"]] = Field(alias="_and", default=None)
    not_: Optional["omics_profile_bool_exp"] = Field(alias="_not", default=None)
    or_: Optional[List["omics_profile_bool_exp"]] = Field(alias="_or", default=None)
    actual_seq_technology: Optional["String_comparison_exp"] = None
    baits: Optional["String_comparison_exp"] = None
    bam_public_sra_path: Optional["String_comparison_exp"] = None
    billing_date: Optional["date_comparison_exp"] = None
    blacklist_expiration_date: Optional["date_comparison_exp"] = None
    blacklist_omics: Optional["Boolean_comparison_exp"] = None
    blacklist_reason: Optional["String_comparison_exp"] = None
    bsp_sample_id_csv: Optional["String_comparison_exp"] = None
    cell_available: Optional["Boolean_comparison_exp"] = None
    cell_pellet_needed: Optional["Boolean_comparison_exp"] = None
    collaborator_sample_id: Optional["String_comparison_exp"] = None
    consortium_release_date: Optional["date_comparison_exp"] = None
    consortium_retracted_date: Optional["date_comparison_exp"] = None
    datatype: Optional["String_comparison_exp"] = None
    deliverables: Optional["String_comparison_exp"] = None
    destination_datasets: Optional["String_comparison_exp"] = None
    drop_reason: Optional["String_comparison_exp"] = None
    eta_for_omics_completion: Optional["date_comparison_exp"] = None
    extraction_needed: Optional["Boolean_comparison_exp"] = None
    ibm_release_date: Optional["date_comparison_exp"] = None
    internal_release_date: Optional["date_comparison_exp"] = None
    internal_retracted_date: Optional["date_comparison_exp"] = None
    issue: Optional["String_comparison_exp"] = None
    kit_id: Optional["String_comparison_exp"] = None
    lcset_protocol: Optional["String_comparison_exp"] = None
    lcsets: Optional["String_comparison_exp"] = None
    line_received_by_gp: Optional["date_comparison_exp"] = None
    line_sent_to_gp: Optional["date_comparison_exp"] = None
    main_sequencing_id: Optional["String_comparison_exp"] = None
    model_condition_id: Optional["String_comparison_exp"] = None
    omics_order_date: Optional["date_comparison_exp"] = None
    omics_profile_flagship: Optional["String_comparison_exp"] = None
    omics_profile_funding_source: Optional["String_comparison_exp"] = None
    omics_return_date: Optional["date_comparison_exp"] = None
    omics_sequencings: Optional["omics_sequencing_bool_exp"] = None
    omics_sequencings_aggregate: Optional["omics_sequencing_aggregate_bool_exp"] = None
    pdo_title: Optional["String_comparison_exp"] = None
    pdoid: Optional["String_comparison_exp"] = None
    pf_bases_bc: Optional["String_comparison_exp"] = None
    prioritized: Optional["Boolean_comparison_exp"] = None
    product: Optional["String_comparison_exp"] = None
    product_goal: Optional["String_comparison_exp"] = None
    profile_id: Optional["String_comparison_exp"] = None
    profile_source: Optional["String_comparison_exp"] = None
    project: Optional["String_comparison_exp"] = None
    proposed_release_date: Optional["date_comparison_exp"] = None
    public_release_date: Optional["date_comparison_exp"] = None
    public_retracted_date: Optional["date_comparison_exp"] = None
    quote_to_bill: Optional["String_comparison_exp"] = None
    registered: Optional["Boolean_comparison_exp"] = None
    resubmit_for_extraction: Optional["Boolean_comparison_exp"] = None
    rna_delivery_date: Optional["date_comparison_exp"] = None
    sample_coverage_normalized: Optional["String_comparison_exp"] = None
    sample_coverage_rounded: Optional["String_comparison_exp"] = None
    sample_is_on_risk: Optional["Boolean_comparison_exp"] = None
    sample_type: Optional["String_comparison_exp"] = None
    sm_id_matched: Optional["String_comparison_exp"] = None
    smid_ordered: Optional["String_comparison_exp"] = None
    smid_returned: Optional["String_comparison_exp"] = None
    status: Optional["String_comparison_exp"] = None
    the_model_condition: Optional["model_condition_bool_exp"] = None
    version: Optional["String_comparison_exp"] = None
    wgs_delivery_date: Optional["date_comparison_exp"] = None
    workspace: Optional["String_comparison_exp"] = None


class omics_profile_insert_input(BaseModel):
    actual_seq_technology: Optional[str] = None
    baits: Optional[str] = None
    bam_public_sra_path: Optional[str] = None
    billing_date: Optional[Any] = None
    blacklist_expiration_date: Optional[Any] = None
    blacklist_omics: Optional[bool] = None
    blacklist_reason: Optional[str] = None
    bsp_sample_id_csv: Optional[str] = None
    cell_available: Optional[bool] = None
    cell_pellet_needed: Optional[bool] = None
    collaborator_sample_id: Optional[str] = None
    consortium_release_date: Optional[Any] = None
    consortium_retracted_date: Optional[Any] = None
    datatype: Optional[str] = None
    deliverables: Optional[str] = None
    destination_datasets: Optional[str] = None
    drop_reason: Optional[str] = None
    eta_for_omics_completion: Optional[Any] = None
    extraction_needed: Optional[bool] = None
    ibm_release_date: Optional[Any] = None
    internal_release_date: Optional[Any] = None
    internal_retracted_date: Optional[Any] = None
    issue: Optional[str] = None
    kit_id: Optional[str] = None
    lcset_protocol: Optional[str] = None
    lcsets: Optional[str] = None
    line_received_by_gp: Optional[Any] = None
    line_sent_to_gp: Optional[Any] = None
    main_sequencing_id: Optional[str] = None
    model_condition_id: Optional[str] = None
    omics_order_date: Optional[Any] = None
    omics_profile_flagship: Optional[str] = None
    omics_profile_funding_source: Optional[str] = None
    omics_return_date: Optional[Any] = None
    omics_sequencings: Optional["omics_sequencing_arr_rel_insert_input"] = None
    pdo_title: Optional[str] = None
    pdoid: Optional[str] = None
    pf_bases_bc: Optional[str] = None
    prioritized: Optional[bool] = None
    product: Optional[str] = None
    product_goal: Optional[str] = None
    profile_id: Optional[str] = None
    profile_source: Optional[str] = None
    project: Optional[str] = None
    proposed_release_date: Optional[Any] = None
    public_release_date: Optional[Any] = None
    public_retracted_date: Optional[Any] = None
    quote_to_bill: Optional[str] = None
    registered: Optional[bool] = None
    resubmit_for_extraction: Optional[bool] = None
    rna_delivery_date: Optional[Any] = None
    sample_coverage_normalized: Optional[str] = None
    sample_coverage_rounded: Optional[str] = None
    sample_is_on_risk: Optional[bool] = None
    sample_type: Optional[str] = None
    sm_id_matched: Optional[str] = None
    smid_ordered: Optional[str] = None
    smid_returned: Optional[str] = None
    status: Optional[str] = None
    the_model_condition: Optional["model_condition_obj_rel_insert_input"] = None
    version: Optional[str] = None
    wgs_delivery_date: Optional[Any] = None
    workspace: Optional[str] = None


class omics_profile_max_order_by(BaseModel):
    actual_seq_technology: Optional[order_by] = None
    baits: Optional[order_by] = None
    bam_public_sra_path: Optional[order_by] = None
    billing_date: Optional[order_by] = None
    blacklist_expiration_date: Optional[order_by] = None
    blacklist_reason: Optional[order_by] = None
    bsp_sample_id_csv: Optional[order_by] = None
    collaborator_sample_id: Optional[order_by] = None
    consortium_release_date: Optional[order_by] = None
    consortium_retracted_date: Optional[order_by] = None
    datatype: Optional[order_by] = None
    deliverables: Optional[order_by] = None
    destination_datasets: Optional[order_by] = None
    drop_reason: Optional[order_by] = None
    eta_for_omics_completion: Optional[order_by] = None
    ibm_release_date: Optional[order_by] = None
    internal_release_date: Optional[order_by] = None
    internal_retracted_date: Optional[order_by] = None
    issue: Optional[order_by] = None
    kit_id: Optional[order_by] = None
    lcset_protocol: Optional[order_by] = None
    lcsets: Optional[order_by] = None
    line_received_by_gp: Optional[order_by] = None
    line_sent_to_gp: Optional[order_by] = None
    main_sequencing_id: Optional[order_by] = None
    model_condition_id: Optional[order_by] = None
    omics_order_date: Optional[order_by] = None
    omics_profile_flagship: Optional[order_by] = None
    omics_profile_funding_source: Optional[order_by] = None
    omics_return_date: Optional[order_by] = None
    pdo_title: Optional[order_by] = None
    pdoid: Optional[order_by] = None
    pf_bases_bc: Optional[order_by] = None
    product: Optional[order_by] = None
    product_goal: Optional[order_by] = None
    profile_id: Optional[order_by] = None
    profile_source: Optional[order_by] = None
    project: Optional[order_by] = None
    proposed_release_date: Optional[order_by] = None
    public_release_date: Optional[order_by] = None
    public_retracted_date: Optional[order_by] = None
    quote_to_bill: Optional[order_by] = None
    rna_delivery_date: Optional[order_by] = None
    sample_coverage_normalized: Optional[order_by] = None
    sample_coverage_rounded: Optional[order_by] = None
    sample_type: Optional[order_by] = None
    sm_id_matched: Optional[order_by] = None
    smid_ordered: Optional[order_by] = None
    smid_returned: Optional[order_by] = None
    status: Optional[order_by] = None
    version: Optional[order_by] = None
    wgs_delivery_date: Optional[order_by] = None
    workspace: Optional[order_by] = None


class omics_profile_min_order_by(BaseModel):
    actual_seq_technology: Optional[order_by] = None
    baits: Optional[order_by] = None
    bam_public_sra_path: Optional[order_by] = None
    billing_date: Optional[order_by] = None
    blacklist_expiration_date: Optional[order_by] = None
    blacklist_reason: Optional[order_by] = None
    bsp_sample_id_csv: Optional[order_by] = None
    collaborator_sample_id: Optional[order_by] = None
    consortium_release_date: Optional[order_by] = None
    consortium_retracted_date: Optional[order_by] = None
    datatype: Optional[order_by] = None
    deliverables: Optional[order_by] = None
    destination_datasets: Optional[order_by] = None
    drop_reason: Optional[order_by] = None
    eta_for_omics_completion: Optional[order_by] = None
    ibm_release_date: Optional[order_by] = None
    internal_release_date: Optional[order_by] = None
    internal_retracted_date: Optional[order_by] = None
    issue: Optional[order_by] = None
    kit_id: Optional[order_by] = None
    lcset_protocol: Optional[order_by] = None
    lcsets: Optional[order_by] = None
    line_received_by_gp: Optional[order_by] = None
    line_sent_to_gp: Optional[order_by] = None
    main_sequencing_id: Optional[order_by] = None
    model_condition_id: Optional[order_by] = None
    omics_order_date: Optional[order_by] = None
    omics_profile_flagship: Optional[order_by] = None
    omics_profile_funding_source: Optional[order_by] = None
    omics_return_date: Optional[order_by] = None
    pdo_title: Optional[order_by] = None
    pdoid: Optional[order_by] = None
    pf_bases_bc: Optional[order_by] = None
    product: Optional[order_by] = None
    product_goal: Optional[order_by] = None
    profile_id: Optional[order_by] = None
    profile_source: Optional[order_by] = None
    project: Optional[order_by] = None
    proposed_release_date: Optional[order_by] = None
    public_release_date: Optional[order_by] = None
    public_retracted_date: Optional[order_by] = None
    quote_to_bill: Optional[order_by] = None
    rna_delivery_date: Optional[order_by] = None
    sample_coverage_normalized: Optional[order_by] = None
    sample_coverage_rounded: Optional[order_by] = None
    sample_type: Optional[order_by] = None
    sm_id_matched: Optional[order_by] = None
    smid_ordered: Optional[order_by] = None
    smid_returned: Optional[order_by] = None
    status: Optional[order_by] = None
    version: Optional[order_by] = None
    wgs_delivery_date: Optional[order_by] = None
    workspace: Optional[order_by] = None


class omics_profile_obj_rel_insert_input(BaseModel):
    data: "omics_profile_insert_input"
    on_conflict: Optional["omics_profile_on_conflict"] = None


class omics_profile_on_conflict(BaseModel):
    constraint: omics_profile_constraint
    update_columns: List[omics_profile_update_column]
    where: Optional["omics_profile_bool_exp"] = None


class omics_profile_order_by(BaseModel):
    actual_seq_technology: Optional[order_by] = None
    baits: Optional[order_by] = None
    bam_public_sra_path: Optional[order_by] = None
    billing_date: Optional[order_by] = None
    blacklist_expiration_date: Optional[order_by] = None
    blacklist_omics: Optional[order_by] = None
    blacklist_reason: Optional[order_by] = None
    bsp_sample_id_csv: Optional[order_by] = None
    cell_available: Optional[order_by] = None
    cell_pellet_needed: Optional[order_by] = None
    collaborator_sample_id: Optional[order_by] = None
    consortium_release_date: Optional[order_by] = None
    consortium_retracted_date: Optional[order_by] = None
    datatype: Optional[order_by] = None
    deliverables: Optional[order_by] = None
    destination_datasets: Optional[order_by] = None
    drop_reason: Optional[order_by] = None
    eta_for_omics_completion: Optional[order_by] = None
    extraction_needed: Optional[order_by] = None
    ibm_release_date: Optional[order_by] = None
    internal_release_date: Optional[order_by] = None
    internal_retracted_date: Optional[order_by] = None
    issue: Optional[order_by] = None
    kit_id: Optional[order_by] = None
    lcset_protocol: Optional[order_by] = None
    lcsets: Optional[order_by] = None
    line_received_by_gp: Optional[order_by] = None
    line_sent_to_gp: Optional[order_by] = None
    main_sequencing_id: Optional[order_by] = None
    model_condition_id: Optional[order_by] = None
    omics_order_date: Optional[order_by] = None
    omics_profile_flagship: Optional[order_by] = None
    omics_profile_funding_source: Optional[order_by] = None
    omics_return_date: Optional[order_by] = None
    omics_sequencings_aggregate: Optional["omics_sequencing_aggregate_order_by"] = None
    pdo_title: Optional[order_by] = None
    pdoid: Optional[order_by] = None
    pf_bases_bc: Optional[order_by] = None
    prioritized: Optional[order_by] = None
    product: Optional[order_by] = None
    product_goal: Optional[order_by] = None
    profile_id: Optional[order_by] = None
    profile_source: Optional[order_by] = None
    project: Optional[order_by] = None
    proposed_release_date: Optional[order_by] = None
    public_release_date: Optional[order_by] = None
    public_retracted_date: Optional[order_by] = None
    quote_to_bill: Optional[order_by] = None
    registered: Optional[order_by] = None
    resubmit_for_extraction: Optional[order_by] = None
    rna_delivery_date: Optional[order_by] = None
    sample_coverage_normalized: Optional[order_by] = None
    sample_coverage_rounded: Optional[order_by] = None
    sample_is_on_risk: Optional[order_by] = None
    sample_type: Optional[order_by] = None
    sm_id_matched: Optional[order_by] = None
    smid_ordered: Optional[order_by] = None
    smid_returned: Optional[order_by] = None
    status: Optional[order_by] = None
    the_model_condition: Optional["model_condition_order_by"] = None
    version: Optional[order_by] = None
    wgs_delivery_date: Optional[order_by] = None
    workspace: Optional[order_by] = None


class omics_profile_pk_columns_input(BaseModel):
    profile_id: str


class omics_profile_set_input(BaseModel):
    actual_seq_technology: Optional[str] = None
    baits: Optional[str] = None
    bam_public_sra_path: Optional[str] = None
    billing_date: Optional[Any] = None
    blacklist_expiration_date: Optional[Any] = None
    blacklist_omics: Optional[bool] = None
    blacklist_reason: Optional[str] = None
    bsp_sample_id_csv: Optional[str] = None
    cell_available: Optional[bool] = None
    cell_pellet_needed: Optional[bool] = None
    collaborator_sample_id: Optional[str] = None
    consortium_release_date: Optional[Any] = None
    consortium_retracted_date: Optional[Any] = None
    datatype: Optional[str] = None
    deliverables: Optional[str] = None
    destination_datasets: Optional[str] = None
    drop_reason: Optional[str] = None
    eta_for_omics_completion: Optional[Any] = None
    extraction_needed: Optional[bool] = None
    ibm_release_date: Optional[Any] = None
    internal_release_date: Optional[Any] = None
    internal_retracted_date: Optional[Any] = None
    issue: Optional[str] = None
    kit_id: Optional[str] = None
    lcset_protocol: Optional[str] = None
    lcsets: Optional[str] = None
    line_received_by_gp: Optional[Any] = None
    line_sent_to_gp: Optional[Any] = None
    main_sequencing_id: Optional[str] = None
    model_condition_id: Optional[str] = None
    omics_order_date: Optional[Any] = None
    omics_profile_flagship: Optional[str] = None
    omics_profile_funding_source: Optional[str] = None
    omics_return_date: Optional[Any] = None
    pdo_title: Optional[str] = None
    pdoid: Optional[str] = None
    pf_bases_bc: Optional[str] = None
    prioritized: Optional[bool] = None
    product: Optional[str] = None
    product_goal: Optional[str] = None
    profile_id: Optional[str] = None
    profile_source: Optional[str] = None
    project: Optional[str] = None
    proposed_release_date: Optional[Any] = None
    public_release_date: Optional[Any] = None
    public_retracted_date: Optional[Any] = None
    quote_to_bill: Optional[str] = None
    registered: Optional[bool] = None
    resubmit_for_extraction: Optional[bool] = None
    rna_delivery_date: Optional[Any] = None
    sample_coverage_normalized: Optional[str] = None
    sample_coverage_rounded: Optional[str] = None
    sample_is_on_risk: Optional[bool] = None
    sample_type: Optional[str] = None
    sm_id_matched: Optional[str] = None
    smid_ordered: Optional[str] = None
    smid_returned: Optional[str] = None
    status: Optional[str] = None
    version: Optional[str] = None
    wgs_delivery_date: Optional[Any] = None
    workspace: Optional[str] = None


class omics_profile_stream_cursor_input(BaseModel):
    initial_value: "omics_profile_stream_cursor_value_input"
    ordering: Optional[cursor_ordering] = None


class omics_profile_stream_cursor_value_input(BaseModel):
    actual_seq_technology: Optional[str] = None
    baits: Optional[str] = None
    bam_public_sra_path: Optional[str] = None
    billing_date: Optional[Any] = None
    blacklist_expiration_date: Optional[Any] = None
    blacklist_omics: Optional[bool] = None
    blacklist_reason: Optional[str] = None
    bsp_sample_id_csv: Optional[str] = None
    cell_available: Optional[bool] = None
    cell_pellet_needed: Optional[bool] = None
    collaborator_sample_id: Optional[str] = None
    consortium_release_date: Optional[Any] = None
    consortium_retracted_date: Optional[Any] = None
    datatype: Optional[str] = None
    deliverables: Optional[str] = None
    destination_datasets: Optional[str] = None
    drop_reason: Optional[str] = None
    eta_for_omics_completion: Optional[Any] = None
    extraction_needed: Optional[bool] = None
    ibm_release_date: Optional[Any] = None
    internal_release_date: Optional[Any] = None
    internal_retracted_date: Optional[Any] = None
    issue: Optional[str] = None
    kit_id: Optional[str] = None
    lcset_protocol: Optional[str] = None
    lcsets: Optional[str] = None
    line_received_by_gp: Optional[Any] = None
    line_sent_to_gp: Optional[Any] = None
    main_sequencing_id: Optional[str] = None
    model_condition_id: Optional[str] = None
    omics_order_date: Optional[Any] = None
    omics_profile_flagship: Optional[str] = None
    omics_profile_funding_source: Optional[str] = None
    omics_return_date: Optional[Any] = None
    pdo_title: Optional[str] = None
    pdoid: Optional[str] = None
    pf_bases_bc: Optional[str] = None
    prioritized: Optional[bool] = None
    product: Optional[str] = None
    product_goal: Optional[str] = None
    profile_id: Optional[str] = None
    profile_source: Optional[str] = None
    project: Optional[str] = None
    proposed_release_date: Optional[Any] = None
    public_release_date: Optional[Any] = None
    public_retracted_date: Optional[Any] = None
    quote_to_bill: Optional[str] = None
    registered: Optional[bool] = None
    resubmit_for_extraction: Optional[bool] = None
    rna_delivery_date: Optional[Any] = None
    sample_coverage_normalized: Optional[str] = None
    sample_coverage_rounded: Optional[str] = None
    sample_is_on_risk: Optional[bool] = None
    sample_type: Optional[str] = None
    sm_id_matched: Optional[str] = None
    smid_ordered: Optional[str] = None
    smid_returned: Optional[str] = None
    status: Optional[str] = None
    version: Optional[str] = None
    wgs_delivery_date: Optional[Any] = None
    workspace: Optional[str] = None


class omics_profile_updates(BaseModel):
    set: Optional["omics_profile_set_input"] = Field(alias="_set", default=None)
    where: "omics_profile_bool_exp"


class omics_sequencing_aggregate_bool_exp(BaseModel):
    bool_and: Optional["omics_sequencing_aggregate_bool_exp_bool_and"] = None
    bool_or: Optional["omics_sequencing_aggregate_bool_exp_bool_or"] = None
    count: Optional["omics_sequencing_aggregate_bool_exp_count"] = None


class omics_sequencing_aggregate_bool_exp_bool_and(BaseModel):
    arguments: omics_sequencing_select_column_omics_sequencing_aggregate_bool_exp_bool_and_arguments_columns
    distinct: Optional[bool] = None
    filter: Optional["omics_sequencing_bool_exp"] = None
    predicate: "Boolean_comparison_exp"


class omics_sequencing_aggregate_bool_exp_bool_or(BaseModel):
    arguments: omics_sequencing_select_column_omics_sequencing_aggregate_bool_exp_bool_or_arguments_columns
    distinct: Optional[bool] = None
    filter: Optional["omics_sequencing_bool_exp"] = None
    predicate: "Boolean_comparison_exp"


class omics_sequencing_aggregate_bool_exp_count(BaseModel):
    arguments: Optional[List[omics_sequencing_select_column]] = None
    distinct: Optional[bool] = None
    filter: Optional["omics_sequencing_bool_exp"] = None
    predicate: "Int_comparison_exp"


class omics_sequencing_aggregate_order_by(BaseModel):
    avg: Optional["omics_sequencing_avg_order_by"] = None
    count: Optional[order_by] = None
    max: Optional["omics_sequencing_max_order_by"] = None
    min: Optional["omics_sequencing_min_order_by"] = None
    stddev: Optional["omics_sequencing_stddev_order_by"] = None
    stddev_pop: Optional["omics_sequencing_stddev_pop_order_by"] = None
    stddev_samp: Optional["omics_sequencing_stddev_samp_order_by"] = None
    sum: Optional["omics_sequencing_sum_order_by"] = None
    var_pop: Optional["omics_sequencing_var_pop_order_by"] = None
    var_samp: Optional["omics_sequencing_var_samp_order_by"] = None
    variance: Optional["omics_sequencing_variance_order_by"] = None


class omics_sequencing_arr_rel_insert_input(BaseModel):
    data: List["omics_sequencing_insert_input"]
    on_conflict: Optional["omics_sequencing_on_conflict"] = None


class omics_sequencing_avg_order_by(BaseModel):
    bam_size: Optional[order_by] = None
    hg_19_bam_size: Optional[order_by] = Field(alias="hg19_bam_size", default=None)
    hg_38_cram_size: Optional[order_by] = Field(alias="hg38_cram_size", default=None)
    month_sequencing_billed: Optional[order_by] = None
    version: Optional[order_by] = None
    year_sequencing_billed: Optional[order_by] = None


class omics_sequencing_bool_exp(BaseModel):
    and_: Optional[List["omics_sequencing_bool_exp"]] = Field(
        alias="_and", default=None
    )
    not_: Optional["omics_sequencing_bool_exp"] = Field(alias="_not", default=None)
    or_: Optional[List["omics_sequencing_bool_exp"]] = Field(alias="_or", default=None)
    bai_filepath: Optional["String_comparison_exp"] = None
    bam_crc_32_c_hash: Optional["String_comparison_exp"] = Field(
        alias="bam_crc32c_hash", default=None
    )
    bam_filepath: Optional["String_comparison_exp"] = None
    bam_qc: Optional["String_comparison_exp"] = None
    bam_size: Optional["bigint_comparison_exp"] = None
    blacklist: Optional["Boolean_comparison_exp"] = None
    expected_type: Optional["String_comparison_exp"] = None
    gp_alignment: Optional["String_comparison_exp"] = None
    hg_19_bai_filepath: Optional["String_comparison_exp"] = Field(
        alias="hg19_bai_filepath", default=None
    )
    hg_19_bam_crc_32_c_hash: Optional["String_comparison_exp"] = Field(
        alias="hg19_bam_crc32c_hash", default=None
    )
    hg_19_bam_filepath: Optional["String_comparison_exp"] = Field(
        alias="hg19_bam_filepath", default=None
    )
    hg_19_bam_size: Optional["bigint_comparison_exp"] = Field(
        alias="hg19_bam_size", default=None
    )
    hg_38_crai_filepath: Optional["String_comparison_exp"] = Field(
        alias="hg38_crai_filepath", default=None
    )
    hg_38_cram_crc_32_c_hash: Optional["String_comparison_exp"] = Field(
        alias="hg38_cram_crc32c_hash", default=None
    )
    hg_38_cram_filepath: Optional["String_comparison_exp"] = Field(
        alias="hg38_cram_filepath", default=None
    )
    hg_38_cram_size: Optional["bigint_comparison_exp"] = Field(
        alias="hg38_cram_size", default=None
    )
    issue: Optional["String_comparison_exp"] = None
    low_quality: Optional["Boolean_comparison_exp"] = None
    month_sequencing_billed: Optional["Int_comparison_exp"] = None
    omics_profile: Optional["omics_profile_bool_exp"] = None
    pdo_id: Optional["String_comparison_exp"] = None
    prioritized: Optional["Boolean_comparison_exp"] = None
    processed_sequence: Optional["Boolean_comparison_exp"] = None
    processing_qc: Optional["String_comparison_exp"] = None
    profile_id: Optional["String_comparison_exp"] = None
    sequencing_date: Optional["date_comparison_exp"] = None
    sequencing_id: Optional["String_comparison_exp"] = None
    sm_id: Optional["String_comparison_exp"] = None
    source: Optional["String_comparison_exp"] = None
    str_profile: Optional["String_comparison_exp"] = None
    stranded: Optional["Boolean_comparison_exp"] = None
    update_time: Optional["date_comparison_exp"] = None
    version: Optional["Int_comparison_exp"] = None
    year_sequencing_billed: Optional["Int_comparison_exp"] = None


class omics_sequencing_inc_input(BaseModel):
    bam_size: Optional[int] = None
    hg_19_bam_size: Optional[int] = Field(alias="hg19_bam_size", default=None)
    hg_38_cram_size: Optional[int] = Field(alias="hg38_cram_size", default=None)
    month_sequencing_billed: Optional[int] = None
    version: Optional[int] = None
    year_sequencing_billed: Optional[int] = None


class omics_sequencing_insert_input(BaseModel):
    bai_filepath: Optional[str] = None
    bam_crc_32_c_hash: Optional[str] = Field(alias="bam_crc32c_hash", default=None)
    bam_filepath: Optional[str] = None
    bam_qc: Optional[str] = None
    bam_size: Optional[int] = None
    blacklist: Optional[bool] = None
    expected_type: Optional[str] = None
    gp_alignment: Optional[str] = None
    hg_19_bai_filepath: Optional[str] = Field(alias="hg19_bai_filepath", default=None)
    hg_19_bam_crc_32_c_hash: Optional[str] = Field(
        alias="hg19_bam_crc32c_hash", default=None
    )
    hg_19_bam_filepath: Optional[str] = Field(alias="hg19_bam_filepath", default=None)
    hg_19_bam_size: Optional[int] = Field(alias="hg19_bam_size", default=None)
    hg_38_crai_filepath: Optional[str] = Field(alias="hg38_crai_filepath", default=None)
    hg_38_cram_crc_32_c_hash: Optional[str] = Field(
        alias="hg38_cram_crc32c_hash", default=None
    )
    hg_38_cram_filepath: Optional[str] = Field(alias="hg38_cram_filepath", default=None)
    hg_38_cram_size: Optional[int] = Field(alias="hg38_cram_size", default=None)
    issue: Optional[str] = None
    low_quality: Optional[bool] = None
    month_sequencing_billed: Optional[int] = None
    omics_profile: Optional["omics_profile_obj_rel_insert_input"] = None
    pdo_id: Optional[str] = None
    prioritized: Optional[bool] = None
    processed_sequence: Optional[bool] = None
    processing_qc: Optional[str] = None
    profile_id: Optional[str] = None
    sequencing_date: Optional[Any] = None
    sequencing_id: Optional[str] = None
    sm_id: Optional[str] = None
    source: Optional[str] = None
    str_profile: Optional[str] = None
    stranded: Optional[bool] = None
    update_time: Optional[Any] = None
    version: Optional[int] = None
    year_sequencing_billed: Optional[int] = None


class omics_sequencing_max_order_by(BaseModel):
    bai_filepath: Optional[order_by] = None
    bam_crc_32_c_hash: Optional[order_by] = Field(alias="bam_crc32c_hash", default=None)
    bam_filepath: Optional[order_by] = None
    bam_qc: Optional[order_by] = None
    bam_size: Optional[order_by] = None
    expected_type: Optional[order_by] = None
    gp_alignment: Optional[order_by] = None
    hg_19_bai_filepath: Optional[order_by] = Field(
        alias="hg19_bai_filepath", default=None
    )
    hg_19_bam_crc_32_c_hash: Optional[order_by] = Field(
        alias="hg19_bam_crc32c_hash", default=None
    )
    hg_19_bam_filepath: Optional[order_by] = Field(
        alias="hg19_bam_filepath", default=None
    )
    hg_19_bam_size: Optional[order_by] = Field(alias="hg19_bam_size", default=None)
    hg_38_crai_filepath: Optional[order_by] = Field(
        alias="hg38_crai_filepath", default=None
    )
    hg_38_cram_crc_32_c_hash: Optional[order_by] = Field(
        alias="hg38_cram_crc32c_hash", default=None
    )
    hg_38_cram_filepath: Optional[order_by] = Field(
        alias="hg38_cram_filepath", default=None
    )
    hg_38_cram_size: Optional[order_by] = Field(alias="hg38_cram_size", default=None)
    issue: Optional[order_by] = None
    month_sequencing_billed: Optional[order_by] = None
    pdo_id: Optional[order_by] = None
    processing_qc: Optional[order_by] = None
    profile_id: Optional[order_by] = None
    sequencing_date: Optional[order_by] = None
    sequencing_id: Optional[order_by] = None
    sm_id: Optional[order_by] = None
    source: Optional[order_by] = None
    str_profile: Optional[order_by] = None
    update_time: Optional[order_by] = None
    version: Optional[order_by] = None
    year_sequencing_billed: Optional[order_by] = None


class omics_sequencing_min_order_by(BaseModel):
    bai_filepath: Optional[order_by] = None
    bam_crc_32_c_hash: Optional[order_by] = Field(alias="bam_crc32c_hash", default=None)
    bam_filepath: Optional[order_by] = None
    bam_qc: Optional[order_by] = None
    bam_size: Optional[order_by] = None
    expected_type: Optional[order_by] = None
    gp_alignment: Optional[order_by] = None
    hg_19_bai_filepath: Optional[order_by] = Field(
        alias="hg19_bai_filepath", default=None
    )
    hg_19_bam_crc_32_c_hash: Optional[order_by] = Field(
        alias="hg19_bam_crc32c_hash", default=None
    )
    hg_19_bam_filepath: Optional[order_by] = Field(
        alias="hg19_bam_filepath", default=None
    )
    hg_19_bam_size: Optional[order_by] = Field(alias="hg19_bam_size", default=None)
    hg_38_crai_filepath: Optional[order_by] = Field(
        alias="hg38_crai_filepath", default=None
    )
    hg_38_cram_crc_32_c_hash: Optional[order_by] = Field(
        alias="hg38_cram_crc32c_hash", default=None
    )
    hg_38_cram_filepath: Optional[order_by] = Field(
        alias="hg38_cram_filepath", default=None
    )
    hg_38_cram_size: Optional[order_by] = Field(alias="hg38_cram_size", default=None)
    issue: Optional[order_by] = None
    month_sequencing_billed: Optional[order_by] = None
    pdo_id: Optional[order_by] = None
    processing_qc: Optional[order_by] = None
    profile_id: Optional[order_by] = None
    sequencing_date: Optional[order_by] = None
    sequencing_id: Optional[order_by] = None
    sm_id: Optional[order_by] = None
    source: Optional[order_by] = None
    str_profile: Optional[order_by] = None
    update_time: Optional[order_by] = None
    version: Optional[order_by] = None
    year_sequencing_billed: Optional[order_by] = None


class omics_sequencing_obj_rel_insert_input(BaseModel):
    data: "omics_sequencing_insert_input"
    on_conflict: Optional["omics_sequencing_on_conflict"] = None


class omics_sequencing_on_conflict(BaseModel):
    constraint: omics_sequencing_constraint
    update_columns: List[omics_sequencing_update_column]
    where: Optional["omics_sequencing_bool_exp"] = None


class omics_sequencing_order_by(BaseModel):
    bai_filepath: Optional[order_by] = None
    bam_crc_32_c_hash: Optional[order_by] = Field(alias="bam_crc32c_hash", default=None)
    bam_filepath: Optional[order_by] = None
    bam_qc: Optional[order_by] = None
    bam_size: Optional[order_by] = None
    blacklist: Optional[order_by] = None
    expected_type: Optional[order_by] = None
    gp_alignment: Optional[order_by] = None
    hg_19_bai_filepath: Optional[order_by] = Field(
        alias="hg19_bai_filepath", default=None
    )
    hg_19_bam_crc_32_c_hash: Optional[order_by] = Field(
        alias="hg19_bam_crc32c_hash", default=None
    )
    hg_19_bam_filepath: Optional[order_by] = Field(
        alias="hg19_bam_filepath", default=None
    )
    hg_19_bam_size: Optional[order_by] = Field(alias="hg19_bam_size", default=None)
    hg_38_crai_filepath: Optional[order_by] = Field(
        alias="hg38_crai_filepath", default=None
    )
    hg_38_cram_crc_32_c_hash: Optional[order_by] = Field(
        alias="hg38_cram_crc32c_hash", default=None
    )
    hg_38_cram_filepath: Optional[order_by] = Field(
        alias="hg38_cram_filepath", default=None
    )
    hg_38_cram_size: Optional[order_by] = Field(alias="hg38_cram_size", default=None)
    issue: Optional[order_by] = None
    low_quality: Optional[order_by] = None
    month_sequencing_billed: Optional[order_by] = None
    omics_profile: Optional["omics_profile_order_by"] = None
    pdo_id: Optional[order_by] = None
    prioritized: Optional[order_by] = None
    processed_sequence: Optional[order_by] = None
    processing_qc: Optional[order_by] = None
    profile_id: Optional[order_by] = None
    sequencing_date: Optional[order_by] = None
    sequencing_id: Optional[order_by] = None
    sm_id: Optional[order_by] = None
    source: Optional[order_by] = None
    str_profile: Optional[order_by] = None
    stranded: Optional[order_by] = None
    update_time: Optional[order_by] = None
    version: Optional[order_by] = None
    year_sequencing_billed: Optional[order_by] = None


class omics_sequencing_pk_columns_input(BaseModel):
    sequencing_id: str


class omics_sequencing_set_input(BaseModel):
    bai_filepath: Optional[str] = None
    bam_crc_32_c_hash: Optional[str] = Field(alias="bam_crc32c_hash", default=None)
    bam_filepath: Optional[str] = None
    bam_qc: Optional[str] = None
    bam_size: Optional[int] = None
    blacklist: Optional[bool] = None
    expected_type: Optional[str] = None
    gp_alignment: Optional[str] = None
    hg_19_bai_filepath: Optional[str] = Field(alias="hg19_bai_filepath", default=None)
    hg_19_bam_crc_32_c_hash: Optional[str] = Field(
        alias="hg19_bam_crc32c_hash", default=None
    )
    hg_19_bam_filepath: Optional[str] = Field(alias="hg19_bam_filepath", default=None)
    hg_19_bam_size: Optional[int] = Field(alias="hg19_bam_size", default=None)
    hg_38_crai_filepath: Optional[str] = Field(alias="hg38_crai_filepath", default=None)
    hg_38_cram_crc_32_c_hash: Optional[str] = Field(
        alias="hg38_cram_crc32c_hash", default=None
    )
    hg_38_cram_filepath: Optional[str] = Field(alias="hg38_cram_filepath", default=None)
    hg_38_cram_size: Optional[int] = Field(alias="hg38_cram_size", default=None)
    issue: Optional[str] = None
    low_quality: Optional[bool] = None
    month_sequencing_billed: Optional[int] = None
    pdo_id: Optional[str] = None
    prioritized: Optional[bool] = None
    processed_sequence: Optional[bool] = None
    processing_qc: Optional[str] = None
    profile_id: Optional[str] = None
    sequencing_date: Optional[Any] = None
    sequencing_id: Optional[str] = None
    sm_id: Optional[str] = None
    source: Optional[str] = None
    str_profile: Optional[str] = None
    stranded: Optional[bool] = None
    update_time: Optional[Any] = None
    version: Optional[int] = None
    year_sequencing_billed: Optional[int] = None


class omics_sequencing_stddev_order_by(BaseModel):
    bam_size: Optional[order_by] = None
    hg_19_bam_size: Optional[order_by] = Field(alias="hg19_bam_size", default=None)
    hg_38_cram_size: Optional[order_by] = Field(alias="hg38_cram_size", default=None)
    month_sequencing_billed: Optional[order_by] = None
    version: Optional[order_by] = None
    year_sequencing_billed: Optional[order_by] = None


class omics_sequencing_stddev_pop_order_by(BaseModel):
    bam_size: Optional[order_by] = None
    hg_19_bam_size: Optional[order_by] = Field(alias="hg19_bam_size", default=None)
    hg_38_cram_size: Optional[order_by] = Field(alias="hg38_cram_size", default=None)
    month_sequencing_billed: Optional[order_by] = None
    version: Optional[order_by] = None
    year_sequencing_billed: Optional[order_by] = None


class omics_sequencing_stddev_samp_order_by(BaseModel):
    bam_size: Optional[order_by] = None
    hg_19_bam_size: Optional[order_by] = Field(alias="hg19_bam_size", default=None)
    hg_38_cram_size: Optional[order_by] = Field(alias="hg38_cram_size", default=None)
    month_sequencing_billed: Optional[order_by] = None
    version: Optional[order_by] = None
    year_sequencing_billed: Optional[order_by] = None


class omics_sequencing_stream_cursor_input(BaseModel):
    initial_value: "omics_sequencing_stream_cursor_value_input"
    ordering: Optional[cursor_ordering] = None


class omics_sequencing_stream_cursor_value_input(BaseModel):
    bai_filepath: Optional[str] = None
    bam_crc_32_c_hash: Optional[str] = Field(alias="bam_crc32c_hash", default=None)
    bam_filepath: Optional[str] = None
    bam_qc: Optional[str] = None
    bam_size: Optional[int] = None
    blacklist: Optional[bool] = None
    expected_type: Optional[str] = None
    gp_alignment: Optional[str] = None
    hg_19_bai_filepath: Optional[str] = Field(alias="hg19_bai_filepath", default=None)
    hg_19_bam_crc_32_c_hash: Optional[str] = Field(
        alias="hg19_bam_crc32c_hash", default=None
    )
    hg_19_bam_filepath: Optional[str] = Field(alias="hg19_bam_filepath", default=None)
    hg_19_bam_size: Optional[int] = Field(alias="hg19_bam_size", default=None)
    hg_38_crai_filepath: Optional[str] = Field(alias="hg38_crai_filepath", default=None)
    hg_38_cram_crc_32_c_hash: Optional[str] = Field(
        alias="hg38_cram_crc32c_hash", default=None
    )
    hg_38_cram_filepath: Optional[str] = Field(alias="hg38_cram_filepath", default=None)
    hg_38_cram_size: Optional[int] = Field(alias="hg38_cram_size", default=None)
    issue: Optional[str] = None
    low_quality: Optional[bool] = None
    month_sequencing_billed: Optional[int] = None
    pdo_id: Optional[str] = None
    prioritized: Optional[bool] = None
    processed_sequence: Optional[bool] = None
    processing_qc: Optional[str] = None
    profile_id: Optional[str] = None
    sequencing_date: Optional[Any] = None
    sequencing_id: Optional[str] = None
    sm_id: Optional[str] = None
    source: Optional[str] = None
    str_profile: Optional[str] = None
    stranded: Optional[bool] = None
    update_time: Optional[Any] = None
    version: Optional[int] = None
    year_sequencing_billed: Optional[int] = None


class omics_sequencing_sum_order_by(BaseModel):
    bam_size: Optional[order_by] = None
    hg_19_bam_size: Optional[order_by] = Field(alias="hg19_bam_size", default=None)
    hg_38_cram_size: Optional[order_by] = Field(alias="hg38_cram_size", default=None)
    month_sequencing_billed: Optional[order_by] = None
    version: Optional[order_by] = None
    year_sequencing_billed: Optional[order_by] = None


class omics_sequencing_updates(BaseModel):
    inc: Optional["omics_sequencing_inc_input"] = Field(alias="_inc", default=None)
    set: Optional["omics_sequencing_set_input"] = Field(alias="_set", default=None)
    where: "omics_sequencing_bool_exp"


class omics_sequencing_var_pop_order_by(BaseModel):
    bam_size: Optional[order_by] = None
    hg_19_bam_size: Optional[order_by] = Field(alias="hg19_bam_size", default=None)
    hg_38_cram_size: Optional[order_by] = Field(alias="hg38_cram_size", default=None)
    month_sequencing_billed: Optional[order_by] = None
    version: Optional[order_by] = None
    year_sequencing_billed: Optional[order_by] = None


class omics_sequencing_var_samp_order_by(BaseModel):
    bam_size: Optional[order_by] = None
    hg_19_bam_size: Optional[order_by] = Field(alias="hg19_bam_size", default=None)
    hg_38_cram_size: Optional[order_by] = Field(alias="hg38_cram_size", default=None)
    month_sequencing_billed: Optional[order_by] = None
    version: Optional[order_by] = None
    year_sequencing_billed: Optional[order_by] = None


class omics_sequencing_variance_order_by(BaseModel):
    bam_size: Optional[order_by] = None
    hg_19_bam_size: Optional[order_by] = Field(alias="hg19_bam_size", default=None)
    hg_38_cram_size: Optional[order_by] = Field(alias="hg38_cram_size", default=None)
    month_sequencing_billed: Optional[order_by] = None
    version: Optional[order_by] = None
    year_sequencing_billed: Optional[order_by] = None


class patient_bool_exp(BaseModel):
    and_: Optional[List["patient_bool_exp"]] = Field(alias="_and", default=None)
    not_: Optional["patient_bool_exp"] = Field(alias="_not", default=None)
    or_: Optional[List["patient_bool_exp"]] = Field(alias="_or", default=None)
    patient_id: Optional["String_comparison_exp"] = None


class patient_insert_input(BaseModel):
    patient_id: Optional[str] = None


class patient_on_conflict(BaseModel):
    constraint: patient_constraint
    update_columns: List[patient_update_column]
    where: Optional["patient_bool_exp"] = None


class patient_order_by(BaseModel):
    patient_id: Optional[order_by] = None


class patient_pk_columns_input(BaseModel):
    patient_id: str


class patient_set_input(BaseModel):
    patient_id: Optional[str] = None


class patient_stream_cursor_input(BaseModel):
    initial_value: "patient_stream_cursor_value_input"
    ordering: Optional[cursor_ordering] = None


class patient_stream_cursor_value_input(BaseModel):
    patient_id: Optional[str] = None


class patient_updates(BaseModel):
    set: Optional["patient_set_input"] = Field(alias="_set", default=None)
    where: "patient_bool_exp"


class set_username_args(BaseModel):
    username: Optional[str] = Field(alias="_username", default=None)


class smallint_comparison_exp(BaseModel):
    eq: Optional[Any] = Field(alias="_eq", default=None)
    gt: Optional[Any] = Field(alias="_gt", default=None)
    gte: Optional[Any] = Field(alias="_gte", default=None)
    in_: Optional[List[Any]] = Field(alias="_in", default=None)
    is_null: Optional[bool] = Field(alias="_is_null", default=None)
    lt: Optional[Any] = Field(alias="_lt", default=None)
    lte: Optional[Any] = Field(alias="_lte", default=None)
    neq: Optional[Any] = Field(alias="_neq", default=None)
    nin: Optional[List[Any]] = Field(alias="_nin", default=None)


class snp_fingerprint_bool_exp(BaseModel):
    and_: Optional[List["snp_fingerprint_bool_exp"]] = Field(alias="_and", default=None)
    not_: Optional["snp_fingerprint_bool_exp"] = Field(alias="_not", default=None)
    or_: Optional[List["snp_fingerprint_bool_exp"]] = Field(alias="_or", default=None)
    comments: Optional["String_comparison_exp"] = None
    created_at: Optional["timestamptz_comparison_exp"] = None
    genotypes: Optional["String_comparison_exp"] = None
    id: Optional["bigint_comparison_exp"] = None
    omics_sequencing_id: Optional["String_comparison_exp"] = None
    updated_at: Optional["timestamptz_comparison_exp"] = None
    vcf_uri: Optional["String_comparison_exp"] = None


class snp_fingerprint_comparison_bool_exp(BaseModel):
    and_: Optional[List["snp_fingerprint_comparison_bool_exp"]] = Field(
        alias="_and", default=None
    )
    not_: Optional["snp_fingerprint_comparison_bool_exp"] = Field(
        alias="_not", default=None
    )
    or_: Optional[List["snp_fingerprint_comparison_bool_exp"]] = Field(
        alias="_or", default=None
    )
    acknowledged: Optional["Boolean_comparison_exp"] = None
    comments: Optional["String_comparison_exp"] = None
    created_at: Optional["timestamptz_comparison_exp"] = None
    id: Optional["bigint_comparison_exp"] = None
    issue: Optional["String_comparison_exp"] = None
    n_common_snps: Optional["smallint_comparison_exp"] = None
    n_matching_genotypes: Optional["smallint_comparison_exp"] = None
    patient_id_1: Optional["String_comparison_exp"] = Field(
        alias="patient_id1", default=None
    )
    patient_id_2: Optional["String_comparison_exp"] = Field(
        alias="patient_id2", default=None
    )
    score: Optional["float8_comparison_exp"] = None
    snp_fingerprint_id_1: Optional["bigint_comparison_exp"] = Field(
        alias="snp_fingerprint_id1", default=None
    )
    snp_fingerprint_id_2: Optional["bigint_comparison_exp"] = Field(
        alias="snp_fingerprint_id2", default=None
    )
    updated_at: Optional["timestamptz_comparison_exp"] = None


class snp_fingerprint_comparison_inc_input(BaseModel):
    id: Optional[int] = None
    n_common_snps: Optional[Any] = None
    n_matching_genotypes: Optional[Any] = None
    score: Optional[Any] = None
    snp_fingerprint_id_1: Optional[int] = Field(
        alias="snp_fingerprint_id1", default=None
    )
    snp_fingerprint_id_2: Optional[int] = Field(
        alias="snp_fingerprint_id2", default=None
    )


class snp_fingerprint_comparison_insert_input(BaseModel):
    acknowledged: Optional[bool] = None
    comments: Optional[str] = None
    created_at: Optional[Any] = None
    id: Optional[int] = None
    issue: Optional[str] = None
    n_common_snps: Optional[Any] = None
    n_matching_genotypes: Optional[Any] = None
    patient_id_1: Optional[str] = Field(alias="patient_id1", default=None)
    patient_id_2: Optional[str] = Field(alias="patient_id2", default=None)
    score: Optional[Any] = None
    snp_fingerprint_id_1: Optional[int] = Field(
        alias="snp_fingerprint_id1", default=None
    )
    snp_fingerprint_id_2: Optional[int] = Field(
        alias="snp_fingerprint_id2", default=None
    )
    updated_at: Optional[Any] = None


class snp_fingerprint_comparison_on_conflict(BaseModel):
    constraint: snp_fingerprint_comparison_constraint
    update_columns: List[snp_fingerprint_comparison_update_column]
    where: Optional["snp_fingerprint_comparison_bool_exp"] = None


class snp_fingerprint_comparison_order_by(BaseModel):
    acknowledged: Optional[order_by] = None
    comments: Optional[order_by] = None
    created_at: Optional[order_by] = None
    id: Optional[order_by] = None
    issue: Optional[order_by] = None
    n_common_snps: Optional[order_by] = None
    n_matching_genotypes: Optional[order_by] = None
    patient_id_1: Optional[order_by] = Field(alias="patient_id1", default=None)
    patient_id_2: Optional[order_by] = Field(alias="patient_id2", default=None)
    score: Optional[order_by] = None
    snp_fingerprint_id_1: Optional[order_by] = Field(
        alias="snp_fingerprint_id1", default=None
    )
    snp_fingerprint_id_2: Optional[order_by] = Field(
        alias="snp_fingerprint_id2", default=None
    )
    updated_at: Optional[order_by] = None


class snp_fingerprint_comparison_pk_columns_input(BaseModel):
    id: int


class snp_fingerprint_comparison_set_input(BaseModel):
    acknowledged: Optional[bool] = None
    comments: Optional[str] = None
    created_at: Optional[Any] = None
    id: Optional[int] = None
    issue: Optional[str] = None
    n_common_snps: Optional[Any] = None
    n_matching_genotypes: Optional[Any] = None
    patient_id_1: Optional[str] = Field(alias="patient_id1", default=None)
    patient_id_2: Optional[str] = Field(alias="patient_id2", default=None)
    score: Optional[Any] = None
    snp_fingerprint_id_1: Optional[int] = Field(
        alias="snp_fingerprint_id1", default=None
    )
    snp_fingerprint_id_2: Optional[int] = Field(
        alias="snp_fingerprint_id2", default=None
    )
    updated_at: Optional[Any] = None


class snp_fingerprint_comparison_stream_cursor_input(BaseModel):
    initial_value: "snp_fingerprint_comparison_stream_cursor_value_input"
    ordering: Optional[cursor_ordering] = None


class snp_fingerprint_comparison_stream_cursor_value_input(BaseModel):
    acknowledged: Optional[bool] = None
    comments: Optional[str] = None
    created_at: Optional[Any] = None
    id: Optional[int] = None
    issue: Optional[str] = None
    n_common_snps: Optional[Any] = None
    n_matching_genotypes: Optional[Any] = None
    patient_id_1: Optional[str] = Field(alias="patient_id1", default=None)
    patient_id_2: Optional[str] = Field(alias="patient_id2", default=None)
    score: Optional[Any] = None
    snp_fingerprint_id_1: Optional[int] = Field(
        alias="snp_fingerprint_id1", default=None
    )
    snp_fingerprint_id_2: Optional[int] = Field(
        alias="snp_fingerprint_id2", default=None
    )
    updated_at: Optional[Any] = None


class snp_fingerprint_comparison_updates(BaseModel):
    inc: Optional["snp_fingerprint_comparison_inc_input"] = Field(
        alias="_inc", default=None
    )
    set: Optional["snp_fingerprint_comparison_set_input"] = Field(
        alias="_set", default=None
    )
    where: "snp_fingerprint_comparison_bool_exp"


class snp_fingerprint_inc_input(BaseModel):
    id: Optional[int] = None


class snp_fingerprint_insert_input(BaseModel):
    comments: Optional[str] = None
    created_at: Optional[Any] = None
    genotypes: Optional[str] = None
    id: Optional[int] = None
    omics_sequencing_id: Optional[str] = None
    updated_at: Optional[Any] = None
    vcf_uri: Optional[str] = None


class snp_fingerprint_on_conflict(BaseModel):
    constraint: snp_fingerprint_constraint
    update_columns: List[snp_fingerprint_update_column]
    where: Optional["snp_fingerprint_bool_exp"] = None


class snp_fingerprint_order_by(BaseModel):
    comments: Optional[order_by] = None
    created_at: Optional[order_by] = None
    genotypes: Optional[order_by] = None
    id: Optional[order_by] = None
    omics_sequencing_id: Optional[order_by] = None
    updated_at: Optional[order_by] = None
    vcf_uri: Optional[order_by] = None


class snp_fingerprint_pk_columns_input(BaseModel):
    id: int


class snp_fingerprint_qc_bool_exp(BaseModel):
    and_: Optional[List["snp_fingerprint_qc_bool_exp"]] = Field(
        alias="_and", default=None
    )
    not_: Optional["snp_fingerprint_qc_bool_exp"] = Field(alias="_not", default=None)
    or_: Optional[List["snp_fingerprint_qc_bool_exp"]] = Field(
        alias="_or", default=None
    )
    acknowledged: Optional["Boolean_comparison_exp"] = None
    cell_line_name_1: Optional["String_comparison_exp"] = Field(
        alias="cell_line_name1", default=None
    )
    cell_line_name_2: Optional["String_comparison_exp"] = Field(
        alias="cell_line_name2", default=None
    )
    comments: Optional["String_comparison_exp"] = None
    compared_patient_id_1: Optional["String_comparison_exp"] = Field(
        alias="compared_patient_id1", default=None
    )
    compared_patient_id_2: Optional["String_comparison_exp"] = Field(
        alias="compared_patient_id2", default=None
    )
    created_at: Optional["timestamptz_comparison_exp"] = None
    expected_type_1: Optional["String_comparison_exp"] = Field(
        alias="expected_type1", default=None
    )
    expected_type_2: Optional["String_comparison_exp"] = Field(
        alias="expected_type2", default=None
    )
    id: Optional["bigint_comparison_exp"] = None
    issue: Optional["String_comparison_exp"] = None
    linked_patient_id_1: Optional["String_comparison_exp"] = Field(
        alias="linked_patient_id1", default=None
    )
    linked_patient_id_2: Optional["String_comparison_exp"] = Field(
        alias="linked_patient_id2", default=None
    )
    low_quality_1: Optional["Boolean_comparison_exp"] = Field(
        alias="low_quality1", default=None
    )
    low_quality_2: Optional["Boolean_comparison_exp"] = Field(
        alias="low_quality2", default=None
    )
    model_condition_1: Optional["String_comparison_exp"] = Field(
        alias="model_condition1", default=None
    )
    model_condition_2: Optional["String_comparison_exp"] = Field(
        alias="model_condition2", default=None
    )
    model_id_1: Optional["String_comparison_exp"] = Field(
        alias="model_id1", default=None
    )
    model_id_2: Optional["String_comparison_exp"] = Field(
        alias="model_id2", default=None
    )
    omics_sequencing_blacklist_1: Optional["Boolean_comparison_exp"] = Field(
        alias="omics_sequencing_blacklist1", default=None
    )
    omics_sequencing_blacklist_2: Optional["Boolean_comparison_exp"] = Field(
        alias="omics_sequencing_blacklist2", default=None
    )
    profile_blacklist_omics_1: Optional["Boolean_comparison_exp"] = Field(
        alias="profile_blacklist_omics1", default=None
    )
    profile_blacklist_omics_2: Optional["Boolean_comparison_exp"] = Field(
        alias="profile_blacklist_omics2", default=None
    )
    profile_id_1: Optional["String_comparison_exp"] = Field(
        alias="profile_id1", default=None
    )
    profile_id_2: Optional["String_comparison_exp"] = Field(
        alias="profile_id2", default=None
    )
    score: Optional["float8_comparison_exp"] = None
    sequencing_id_1: Optional["String_comparison_exp"] = Field(
        alias="sequencing_id1", default=None
    )
    sequencing_id_2: Optional["String_comparison_exp"] = Field(
        alias="sequencing_id2", default=None
    )
    snp_fingerprint_comparison_id: Optional["bigint_comparison_exp"] = None


class snp_fingerprint_qc_order_by(BaseModel):
    acknowledged: Optional[order_by] = None
    cell_line_name_1: Optional[order_by] = Field(alias="cell_line_name1", default=None)
    cell_line_name_2: Optional[order_by] = Field(alias="cell_line_name2", default=None)
    comments: Optional[order_by] = None
    compared_patient_id_1: Optional[order_by] = Field(
        alias="compared_patient_id1", default=None
    )
    compared_patient_id_2: Optional[order_by] = Field(
        alias="compared_patient_id2", default=None
    )
    created_at: Optional[order_by] = None
    expected_type_1: Optional[order_by] = Field(alias="expected_type1", default=None)
    expected_type_2: Optional[order_by] = Field(alias="expected_type2", default=None)
    id: Optional[order_by] = None
    issue: Optional[order_by] = None
    linked_patient_id_1: Optional[order_by] = Field(
        alias="linked_patient_id1", default=None
    )
    linked_patient_id_2: Optional[order_by] = Field(
        alias="linked_patient_id2", default=None
    )
    low_quality_1: Optional[order_by] = Field(alias="low_quality1", default=None)
    low_quality_2: Optional[order_by] = Field(alias="low_quality2", default=None)
    model_condition_1: Optional[order_by] = Field(
        alias="model_condition1", default=None
    )
    model_condition_2: Optional[order_by] = Field(
        alias="model_condition2", default=None
    )
    model_id_1: Optional[order_by] = Field(alias="model_id1", default=None)
    model_id_2: Optional[order_by] = Field(alias="model_id2", default=None)
    omics_sequencing_blacklist_1: Optional[order_by] = Field(
        alias="omics_sequencing_blacklist1", default=None
    )
    omics_sequencing_blacklist_2: Optional[order_by] = Field(
        alias="omics_sequencing_blacklist2", default=None
    )
    profile_blacklist_omics_1: Optional[order_by] = Field(
        alias="profile_blacklist_omics1", default=None
    )
    profile_blacklist_omics_2: Optional[order_by] = Field(
        alias="profile_blacklist_omics2", default=None
    )
    profile_id_1: Optional[order_by] = Field(alias="profile_id1", default=None)
    profile_id_2: Optional[order_by] = Field(alias="profile_id2", default=None)
    score: Optional[order_by] = None
    sequencing_id_1: Optional[order_by] = Field(alias="sequencing_id1", default=None)
    sequencing_id_2: Optional[order_by] = Field(alias="sequencing_id2", default=None)
    snp_fingerprint_comparison_id: Optional[order_by] = None


class snp_fingerprint_qc_stream_cursor_input(BaseModel):
    initial_value: "snp_fingerprint_qc_stream_cursor_value_input"
    ordering: Optional[cursor_ordering] = None


class snp_fingerprint_qc_stream_cursor_value_input(BaseModel):
    acknowledged: Optional[bool] = None
    cell_line_name_1: Optional[str] = Field(alias="cell_line_name1", default=None)
    cell_line_name_2: Optional[str] = Field(alias="cell_line_name2", default=None)
    comments: Optional[str] = None
    compared_patient_id_1: Optional[str] = Field(
        alias="compared_patient_id1", default=None
    )
    compared_patient_id_2: Optional[str] = Field(
        alias="compared_patient_id2", default=None
    )
    created_at: Optional[Any] = None
    expected_type_1: Optional[str] = Field(alias="expected_type1", default=None)
    expected_type_2: Optional[str] = Field(alias="expected_type2", default=None)
    id: Optional[int] = None
    issue: Optional[str] = None
    linked_patient_id_1: Optional[str] = Field(alias="linked_patient_id1", default=None)
    linked_patient_id_2: Optional[str] = Field(alias="linked_patient_id2", default=None)
    low_quality_1: Optional[bool] = Field(alias="low_quality1", default=None)
    low_quality_2: Optional[bool] = Field(alias="low_quality2", default=None)
    model_condition_1: Optional[str] = Field(alias="model_condition1", default=None)
    model_condition_2: Optional[str] = Field(alias="model_condition2", default=None)
    model_id_1: Optional[str] = Field(alias="model_id1", default=None)
    model_id_2: Optional[str] = Field(alias="model_id2", default=None)
    omics_sequencing_blacklist_1: Optional[bool] = Field(
        alias="omics_sequencing_blacklist1", default=None
    )
    omics_sequencing_blacklist_2: Optional[bool] = Field(
        alias="omics_sequencing_blacklist2", default=None
    )
    profile_blacklist_omics_1: Optional[bool] = Field(
        alias="profile_blacklist_omics1", default=None
    )
    profile_blacklist_omics_2: Optional[bool] = Field(
        alias="profile_blacklist_omics2", default=None
    )
    profile_id_1: Optional[str] = Field(alias="profile_id1", default=None)
    profile_id_2: Optional[str] = Field(alias="profile_id2", default=None)
    score: Optional[Any] = None
    sequencing_id_1: Optional[str] = Field(alias="sequencing_id1", default=None)
    sequencing_id_2: Optional[str] = Field(alias="sequencing_id2", default=None)
    snp_fingerprint_comparison_id: Optional[int] = None


class snp_fingerprint_set_input(BaseModel):
    comments: Optional[str] = None
    created_at: Optional[Any] = None
    genotypes: Optional[str] = None
    id: Optional[int] = None
    omics_sequencing_id: Optional[str] = None
    updated_at: Optional[Any] = None
    vcf_uri: Optional[str] = None


class snp_fingerprint_stream_cursor_input(BaseModel):
    initial_value: "snp_fingerprint_stream_cursor_value_input"
    ordering: Optional[cursor_ordering] = None


class snp_fingerprint_stream_cursor_value_input(BaseModel):
    comments: Optional[str] = None
    created_at: Optional[Any] = None
    genotypes: Optional[str] = None
    id: Optional[int] = None
    omics_sequencing_id: Optional[str] = None
    updated_at: Optional[Any] = None
    vcf_uri: Optional[str] = None


class snp_fingerprint_updates(BaseModel):
    inc: Optional["snp_fingerprint_inc_input"] = Field(alias="_inc", default=None)
    set: Optional["snp_fingerprint_set_input"] = Field(alias="_set", default=None)
    where: "snp_fingerprint_bool_exp"


class str_profile_bool_exp(BaseModel):
    and_: Optional[List["str_profile_bool_exp"]] = Field(alias="_and", default=None)
    not_: Optional["str_profile_bool_exp"] = Field(alias="_not", default=None)
    or_: Optional[List["str_profile_bool_exp"]] = Field(alias="_or", default=None)
    amelogenin: Optional["String_comparison_exp"] = None
    comments: Optional["String_comparison_exp"] = None
    created_at: Optional["timestamptz_comparison_exp"] = None
    csf_1_po: Optional["String_comparison_exp"] = Field(alias="csf1po", default=None)
    d_13_s_317: Optional["String_comparison_exp"] = Field(alias="d13s317", default=None)
    d_16_s_539: Optional["String_comparison_exp"] = Field(alias="d16s539", default=None)
    d_18_s_51: Optional["String_comparison_exp"] = Field(alias="d18s51", default=None)
    d_21_s_11: Optional["String_comparison_exp"] = Field(alias="d21s11", default=None)
    d_3_s_1358: Optional["String_comparison_exp"] = Field(alias="d3s1358", default=None)
    d_5_s_818: Optional["String_comparison_exp"] = Field(alias="d5s818", default=None)
    d_7_s_820: Optional["String_comparison_exp"] = Field(alias="d7s820", default=None)
    d_8_s_1179: Optional["String_comparison_exp"] = Field(alias="d8s1179", default=None)
    fga: Optional["String_comparison_exp"] = None
    id: Optional["String_comparison_exp"] = None
    is_reference: Optional["Boolean_comparison_exp"] = None
    lab_corp_case_nbr: Optional["String_comparison_exp"] = None
    lab_corp_spec_nbr: Optional["String_comparison_exp"] = None
    model_condition_id: Optional["String_comparison_exp"] = None
    mouse: Optional["String_comparison_exp"] = None
    mycoplasma: Optional["String_comparison_exp"] = None
    patient_id: Optional["String_comparison_exp"] = None
    pellet_creation_date: Optional["date_comparison_exp"] = None
    pellet_submitted_date: Optional["date_comparison_exp"] = None
    penta_d: Optional["String_comparison_exp"] = None
    penta_e: Optional["String_comparison_exp"] = None
    percentage_match_to_parental: Optional["float8_comparison_exp"] = None
    sample_reference: Optional["String_comparison_exp"] = None
    source: Optional["String_comparison_exp"] = None
    source_group: Optional["String_comparison_exp"] = None
    th_01: Optional["String_comparison_exp"] = Field(alias="th01", default=None)
    tpox: Optional["String_comparison_exp"] = None
    updated_at: Optional["timestamptz_comparison_exp"] = None
    vwa: Optional["String_comparison_exp"] = None


class str_profile_inc_input(BaseModel):
    percentage_match_to_parental: Optional[Any] = None


class str_profile_insert_input(BaseModel):
    amelogenin: Optional[str] = None
    comments: Optional[str] = None
    created_at: Optional[Any] = None
    csf_1_po: Optional[str] = Field(alias="csf1po", default=None)
    d_13_s_317: Optional[str] = Field(alias="d13s317", default=None)
    d_16_s_539: Optional[str] = Field(alias="d16s539", default=None)
    d_18_s_51: Optional[str] = Field(alias="d18s51", default=None)
    d_21_s_11: Optional[str] = Field(alias="d21s11", default=None)
    d_3_s_1358: Optional[str] = Field(alias="d3s1358", default=None)
    d_5_s_818: Optional[str] = Field(alias="d5s818", default=None)
    d_7_s_820: Optional[str] = Field(alias="d7s820", default=None)
    d_8_s_1179: Optional[str] = Field(alias="d8s1179", default=None)
    fga: Optional[str] = None
    id: Optional[str] = None
    is_reference: Optional[bool] = None
    lab_corp_case_nbr: Optional[str] = None
    lab_corp_spec_nbr: Optional[str] = None
    model_condition_id: Optional[str] = None
    mouse: Optional[str] = None
    mycoplasma: Optional[str] = None
    patient_id: Optional[str] = None
    pellet_creation_date: Optional[Any] = None
    pellet_submitted_date: Optional[Any] = None
    penta_d: Optional[str] = None
    penta_e: Optional[str] = None
    percentage_match_to_parental: Optional[Any] = None
    sample_reference: Optional[str] = None
    source: Optional[str] = None
    source_group: Optional[str] = None
    th_01: Optional[str] = Field(alias="th01", default=None)
    tpox: Optional[str] = None
    updated_at: Optional[Any] = None
    vwa: Optional[str] = None


class str_profile_on_conflict(BaseModel):
    constraint: str_profile_constraint
    update_columns: List[str_profile_update_column]
    where: Optional["str_profile_bool_exp"] = None


class str_profile_order_by(BaseModel):
    amelogenin: Optional[order_by] = None
    comments: Optional[order_by] = None
    created_at: Optional[order_by] = None
    csf_1_po: Optional[order_by] = Field(alias="csf1po", default=None)
    d_13_s_317: Optional[order_by] = Field(alias="d13s317", default=None)
    d_16_s_539: Optional[order_by] = Field(alias="d16s539", default=None)
    d_18_s_51: Optional[order_by] = Field(alias="d18s51", default=None)
    d_21_s_11: Optional[order_by] = Field(alias="d21s11", default=None)
    d_3_s_1358: Optional[order_by] = Field(alias="d3s1358", default=None)
    d_5_s_818: Optional[order_by] = Field(alias="d5s818", default=None)
    d_7_s_820: Optional[order_by] = Field(alias="d7s820", default=None)
    d_8_s_1179: Optional[order_by] = Field(alias="d8s1179", default=None)
    fga: Optional[order_by] = None
    id: Optional[order_by] = None
    is_reference: Optional[order_by] = None
    lab_corp_case_nbr: Optional[order_by] = None
    lab_corp_spec_nbr: Optional[order_by] = None
    model_condition_id: Optional[order_by] = None
    mouse: Optional[order_by] = None
    mycoplasma: Optional[order_by] = None
    patient_id: Optional[order_by] = None
    pellet_creation_date: Optional[order_by] = None
    pellet_submitted_date: Optional[order_by] = None
    penta_d: Optional[order_by] = None
    penta_e: Optional[order_by] = None
    percentage_match_to_parental: Optional[order_by] = None
    sample_reference: Optional[order_by] = None
    source: Optional[order_by] = None
    source_group: Optional[order_by] = None
    th_01: Optional[order_by] = Field(alias="th01", default=None)
    tpox: Optional[order_by] = None
    updated_at: Optional[order_by] = None
    vwa: Optional[order_by] = None


class str_profile_pk_columns_input(BaseModel):
    id: str


class str_profile_set_input(BaseModel):
    amelogenin: Optional[str] = None
    comments: Optional[str] = None
    created_at: Optional[Any] = None
    csf_1_po: Optional[str] = Field(alias="csf1po", default=None)
    d_13_s_317: Optional[str] = Field(alias="d13s317", default=None)
    d_16_s_539: Optional[str] = Field(alias="d16s539", default=None)
    d_18_s_51: Optional[str] = Field(alias="d18s51", default=None)
    d_21_s_11: Optional[str] = Field(alias="d21s11", default=None)
    d_3_s_1358: Optional[str] = Field(alias="d3s1358", default=None)
    d_5_s_818: Optional[str] = Field(alias="d5s818", default=None)
    d_7_s_820: Optional[str] = Field(alias="d7s820", default=None)
    d_8_s_1179: Optional[str] = Field(alias="d8s1179", default=None)
    fga: Optional[str] = None
    id: Optional[str] = None
    is_reference: Optional[bool] = None
    lab_corp_case_nbr: Optional[str] = None
    lab_corp_spec_nbr: Optional[str] = None
    model_condition_id: Optional[str] = None
    mouse: Optional[str] = None
    mycoplasma: Optional[str] = None
    patient_id: Optional[str] = None
    pellet_creation_date: Optional[Any] = None
    pellet_submitted_date: Optional[Any] = None
    penta_d: Optional[str] = None
    penta_e: Optional[str] = None
    percentage_match_to_parental: Optional[Any] = None
    sample_reference: Optional[str] = None
    source: Optional[str] = None
    source_group: Optional[str] = None
    th_01: Optional[str] = Field(alias="th01", default=None)
    tpox: Optional[str] = None
    updated_at: Optional[Any] = None
    vwa: Optional[str] = None


class str_profile_stream_cursor_input(BaseModel):
    initial_value: "str_profile_stream_cursor_value_input"
    ordering: Optional[cursor_ordering] = None


class str_profile_stream_cursor_value_input(BaseModel):
    amelogenin: Optional[str] = None
    comments: Optional[str] = None
    created_at: Optional[Any] = None
    csf_1_po: Optional[str] = Field(alias="csf1po", default=None)
    d_13_s_317: Optional[str] = Field(alias="d13s317", default=None)
    d_16_s_539: Optional[str] = Field(alias="d16s539", default=None)
    d_18_s_51: Optional[str] = Field(alias="d18s51", default=None)
    d_21_s_11: Optional[str] = Field(alias="d21s11", default=None)
    d_3_s_1358: Optional[str] = Field(alias="d3s1358", default=None)
    d_5_s_818: Optional[str] = Field(alias="d5s818", default=None)
    d_7_s_820: Optional[str] = Field(alias="d7s820", default=None)
    d_8_s_1179: Optional[str] = Field(alias="d8s1179", default=None)
    fga: Optional[str] = None
    id: Optional[str] = None
    is_reference: Optional[bool] = None
    lab_corp_case_nbr: Optional[str] = None
    lab_corp_spec_nbr: Optional[str] = None
    model_condition_id: Optional[str] = None
    mouse: Optional[str] = None
    mycoplasma: Optional[str] = None
    patient_id: Optional[str] = None
    pellet_creation_date: Optional[Any] = None
    pellet_submitted_date: Optional[Any] = None
    penta_d: Optional[str] = None
    penta_e: Optional[str] = None
    percentage_match_to_parental: Optional[Any] = None
    sample_reference: Optional[str] = None
    source: Optional[str] = None
    source_group: Optional[str] = None
    th_01: Optional[str] = Field(alias="th01", default=None)
    tpox: Optional[str] = None
    updated_at: Optional[Any] = None
    vwa: Optional[str] = None


class str_profile_updates(BaseModel):
    inc: Optional["str_profile_inc_input"] = Field(alias="_inc", default=None)
    set: Optional["str_profile_set_input"] = Field(alias="_set", default=None)
    where: "str_profile_bool_exp"


class task_entity_bool_exp(BaseModel):
    and_: Optional[List["task_entity_bool_exp"]] = Field(alias="_and", default=None)
    not_: Optional["task_entity_bool_exp"] = Field(alias="_not", default=None)
    or_: Optional[List["task_entity_bool_exp"]] = Field(alias="_or", default=None)
    id: Optional["bigint_comparison_exp"] = None
    omics_sequencing: Optional["omics_sequencing_bool_exp"] = None
    sequencing_id: Optional["String_comparison_exp"] = None
    task_results: Optional["task_result_bool_exp"] = None
    task_results_aggregate: Optional["task_result_aggregate_bool_exp"] = None


class task_entity_inc_input(BaseModel):
    id: Optional[int] = None


class task_entity_insert_input(BaseModel):
    id: Optional[int] = None
    omics_sequencing: Optional["omics_sequencing_obj_rel_insert_input"] = None
    sequencing_id: Optional[str] = None
    task_results: Optional["task_result_arr_rel_insert_input"] = None


class task_entity_obj_rel_insert_input(BaseModel):
    data: "task_entity_insert_input"
    on_conflict: Optional["task_entity_on_conflict"] = None


class task_entity_on_conflict(BaseModel):
    constraint: task_entity_constraint
    update_columns: List[task_entity_update_column]
    where: Optional["task_entity_bool_exp"] = None


class task_entity_order_by(BaseModel):
    id: Optional[order_by] = None
    omics_sequencing: Optional["omics_sequencing_order_by"] = None
    sequencing_id: Optional[order_by] = None
    task_results_aggregate: Optional["task_result_aggregate_order_by"] = None


class task_entity_pk_columns_input(BaseModel):
    id: int


class task_entity_set_input(BaseModel):
    id: Optional[int] = None
    sequencing_id: Optional[str] = None


class task_entity_stream_cursor_input(BaseModel):
    initial_value: "task_entity_stream_cursor_value_input"
    ordering: Optional[cursor_ordering] = None


class task_entity_stream_cursor_value_input(BaseModel):
    id: Optional[int] = None
    sequencing_id: Optional[str] = None


class task_entity_updates(BaseModel):
    inc: Optional["task_entity_inc_input"] = Field(alias="_inc", default=None)
    set: Optional["task_entity_set_input"] = Field(alias="_set", default=None)
    where: "task_entity_bool_exp"


class task_result_aggregate_bool_exp(BaseModel):
    count: Optional["task_result_aggregate_bool_exp_count"] = None


class task_result_aggregate_bool_exp_count(BaseModel):
    arguments: Optional[List[task_result_select_column]] = None
    distinct: Optional[bool] = None
    filter: Optional["task_result_bool_exp"] = None
    predicate: "Int_comparison_exp"


class task_result_aggregate_order_by(BaseModel):
    avg: Optional["task_result_avg_order_by"] = None
    count: Optional[order_by] = None
    max: Optional["task_result_max_order_by"] = None
    min: Optional["task_result_min_order_by"] = None
    stddev: Optional["task_result_stddev_order_by"] = None
    stddev_pop: Optional["task_result_stddev_pop_order_by"] = None
    stddev_samp: Optional["task_result_stddev_samp_order_by"] = None
    sum: Optional["task_result_sum_order_by"] = None
    var_pop: Optional["task_result_var_pop_order_by"] = None
    var_samp: Optional["task_result_var_samp_order_by"] = None
    variance: Optional["task_result_variance_order_by"] = None


class task_result_append_input(BaseModel):
    terra_workflow_inputs: Optional[Any] = None
    value: Optional[Any] = None


class task_result_arr_rel_insert_input(BaseModel):
    data: List["task_result_insert_input"]
    on_conflict: Optional["task_result_on_conflict"] = None


class task_result_avg_order_by(BaseModel):
    size: Optional[order_by] = None
    task_entity_id: Optional[order_by] = None
    terra_sync_id: Optional[order_by] = None


class task_result_bool_exp(BaseModel):
    and_: Optional[List["task_result_bool_exp"]] = Field(alias="_and", default=None)
    not_: Optional["task_result_bool_exp"] = Field(alias="_not", default=None)
    or_: Optional[List["task_result_bool_exp"]] = Field(alias="_or", default=None)
    completed_at: Optional["timestamptz_comparison_exp"] = None
    crc_32_c_hash: Optional["String_comparison_exp"] = Field(
        alias="crc32c_hash", default=None
    )
    created_at: Optional["timestamptz_comparison_exp"] = None
    format: Optional["String_comparison_exp"] = None
    id: Optional["uuid_comparison_exp"] = None
    label: Optional["String_comparison_exp"] = None
    size: Optional["bigint_comparison_exp"] = None
    task_entity: Optional["task_entity_bool_exp"] = None
    task_entity_id: Optional["bigint_comparison_exp"] = None
    terra_entity_name: Optional["String_comparison_exp"] = None
    terra_entity_type: Optional["String_comparison_exp"] = None
    terra_method_config_name: Optional["String_comparison_exp"] = None
    terra_method_config_namespace: Optional["String_comparison_exp"] = None
    terra_submission_id: Optional["String_comparison_exp"] = None
    terra_sync: Optional["terra_sync_bool_exp"] = None
    terra_sync_id: Optional["bigint_comparison_exp"] = None
    terra_workflow_id: Optional["String_comparison_exp"] = None
    terra_workflow_inputs: Optional["jsonb_comparison_exp"] = None
    terra_workflow_root_dir: Optional["String_comparison_exp"] = None
    terra_workspace_id: Optional["String_comparison_exp"] = None
    terra_workspace_name: Optional["String_comparison_exp"] = None
    terra_workspace_namespace: Optional["String_comparison_exp"] = None
    url: Optional["String_comparison_exp"] = None
    value: Optional["jsonb_comparison_exp"] = None
    workflow_name: Optional["String_comparison_exp"] = None
    workflow_source_url: Optional["String_comparison_exp"] = None
    workflow_version: Optional["String_comparison_exp"] = None


class task_result_delete_at_path_input(BaseModel):
    terra_workflow_inputs: Optional[List[str]] = None
    value: Optional[List[str]] = None


class task_result_delete_elem_input(BaseModel):
    terra_workflow_inputs: Optional[int] = None
    value: Optional[int] = None


class task_result_delete_key_input(BaseModel):
    terra_workflow_inputs: Optional[str] = None
    value: Optional[str] = None


class task_result_inc_input(BaseModel):
    size: Optional[int] = None
    task_entity_id: Optional[int] = None
    terra_sync_id: Optional[int] = None


class task_result_insert_input(BaseModel):
    completed_at: Optional[Any] = None
    crc_32_c_hash: Optional[str] = Field(alias="crc32c_hash", default=None)
    created_at: Optional[Any] = None
    format: Optional[str] = None
    id: Optional[Any] = None
    label: Optional[str] = None
    size: Optional[int] = None
    task_entity: Optional["task_entity_obj_rel_insert_input"] = None
    task_entity_id: Optional[int] = None
    terra_entity_name: Optional[str] = None
    terra_entity_type: Optional[str] = None
    terra_method_config_name: Optional[str] = None
    terra_method_config_namespace: Optional[str] = None
    terra_submission_id: Optional[str] = None
    terra_sync: Optional["terra_sync_obj_rel_insert_input"] = None
    terra_sync_id: Optional[int] = None
    terra_workflow_id: Optional[str] = None
    terra_workflow_inputs: Optional[Any] = None
    terra_workflow_root_dir: Optional[str] = None
    terra_workspace_id: Optional[str] = None
    terra_workspace_name: Optional[str] = None
    terra_workspace_namespace: Optional[str] = None
    url: Optional[str] = None
    value: Optional[Any] = None
    workflow_name: Optional[str] = None
    workflow_source_url: Optional[str] = None
    workflow_version: Optional[str] = None


class task_result_max_order_by(BaseModel):
    completed_at: Optional[order_by] = None
    crc_32_c_hash: Optional[order_by] = Field(alias="crc32c_hash", default=None)
    created_at: Optional[order_by] = None
    format: Optional[order_by] = None
    id: Optional[order_by] = None
    label: Optional[order_by] = None
    size: Optional[order_by] = None
    task_entity_id: Optional[order_by] = None
    terra_entity_name: Optional[order_by] = None
    terra_entity_type: Optional[order_by] = None
    terra_method_config_name: Optional[order_by] = None
    terra_method_config_namespace: Optional[order_by] = None
    terra_submission_id: Optional[order_by] = None
    terra_sync_id: Optional[order_by] = None
    terra_workflow_id: Optional[order_by] = None
    terra_workflow_root_dir: Optional[order_by] = None
    terra_workspace_id: Optional[order_by] = None
    terra_workspace_name: Optional[order_by] = None
    terra_workspace_namespace: Optional[order_by] = None
    url: Optional[order_by] = None
    workflow_name: Optional[order_by] = None
    workflow_source_url: Optional[order_by] = None
    workflow_version: Optional[order_by] = None


class task_result_min_order_by(BaseModel):
    completed_at: Optional[order_by] = None
    crc_32_c_hash: Optional[order_by] = Field(alias="crc32c_hash", default=None)
    created_at: Optional[order_by] = None
    format: Optional[order_by] = None
    id: Optional[order_by] = None
    label: Optional[order_by] = None
    size: Optional[order_by] = None
    task_entity_id: Optional[order_by] = None
    terra_entity_name: Optional[order_by] = None
    terra_entity_type: Optional[order_by] = None
    terra_method_config_name: Optional[order_by] = None
    terra_method_config_namespace: Optional[order_by] = None
    terra_submission_id: Optional[order_by] = None
    terra_sync_id: Optional[order_by] = None
    terra_workflow_id: Optional[order_by] = None
    terra_workflow_root_dir: Optional[order_by] = None
    terra_workspace_id: Optional[order_by] = None
    terra_workspace_name: Optional[order_by] = None
    terra_workspace_namespace: Optional[order_by] = None
    url: Optional[order_by] = None
    workflow_name: Optional[order_by] = None
    workflow_source_url: Optional[order_by] = None
    workflow_version: Optional[order_by] = None


class task_result_on_conflict(BaseModel):
    constraint: task_result_constraint
    update_columns: List[task_result_update_column]
    where: Optional["task_result_bool_exp"] = None


class task_result_order_by(BaseModel):
    completed_at: Optional[order_by] = None
    crc_32_c_hash: Optional[order_by] = Field(alias="crc32c_hash", default=None)
    created_at: Optional[order_by] = None
    format: Optional[order_by] = None
    id: Optional[order_by] = None
    label: Optional[order_by] = None
    size: Optional[order_by] = None
    task_entity: Optional["task_entity_order_by"] = None
    task_entity_id: Optional[order_by] = None
    terra_entity_name: Optional[order_by] = None
    terra_entity_type: Optional[order_by] = None
    terra_method_config_name: Optional[order_by] = None
    terra_method_config_namespace: Optional[order_by] = None
    terra_submission_id: Optional[order_by] = None
    terra_sync: Optional["terra_sync_order_by"] = None
    terra_sync_id: Optional[order_by] = None
    terra_workflow_id: Optional[order_by] = None
    terra_workflow_inputs: Optional[order_by] = None
    terra_workflow_root_dir: Optional[order_by] = None
    terra_workspace_id: Optional[order_by] = None
    terra_workspace_name: Optional[order_by] = None
    terra_workspace_namespace: Optional[order_by] = None
    url: Optional[order_by] = None
    value: Optional[order_by] = None
    workflow_name: Optional[order_by] = None
    workflow_source_url: Optional[order_by] = None
    workflow_version: Optional[order_by] = None


class task_result_pk_columns_input(BaseModel):
    id: Any


class task_result_prepend_input(BaseModel):
    terra_workflow_inputs: Optional[Any] = None
    value: Optional[Any] = None


class task_result_set_input(BaseModel):
    completed_at: Optional[Any] = None
    crc_32_c_hash: Optional[str] = Field(alias="crc32c_hash", default=None)
    created_at: Optional[Any] = None
    format: Optional[str] = None
    id: Optional[Any] = None
    label: Optional[str] = None
    size: Optional[int] = None
    task_entity_id: Optional[int] = None
    terra_entity_name: Optional[str] = None
    terra_entity_type: Optional[str] = None
    terra_method_config_name: Optional[str] = None
    terra_method_config_namespace: Optional[str] = None
    terra_submission_id: Optional[str] = None
    terra_sync_id: Optional[int] = None
    terra_workflow_id: Optional[str] = None
    terra_workflow_inputs: Optional[Any] = None
    terra_workflow_root_dir: Optional[str] = None
    terra_workspace_id: Optional[str] = None
    terra_workspace_name: Optional[str] = None
    terra_workspace_namespace: Optional[str] = None
    url: Optional[str] = None
    value: Optional[Any] = None
    workflow_name: Optional[str] = None
    workflow_source_url: Optional[str] = None
    workflow_version: Optional[str] = None


class task_result_stddev_order_by(BaseModel):
    size: Optional[order_by] = None
    task_entity_id: Optional[order_by] = None
    terra_sync_id: Optional[order_by] = None


class task_result_stddev_pop_order_by(BaseModel):
    size: Optional[order_by] = None
    task_entity_id: Optional[order_by] = None
    terra_sync_id: Optional[order_by] = None


class task_result_stddev_samp_order_by(BaseModel):
    size: Optional[order_by] = None
    task_entity_id: Optional[order_by] = None
    terra_sync_id: Optional[order_by] = None


class task_result_stream_cursor_input(BaseModel):
    initial_value: "task_result_stream_cursor_value_input"
    ordering: Optional[cursor_ordering] = None


class task_result_stream_cursor_value_input(BaseModel):
    completed_at: Optional[Any] = None
    crc_32_c_hash: Optional[str] = Field(alias="crc32c_hash", default=None)
    created_at: Optional[Any] = None
    format: Optional[str] = None
    id: Optional[Any] = None
    label: Optional[str] = None
    size: Optional[int] = None
    task_entity_id: Optional[int] = None
    terra_entity_name: Optional[str] = None
    terra_entity_type: Optional[str] = None
    terra_method_config_name: Optional[str] = None
    terra_method_config_namespace: Optional[str] = None
    terra_submission_id: Optional[str] = None
    terra_sync_id: Optional[int] = None
    terra_workflow_id: Optional[str] = None
    terra_workflow_inputs: Optional[Any] = None
    terra_workflow_root_dir: Optional[str] = None
    terra_workspace_id: Optional[str] = None
    terra_workspace_name: Optional[str] = None
    terra_workspace_namespace: Optional[str] = None
    url: Optional[str] = None
    value: Optional[Any] = None
    workflow_name: Optional[str] = None
    workflow_source_url: Optional[str] = None
    workflow_version: Optional[str] = None


class task_result_sum_order_by(BaseModel):
    size: Optional[order_by] = None
    task_entity_id: Optional[order_by] = None
    terra_sync_id: Optional[order_by] = None


class task_result_updates(BaseModel):
    append: Optional["task_result_append_input"] = Field(alias="_append", default=None)
    delete_at_path: Optional["task_result_delete_at_path_input"] = Field(
        alias="_delete_at_path", default=None
    )
    delete_elem: Optional["task_result_delete_elem_input"] = Field(
        alias="_delete_elem", default=None
    )
    delete_key: Optional["task_result_delete_key_input"] = Field(
        alias="_delete_key", default=None
    )
    inc: Optional["task_result_inc_input"] = Field(alias="_inc", default=None)
    prepend: Optional["task_result_prepend_input"] = Field(
        alias="_prepend", default=None
    )
    set: Optional["task_result_set_input"] = Field(alias="_set", default=None)
    where: "task_result_bool_exp"


class task_result_var_pop_order_by(BaseModel):
    size: Optional[order_by] = None
    task_entity_id: Optional[order_by] = None
    terra_sync_id: Optional[order_by] = None


class task_result_var_samp_order_by(BaseModel):
    size: Optional[order_by] = None
    task_entity_id: Optional[order_by] = None
    terra_sync_id: Optional[order_by] = None


class task_result_variance_order_by(BaseModel):
    size: Optional[order_by] = None
    task_entity_id: Optional[order_by] = None
    terra_sync_id: Optional[order_by] = None


class terra_sync_bool_exp(BaseModel):
    and_: Optional[List["terra_sync_bool_exp"]] = Field(alias="_and", default=None)
    not_: Optional["terra_sync_bool_exp"] = Field(alias="_not", default=None)
    or_: Optional[List["terra_sync_bool_exp"]] = Field(alias="_or", default=None)
    created_at: Optional["timestamptz_comparison_exp"] = None
    id: Optional["bigint_comparison_exp"] = None
    task_results: Optional["task_result_bool_exp"] = None
    task_results_aggregate: Optional["task_result_aggregate_bool_exp"] = None
    terra_workspace_name: Optional["String_comparison_exp"] = None
    terra_workspace_namespace: Optional["String_comparison_exp"] = None


class terra_sync_inc_input(BaseModel):
    id: Optional[int] = None


class terra_sync_insert_input(BaseModel):
    created_at: Optional[Any] = None
    id: Optional[int] = None
    task_results: Optional["task_result_arr_rel_insert_input"] = None
    terra_workspace_name: Optional[str] = None
    terra_workspace_namespace: Optional[str] = None


class terra_sync_obj_rel_insert_input(BaseModel):
    data: "terra_sync_insert_input"
    on_conflict: Optional["terra_sync_on_conflict"] = None


class terra_sync_on_conflict(BaseModel):
    constraint: terra_sync_constraint
    update_columns: List[terra_sync_update_column]
    where: Optional["terra_sync_bool_exp"] = None


class terra_sync_order_by(BaseModel):
    created_at: Optional[order_by] = None
    id: Optional[order_by] = None
    task_results_aggregate: Optional["task_result_aggregate_order_by"] = None
    terra_workspace_name: Optional[order_by] = None
    terra_workspace_namespace: Optional[order_by] = None


class terra_sync_pk_columns_input(BaseModel):
    id: int


class terra_sync_set_input(BaseModel):
    created_at: Optional[Any] = None
    id: Optional[int] = None
    terra_workspace_name: Optional[str] = None
    terra_workspace_namespace: Optional[str] = None


class terra_sync_stream_cursor_input(BaseModel):
    initial_value: "terra_sync_stream_cursor_value_input"
    ordering: Optional[cursor_ordering] = None


class terra_sync_stream_cursor_value_input(BaseModel):
    created_at: Optional[Any] = None
    id: Optional[int] = None
    terra_workspace_name: Optional[str] = None
    terra_workspace_namespace: Optional[str] = None


class terra_sync_updates(BaseModel):
    inc: Optional["terra_sync_inc_input"] = Field(alias="_inc", default=None)
    set: Optional["terra_sync_set_input"] = Field(alias="_set", default=None)
    where: "terra_sync_bool_exp"


class timestamptz_comparison_exp(BaseModel):
    eq: Optional[Any] = Field(alias="_eq", default=None)
    gt: Optional[Any] = Field(alias="_gt", default=None)
    gte: Optional[Any] = Field(alias="_gte", default=None)
    in_: Optional[List[Any]] = Field(alias="_in", default=None)
    is_null: Optional[bool] = Field(alias="_is_null", default=None)
    lt: Optional[Any] = Field(alias="_lt", default=None)
    lte: Optional[Any] = Field(alias="_lte", default=None)
    neq: Optional[Any] = Field(alias="_neq", default=None)
    nin: Optional[List[Any]] = Field(alias="_nin", default=None)


class uuid_comparison_exp(BaseModel):
    eq: Optional[Any] = Field(alias="_eq", default=None)
    gt: Optional[Any] = Field(alias="_gt", default=None)
    gte: Optional[Any] = Field(alias="_gte", default=None)
    in_: Optional[List[Any]] = Field(alias="_in", default=None)
    is_null: Optional[bool] = Field(alias="_is_null", default=None)
    lt: Optional[Any] = Field(alias="_lt", default=None)
    lte: Optional[Any] = Field(alias="_lte", default=None)
    neq: Optional[Any] = Field(alias="_neq", default=None)
    nin: Optional[List[Any]] = Field(alias="_nin", default=None)
