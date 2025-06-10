from enum import Enum


class audit_user_constraint(str, Enum):
    audit_user_pkey = "audit_user_pkey"


class audit_user_select_column(str, Enum):
    username = "username"


class audit_user_update_column(str, Enum):
    username = "username"


class cursor_ordering(str, Enum):
    ASC = "ASC"
    DESC = "DESC"


class depmap_model_type_constraint(str, Enum):
    depmap_model_type_pkey = "depmap_model_type_pkey"


class depmap_model_type_select_column(str, Enum):
    depmap_code = "depmap_code"
    lineage = "lineage"
    oncotree_code = "oncotree_code"
    primary_disease = "primary_disease"
    subtype = "subtype"


class depmap_model_type_update_column(str, Enum):
    depmap_code = "depmap_code"
    lineage = "lineage"
    oncotree_code = "oncotree_code"
    primary_disease = "primary_disease"
    subtype = "subtype"


class genomic_fingerprint_comparison_constraint(str, Enum):
    genomic_fingerprint_comparison_pkey = "genomic_fingerprint_comparison_pkey"


class genomic_fingerprint_comparison_select_column(str, Enum):
    created_at = "created_at"
    genomic_fingerprint_id1 = "genomic_fingerprint_id1"
    genomic_fingerprint_id2 = "genomic_fingerprint_id2"
    id = "id"
    n_common_snps = "n_common_snps"
    n_matching_genotypes = "n_matching_genotypes"
    patient_id1 = "patient_id1"
    patient_id2 = "patient_id2"
    score = "score"


class genomic_fingerprint_comparison_select_column_genomic_fingerprint_comparison_aggregate_bool_exp_avg_arguments_columns(
    str, Enum
):
    score = "score"


class genomic_fingerprint_comparison_select_column_genomic_fingerprint_comparison_aggregate_bool_exp_corr_arguments_columns(
    str, Enum
):
    score = "score"


class genomic_fingerprint_comparison_select_column_genomic_fingerprint_comparison_aggregate_bool_exp_covar_samp_arguments_columns(
    str, Enum
):
    score = "score"


class genomic_fingerprint_comparison_select_column_genomic_fingerprint_comparison_aggregate_bool_exp_max_arguments_columns(
    str, Enum
):
    score = "score"


class genomic_fingerprint_comparison_select_column_genomic_fingerprint_comparison_aggregate_bool_exp_min_arguments_columns(
    str, Enum
):
    score = "score"


class genomic_fingerprint_comparison_select_column_genomic_fingerprint_comparison_aggregate_bool_exp_stddev_samp_arguments_columns(
    str, Enum
):
    score = "score"


class genomic_fingerprint_comparison_select_column_genomic_fingerprint_comparison_aggregate_bool_exp_sum_arguments_columns(
    str, Enum
):
    score = "score"


class genomic_fingerprint_comparison_select_column_genomic_fingerprint_comparison_aggregate_bool_exp_var_samp_arguments_columns(
    str, Enum
):
    score = "score"


class genomic_fingerprint_comparison_update_column(str, Enum):
    created_at = "created_at"
    genomic_fingerprint_id1 = "genomic_fingerprint_id1"
    genomic_fingerprint_id2 = "genomic_fingerprint_id2"
    id = "id"
    n_common_snps = "n_common_snps"
    n_matching_genotypes = "n_matching_genotypes"
    patient_id1 = "patient_id1"
    patient_id2 = "patient_id2"
    score = "score"


class genomic_fingerprint_constraint(str, Enum):
    genomic_fingerprint_pkey = "genomic_fingerprint_pkey"


class genomic_fingerprint_failure_constraint(str, Enum):
    genomic_fingerprint_failure_pkey = "genomic_fingerprint_failure_pkey"


class genomic_fingerprint_failure_genomic_fingerprint_comparison_constraint(str, Enum):
    genomic_fingerprint_fail_genomicfingerprintfailur_151d88a5_uniq = (
        "genomic_fingerprint_fail_genomicfingerprintfailur_151d88a5_uniq"
    )
    genomic_fingerprint_failure_genomic_fingerprint_comparison_pkey = (
        "genomic_fingerprint_failure_genomic_fingerprint_comparison_pkey"
    )


class genomic_fingerprint_failure_genomic_fingerprint_comparison_select_column(
    str, Enum
):
    genomicfingerprintcomparison_id = "genomicfingerprintcomparison_id"
    genomicfingerprintfailure_id = "genomicfingerprintfailure_id"
    id = "id"


class genomic_fingerprint_failure_genomic_fingerprint_comparison_update_column(
    str, Enum
):
    genomicfingerprintcomparison_id = "genomicfingerprintcomparison_id"
    genomicfingerprintfailure_id = "genomicfingerprintfailure_id"
    id = "id"


class genomic_fingerprint_failure_select_column(str, Enum):
    acknowledged = "acknowledged"
    comments = "comments"
    created_at = "created_at"
    id = "id"
    updated_at = "updated_at"


class genomic_fingerprint_failure_update_column(str, Enum):
    acknowledged = "acknowledged"
    comments = "comments"
    created_at = "created_at"
    id = "id"
    updated_at = "updated_at"


class genomic_fingerprint_select_column(str, Enum):
    created_at = "created_at"
    genotypes = "genotypes"
    id = "id"
    sequencing_alignment_id = "sequencing_alignment_id"
    vcf_url = "vcf_url"


class genomic_fingerprint_update_column(str, Enum):
    created_at = "created_at"
    genotypes = "genotypes"
    id = "id"
    sequencing_alignment_id = "sequencing_alignment_id"
    vcf_url = "vcf_url"


class media_constraint(str, Enum):
    media_pkey = "media_pkey"


class media_select_column(str, Enum):
    formulation = "formulation"
    media_id = "media_id"
    serum_free = "serum_free"


class media_update_column(str, Enum):
    formulation = "formulation"
    media_id = "media_id"
    serum_free = "serum_free"


class model_condition_constraint(str, Enum):
    model_condition_condition_only_key = "model_condition_condition_only_key"
    model_condition_pkey = "model_condition_pkey"


class model_condition_select_column(str, Enum):
    batch_doubling_time = "batch_doubling_time"
    cell_characteristics = "cell_characteristics"
    cell_format = "cell_format"
    cell_grouping = "cell_grouping"
    cell_morphology = "cell_morphology"
    cell_shape = "cell_shape"
    cell_size = "cell_size"
    comments = "comments"
    condition_only = "condition_only"
    contaminated = "contaminated"
    contamination_details = "contamination_details"
    days_with_drug = "days_with_drug"
    drug = "drug"
    drug_concentration = "drug_concentration"
    expansion_team = "expansion_team"
    growth_media = "growth_media"
    model_condition_id = "model_condition_id"
    model_id = "model_id"
    parent_model_condition_id = "parent_model_condition_id"
    passage_number = "passage_number"
    plate_coating = "plate_coating"
    prescreen_treatment_days = "prescreen_treatment_days"
    prescreen_treatment_drug = "prescreen_treatment_drug"
    project = "project"
    resistance_mechanism = "resistance_mechanism"
    source = "source"
    source_doubling_time = "source_doubling_time"
    supplements = "supplements"


class model_condition_select_column_model_condition_aggregate_bool_exp_bool_and_arguments_columns(
    str, Enum
):
    contaminated = "contaminated"


class model_condition_select_column_model_condition_aggregate_bool_exp_bool_or_arguments_columns(
    str, Enum
):
    contaminated = "contaminated"


class model_condition_update_column(str, Enum):
    batch_doubling_time = "batch_doubling_time"
    cell_characteristics = "cell_characteristics"
    cell_format = "cell_format"
    cell_grouping = "cell_grouping"
    cell_morphology = "cell_morphology"
    cell_shape = "cell_shape"
    cell_size = "cell_size"
    comments = "comments"
    condition_only = "condition_only"
    contaminated = "contaminated"
    contamination_details = "contamination_details"
    days_with_drug = "days_with_drug"
    drug = "drug"
    drug_concentration = "drug_concentration"
    expansion_team = "expansion_team"
    growth_media = "growth_media"
    model_condition_id = "model_condition_id"
    model_id = "model_id"
    parent_model_condition_id = "parent_model_condition_id"
    passage_number = "passage_number"
    plate_coating = "plate_coating"
    prescreen_treatment_days = "prescreen_treatment_days"
    prescreen_treatment_drug = "prescreen_treatment_drug"
    project = "project"
    resistance_mechanism = "resistance_mechanism"
    source = "source"
    source_doubling_time = "source_doubling_time"
    supplements = "supplements"


class model_constraint(str, Enum):
    model_pkey = "model_pkey"


class model_select_column(str, Enum):
    age = "age"
    age_category = "age_category"
    ancestry = "ancestry"
    catalog_number = "catalog_number"
    ccle_line = "ccle_line"
    ccle_name = "ccle_name"
    cell_line_aliases = "cell_line_aliases"
    cell_line_in_stock = "cell_line_in_stock"
    cell_line_name = "cell_line_name"
    cell_line_ordered_date = "cell_line_ordered_date"
    cell_line_received = "cell_line_received"
    comments = "comments"
    consent_2015 = "consent_2015"
    converge_id = "converge_id"
    cosmic_id = "cosmic_id"
    cultured_drug_resistance = "cultured_drug_resistance"
    date_cell_line_received = "date_cell_line_received"
    date_first_publication = "date_first_publication"
    date_model_derived = "date_model_derived"
    date_shared_in_dbgap = "date_shared_in_dbgap"
    dbgap = "dbgap"
    depmap_model_type_id = "depmap_model_type_id"
    derived_outside_us = "derived_outside_us"
    do_not_screen = "do_not_screen"
    engineered_model = "engineered_model"
    engineered_model_details = "engineered_model_details"
    first_publication_link = "first_publication_link"
    geo_loc = "geo_loc"
    growth_pattern = "growth_pattern"
    hcmi_id = "hcmi_id"
    inferred_ethnicity = "inferred_ethnicity"
    model_data_sharing = "model_data_sharing"
    model_data_sharing_comments = "model_data_sharing_comments"
    model_derivation_material = "model_derivation_material"
    model_id = "model_id"
    model_id_alias = "model_id_alias"
    model_subtype_features = "model_subtype_features"
    model_transfer = "model_transfer"
    model_transfer_comments = "model_transfer_comments"
    model_transferred_to_stjude = "model_transferred_to_stjude"
    model_type = "model_type"
    new_histological_subtype = "new_histological_subtype"
    onboarded_doubling_time = "onboarded_doubling_time"
    onboarded_media = "onboarded_media"
    orspid = "orspid"
    patient_id = "patient_id"
    patient_resistance = "patient_resistance"
    patient_subtype_features = "patient_subtype_features"
    patient_treatment_type = "patient_treatment_type"
    patient_tumor_grade = "patient_tumor_grade"
    peddep_line = "peddep_line"
    peddep_nominated = "peddep_nominated"
    peddep_subgroup = "peddep_subgroup"
    permission_to_release = "permission_to_release"
    plate_coating = "plate_coating"
    primary_disease = "primary_disease"
    primary_or_metastasis = "primary_or_metastasis"
    proposed_deliverable = "proposed_deliverable"
    proposed_release_date = "proposed_release_date"
    public_comments = "public_comments"
    rrid = "rrid"
    sample_collection_site = "sample_collection_site"
    sanger_model_id = "sanger_model_id"
    screen_comments = "screen_comments"
    sex = "sex"
    sj_compbio_id = "sj_compbio_id"
    source_detail = "source_detail"
    source_type = "source_type"
    stage = "stage"
    staging_system = "staging_system"
    stated_race = "stated_race"
    stjude_derived = "stjude_derived"
    stripped_cell_line_name = "stripped_cell_line_name"
    tissue_origin = "tissue_origin"
    transformed_type = "transformed_type"
    treatment_details = "treatment_details"
    treatment_status = "treatment_status"
    wtsi_master_cell_id = "wtsi_master_cell_id"


class model_update_column(str, Enum):
    age = "age"
    age_category = "age_category"
    ancestry = "ancestry"
    catalog_number = "catalog_number"
    ccle_line = "ccle_line"
    ccle_name = "ccle_name"
    cell_line_aliases = "cell_line_aliases"
    cell_line_in_stock = "cell_line_in_stock"
    cell_line_name = "cell_line_name"
    cell_line_ordered_date = "cell_line_ordered_date"
    cell_line_received = "cell_line_received"
    comments = "comments"
    consent_2015 = "consent_2015"
    converge_id = "converge_id"
    cosmic_id = "cosmic_id"
    cultured_drug_resistance = "cultured_drug_resistance"
    date_cell_line_received = "date_cell_line_received"
    date_first_publication = "date_first_publication"
    date_model_derived = "date_model_derived"
    date_shared_in_dbgap = "date_shared_in_dbgap"
    dbgap = "dbgap"
    depmap_model_type_id = "depmap_model_type_id"
    derived_outside_us = "derived_outside_us"
    do_not_screen = "do_not_screen"
    engineered_model = "engineered_model"
    engineered_model_details = "engineered_model_details"
    first_publication_link = "first_publication_link"
    geo_loc = "geo_loc"
    growth_pattern = "growth_pattern"
    hcmi_id = "hcmi_id"
    inferred_ethnicity = "inferred_ethnicity"
    model_data_sharing = "model_data_sharing"
    model_data_sharing_comments = "model_data_sharing_comments"
    model_derivation_material = "model_derivation_material"
    model_id = "model_id"
    model_id_alias = "model_id_alias"
    model_subtype_features = "model_subtype_features"
    model_transfer = "model_transfer"
    model_transfer_comments = "model_transfer_comments"
    model_transferred_to_stjude = "model_transferred_to_stjude"
    model_type = "model_type"
    new_histological_subtype = "new_histological_subtype"
    onboarded_doubling_time = "onboarded_doubling_time"
    onboarded_media = "onboarded_media"
    orspid = "orspid"
    patient_id = "patient_id"
    patient_resistance = "patient_resistance"
    patient_subtype_features = "patient_subtype_features"
    patient_treatment_type = "patient_treatment_type"
    patient_tumor_grade = "patient_tumor_grade"
    peddep_line = "peddep_line"
    peddep_nominated = "peddep_nominated"
    peddep_subgroup = "peddep_subgroup"
    permission_to_release = "permission_to_release"
    plate_coating = "plate_coating"
    primary_disease = "primary_disease"
    primary_or_metastasis = "primary_or_metastasis"
    proposed_deliverable = "proposed_deliverable"
    proposed_release_date = "proposed_release_date"
    public_comments = "public_comments"
    rrid = "rrid"
    sample_collection_site = "sample_collection_site"
    sanger_model_id = "sanger_model_id"
    screen_comments = "screen_comments"
    sex = "sex"
    sj_compbio_id = "sj_compbio_id"
    source_detail = "source_detail"
    source_type = "source_type"
    stage = "stage"
    staging_system = "staging_system"
    stated_race = "stated_race"
    stjude_derived = "stjude_derived"
    stripped_cell_line_name = "stripped_cell_line_name"
    tissue_origin = "tissue_origin"
    transformed_type = "transformed_type"
    treatment_details = "treatment_details"
    treatment_status = "treatment_status"
    wtsi_master_cell_id = "wtsi_master_cell_id"


class omics_mapping_select_column(str, Enum):
    datatype = "datatype"
    id = "id"
    model_condition_id = "model_condition_id"
    model_id = "model_id"
    omics_profile_id = "omics_profile_id"
    omics_sequencing_id = "omics_sequencing_id"
    priority = "priority"


class omics_profile_constraint(str, Enum):
    omics_profile_pkey = "omics_profile_pkey"


class omics_profile_select_column(str, Enum):
    actual_seq_technology = "actual_seq_technology"
    baits = "baits"
    bam_public_sra_path = "bam_public_sra_path"
    billing_date = "billing_date"
    blacklist_expiration_date = "blacklist_expiration_date"
    blacklist_omics = "blacklist_omics"
    blacklist_reason = "blacklist_reason"
    bsp_sample_id_csv = "bsp_sample_id_csv"
    cell_available = "cell_available"
    cell_pellet_needed = "cell_pellet_needed"
    collaborator_sample_id = "collaborator_sample_id"
    consortium_release_date = "consortium_release_date"
    consortium_retracted_date = "consortium_retracted_date"
    datatype = "datatype"
    deliverables = "deliverables"
    destination_datasets = "destination_datasets"
    drop_reason = "drop_reason"
    eta_for_omics_completion = "eta_for_omics_completion"
    extraction_needed = "extraction_needed"
    ibm_release_date = "ibm_release_date"
    internal_release_date = "internal_release_date"
    internal_retracted_date = "internal_retracted_date"
    issue = "issue"
    kit_id = "kit_id"
    lcset_protocol = "lcset_protocol"
    lcsets = "lcsets"
    line_received_by_gp = "line_received_by_gp"
    line_sent_to_gp = "line_sent_to_gp"
    main_sequencing_id = "main_sequencing_id"
    model_condition_id = "model_condition_id"
    omics_order_date = "omics_order_date"
    omics_profile_flagship = "omics_profile_flagship"
    omics_profile_funding_source = "omics_profile_funding_source"
    omics_return_date = "omics_return_date"
    pdo_title = "pdo_title"
    pdoid = "pdoid"
    pf_bases_bc = "pf_bases_bc"
    prioritized = "prioritized"
    product = "product"
    product_goal = "product_goal"
    profile_id = "profile_id"
    profile_source = "profile_source"
    project = "project"
    proposed_release_date = "proposed_release_date"
    public_release_date = "public_release_date"
    public_retracted_date = "public_retracted_date"
    quote_to_bill = "quote_to_bill"
    registered = "registered"
    resubmit_for_extraction = "resubmit_for_extraction"
    rna_delivery_date = "rna_delivery_date"
    sample_coverage_normalized = "sample_coverage_normalized"
    sample_coverage_rounded = "sample_coverage_rounded"
    sample_is_on_risk = "sample_is_on_risk"
    sample_type = "sample_type"
    shared_to_dbgap = "shared_to_dbgap"
    sm_id_matched = "sm_id_matched"
    smid_ordered = "smid_ordered"
    smid_returned = "smid_returned"
    status = "status"
    version = "version"
    wgs_delivery_date = "wgs_delivery_date"
    workspace = "workspace"


class omics_profile_select_column_omics_profile_aggregate_bool_exp_bool_and_arguments_columns(
    str, Enum
):
    blacklist_omics = "blacklist_omics"
    cell_available = "cell_available"
    cell_pellet_needed = "cell_pellet_needed"
    extraction_needed = "extraction_needed"
    prioritized = "prioritized"
    registered = "registered"
    resubmit_for_extraction = "resubmit_for_extraction"
    sample_is_on_risk = "sample_is_on_risk"
    shared_to_dbgap = "shared_to_dbgap"


class omics_profile_select_column_omics_profile_aggregate_bool_exp_bool_or_arguments_columns(
    str, Enum
):
    blacklist_omics = "blacklist_omics"
    cell_available = "cell_available"
    cell_pellet_needed = "cell_pellet_needed"
    extraction_needed = "extraction_needed"
    prioritized = "prioritized"
    registered = "registered"
    resubmit_for_extraction = "resubmit_for_extraction"
    sample_is_on_risk = "sample_is_on_risk"
    shared_to_dbgap = "shared_to_dbgap"


class omics_profile_update_column(str, Enum):
    actual_seq_technology = "actual_seq_technology"
    baits = "baits"
    bam_public_sra_path = "bam_public_sra_path"
    billing_date = "billing_date"
    blacklist_expiration_date = "blacklist_expiration_date"
    blacklist_omics = "blacklist_omics"
    blacklist_reason = "blacklist_reason"
    bsp_sample_id_csv = "bsp_sample_id_csv"
    cell_available = "cell_available"
    cell_pellet_needed = "cell_pellet_needed"
    collaborator_sample_id = "collaborator_sample_id"
    consortium_release_date = "consortium_release_date"
    consortium_retracted_date = "consortium_retracted_date"
    datatype = "datatype"
    deliverables = "deliverables"
    destination_datasets = "destination_datasets"
    drop_reason = "drop_reason"
    eta_for_omics_completion = "eta_for_omics_completion"
    extraction_needed = "extraction_needed"
    ibm_release_date = "ibm_release_date"
    internal_release_date = "internal_release_date"
    internal_retracted_date = "internal_retracted_date"
    issue = "issue"
    kit_id = "kit_id"
    lcset_protocol = "lcset_protocol"
    lcsets = "lcsets"
    line_received_by_gp = "line_received_by_gp"
    line_sent_to_gp = "line_sent_to_gp"
    main_sequencing_id = "main_sequencing_id"
    model_condition_id = "model_condition_id"
    omics_order_date = "omics_order_date"
    omics_profile_flagship = "omics_profile_flagship"
    omics_profile_funding_source = "omics_profile_funding_source"
    omics_return_date = "omics_return_date"
    pdo_title = "pdo_title"
    pdoid = "pdoid"
    pf_bases_bc = "pf_bases_bc"
    prioritized = "prioritized"
    product = "product"
    product_goal = "product_goal"
    profile_id = "profile_id"
    profile_source = "profile_source"
    project = "project"
    proposed_release_date = "proposed_release_date"
    public_release_date = "public_release_date"
    public_retracted_date = "public_retracted_date"
    quote_to_bill = "quote_to_bill"
    registered = "registered"
    resubmit_for_extraction = "resubmit_for_extraction"
    rna_delivery_date = "rna_delivery_date"
    sample_coverage_normalized = "sample_coverage_normalized"
    sample_coverage_rounded = "sample_coverage_rounded"
    sample_is_on_risk = "sample_is_on_risk"
    sample_type = "sample_type"
    shared_to_dbgap = "shared_to_dbgap"
    sm_id_matched = "sm_id_matched"
    smid_ordered = "smid_ordered"
    smid_returned = "smid_returned"
    status = "status"
    version = "version"
    wgs_delivery_date = "wgs_delivery_date"
    workspace = "workspace"


class omics_sequencing_constraint(str, Enum):
    omics_sequencing_pkey = "omics_sequencing_pkey"


class omics_sequencing_select_column(str, Enum):
    bam_qc = "bam_qc"
    blacklist = "blacklist"
    expected_type = "expected_type"
    gp_alignment = "gp_alignment"
    issue = "issue"
    month_sequencing_billed = "month_sequencing_billed"
    pdo_id = "pdo_id"
    prioritized = "prioritized"
    processed_sequence = "processed_sequence"
    processing_qc = "processing_qc"
    profile_id = "profile_id"
    sequencing_date = "sequencing_date"
    sequencing_id = "sequencing_id"
    sm_id = "sm_id"
    source = "source"
    str_profile = "str_profile"
    stranded = "stranded"
    update_time = "update_time"
    version = "version"
    year_sequencing_billed = "year_sequencing_billed"


class omics_sequencing_select_column_omics_sequencing_aggregate_bool_exp_bool_and_arguments_columns(
    str, Enum
):
    blacklist = "blacklist"
    prioritized = "prioritized"
    processed_sequence = "processed_sequence"
    stranded = "stranded"


class omics_sequencing_select_column_omics_sequencing_aggregate_bool_exp_bool_or_arguments_columns(
    str, Enum
):
    blacklist = "blacklist"
    prioritized = "prioritized"
    processed_sequence = "processed_sequence"
    stranded = "stranded"


class omics_sequencing_update_column(str, Enum):
    bam_qc = "bam_qc"
    blacklist = "blacklist"
    expected_type = "expected_type"
    gp_alignment = "gp_alignment"
    issue = "issue"
    month_sequencing_billed = "month_sequencing_billed"
    pdo_id = "pdo_id"
    prioritized = "prioritized"
    processed_sequence = "processed_sequence"
    processing_qc = "processing_qc"
    profile_id = "profile_id"
    sequencing_date = "sequencing_date"
    sequencing_id = "sequencing_id"
    sm_id = "sm_id"
    source = "source"
    str_profile = "str_profile"
    stranded = "stranded"
    update_time = "update_time"
    version = "version"
    year_sequencing_billed = "year_sequencing_billed"


class onboarding_job_constraint(str, Enum):
    onboarding_job_pkey = "onboarding_job_pkey"


class onboarding_job_select_column(str, Enum):
    created_at = "created_at"
    id = "id"
    n_samples = "n_samples"
    n_samples_excluded = "n_samples_excluded"
    n_samples_failed = "n_samples_failed"
    n_samples_new = "n_samples_new"
    n_samples_succeeded = "n_samples_succeeded"
    onboarding_workspace_id = "onboarding_workspace_id"
    succeeded = "succeeded"


class onboarding_job_select_column_onboarding_job_aggregate_bool_exp_bool_and_arguments_columns(
    str, Enum
):
    succeeded = "succeeded"


class onboarding_job_select_column_onboarding_job_aggregate_bool_exp_bool_or_arguments_columns(
    str, Enum
):
    succeeded = "succeeded"


class onboarding_job_update_column(str, Enum):
    created_at = "created_at"
    id = "id"
    n_samples = "n_samples"
    n_samples_excluded = "n_samples_excluded"
    n_samples_failed = "n_samples_failed"
    n_samples_new = "n_samples_new"
    n_samples_succeeded = "n_samples_succeeded"
    onboarding_workspace_id = "onboarding_workspace_id"
    succeeded = "succeeded"


class onboarding_sample_constraint(str, Enum):
    onboarding_sample_pkey = "onboarding_sample_pkey"


class onboarding_sample_select_column(str, Enum):
    created_at = "created_at"
    id = "id"
    issue = "issue"
    omics_profile_id = "omics_profile_id"
    onboarding_job_id = "onboarding_job_id"
    sequencing_alignment_id = "sequencing_alignment_id"
    sm_id = "sm_id"
    terra_sample_id = "terra_sample_id"


class onboarding_sample_update_column(str, Enum):
    created_at = "created_at"
    id = "id"
    issue = "issue"
    omics_profile_id = "omics_profile_id"
    onboarding_job_id = "onboarding_job_id"
    sequencing_alignment_id = "sequencing_alignment_id"
    sm_id = "sm_id"
    terra_sample_id = "terra_sample_id"


class onboarding_workspace_constraint(str, Enum):
    onboarding_workspace_pkey = "onboarding_workspace_pkey"


class onboarding_workspace_select_column(str, Enum):
    active = "active"
    created_at = "created_at"
    custom_values = "custom_values"
    excluded_terra_sample_ids = "excluded_terra_sample_ids"
    expected_type = "expected_type"
    gcs_destination_bucket = "gcs_destination_bucket"
    gcs_destination_prefix = "gcs_destination_prefix"
    id = "id"
    max_file_size = "max_file_size"
    min_file_size = "min_file_size"
    reference_genome = "reference_genome"
    source = "source"
    terra_col_names = "terra_col_names"
    updated_at = "updated_at"
    workspace_name = "workspace_name"
    workspace_namespace = "workspace_namespace"


class onboarding_workspace_update_column(str, Enum):
    active = "active"
    created_at = "created_at"
    custom_values = "custom_values"
    excluded_terra_sample_ids = "excluded_terra_sample_ids"
    expected_type = "expected_type"
    gcs_destination_bucket = "gcs_destination_bucket"
    gcs_destination_prefix = "gcs_destination_prefix"
    id = "id"
    max_file_size = "max_file_size"
    min_file_size = "min_file_size"
    reference_genome = "reference_genome"
    source = "source"
    terra_col_names = "terra_col_names"
    updated_at = "updated_at"
    workspace_name = "workspace_name"
    workspace_namespace = "workspace_namespace"


class order_by(str, Enum):
    asc = "asc"
    asc_nulls_first = "asc_nulls_first"
    asc_nulls_last = "asc_nulls_last"
    desc = "desc"
    desc_nulls_first = "desc_nulls_first"
    desc_nulls_last = "desc_nulls_last"


class patient_constraint(str, Enum):
    patient_pkey = "patient_pkey"


class patient_select_column(str, Enum):
    patient_id = "patient_id"


class patient_update_column(str, Enum):
    patient_id = "patient_id"


class sequencing_alignment_constraint(str, Enum):
    sequencing_alignment_pkey = "sequencing_alignment_pkey"


class sequencing_alignment_select_column(str, Enum):
    crc32c_hash = "crc32c_hash"
    created_at = "created_at"
    id = "id"
    index_url = "index_url"
    omics_sequencing_id = "omics_sequencing_id"
    reference_genome = "reference_genome"
    sequencing_alignment_source = "sequencing_alignment_source"
    size = "size"
    str_profile_id = "str_profile_id"
    updated_at = "updated_at"
    url = "url"


class sequencing_alignment_update_column(str, Enum):
    crc32c_hash = "crc32c_hash"
    created_at = "created_at"
    id = "id"
    index_url = "index_url"
    omics_sequencing_id = "omics_sequencing_id"
    reference_genome = "reference_genome"
    sequencing_alignment_source = "sequencing_alignment_source"
    size = "size"
    str_profile_id = "str_profile_id"
    updated_at = "updated_at"
    url = "url"


class str_profile_constraint(str, Enum):
    single_reference_per_patient = "single_reference_per_patient"
    str_profile_pkey = "str_profile_pkey"


class str_profile_select_column(str, Enum):
    amelogenin = "amelogenin"
    comments = "comments"
    created_at = "created_at"
    csf1po = "csf1po"
    d13s317 = "d13s317"
    d16s539 = "d16s539"
    d18s51 = "d18s51"
    d21s11 = "d21s11"
    d3s1358 = "d3s1358"
    d5s818 = "d5s818"
    d7s820 = "d7s820"
    d8s1179 = "d8s1179"
    fga = "fga"
    id = "id"
    is_reference = "is_reference"
    lab_corp_case_nbr = "lab_corp_case_nbr"
    lab_corp_spec_nbr = "lab_corp_spec_nbr"
    model_condition_id = "model_condition_id"
    mouse = "mouse"
    mycoplasma = "mycoplasma"
    patient_id = "patient_id"
    pellet_creation_date = "pellet_creation_date"
    pellet_submitted_date = "pellet_submitted_date"
    penta_d = "penta_d"
    penta_e = "penta_e"
    percentage_match_to_parental = "percentage_match_to_parental"
    sample_reference = "sample_reference"
    source = "source"
    source_group = "source_group"
    th01 = "th01"
    tpox = "tpox"
    updated_at = "updated_at"
    vwa = "vwa"


class str_profile_update_column(str, Enum):
    amelogenin = "amelogenin"
    comments = "comments"
    created_at = "created_at"
    csf1po = "csf1po"
    d13s317 = "d13s317"
    d16s539 = "d16s539"
    d18s51 = "d18s51"
    d21s11 = "d21s11"
    d3s1358 = "d3s1358"
    d5s818 = "d5s818"
    d7s820 = "d7s820"
    d8s1179 = "d8s1179"
    fga = "fga"
    id = "id"
    is_reference = "is_reference"
    lab_corp_case_nbr = "lab_corp_case_nbr"
    lab_corp_spec_nbr = "lab_corp_spec_nbr"
    model_condition_id = "model_condition_id"
    mouse = "mouse"
    mycoplasma = "mycoplasma"
    patient_id = "patient_id"
    pellet_creation_date = "pellet_creation_date"
    pellet_submitted_date = "pellet_submitted_date"
    penta_d = "penta_d"
    penta_e = "penta_e"
    percentage_match_to_parental = "percentage_match_to_parental"
    sample_reference = "sample_reference"
    source = "source"
    source_group = "source_group"
    th01 = "th01"
    tpox = "tpox"
    updated_at = "updated_at"
    vwa = "vwa"


class task_entity_constraint(str, Enum):
    task_entity_pkey = "task_entity_pkey"


class task_entity_select_column(str, Enum):
    id = "id"
    sequencing_id = "sequencing_id"


class task_entity_update_column(str, Enum):
    id = "id"
    sequencing_id = "sequencing_id"


class task_result_constraint(str, Enum):
    task_result_pkey = "task_result_pkey"


class task_result_select_column(str, Enum):
    completed_at = "completed_at"
    crc32c_hash = "crc32c_hash"
    created_at = "created_at"
    format = "format"
    id = "id"
    label = "label"
    size = "size"
    task_entity_id = "task_entity_id"
    terra_entity_name = "terra_entity_name"
    terra_entity_type = "terra_entity_type"
    terra_method_config_name = "terra_method_config_name"
    terra_method_config_namespace = "terra_method_config_namespace"
    terra_submission_id = "terra_submission_id"
    terra_sync_id = "terra_sync_id"
    terra_workflow_id = "terra_workflow_id"
    terra_workflow_inputs = "terra_workflow_inputs"
    terra_workflow_root_dir = "terra_workflow_root_dir"
    terra_workspace_id = "terra_workspace_id"
    terra_workspace_name = "terra_workspace_name"
    terra_workspace_namespace = "terra_workspace_namespace"
    url = "url"
    value = "value"
    workflow_name = "workflow_name"
    workflow_source_url = "workflow_source_url"
    workflow_version = "workflow_version"


class task_result_update_column(str, Enum):
    completed_at = "completed_at"
    crc32c_hash = "crc32c_hash"
    created_at = "created_at"
    format = "format"
    id = "id"
    label = "label"
    size = "size"
    task_entity_id = "task_entity_id"
    terra_entity_name = "terra_entity_name"
    terra_entity_type = "terra_entity_type"
    terra_method_config_name = "terra_method_config_name"
    terra_method_config_namespace = "terra_method_config_namespace"
    terra_submission_id = "terra_submission_id"
    terra_sync_id = "terra_sync_id"
    terra_workflow_id = "terra_workflow_id"
    terra_workflow_inputs = "terra_workflow_inputs"
    terra_workflow_root_dir = "terra_workflow_root_dir"
    terra_workspace_id = "terra_workspace_id"
    terra_workspace_name = "terra_workspace_name"
    terra_workspace_namespace = "terra_workspace_namespace"
    url = "url"
    value = "value"
    workflow_name = "workflow_name"
    workflow_source_url = "workflow_source_url"
    workflow_version = "workflow_version"


class terra_sync_constraint(str, Enum):
    terra_sync_pkey = "terra_sync_pkey"


class terra_sync_select_column(str, Enum):
    created_at = "created_at"
    id = "id"
    terra_workspace_name = "terra_workspace_name"
    terra_workspace_namespace = "terra_workspace_namespace"


class terra_sync_update_column(str, Enum):
    created_at = "created_at"
    id = "id"
    terra_workspace_name = "terra_workspace_name"
    terra_workspace_namespace = "terra_workspace_namespace"
