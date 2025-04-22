version 1.0

workflow quantify_lr_rna {
    input {
        String workflow_version = "1.0" # internal semver
        String workflow_source_url # populated automatically with URL of this script

        String sample_id
        File input_bam
        File input_bai
        File ref_fasta
        File ref_annotation_gtf
        File ref_annotation_db
        File? star_junctions
    }

    call run_isoquant {
        input:
            sample_id = sample_id,
            input_bam = input_bam,
            input_bai = input_bai,
            ref_annotation_db = ref_annotation_db,
            ref_fasta = ref_fasta
    }

    call run_sqanti3 {
        input:
            sample_id = sample_id,
            isoquant_gtf = run_isoquant.extended_annotation,
            star_junctions = star_junctions,
            ref_annotation_gtf = ref_annotation_gtf,
            ref_fasta = ref_fasta
    }

    output {
        File transcript_counts = run_isoquant.transcript_counts
        File model_counts = run_isoquant.model_counts
        File transcript_tpm = run_isoquant.transcript_tpm
        File transcript_model_tpm = run_isoquant.transcript_model_tpm
        File gene_tpm = run_isoquant.gene_tpm
        File extended_annotation = run_isoquant.extended_annotation
        File gene_counts = run_isoquant.gene_counts
        File sq_junctions = run_sqanti3.sq_junctions
        File sq_class = run_sqanti3.sq_class
    }
}

task run_isoquant {
    input {
        String sample_id
        File input_bam
        File input_bai
        File ref_fasta
        File ref_annotation_db

        String docker_image
        String docker_image_hash_or_tag
        Int cpu = 8
        Int mem_gb = 32
        Int preemptible = 1
        Int max_retries = 1
        Int additional_disk_gb = 0
    }

    Int disk_space = (
        ceil(
            size(input_bam, "GiB") + size(ref_annotation_db, "GiB")
            + size(ref_fasta, "GiB")
        ) + 20 + additional_disk_gb
    )

    command <<<
        set -euo pipefail

        touch "~{input_bai}" # BAI must be newer than BAM to avoid warning

        /usr/local/bin/isoquant.py \
            --reference ~{ref_fasta} \
            --genedb ~{ref_annotation_db} \
            --complete_genedb \
            --bam ~{input_bam} \
            --data_type pacbio_ccs \
            --stranded forward \
            --transcript_quantification unique_splicing_consistent \
            --model_construction_strategy fl_pacbio \
            --gene_quantification unique_splicing_consistent \
            --report_canonical all \
            --report_novel_unspliced false \
            --threads ~{cpu} \
            --labels ~{sample_id} \
            --prefix ~{sample_id} \
            --counts_format both

        find isoquant_output/~{sample_id}/ -maxdepth 1 -type f -not -name '*.gz' -exec gzip {} +
    >>>

    output {
        File transcript_counts = "isoquant_output/~{sample_id}/~{sample_id}.transcript_counts.tsv.gz"
        File model_counts = "isoquant_output/~{sample_id}/~{sample_id}.transcript_model_counts.tsv.gz"
        File transcript_tpm = "isoquant_output/~{sample_id}/~{sample_id}.transcript_tpm.tsv.gz"
        File transcript_model_tpm = "isoquant_output/~{sample_id}/~{sample_id}.transcript_model_tpm.tsv.gz"
        File gene_tpm = "isoquant_output/~{sample_id}/~{sample_id}.gene_tpm.tsv.gz"
        File gene_counts = "isoquant_output/~{sample_id}/~{sample_id}.gene_counts.tsv.gz"
        File extended_annotation = "isoquant_output/~{sample_id}/~{sample_id}.extended_annotation.gtf.gz"
    }

    runtime {
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_space + " SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }

    meta {
        allowNestedInputs: true
    }
}

task run_sqanti3 {
	input {
		String sample_id
        File isoquant_gtf
		File? star_junctions
        File ref_annotation_gtf
        File ref_fasta

        String docker_image
        String docker_image_hash_or_tag
        Int cpu = 2
        Int mem_gb = 8
        Int preemptible = 2
        Int max_retries = 1
        Int additional_disk_gb = 0
    }

    Int disk_space = (
        ceil(
            size(isoquant_gtf, "GiB") * 3 + size(star_junctions, "GiB") * 3
            + size(ref_annotation_gtf, "GiB") + size(ref_fasta, "GiB")
        ) + 20 + additional_disk_gb
    )

    command <<<
        set -euo pipefail

        zcat "~{isoquant_gtf}" \
            | awk '{ if ($7 != ".") print }' \
            > "~{isoquant_gtf}.unzipped"

        if [[ -f "~{star_junctions}" ]];
        then
            zcat "~{star_junctions}" > "~{star_junctions}.unzipped"
            coverage_opt="--coverage ~{star_junctions}.unzipped"
        else
            coverage_opt=""
        fi

        python /usr/local/src/SQANTI3-5.1.2/sqanti3_qc.py \
            --report skip \
            --skipORF \
            $coverage_opt \
            --output "~{sample_id}" \
            "~{isoquant_gtf}.unzipped" \
            "~{ref_annotation_gtf}" \
            "~{ref_fasta}"
    >>>

    output {
    	File sq_junctions = "~{sample_id}_junctions.txt"
        File sq_class = "~{sample_id}_classification.txt"
    }

    runtime {
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_space + " SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }

    meta {
        allowNestedInputs: true
    }
}
