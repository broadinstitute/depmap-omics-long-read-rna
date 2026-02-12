version 1.0

workflow requantify_lr_rna {
    meta {
        description: "TODO"
    }

    parameter_meta {
        # inputs
        sample_id: "ID of this sample"
        input_bam: "aligned BAM file"
        input_bai: "index of BAM file"
        transcript_counts: "transcript counts from first-pass of IsoQuant"
        combined_sorted_gtf: "combined+sorted GTF from latest run of `combine_gtfs`"
        updated_tracking_sq_filtered: "TODO from latest run of `combine_gtfs`"
        data_type: "IsoQuant option"
        model_construction_strategy: "IsoQuant option"
        stranded: "IsoQuant option"
        transcript_quantification: "IsoQuant option"
        gene_quantification: "IsoQuant option"
        report_novel_unspliced: "IsoQuant option"
        report_canonical: "IsoQuant option"
        ref_fasta: "reference sequence FASTA"

        # outputs
        requantified_exon_counts: "TODO"
        requantified_intron_counts: "TODO"
        requantified_read_assignments: "TODO"
        requantified_transcript_counts: "TODO"
        requantified_transcript_tpm: "TODO"
    }

    input {
        # per-sample inputs
        String sample_id
        File input_bam
        File input_bai
        File transcript_counts

        # outputs from `combine_gtfs`
        File combined_sorted_gtf
        File updated_tracking_sq_filtered

        # isoquant options
        String data_type
        String model_construction_strategy
        String stranded
        String transcript_quantification
        String gene_quantification
        String report_novel_unspliced
        String report_canonical

        File ref_fasta = "gs://ccleparams/hg38ref_no_alt/GRCh38_no_alt.fa"
    }

    call filter_gtf_per_sample {
        input:
            sample_id = sample_id,
            transcript_counts = transcript_counts,
            combined_sorted_gtf = combined_sorted_gtf,
            updated_tracking_sq_filtered = updated_tracking_sq_filtered
    }

    call gtf_to_db {
        input:
            gtf = filter_gtf_per_sample.filtered_gtf
    }

    call run_isoquant {
        input:
            sample_id = sample_id,
            input_bam = input_bam,
            input_bai = input_bai,
            ref_annotation_db = gtf_to_db.db,
            ref_fasta = ref_fasta,
            data_type = data_type,
            model_construction_strategy = model_construction_strategy,
            stranded = stranded,
            transcript_quantification = transcript_quantification,
            gene_quantification = gene_quantification,
            report_novel_unspliced = report_novel_unspliced,
            report_canonical = report_canonical
    }

    output {
        File requantified_exon_counts = run_isoquant.exon_counts
        File requantified_intron_counts = run_isoquant.intron_counts
        File requantified_read_assignments = run_isoquant.read_assignments
        File requantified_transcript_counts = run_isoquant.transcript_counts
        File requantified_transcript_tpm = run_isoquant.transcript_tpm
    }
}

task filter_gtf_per_sample {
    meta {
        description: "TODO"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        sample_id: "ID of this sample"
        transcript_counts: "transcript counts from first-pass of IsoQuant"
        combined_sorted_gtf: "combined+sorted GTF from `combine_gtfs`"
        updated_tracking_sq_filtered: "TODO"

        # outputs
        filtered_gtf: "TODO"
    }

    input {
        String sample_id
        File transcript_counts
        File combined_sorted_gtf
        File updated_tracking_sq_filtered

        String docker_image  = "us-central1-docker.pkg.dev/depmap-omics/terra-images/combine-requantify-tools"
        String docker_image_hash_or_tag = ":production"
        Int cpu = 2
        Int mem_gb = 8
        Int preemptible = 1
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = (
        ceil(
            size(transcript_counts, "GiB")
            + 2 * size(combined_sorted_gtf, "GiB")
            + size(updated_tracking_sq_filtered, "GiB")
        ) + 20 + additional_disk_gb
    )

    command <<<
        set -euo pipefail

        python -m combine_requantify_tools \
            filter-sample-gtf \
            --gtf-in="~{combined_sorted_gtf}" \
            --tracking="~{updated_tracking_sq_filtered}" \
            --sample-id="~{sample_id}" \
            --transcript-counts="~{transcript_counts}" \
            --gtf-out="~{sample_id}_filtered.gtf"
    >>>

    output {
        File filtered_gtf = "~{sample_id}_filtered.gtf"
    }

    runtime {
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_space + " SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}

task gtf_to_db {
    meta {
        description: "TODO"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        gtf: "GTF convert to SQLite"

        # outputs
        db: "SQLite version of the input GTF"
    }

    input {
        File gtf

        String docker_image  = "us-central1-docker.pkg.dev/depmap-omics/terra-images/combine-requantify-tools"
        String docker_image_hash_or_tag = ":production"
        Int cpu = 2
        Int mem_gb = 16
        Int preemptible = 1
        Int max_retries = 0
        Int additional_disk_gb = 0
   }

    String db_name = basename(gtf, ".gtf")

    Int disk_space = (
        ceil(size(gtf, "GiB")) * 3 + 20 + additional_disk_gb
    )

    command <<<
        set -euo pipefail

        python -m combine_requantify_tools \
            gtf-to-db \
            --gtf-in="~{gtf}" \
            --db-out="~{db_name}.db"
    >>>

    output {
        File db = "~{db_name}.db"
    }

    runtime {
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_space + " SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}

task run_isoquant {
    meta {
        description: "TODO"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        sample_id: "ID of this sample"
        input_bam: "aligned BAM file"
        input_bai: "Iindex of BAM file"
        ref_annotation_db: "reference annotation database from gtf_to_db"
        data_type: "IsoQuant option"
        model_construction_strategy: "IsoQuant option"
        stranded: "IsoQuant option"
        transcript_quantification: "IsoQuant option"
        gene_quantification: "IsoQuant option"
        report_novel_unspliced: "IsoQuant option"
        report_canonical: "IsoQuant option"
        ref_fasta: "reference sequence FASTA"

        # outputs
        exon_counts: "TODO"
        intron_counts: "TODO"
        read_assignments: "TODO"
        transcript_counts: "TODO"
        transcript_tpm: "TODO"
    }

    input {
        String sample_id
        File input_bam
        File input_bai
        File ref_annotation_db
        String data_type
        String model_construction_strategy
        String stranded
        String transcript_quantification
        String gene_quantification
        String report_novel_unspliced
        String report_canonical
        File ref_fasta

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/isoquant"
        String docker_image_hash_or_tag = ":3.7.0--hdfd78af_0"
        Int cpu = 8
        Int mem_gb = 32
        Int preemptible = 1
        Int max_retries = 0
        Int additional_disk_gb = 0
   }

    Int disk_space = (
        ceil(
            size(input_bam, "GiB")
            + size(ref_annotation_db, "GiB")
            + size(ref_fasta, "GiB")
        ) + 20 + additional_disk_gb
    )

    command <<<
        set -euo pipefail

        # Ensure BAI is newer than BAM to avoid warnings
        touch "~{input_bai}"

        # Run IsoQuant with custom annotation database
        /usr/local/bin/isoquant.py \
            --bam "~{input_bam}" \
            --count_exons \
            --no_model_construction \
            --data_type "~{data_type}" \
            --gene_quantification "~{gene_quantification}" \
            --genedb "~{ref_annotation_db}" \
            --labels "~{sample_id}" \
            --prefix "~{sample_id}" \
            --reference "~{ref_fasta}" \
            --report_canonical "~{report_canonical}" \
            --stranded "~{stranded}" \
            --threads ~{cpu} \
            --transcript_quantification "~{transcript_quantification}"

        echo "Compressing output TSVS"
        find isoquant_output/~{sample_id}/ -maxdepth 1 -type f -not -name '*.gz' -exec gzip {} +
    >>>

    output {
        File exon_counts = "isoquant_output/~{sample_id}/~{sample_id}.exon_counts.tsv.gz"
        File intron_counts = "isoquant_output/~{sample_id}/~{sample_id}.intron_counts.tsv.gz"
        File read_assignments = "isoquant_output/~{sample_id}/~{sample_id}.read_assignments.tsv.gz"
        File transcript_counts = "isoquant_output/~{sample_id}/~{sample_id}.transcript_counts.tsv.gz"
        File transcript_tpm = "isoquant_output/~{sample_id}/~{sample_id}.transcript_tpm.tsv.gz"
    }

    runtime {
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_space + " SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}
