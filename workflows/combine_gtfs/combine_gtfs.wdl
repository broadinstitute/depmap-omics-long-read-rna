version 1.0

workflow combine_gtfs {
    meta {
        description: "TODO"
    }

    parameter_meta {
        # inputs
        sample_set_id: "identifier for this set of samples"
        sample_ids: "array of the sample IDs in this sample set"
        extended_annotation: "array of extended_annotation output files from quantify_lr_rna"
        discovered_transcript_counts: "array of discovered_transcript_counts output files from quantify_lr_rna"
        gencode_gtf: "Gencode GTF file (if this is the initial run), or a combined_sorted_gtf (from a previous run)"
        prefix: "annotation for transcripts identified in this run (mostly for internal consistency between tasks)"
        ref_fasta: "reference sequence FASTA"

        # outputs
        combined_sorted_gtf: "TODO"
        sq_class: "TODO"
        transcriptome_fasta: "TODO"
        updated_tracking_sq_filtered: "TODO"
    }

    input {
        String sample_set_id
        Array[String] sample_ids
        Array[File] extended_annotation
        Array[File] discovered_transcript_counts
        File gencode_gtf
        String prefix = "TCONS"
        File ref_fasta = "gs://ccleparams/hg38ref_no_alt/GRCh38_no_alt.fa"
    }

    Array[Pair[String, File]] sample_ids_annots = zip(sample_ids, extended_annotation)

    scatter(pair in sample_ids_annots) {
        call filter_isoquant {
            input:
                sample_id = pair.left,
                extended_annotation = pair.right
        }
    }

    call replace_transcript_ids {
        input:
            sample_set_id = sample_set_id,
            gtf_list = filter_isoquant.filtered_gtf,
            sample_ids = sample_ids,
            gencode_gtf = gencode_gtf,
            prefix = prefix
    }

    call run_sqanti3 {
        input:
            sample_set_id = sample_set_id,
            isoquant_gtf = replace_transcript_ids.gtf_out,
            gencode_gtf = gencode_gtf,
            ref_fasta = ref_fasta
    }

    call process_tracking_file {
        input:
            sample_set_id = sample_set_id,
            tracking_file = replace_transcript_ids.tracking_file,
            sample_ids = sample_ids,
            discovered_transcript_counts = discovered_transcript_counts
   }

    call filter_gtf_and_tracking {
        input:
            sample_set_id = sample_set_id,
            combined_gtf = replace_transcript_ids.gtf_out,
            gencode_gtf = gencode_gtf,
            squanti_classification = run_sqanti3.sq_class,
            updated_tracking = process_tracking_file.updated_tracking,
            prefix = prefix,
            ref_fasta = ref_fasta
   }

    call process_gtf {
        input:
            sample_set_id = sample_set_id,
            filtered_gtf = filter_gtf_and_tracking.filtered_gtf,
            gencode_gtf = gencode_gtf,
            ref_fasta = ref_fasta
   }

    output {
        File combined_sorted_gtf = process_gtf.sorted_gtf
        File sq_class = run_sqanti3.sq_class
        File transcriptome_fasta = filter_gtf_and_tracking.transcriptome_fasta
        File updated_tracking_sq_filtered = filter_gtf_and_tracking.updated_tracking_sq_filtered
    }
}

task filter_isoquant {
    meta {
        description: "TODO"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        sample_id: "identifier for this sample"
        extended_annotation: "quantify_lr_rna.extended_annotation output file for this sample"

        # outputs
        filtered_gtf: "GTF file containing only IsoQuant entries"
    }

    input {
        String sample_id
        File extended_annotation

        String docker_image  = "us-central1-docker.pkg.dev/depmap-omics/terra-images/combine-requantify-tools"
        String docker_image_hash_or_tag = ":production"
        Int cpu = 1
        Int mem_gb = 2
        Int preemptible = 2
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(size(extended_annotation, "GiB")) + 20 + additional_disk_gb

    command <<<
        set -euo pipefail

        # Decompress GTF file if needed
        if [[ "~{extended_annotation}" == *.gz ]]; then
            gunzip -c "~{extended_annotation}" > extended_annotation.gtf
        else
            cp "~{extended_annotation}" extended_annotation.gtf
        fi

        # Error handling for decompression issues
        if [[ ! -f extended_annotation.gtf ]]; then
            echo "ERROR: extended_annotation.gtf not found after decompression." >&2
            touch filtered.gtf
        else
            # Filter for lines where the source column (column 2) is "IsoQuant"
                awk '$2 == "IsoQuant"' extended_annotation.gtf > filtered.gtf
            if [[ ! -s filtered.gtf ]]; then
                echo "WARNING: No IsoQuant entries found for ~{sample_id}" >&2
                touch filtered.gtf
            fi
        fi

        mv filtered.gtf "~{sample_id}_filtered.gtf"
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

task replace_transcript_ids {
    meta {
        description: "TODO"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        sample_set_id: "identifier for this set of samples"
        gtf_list: "array of the filtered GTFs from filter_isoquant"
        sample_ids: "array of the sample IDs in this sample set"
        gencode_gtf: "Gencode GTF file (if this is the initial run), or a combined_sorted_gtf (from a previous run)"
        prefix: "annotation for transcripts identified in this run (mostly for internal consistency between tasks)"

        # outputs
        tracking_file: "tracking file showing transcript relationships"
        gtf_out: "combined GTF from all samples"
    }

    input {
        String sample_set_id
        Array[File] gtf_list
        Array[String] sample_ids
        File gencode_gtf
        String prefix

        String docker_image  = "us-central1-docker.pkg.dev/depmap-omics/terra-images/combine-requantify-tools"
        String docker_image_hash_or_tag = ":production"
        Int cpu = 8
        Int mem_gb = 64
        Int preemptible = 1
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = (
        ceil(size(gtf_list, "GiB") + size(gencode_gtf, "GiB")) + 20 + additional_disk_gb
    )

    command <<<
        set -euo pipefail

        echo "Running gffcompare"
        gffcompare -r "~{gencode_gtf}" -X -p "~{prefix}" -V -S -o gffcomp_out ~{sep=' ' gtf_list}

        # Write sample_ids to newline-delimited text file
        printf '%s\n' ~{sep=' ' sample_ids} > sample_ids.txt

        python -m combine_requantify_tools \
            map-transcript-ids \
            --tracking-in="gffcomp_out.tracking" \
            --gtf-in="gffcomp_out.combined.gtf" \
            --sample-ids-list="sample_ids.txt" \
            --tracking-out="~{sample_set_id}_gffcompared.tracking.parquet" \
            --gtf-out="~{sample_set_id}_gffcompared.gtf"
    >>>

    output {
        File tracking_file = "~{sample_set_id}_gffcompared.tracking.parquet"
        File gtf_out = "~{sample_set_id}_gffcompared.gtf"
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

task run_sqanti3 {
    meta {
        description: "TODO"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        sample_set_id: "identifier for this set of samples"
        isoquant_gtf: "combined GTF from replace_transcript_ids"
        gencode_gtf: "Gencode GTF file (if this is the initial run), or a combined_sorted_gtf (from a previous run)"
        ref_fasta: "reference sequence FASTA"

        # outputs
        sq_class: "TODO"
        sq_junctions: "TODO"
    }

    input {
        String sample_set_id
        File isoquant_gtf
        File gencode_gtf
        File ref_fasta

        String docker_image = "us-central1-docker.pkg.dev/methods-dev-lab/lrtools-sqanti3/lrtools-sqanti3-plus"
        String docker_image_hash_or_tag = "@sha256:796ba14856e0e2bc55b3e4770fdc8d2b18af48251fba2a08747144501041437b"
        Int cpu = 2
        Int mem_gb = 64
        Int preemptible = 2
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = (
        ceil(
            size(isoquant_gtf, "GiB")
            + size(gencode_gtf, "GiB")
            + size(ref_fasta, "GiB")
        )
        + 20 + additional_disk_gb
    )

    command <<<
        set -euo pipefail

        awk '{ if ($7 != ".") print }' "~{isoquant_gtf}" \
            > "~{sample_set_id}.in.gtf"

        python /usr/local/src/SQANTI3-5.3.6/sqanti3_qc.py \
            --report both \
            --output "~{sample_set_id}" \
            "~{sample_set_id}.in.gtf" \
            "~{gencode_gtf}" \
            "~{ref_fasta}"
    >>>

    output {
        File sq_class = "~{sample_set_id}_classification.txt"
        File sq_junctions = "~{sample_set_id}_junctions.txt"
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

task process_tracking_file {
    meta {
        description: "TODO"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        sample_set_id: "identifier for this set of samples"
        tracking_file: "output tracking_file from replace_transcript_ids"
        sample_ids: "array of the sample IDs in this sample set"
        discovered_transcript_counts: "array of discovered_transcript_counts output files from quantify_lr_rna"

        # outputs
        updated_tracking: "TODO"
    }

    input {
        String sample_set_id
        File tracking_file
        Int min_count = 5
        Array[String] sample_ids
        Array[File] discovered_transcript_counts

        String docker_image  = "us-central1-docker.pkg.dev/depmap-omics/terra-images/combine-requantify-tools"
        String docker_image_hash_or_tag = ":production"
        Int cpu = 8
        Int mem_gb = 32
        Int preemptible = 2
        Int additional_disk_gb = 0
        Int max_retries = 0
    }

    Int disk_space = (
        ceil(size(tracking_file, "GiB") * 3 + size(discovered_transcript_counts, "GiB"))
        + 20 + additional_disk_gb
    )

    command <<<
        set -euo pipefail

        # Write sample_ids to newline-delimited text file
        printf '%s\n' ~{sep=' ' sample_ids} > sample_ids.txt

        # Write discovered_transcript_counts file paths to newline-delimited text file
        printf '%s\n' ~{sep=' ' discovered_transcript_counts} \
            > discovered_transcript_counts_files.txt

        python -m combine_requantify_tools \
            process-tracking-file \
            --tracking-in="~{tracking_file}" \
            --sample-ids-list="sample_ids.txt" \
            --discovered-transcript-counts-file-list="discovered_transcript_counts_files.txt" \
            --min-count=~{min_count} \
            --tracking-out="~{sample_set_id}_updated_tracking.parquet"
    >>>

    output {
        File updated_tracking = "~{sample_set_id}_updated_tracking.parquet"
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

task filter_gtf_and_tracking {
    meta {
        description: "TODO"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        sample_set_id: "identifier for this set of samples"
        gencode_gtf: "Gencode GTF file (if this is the initial run), or a combined_sorted_gtf (from a previous run)"
        combined_gtf: "Combined GTF from replace_transcript_ids"
        squanti_classification: "sq_class file from run_sqanti3"
        updated_tracking: "updated_tracking file from process_tracking_file"
        prefix: "annotation for transcripts identified in this run (mostly for internal consistency between tasks)"
        ref_fasta: "reference sequence FASTA"

        # outputs
        transcriptome_fasta: "TODO"
        filtered_gtf: "TODO"
        updated_tracking_sq_filtered: "TODO"
    }

    input {
        String sample_set_id
        File gencode_gtf
        File combined_gtf
        File squanti_classification
        File updated_tracking
        String prefix
        File ref_fasta

        String docker_image  = "us-central1-docker.pkg.dev/depmap-omics/terra-images/combine-requantify-tools"
        String docker_image_hash_or_tag = ":production"
        Int cpu = 2
        Int mem_gb = 8
        Int preemptible = 2
        Int max_retries = 0
        Int additional_disk_gb = 0
   }

    Int disk_space = (
        ceil(
            size(gencode_gtf, "GiB")
            + size(combined_gtf, "GiB")
            + size(squanti_classification, "GiB")
            + size(updated_tracking, "GiB")
            + size(ref_fasta, "GiB")
        )
        + 20 + additional_disk_gb
    )

    command <<<
        set -euo pipefail

        python -m combine_requantify_tools \
            filter-gtf-and-tracking \
            --tracking-in="~{updated_tracking}" \
            --squanti-classification="~{squanti_classification}" \
            --gtf-in="~{combined_gtf}" \
            --prefix="~{prefix}" \
            --tracking-out="~{sample_set_id}_updated_tracking_sq_filtered.parquet" \
            --gtf-out="~{sample_set_id}_filtered.gtf"

        echo "Recombining GTFs"
        cat "~{sample_set_id}_filtered.gtf" "~{gencode_gtf}" > recombined.gtf

        # restrict to contigs present in the reference genome
        echo "Filtering GTF contigs"
        grep '^>' ~{ref_fasta} | cut -d ' ' -f1 | sed 's/^>//' > contigs.txt
        awk 'NR==FNR {contigs[$1]; next} $1 in contigs' contigs.txt recombined.gtf \
            > "~{sample_set_id}_filtered.gtf"

        echo "Creating transcriptome FASTA"
        gffread "~{sample_set_id}_filtered.gtf" \
            -g "~{ref_fasta}" \
            -w "~{sample_set_id}_transcriptome.fa"
    >>>

    output {
        File transcriptome_fasta = "~{sample_set_id}_transcriptome.fa"
        File filtered_gtf = "~{sample_set_id}_filtered.gtf"
        File updated_tracking_sq_filtered = "~{sample_set_id}_updated_tracking_sq_filtered.parquet"
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

task process_gtf {
    meta {
        description: "TODO"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        sample_set_id: "identifier for this set of samples"
        filtered_gtf: "filtered_gtf file from filter_gtf_and_tracking"
        gencode_gtf: "Gencode GTF file (if this is the initial run), or a combined_sorted_gtf (from a previous run)"
        ref_fasta: "reference sequence FASTA"

        # outputs
        sorted_gtf: "TODO"
    }

    input {
        String sample_set_id
        File filtered_gtf
        File gencode_gtf
        File ref_fasta

        String docker_image  = "us-central1-docker.pkg.dev/depmap-omics/terra-images/combine-requantify-tools"
        String docker_image_hash_or_tag = ":production"
        Int cpu = 2
        Int mem_gb = 32
        Int preemptible = 2
        Int max_retries = 0
        Int additional_disk_gb = 0
   }

    Int disk_space = (
        ceil(
            size(gencode_gtf, "GiB")
            + size(filtered_gtf, "GiB")
            + size(gencode_gtf, "GiB")
    )  + 20 + additional_disk_gb)

    command <<<
        set -euo pipefail

        python -m combine_requantify_tools \
            process-gtf \
            --gtf-in="~{filtered_gtf}" \
            --gtf-out="~{sample_set_id}_processed.gtf"
    >>>

    output {
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
