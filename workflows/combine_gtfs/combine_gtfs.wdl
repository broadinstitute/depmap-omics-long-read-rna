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

    call run_gffcompare {
        input:
            sample_set_id = sample_set_id,
            gtf_list = filter_isoquant.filtered_gtf,
            gencode_gtf = gencode_gtf,
            prefix = prefix
    }

    call run_sqanti3 {
        input:
            sample_set_id = sample_set_id,
            isoquant_gtf = run_gffcompare.combined_gtf,
            gencode_gtf = gencode_gtf,
            ref_fasta = ref_fasta
    }

    call process_tracking_file {
        input:
            sample_set_id = sample_set_id,
            tracking_file = run_gffcompare.tracking_file,
            sample_ids = sample_ids,
            discovered_transcript_counts = discovered_transcript_counts
   }

    call filter_gtf_and_tracking {
        input:
            sample_set_id = sample_set_id,
            combined_gtf = run_gffcompare.combined_gtf,
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

task run_gffcompare {
    meta {
        description: "TODO"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        sample_set_id: "identifier for this set of samples"
        gtf_list: "array of the filtered GTFs from filter_isoquant"
        gencode_gtf: "Gencode GTF file (if this is the initial run), or a combined_sorted_gtf (from a previous run)"
        prefix: "annotation for transcripts identified in this run (mostly for internal consistency between tasks)"

        # outputs
        combined_gtf: "combined GTF from all samples"
        tracking_file: "tracking file showing transcript relationships"
    }

    input {
        String sample_set_id
        Array[File] gtf_list
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
        echo "Running gffcompare with the following GTFs:"

        gffcompare -r "~{gencode_gtf}" -X -p "~{prefix}" -V -S -o gffcomp_out ~{sep=' ' gtf_list}

        awk -F'\t' '
            FNR==NR && $4 == "=" {
                split($3,a,"|"); new_id=a[2];
                split($1,b,"|"); old_id=b[1];
                map[old_id]=new_id;
                next
            }
            FNR != NR {
                split($1,c,"|"); prefix=c[1]; rest=(index($1,"|")?substr($1,index($1,"|")+1):"");
                if(prefix in map) $1=map[prefix]"|"rest;
                print $0
            }
        ' OFS='\t' gffcomp_out.tracking gffcomp_out.tracking > gffcomp_out.newtracking.tsv

        awk -F'\t' '
            FNR==NR {
                if ($4=="=") {
                    split($3,a,"|"); new_id=a[2];
                    split($1,b,"|"); old_id=b[1];
                    map[old_id]=new_id;
                    print "[DEBUG] Mapping added:", old_id, "->", new_id > "/dev/stderr"
                }
                next; # skip printing tracking lines
            }

            FNR!=NR {
                line=$0
                while (match(line, /transcript_id "[^"]+"/)) {
                    tid = substr(line, RSTART+15, RLENGTH-16)  # extract transcript_id value
                    if (tid in map) {
                        replacement = "transcript_id \"" map[tid] "\""
                        line = substr(line, 1, RSTART-1) replacement substr(line, RSTART+RLENGTH)
                        print "[DEBUG] Replacing transcript_id:", tid, "->", map[tid] > "/dev/stderr"
                    } else {
                        print "[DEBUG] transcript_id not in map:", tid > "/dev/stderr"
                        break
                    }
                }
                print line
            }
        ' gffcomp_out.tracking gffcomp_out.combined.gtf > renamed.gtf

        mv renamed.gtf "~{sample_set_id}_gffcompared.gtf"
        mv gffcomp_out.newtracking.tsv "~{sample_set_id}_gffcomp_out.newtracking.tsv"
    >>>

    output {
        File combined_gtf = "~{sample_set_id}_gffcompared.gtf"
        File tracking_file = "~{sample_set_id}_gffcomp_out.newtracking.tsv"
        File stats = "~{sample_set_id}_gffcomp_out.stats"
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
        isoquant_gtf: "combined GTF from run_gffcompare"
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
        tracking_file: "output tracking_file from run_gffcompare"
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
        combined_gtf: "Combined GTF from run_gffcompare"
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
            --tracking-out="~{sample_set_id}_updated_tracking_sq_filtered.tsv" \
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
        File updated_tracking_sq_filtered = "~{sample_set_id}_updated_tracking_sq_filtered.tsv"
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

        grep -E 'HAVANA|ENSEMBL' "~{filtered_gtf}" | \
            awk -F '\t' '
            {
                match($9, /gene_id "([^"]+)"/, gid);
                match($9, /gene_name "([^"]+)"/, gname);
                if (gid[1] && gname[1]) print gname[1] "\t" gid[1];
            }' | sort -u > gene_name_to_id.tsv

        awk -F '\t' -v OFS='\t' '
            BEGIN {
                # Build gene_name → gene_id map
                while ((getline < "gene_name_to_id.tsv") > 0) {
                    map[$1] = $2;
                }
            }

            {
                if ($0 ~ /^#/ || NF < 9) {
                    print; next;
                }

                # Track last gene_name for multi-line transcript models
                has_gid = match($9, /gene_id "([^"]+)"/, gid);
                has_gname = match($9, /gene_name "([^"]+)"/, gname);

                if (has_gname) {
                    current_gname = gname[1];
                }

                if ($2 == "IsoQuant" && has_gid) {
                    if (!has_gname && current_gname in map) {
                        # No gene_name in this line, but we have a remembered one
                        gsub("gene_id \"" gid[1] "\"", "gene_id \"" map[current_gname] "\"", $9);
                    } else if (has_gname && gname[1] in map) {
                        # Replace based on gene_name in this line
                        gsub("gene_id \"" gid[1] "\"", "gene_id \"" map[gname[1]] "\"", $9);
                    }
                }

                print;
            }' "~{filtered_gtf}" > "~{sample_set_id}_updated.gtf"
    >>>

    output {
        File sorted_gtf = "~{sample_set_id}_updated.gtf"
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
