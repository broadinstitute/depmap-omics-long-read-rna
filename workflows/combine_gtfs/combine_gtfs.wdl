version 1.0

workflow combine_gtfs {
    parameter_meta {
        # inputs
        sample_set_id: "identifier for this set of samples"
        sample_ids: "array of the sample IDs in this sample set"
        extended_annotation: "array of extended_annotation output files from quantify_lr_rna"
        model_counts: "array of model_counts output files from quantify_lr_rna"
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
        Array[File] model_counts
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
            model_counts = model_counts
   }

    call gffread {
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
            filtered_gtf = gffread.filtered_gtf,
            gencode_gtf = gencode_gtf,
            ref_fasta = ref_fasta
   }

    output {
        File combined_sorted_gtf = process_gtf.sorted_gtf
        File sq_class = run_sqanti3.sq_class
        File transcriptome_fasta = gffread.transcriptome_fasta
        File updated_tracking_sq_filtered = gffread.updated_tracking_sq_filtered
    }
}

task filter_isoquant {
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

        String docker_image  = "us-central1-docker.pkg.dev/depmap-omics/terra-images/python-pandas-gffcompare-gffutils-gawk"
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

    meta {
        allowNestedInputs: true
    }
}


task run_gffcompare {
    parameter_meta {
        # inputs
        sample_set_id: "identifier for this set of samples"
        gtf_list: "array of the filtered GTFs from filter_isoquant"
        gencode_gtf: "Gencode GTF file (if this is the initial run), or a combined_sorted_gtf (from a previous run)"
        prefix: "annotation for transcripts identified in this run (mostly for internal consistency between tasks)"

        # outputs
        combined_gtf: "combined GTF from all samples"
        tracking_file: "tracking file showing transcript relationships"
        stats: "statistics about the comparison"
    }

    input {
        String sample_set_id
        Array[File] gtf_list
        File gencode_gtf
        String prefix

        String docker_image  = "us-central1-docker.pkg.dev/depmap-omics/terra-images/python-pandas-gffcompare-gffutils-gawk"
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
        ' OFS='\t' gffcomp_out.tracking gffcomp_out.tracking > gffcomp_out.newtracking

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

        mv renamed.gtf "~{sample_set_id}_renamed.gtf"
        mv gffcomp_out.newtracking "~{sample_set_id}_gffcomp_out.newtracking"
        mv gffcomp_out.stats "~{sample_set_id}_gffcomp_out.stats"
    >>>

    output {
        File combined_gtf = "~{sample_set_id}_renamed.gtf"  # Combined GTF from all samples
        File tracking_file = "~{sample_set_id}_gffcomp_out.newtracking"     # Tracking file showing transcript relationships
        File stats = "~{sample_set_id}_gffcomp_out.stats"           # Statistics about the comparison
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
        Int mem_gb = 16
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

    meta {
        allowNestedInputs: true
    }
}

task process_tracking_file {
    parameter_meta {
        # inputs
        sample_set_id: "identifier for this set of samples"
        tracking_file: "output tracking_file from run_gffcompare"
        sample_ids: "array of the sample IDs in this sample set"
        model_counts: "array of model_counts output files from quantify_lr_rna"

        # outputs
        updated_tracking: "TODO"
    }

    input {
        String sample_set_id
        File tracking_file
        Array[String] sample_ids
        Array[File] model_counts

        String docker_image  = "us-central1-docker.pkg.dev/depmap-omics/terra-images/python-pandas-gffcompare-gffutils-gawk"
        String docker_image_hash_or_tag = ":production"
        Int cpu = 2
        Int mem_gb = 8
        Int preemptible = 2
        Int additional_disk_gb = 0
        Int max_retries = 0
    }

    Int disk_space = (
        ceil(size(tracking_file, "GiB") * 3 + size(model_counts, "GiB"))
        + 20 + additional_disk_gb
    )

    command <<<
        set -euo pipefail

        python3 <<EOF
        import pandas as pd

        tracking = pd.read_csv("~{tracking_file}", sep='\t', header=None)

        sample_cols = "~{sep=',' sample_ids}".split(',')

        fixed_cols = ['transcript_id', 'loc', 'gene_id', 'val']
        all_cols = fixed_cols + sample_cols
        tracking.columns = all_cols

        tracking_m = tracking.melt(id_vars=fixed_cols, var_name='sample', value_name='id')
        tracking_m = tracking_m[tracking_m['id'] != "-"]
        tracking_m['id1'] = tracking_m['id'].str.split('|').str[1]

        updated_tracking = pd.DataFrame()

        sample_ids = "~{sep=',' sample_ids}".split(',')

        def extract_sample_id(path):
            filename = path.split('/')[-1]
            #take the part before '.discovered_transcript_tpm.tsv.gz'
            sample_id = filename.split(".")[0]
            return sample_id

        transcript_model_tpms = "~{sep=',' model_counts}".split(',')

        sample_to_tpm = {extract_sample_id(path): path for path in transcript_model_tpms}

        # Process all samples, not just the first one
        for sample in sample_ids:
            if sample not in sample_to_tpm:
                # No count file for this sample, skip
                print(f"Warning: No count file found for sample {sample}")
                continue

            model_counts_path = sample_to_tpm[sample]
            print(f"Processing sample {sample} with TPM file {model_counts_path}")

            model_counts = pd.read_csv(model_counts_path, sep='\t', comment='#', header=None, compression='gzip')
            model_counts = model_counts.rename(columns={0: 'id1', 1: 'count'})

            tracking_m_mask = tracking_m[tracking_m['sample'] == sample]
            tracking_m_mask = tracking_m_mask.merge(model_counts, on='id1', how='left')
            tracking_m_mask = tracking_m_mask[tracking_m_mask['count'] >= 5]
            tracking_m_mask['sample'] = sample
            updated_tracking = pd.concat([updated_tracking, tracking_m_mask], ignore_index=True)

        sample_counts = updated_tracking['sample'].value_counts()
        print("Transcript counts per sample:")
        for sample, count in sample_counts.items():
            print(f"{sample}: {count} transcripts")

        updated_tracking.to_csv("~{sample_set_id}_updated_tracking.tsv", sep='\t', index=False)
        EOF
    >>>

    output {
        File updated_tracking = "~{sample_set_id}_updated_tracking.tsv"
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

task gffread {
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

        String docker_image  = "us-central1-docker.pkg.dev/depmap-omics/terra-images/python-pandas-gffcompare-gffutils-gawk"
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

        # Step 1: Filter out features without a strand (column 7 = '.')
        if file "~{combined_gtf}" | grep -q 'gzip compressed'; then
            zcat "~{combined_gtf}" | awk '{ if ($7 != ".") print }' > annotation_filtered.gtf
        else
            awk '{ if ($7 != ".") print }' "~{combined_gtf}" > annotation_filtered.gtf
        fi

        python3 <<EOF
        import pandas as pd

        updated_tracking = pd.read_csv("~{updated_tracking}", sep="\\t")
        updated_tracking_nodups = updated_tracking[['transcript_id']].drop_duplicates()
        updated_tracking_nodups['transcript_id'] = updated_tracking_nodups['transcript_id'].str.split('|').str[0]

        # Load classification
        sq = pd.read_csv("~{squanti_classification}", sep="\\t")
        sq_annotated = sq[sq['isoform'].str.startswith('ENST')]

        # Filter TCONS + novel_*_catalog + coding
        sq_filtered = sq[
            sq['isoform'].str.startswith("~{prefix}") &
            sq['structural_category'].isin(['novel_not_in_catalog', 'novel_in_catalog']) &
            (sq['coding'] == 'coding')
        ]
        sq_annotated_sm = sq_annotated[['ORF_seq', 'isoform']]
        merged_sq = sq_filtered.merge(sq_annotated_sm, on='ORF_seq', how='left', suffixes=('_tcons', '_enst'), indicator=True)
        merged_sq = merged_sq[merged_sq['_merge'] == 'left_only'] #orf in novel and not in annotated

        sq_filtered = sq_filtered[sq_filtered['isoform'].isin(merged_sq['isoform_tcons'])]
        sq_filtered = sq_filtered[sq_filtered['RTS_stage'] == False]  # Filter for RTS_stage
        sq_filtered = sq_filtered[sq_filtered['isoform'].isin(updated_tracking_nodups['transcript_id'])]

        sq_filtered_previouslyfound = sq[
            sq['structural_category'].isin(['full-splice_match']) &
            ~((sq['associated_transcript'].str.startswith('ENST')) | (sq['associated_transcript'] == "novel")) &
            (sq['coding'] == 'coding')
        ]

        sq_filtered_tracking = pd.concat([sq_filtered, sq_filtered_previouslyfound])

        updated_tracking['transcript_id'] = updated_tracking['transcript_id'].str.split('|').str[0]
        updated_tracking = updated_tracking[updated_tracking['transcript_id'].isin(sq_filtered_tracking['isoform'])]


        updated_tracking.to_csv("~{sample_set_id}_updated_tracking_sq_filtered.tsv", sep='\t', index=False)

        # Load GTF
        gtf = pd.read_csv("annotation_filtered.gtf", sep="\\t", comment='#', header=None)

        # Extract transcript_id from attributes column
        gtf['transcript_id'] = gtf[8].str.extract(r'transcript_id "([^"]+)"')

        # Filter GTF
        gtf_filtered = gtf[gtf['transcript_id'].isin(sq_filtered['isoform'])].copy()
        gtf_filtered.drop(columns=['transcript_id'], inplace=True)

        # Write output
        gtf_filtered.to_csv(
            "filtered1.gtf", sep='\\t',
            index=False, header=False, quoting=3
        )

        # Add GTF header
        header = """##gff-version 3
        ##description: evidence-based annotation of the human genome (GRCh38), version 38 (Ensembl 104)
        ##provider: GENCODE
        ##contact: gencode-help@ebi.ac.uk
        ##format: gtf
        ##date: 2021-03-12"""

        with open("filtered1.gtf", "r") as f:
            lines = f.readlines()
        with open("filtered1.gtf", "w") as f:
            f.write(header + "\\n")
            f.writelines(lines)
        EOF

        cat filtered1.gtf "~{gencode_gtf}" > combined.gtf

        cp combined.gtf sorted.gtf
        # Step 3: Restrict to contigs present in the reference genome
        grep '^>' ~{ref_fasta} | cut -d ' ' -f1 | sed 's/^>//' > contigs.txt
        awk 'NR==FNR {contigs[$1]; next} $1 in contigs' contigs.txt sorted.gtf \
            > "~{sample_set_id}_filtered.gtf"

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

    meta {
        allowNestedInputs: true
    }
}

task process_gtf {
    parameter_meta {
        # inputs
        sample_set_id: "identifier for this set of samples"
        filtered_gtf: "filtered_gtf file from gffread"
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

        String docker_image  = "us-central1-docker.pkg.dev/depmap-omics/terra-images/python-pandas-gffcompare-gffutils-gawk"
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

    meta {
        allowNestedInputs: true
    }
}
