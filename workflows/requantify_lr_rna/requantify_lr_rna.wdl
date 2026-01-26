version 1.0

workflow requantify_lr_rna {
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
            gtf_file = filter_gtf_per_sample.filtered_gtf
    }

    call run_isoquant {
        input:
            sample_id = sample_id,
            input_bam = input_bam,
            input_bai = input_bai,
            ref_annotation_db = gtf_to_db.db_file,
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
        File requantified_transcript_counts = run_isoquant.transcript_counts
        File requantified_transcript_tpm = run_isoquant.transcript_tpm
        File requantified_read_assignments_tsv = run_isoquant.read_assignments_tsv
        File requantified_exon_counts = run_isoquant.exon_counts
        File requantified_intron_counts = run_isoquant.intron_counts
    }
}


task filter_gtf_per_sample {
    input {
        String sample_id
        File transcript_counts
        File combined_sorted_gtf
        File updated_tracking_sq_filtered

        String docker_image  = "us-central1-docker.pkg.dev/depmap-omics/terra-images/python-pandas-gffcompare-gffutils-gawk"
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
            + size(combined_sorted_gtf, "GiB")
            + size(updated_tracking_sq_filtered, "GiB")
        ) + 20 + additional_disk_gb
    )

    command <<<
        set -euo pipefail

        python3 <<EOF
        import pandas as pd

        # Filter tracking file for this sample
        updated_tracking = pd.read_csv("~{updated_tracking_sq_filtered}", sep="\\t")
        updated_tracking = updated_tracking[updated_tracking['sample'] == "~{sample_id}"]
        updated_tracking = updated_tracking['transcript_id']

        # Get transcripts with non-zero counts
        transcript_counts = pd.read_csv("~{transcript_counts}", sep='\t', comment='#', header=None, compression='gzip')
        transcript_counts = transcript_counts[transcript_counts[1] > 0] # Filter for counts > 0
        transcript_counts = transcript_counts[0]
        transcript_counts = transcript_counts.rename('transcript_id')

        # Combine transcripts from tracking and counts
        concat_tx = pd.concat([transcript_counts, updated_tracking], axis=0)
        transcript_list = concat_tx.drop_duplicates().tolist()

        # Parse GTF file
        gtf = pd.read_csv("~{combined_sorted_gtf}", sep='\t', comment='#', header=None)
        gtf['transcript_id'] = gtf[8].str.extract(r'transcript_id "([^"]+)"')
        gtf['gene_id'] = gtf[8].str.extract('gene_id "([^"]+)"')

        # Get gene IDs for the transcripts we want to keep
        match_tx_gene_id = gtf[gtf['transcript_id'].isin(transcript_list)]
        match_tx_gene_id = match_tx_gene_id[['transcript_id', 'gene_id']].drop_duplicates()

        # Get gene entries for these transcripts
        filtered_gene = gtf[(gtf[2] == 'gene') & (gtf['gene_id'].isin(match_tx_gene_id['gene_id']))]

        # Get transcript entries
        gtf_filtered = gtf[gtf['transcript_id'].isin(transcript_list)]

        # Combine gene and transcript entries
        final_gtf = pd.concat([filtered_gene, gtf_filtered])

        # Write filtered GTF
        final_gtf.drop(columns=['transcript_id']).drop(columns=['gene_id']).to_csv(
            "~{sample_id}_filtered_sample.gtf", sep='\\t',
            index=False, header=False, quoting=3
        )
        EOF
    >>>

    output {
        File filtered_gtf = "~{sample_id}_filtered_sample.gtf"  # Sample-specific filtered GTF
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

task gtf_to_db {
    input {
        File gtf_file              # Filtered GTF file from previous step

        String docker_image  = "us-central1-docker.pkg.dev/depmap-omics/terra-images/python-pandas-gffcompare-gffutils-gawk"
        String docker_image_hash_or_tag = ":production"
        Int cpu = 2
        Int mem_gb = 16
        Int preemptible = 1
        Int max_retries = 0
        Int additional_disk_gb = 0
   }

    String db_name = basename(gtf_file, ".gtf")

    Int disk_space = (
        ceil(size(gtf_file, "GiB")) * 3 + 20 + additional_disk_gb
    )

    command <<<
        set -euo pipefail

        python3 <<EOF
        import gffutils

        # Create a database from the GTF file
        db = gffutils.create_db(
            "~{gtf_file}",
            dbfn="~{db_name}.db",
            force=True,
            keep_order=False,
            merge_strategy="merge",
            disable_infer_transcripts=True,
            disable_infer_genes=True,
            id_spec={"transcript": "transcript_id", "gene": "gene_id"}
        )
        EOF
    >>>

    output {
        File db_file = "~{db_name}.db"  # SQLite database of GTF features
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

task run_isoquant {
    input {
        String sample_id                # Sample identifier
        File input_bam                  # Input BAM file
        File input_bai                  # BAM index file
        File ref_annotation_db          # Reference annotation database (from GTFtoDB)
        File ref_fasta                  # Reference genome FASTA
        String data_type
        String model_construction_strategy
        String stranded
        String transcript_quantification
        String gene_quantification
        String report_novel_unspliced
        String report_canonical

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
        File read_assignments_tsv = "isoquant_output/~{sample_id}/~{sample_id}.read_assignments.tsv.gz"
        File transcript_counts = "isoquant_output/~{sample_id}/~{sample_id}.transcript_counts.tsv.gz"
        File transcript_tpm = "isoquant_output/~{sample_id}/~{sample_id}.transcript_tpm.tsv.gz"
        File exon_counts = "isoquant_output/~{sample_id}/~{sample_id}.exon_counts.tsv.gz"
        File intron_counts = "isoquant_output/~{sample_id}/~{sample_id}.intron_counts.tsv.gz"
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
