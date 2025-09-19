version 1.0

workflow call_lr_rna_fusions {
    input {
        String sample_id
        File input_bam
        Int? max_n_reads

        # optional matched short read inputs
        String? sr_sample_id
        String? sr_cram_or_bam
        File? sr_cram_bam
        File? sr_crai_bai
        File? ref_fasta = "gs://ccleparams/hg38ref_no_alt/GRCh38_no_alt.fa"
        File? ref_fasta_index = "gs://ccleparams/hg38ref_no_alt/GRCh38_no_alt.fa.fai"

        # STAR inputs
        File genome_lib_tar_with_star_idx = "gs://mdl-ctat-genome-libs/__genome_libs_StarFv1.10/GRCh38_gencode_v22_CTAT_lib_Mar012021.for_mm2fusion.withSTARidx.tar"
        Int min_per_id = 70
        Int min_j = 1
        Int min_sum_js = 1
        Int min_novel_junction_support = 1
        String? fi_extra_params
    }

    if (defined(sr_cram_or_bam)) {
        call sam_to_fastq as sr_sam_to_fastq {
            input:
                sample_id = select_first([sr_sample_id]),
                cram_or_bam = select_first([sr_cram_or_bam]),
                cram_bam = select_first([sr_cram_bam]),
                crai_bai = sr_crai_bai,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                max_n_reads = max_n_reads
        }

        call ctat_lr_fusion as ctat_lr_fusion_with_sr {
            input:
                sample_id = sample_id,
                input_bam = input_bam,
                genome_lib_tar = genome_lib_tar_with_star_idx,
                sr_fastq1 = sr_sam_to_fastq.fastq1,
                sr_fastq2 = sr_sam_to_fastq.fastq2,
                sr_cram_bam_size_gb = size(sr_cram_bam, "GiB"),
                min_per_id = min_per_id,
                min_j = min_j,
                min_sum_js = min_sum_js,
                min_novel_junction_support = min_novel_junction_support,
                fi_extra_params = fi_extra_params
        }
    }

    call ctat_lr_fusion as ctat_lr_fusion_without_sr {
        input:
            sample_id = sample_id,
            input_bam = input_bam,
            genome_lib_tar = genome_lib_tar_with_star_idx,
            min_per_id = min_per_id,
            min_j = min_j,
            min_sum_js = min_sum_js,
            min_novel_junction_support = min_novel_junction_support,
            fi_extra_params = fi_extra_params
    }

    output {
        File fusion_report = select_first([
            ctat_lr_fusion_with_sr.fusion_report,
            ctat_lr_fusion_without_sr.fusion_report
        ])
        File fusion_report_abridged = select_first([
            ctat_lr_fusion_with_sr.fusion_report_abridged,
            ctat_lr_fusion_without_sr.fusion_report_abridged
        ])
        File prelim_fusion_report = select_first([
            ctat_lr_fusion_with_sr.prelim_fusion_report,
            ctat_lr_fusion_without_sr.prelim_fusion_report
        ])
        File prelim_fusion_report_abridged = select_first([
            ctat_lr_fusion_with_sr.prelim_fusion_report_abridged,
            ctat_lr_fusion_without_sr.prelim_fusion_report_abridged
        ])
        File? fusion_report_html = select_first([
            ctat_lr_fusion_with_sr.fusion_report_html,
            ctat_lr_fusion_without_sr.fusion_report_html
        ])
        File? igv_tar = select_first([
            ctat_lr_fusion_with_sr.igv_tar,
            ctat_lr_fusion_without_sr.igv_tar
        ])
        Boolean used_sr_evidence = defined(sr_cram_bam)
        String? sr_cram_bam_used = sr_cram_bam
    }
}

task sam_to_fastq {
    input {
        String sample_id
        String cram_or_bam
        File cram_bam
        File? crai_bai
        File? ref_fasta
        File? ref_fasta_index
        Int? max_n_reads

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/samtools_seqtk"
        String docker_image_hash_or_tag = ":production"
        Int cpu = 4
        Int mem_gb = 8
        Int preemptible = 2
        Int max_retries = 1
        Int additional_disk_gb = 0
    }

    Float bam_size_gb = if cram_or_bam == "BAM" then (
        size(cram_bam, "GiB")
    ) else 2 * size(cram_bam, "GiB")

    Int disk_space = ceil(bam_size_gb * 10) + 20 + additional_disk_gb
    Int n_threads = cpu - 1

    command <<<
        set -euo pipefail

        if [[ "~{cram_or_bam}" == "CRAM" ]];
        then
            samtools fastq \
                -@ ~{n_threads} \
                --reference "~{ref_fasta}" \
                -1 "~{sample_id}.1.fastq" \
                -2 "~{sample_id}.2.fastq" \
                "~{cram_bam}"
        else
            samtools fastq \
                -@ ~{n_threads} \
                -1 "~{sample_id}.1.fastq" \
                -2 "~{sample_id}.2.fastq" \
                "~{cram_bam}"
        fi

        if [[ -n "${max_n_reads:-}" ]]; then
            N_READS=$(awk 'END {print NR/4}' "~{sample_id}.1.fastq")

            if (( $(echo "$N_READS > ~{max_n_reads}" | bc -l) )); then
                echo "Downsampling $N_READS reads to ~{max_n_reads}"

                seqtk sample -2 -s100 "~{sample_id}.1.fastq" ~{max_n_reads} \
                    > "~{sample_id}.1.less.fastq" \
                    && rm "~{sample_id}.1.fastq" \
                    && mv "~{sample_id}.1.less.fastq" "~{sample_id}.1.fastq"

                seqtk sample -2 -s100 "~{sample_id}.2.fastq" ~{max_n_reads} \
                    > "~{sample_id}.2.less.fastq" \
                    && rm "~{sample_id}.2.fastq" \
                    && mv "~{sample_id}.2.less.fastq" "~{sample_id}.2.fastq"
            else
                echo "$N_READS <= ~{max_n_reads} reads (no downsampling)"
            fi
        fi
    >>>

    output {
        File fastq1 = "~{sample_id}.1.fastq"
        File fastq2 = "~{sample_id}.2.fastq"
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

task ctat_lr_fusion {
    input {
        String sample_id
        File input_bam
        File genome_lib_tar
        File? sr_fastq1
        File? sr_fastq2
        Float? sr_cram_bam_size_gb
        Int min_per_id
        Int min_j
        Int min_sum_js
        Int min_novel_junction_support
        String? fi_extra_params
        Boolean no_ctat_mm2 = false
        Boolean vis = true

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/ctat_lr_fusion"
        String docker_image_hash_or_tag = ":1.1.0"
        Int cpu = 16
        Int mem_gb = 64
        Int preemptible = 1
        Int max_retries = 1
        Int additional_disk_gb = 0
    }

    Int disk_space = (
        ceil(
            size(input_bam, "GiB")
            + 4 * size(genome_lib_tar, "GiB")
            + size(sr_fastq1, "GiB")
            + size(sr_fastq2, "GiB")
        ) + 20 + additional_disk_gb
    )

    Int bam_sort_ram_gb = ceil(mem_gb * 0.85)

    String no_ctat_mm2_flag = if (no_ctat_mm2) then "--no_ctat_mm2" else ""
    String vis_flag = if (vis) then "--vis" else ""

    command <<<
        set -euo pipefail

        # untar the genome lib
        tar xf "~{genome_lib_tar}"
        rm "~{genome_lib_tar}"

        # ctat-LR-fusion
        ctat-LR-fusion \
            --LR_bam "~{input_bam}" \
            --genome_lib_dir ctat_genome_lib_build_dir \
            --min_J ~{min_j} \
            --min_sumJS ~{min_sum_js} \
            --min_novel_junction_support ~{min_novel_junction_support} \
            --min_per_id ~{min_per_id} \
            --CPU ~{cpu} \
            --STAR_xtra_params '--limitBAMsortRAM ~{bam_sort_ram_gb}000000000' \
            ~{"--left_fq " + sr_fastq1} \
            ~{"--right_fq " + sr_fastq2 } \
            ~{no_ctat_mm2_flag} \
            ~{vis_flag} \
            ~{"--FI_extra_params " + fi_extra_params } \
            -o out

        mv out/ctat-LR-fusion.fusion_predictions.preliminary.tsv \
            "~{sample_id}.ctat-LR-fusion.fusion_predictions.preliminary.tsv"
        mv out/ctat-LR-fusion.fusion_predictions.preliminary.abridged.tsv \
            "~{sample_id}.ctat-LR-fusion.fusion_predictions.preliminary.abridged.tsv"
        mv out/ctat-LR-fusion.fusion_predictions.tsv \
            "~{sample_id}.ctat-LR-fusion.fusion_predictions.tsv"
        mv out/ctat-LR-fusion.fusion_predictions.abridged.tsv \
            "~{sample_id}.ctat-LR-fusion.fusion_predictions.abridged.tsv"

        if [ "~{vis}" = "true" ]; then
            mv out/fusion_intermediates_dir/IGV_prep/igv.genome.fa \
                "~{sample_id}.ctat-LR-fusion.igv.genome.fa"
            mv out/fusion_intermediates_dir/IGV_prep/igv.genome.fa.fai \
                "~{sample_id}.ctat-LR-fusion.igv.genome.fa.fai"
            mv out/fusion_intermediates_dir/IGV_prep/igv.annot.gtf \
                "~{sample_id}.ctat-LR-fusion.igv.annot.gtf"
            mv out/fusion_intermediates_dir/IGV_prep/igv.annot.bed \
                "~{sample_id}.ctat-LR-fusion.igv.annot.bed"
            mv out/fusion_intermediates_dir/IGV_prep/igv.LR.sorted.bam \
                "~{sample_id}.ctat-LR-fusion.igv.LR.sorted.bam"
            mv out/fusion_intermediates_dir/IGV_prep/igv.LR.sorted.bam.bai \
                "~{sample_id}.ctat-LR-fusion.igv.LR.sorted.bam.bai"
            mv out/fusion_intermediates_dir/IGV_prep/igv.pfam.bed \
                "~{sample_id}.ctat-LR-fusion.igv.pfam.bed"
            mv out/fusion_intermediates_dir/IGV_prep/igv.seqsimilar.bed \
                "~{sample_id}.ctat-LR-fusion.igv.seqsimilar.bed"
            mv out/fusion_intermediates_dir/IGV_prep/igv.LR.breakoint.roi.bed \
                "~{sample_id}.ctat-LR-fusion.igv.LR.breakoint.roi.bed"
            mv out/ctat-LR-fusion.fusion_inspector_web.html \
                "~{sample_id}.ctat-LR-fusion.fusion_inspector_web.html"

            tar -zchf ~{sample_id}.ctat-LR-fusion.igv.tar.gz ~{sample_id}.ctat-LR-fusion.igv.*
        fi
    >>>

    output {
        File fusion_report = "~{sample_id}.ctat-LR-fusion.fusion_predictions.tsv"
        File fusion_report_abridged = "~{sample_id}.ctat-LR-fusion.fusion_predictions.abridged.tsv"
        File prelim_fusion_report = "~{sample_id}.ctat-LR-fusion.fusion_predictions.preliminary.tsv"
        File prelim_fusion_report_abridged = "~{sample_id}.ctat-LR-fusion.fusion_predictions.preliminary.abridged.tsv"
        File? fusion_report_html = "~{sample_id}.ctat-LR-fusion.fusion_inspector_web.html"
        File? igv_tar = "~{sample_id}.ctat-LR-fusion.igv.tar.gz"
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
