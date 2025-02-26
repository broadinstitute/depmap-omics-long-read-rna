version 1.0

workflow call_fusions {
    input {
        String workflow_version = "1.0" # internal semver
        String workflow_source_url # populated automatically with URL of this script

        String sample_id
        File sr_bam
        File ref_fasta
        File input_bam
        File genome_lib_tar_with_star_idx
        Int min_per_id = 70
        Int min_j = 1
        Int min_sum_js = 1
        Int min_novel_junction_support = 1
        String? fi_extra_params
    }

    call sam_to_fastq as sr_sam_to_fastq {
        input:
            sample_id = sample_id,
            input_bam = sr_bam,
            ref_fasta = ref_fasta
    }

    call ctat_lr_fusion {
        input:
            sample_id = sample_id,
            input_bam = input_bam,
            genome_lib_tar = genome_lib_tar_with_star_idx,
            sr_fastq1 = sr_sam_to_fastq.fastq1,
            sr_fastq2 = sr_sam_to_fastq.fastq2,
            min_per_id = min_per_id,
            min_j = min_j,
            min_sum_js = min_sum_js,
            min_novel_junction_support = min_novel_junction_support,
            fi_extra_params = fi_extra_params
    }

    output {
        File fusion_report = ctat_lr_fusion.fusion_report
        File fusion_report_abridged = ctat_lr_fusion.fusion_report_abridged
        File prelim_fusion_report = ctat_lr_fusion.prelim_fusion_report
        File prelim_fusion_report_abridged = ctat_lr_fusion.prelim_fusion_report_abridged
        File fusion_report_html = ctat_lr_fusion.fusion_report_html
        File igv_tar = ctat_lr_fusion.igv_tar
    }
}

task sam_to_fastq {
    input {
        String sample_id
        File input_bam
        File ref_fasta

        String docker_image
        String docker_image_hash_or_tag
        Int cpu = 4
        Int mem_gb = 16
        Int preemptible = 2
        Int max_retries = 1
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(size(input_bam, "GiB") * 6) + additional_disk_gb

    command <<<
        set -euo pipefail

        samtools fastq \
            -@ ~{cpu} \
            -1 ~{sample_id}.1.fastq \
            -2 ~{sample_id}.2.fastq \
            ~{input_bam}
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
        Int min_per_id
        Int min_j
        Int min_sum_js
        Int min_novel_junction_support
        String? fi_extra_params
        Boolean no_ctat_mm2 = false

        String docker_image
        String docker_image_hash_or_tag
        Int cpu = 12
        Int mem_gb = 64
        Int preemptible = 2
        Int max_retries = 1
        Int additional_disk_gb = 0
    }

    Int disk_space = (
        ceil(
            size(input_bam, "GiB")
            + 3 * size(genome_lib_tar, "GiB")
            + size(sr_fastq1, "GiB")
            + size(sr_fastq2, "GiB")
        ) + 20 + additional_disk_gb
    )

    Int bam_sort_ram_bytes = mem_gb * 1000 * 1000 * 1000 * 0.75

    String no_ctat_mm2_flag = if (no_ctat_mm2) then "--no_ctat_mm2" else ""

    command <<<
        set -euo pipefail

        # untar the genome lib
        tar xvf ~{genome_lib_tar}
        rm ~{genome_lib_tar}

        # ctat-LR-fusion
        ctat-LR-fusion \
            --LR_bam ~{input_bam} \
            --genome_lib_dir ctat_genome_lib_build_dir \
            --min_J ~{min_j} \
            --min_sumJS ~{min_sum_js} \
            --min_novel_junction_support ~{min_novel_junction_support} \
            --min_per_id ~{min_per_id} \
            --CPU ~{cpu} \
            --STAR_xtra_params '--limitBAMsortRAM ~{bam_sort_ram_bytes}' \
            --vis \
            ~{"--left_fq " + sr_fastq1} \
            ~{"--right_fq " + sr_fastq2 } \
            -o out \
            ~{no_ctat_mm2_flag} \
            ~{"--FI_extra_params " + fi_extra_params }

        mv out/ctat-LR-fusion.fusion_predictions.preliminary.tsv ~{sample_id}.ctat-LR-fusion.fusion_predictions.preliminary.tsv
        mv out/ctat-LR-fusion.fusion_predictions.preliminary.abridged.tsv ~{sample_id}.ctat-LR-fusion.fusion_predictions.preliminary.abridged.tsv
        mv out/ctat-LR-fusion.fusion_predictions.tsv ~{sample_id}.ctat-LR-fusion.fusion_predictions.tsv
        mv out/ctat-LR-fusion.fusion_predictions.abridged.tsv ~{sample_id}.ctat-LR-fusion.fusion_predictions.abridged.tsv
        mv out/ctat-LR-fusion.fusion_inspector_web.html ~{sample_id}.ctat-LR-fusion.fusion_inspector_web.html
        mv out/fusion_intermediates_dir/IGV_prep/igv.genome.fa ~{sample_id}.ctat-LR-fusion.igv.genome.fa
        mv out/fusion_intermediates_dir/IGV_prep/igv.genome.fa.fai ~{sample_id}.ctat-LR-fusion.igv.genome.fa.fai
        mv out/fusion_intermediates_dir/IGV_prep/igv.annot.gtf ~{sample_id}.ctat-LR-fusion.igv.annot.gtf
        mv out/fusion_intermediates_dir/IGV_prep/igv.annot.bed ~{sample_id}.ctat-LR-fusion.igv.annot.bed
        mv out/fusion_intermediates_dir/IGV_prep/igv.LR.sorted.bam ~{sample_id}.ctat-LR-fusion.igv.LR.sorted.bam
        mv out/fusion_intermediates_dir/IGV_prep/igv.LR.sorted.bam.bai ~{sample_id}.ctat-LR-fusion.igv.LR.sorted.bam.bai
        mv out/fusion_intermediates_dir/IGV_prep/igv.pfam.bed ~{sample_id}.ctat-LR-fusion.igv.pfam.bed
        mv out/fusion_intermediates_dir/IGV_prep/igv.seqsimilar.bed ~{sample_id}.ctat-LR-fusion.igv.seqsimilar.bed
        mv out/fusion_intermediates_dir/IGV_prep/igv.LR.breakoint.roi.bed ~{sample_id}.ctat-LR-fusion.igv.LR.breakoint.roi.bed

        tar -zcvhf ~{sample_id}.ctat-LR-fusion.igv.tar.gz ~{sample_id}.ctat-LR-fusion.igv.*
    >>>

    output {
        File fusion_report = "~{sample_id}.ctat-LR-fusion.fusion_predictions.tsv"
        File fusion_report_abridged = "~{sample_id}.ctat-LR-fusion.fusion_predictions.abridged.tsv"
        File prelim_fusion_report = "~{sample_id}.ctat-LR-fusion.fusion_predictions.preliminary.tsv"
        File prelim_fusion_report_abridged = "~{sample_id}.ctat-LR-fusion.fusion_predictions.preliminary.abridged.tsv"
        File fusion_report_html = "~{sample_id}.ctat-LR-fusion.fusion_inspector_web.html"
        File igv_tar = "~{sample_id}.ctat-LR-fusion.igv.tar.gz"
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
