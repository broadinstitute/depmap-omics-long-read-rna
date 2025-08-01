version 1.0

workflow align_lr_rna {
    input {
        String sample_id
        File input_bam
        File? junc_bed = "gs://ccleparams/long-read/gencode.v38.primary_assembly.annotation.bed"
        File ref_fasta = "gs://ccleparams/hg38ref_no_alt/GRCh38_no_alt.fa"
    }

    call minimap2 {
        input:
            sample_id = sample_id,
            input_bam = input_bam,
            junc_bed = junc_bed,
            ref_fasta = ref_fasta
    }

    output {
        File aligned_bam = minimap2.aligned_bam
        File aligned_bai = minimap2.aligned_bai
    }
}

task minimap2 {
    input {
        String sample_id
        File input_bam
        File? junc_bed
        File ref_fasta

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/minimap2_lr"
        String docker_image_hash_or_tag = ":production"
        Int cpu = 4
        Int mem_gb = 32
        Int preemptible = 1
        Int max_retries = 1
        Int additional_disk_gb = 0
    }

    Int disk_space = (
        ceil(size(input_bam, "GiB") * 20 + size(ref_fasta, "GiB")) +
        10 + additional_disk_gb
    )
    Int index_threads = cpu - 1

    command <<<
        set -euo pipefail

        juncbed_arg=~{if defined(junc_bed) then '"--junc-bed ${junc_bed}"' else '""'}

        samtools fastq \
            -@ ~{cpu} \
            "~{input_bam}" \
            > temp.fastq && rm "~{input_bam}"

        minimap2 \
            -y \
            -ax splice:hq \
            -uf ${juncbed_arg} \
            -t ~{cpu} \
            "~{ref_fasta}" \
            temp.fastq \
            > temp.sam && rm temp.fastq

        samtools sort \
            -@ ~{cpu} \
            temp.sam \
            > "~{sample_id}.bam" && rm temp.sam

        samtools index \
            -b \
            -@ ~{index_threads} \
            "~{sample_id}.bam"

        mv "~{sample_id}.bam.bai" "~{sample_id}.bai"
    >>>

    output {
        File aligned_bam = "~{sample_id}.bam"
        File aligned_bai = "~{sample_id}.bai"
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
