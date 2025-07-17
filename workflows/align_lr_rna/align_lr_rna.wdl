version 1.0

workflow align_lr_rna {
    input {
        String sample_id
        File input_bam
        File? junc_bed
        File ref_fasta
    }

    call minimap2 {
        input:
            output_basename = sample_id,
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
        String output_basename
        File input_bam
        File? junc_bed
        File ref_fasta

        String docker_image
        String docker_image_hash_or_tag
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
            > temp.fastq

        minimap2 \
            -y \
            -ax splice:hq \
            -uf ${juncbed_arg} \
            -t ~{cpu} \
            "~{ref_fasta}" \
            temp.fastq \
            > temp.sam

        samtools sort \
            -@ ~{cpu} \
            temp.sam \
            > "~{output_basename}.bam"

        samtools index \
            -b \
            -@ ~{index_threads} \
            "~{output_basename}.bam"

        mv "~{output_basename}.bam.bai" "~{output_basename}.bai"
    >>>

    output {
        File aligned_bam = "~{output_basename}.bam"
        File aligned_bai = "~{output_basename}.bai"
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
