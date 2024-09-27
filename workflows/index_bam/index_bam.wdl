version 1.0

workflow index_bam {
    input {
        String workflow_version = "1.0" # internal semver
        String workflow_source_url # populated automatically with URL of this script

        File input_bam
    }

    call do_index_bam {
        input:
            input_bam = input_bam
    }

    output {
        File output_bai = do_index_bam.output_bai
    }
}

task do_index_bam {
    input {
        File input_bam

        String docker_image
        String docker_image_hash_or_tag
        Int cpu = 8
        Int preemptible = 2
        Int max_retries = 1
        Int additional_memory_mb = 0
        Int additional_disk_gb = 0
    }

    Int mem = ceil(size(input_bam, "MiB")) + 10000 + additional_memory_mb
    Int disk_space = ceil(size(input_bam, "GiB") * 1.5) + additional_disk_gb
    Int index_threads = cpu - 1
    String output_bai = basename(input_bam, ".bam") + ".bai"

    command {
        set -euo pipefail

        samtools index \
            -b \
            -@ ~{index_threads} \
            -o ~{output_bai} \
            ~{input_bam}
    }

    output {
        File output_bai = "~{output_bai}"
    }

    runtime {
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        memory: mem + " MiB"
        disks: "local-disk " + disk_space + " SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }

    meta {
        allowNestedInputs: true
    }
}
