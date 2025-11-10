version 1.0

workflow PostProcessingPipeline {
  input {
    String sample_id                
    File all_sample_gtf             
    File tracking_file              
    File transcript_counts          
    String output_basename = "new_annotation"  
    File input_bam                  
    File input_bai                  
    File ref_fasta                  
    File? ref_annotation_db          
    String isoquant_docker_image    
    String isoquant_docker_image_hash_or_tag  
    String filter_gtf_per_sample_docker = "docker.io/hannharris/salmon-with-pandas@sha256:29adeea24ffb7d396e06c0bd7fd3c7ad258d65bd8fb9188f59e2c536023446f3"
    String GTFtoDB_docker = "quay.io/biocontainers/gffutils:0.13--pyh7cba7a3_0"

  }

  call filter_gtf_per_sample {
    input:
      sample_id = sample_id,
      all_sample_gtf = all_sample_gtf,
      tracking_file = tracking_file,
      transcript_counts = transcript_counts,
      docker_image = filter_gtf_per_sample_docker
  }

  call GTFtoDB {
    input:
      gtf_file = filter_gtf_per_sample.filtered_gtf,  
      output_basename = output_basename, 
      docker_image = GTFtoDB_docker
  }

  call run_isoquant {
    input:
      sample_id = sample_id,
      input_bam = input_bam,
      input_bai = input_bai,
      ref_annotation_db = GTFtoDB.db_file,  
      ref_fasta = ref_fasta,
      docker_image = isoquant_docker_image,
      docker_image_hash_or_tag = isoquant_docker_image_hash_or_tag
  }

  output {

    File secondpass_transcript_counts = run_isoquant.transcript_counts
    File secondpass_transcript_tpm = run_isoquant.transcript_tpm
    File secondpass_read_assignments_tsv = run_isoquant.read_assignments_tsv
    File? secondpass_exon_counts = run_isoquant.exon_counts
    File? secondpass_intron_counts = run_isoquant.intron_counts
  }
}

# Step 4 Task - Filter GTF per sample
task filter_gtf_per_sample {
  input {
    String sample_id           
    File transcript_counts     
    File all_sample_gtf        
    File tracking_file         
    Int cpu = 2               
    Int memoryGB = 8           
    Int diskSizeGB = 20        
    String docker_image #= "docker.io/hannharris/salmon-with-pandas@sha256:29adeea24ffb7d396e06c0bd7fd3c7ad258d65bd8fb9188f59e2c536023446f3"
  }

  command <<<
    set -euo pipefail

    python3 <<EOF
    import pandas as pd

    # Filter tracking file for this sample
    updated_tracking = pd.read_csv("~{tracking_file}", sep="\\t")
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
    gtf = pd.read_csv("~{all_sample_gtf}", sep='\t', comment='#', header=None)
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
        "filtered_sample.gtf", sep='\\t',
        index=False, header=False, quoting=3
    )
    EOF
  >>>

  output {
    File filtered_gtf = "filtered_sample.gtf"  # Sample-specific filtered GTF
  }

  runtime {
    cpu: "~{cpu}"
    memory: "~{memoryGB} GiB"
    disks: "local-disk ~{diskSizeGB} HDD"
    docker: "~{docker_image}"
    errorStrategy: "Continue"
  }
}

# Step 5 Task - Convert GTF to Database
task GTFtoDB {
  input {
    File gtf_file              # Filtered GTF file from previous step
    String output_basename     # Base name for output database
    String docker_image
  }

  command <<<
    set -euo pipefail

    python3 <<EOF
    import gffutils

    # Create a database from the GTF file
    db = gffutils.create_db(
        "~{gtf_file}",
        dbfn="~{output_basename}.db",
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
    File db_file = "~{output_basename}.db"  # SQLite database of GTF features
  }

  runtime {
    docker: "~{docker_image}" #"quay.io/biocontainers/gffutils:0.13--pyh7cba7a3_0"
    memory: "16G"
    cpu: 1
  }
}

task run_isoquant {
  input {
    String sample_id                # Sample identifier
    File input_bam                  # Input BAM file
    File input_bai                  # BAM index file
    File ref_annotation_db          # Reference annotation database (from GTFtoDB)
    File ref_fasta                  # Reference genome FASTA
    
    String data_type = "pacbio_ccs"  
    String model_construction_strategy = "fl_pacbio"  
    String stranded = "forward"     
    String transcript_quantification = "unique_only"  
    String gene_quantification = "unique_splicing_consistent"  
    String report_novel_unspliced = "true"  
    String report_canonical = "auto"  

    # Runtime parameters
    String docker_image              
    String docker_image_hash_or_tag  
    Int cpu = 8                      
    Int mem_gb = 32                  
    Int preemptible = 1              
    Int max_retries = 1              
    Int additional_disk_gb = 0       
  }

  # Calculate required disk space based on input file sizes
  Int disk_space = (
    ceil(
      size(input_bam, "GiB") + size(ref_annotation_db, "GiB")
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

    # Compress output files to save space
    find isoquant_output/~{sample_id}/ -maxdepth 1 -type f -not -name '*.gz' -exec gzip {} +
  >>>

  output {
    File transcript_counts = "isoquant_output/~{sample_id}/~{sample_id}.transcript_counts.tsv.gz"  
    File transcript_tpm = "isoquant_output/~{sample_id}/~{sample_id}.transcript_tpm.tsv.gz"        
    File read_assignments_tsv = "isoquant_output/~{sample_id}/~{sample_id}.read_assignments.tsv.gz"  
    File? exon_counts = "isoquant_output/~{sample_id}/~{sample_id}.exon_counts.tsv.gz"             
    File? intron_counts = "isoquant_output/~{sample_id}/~{sample_id}.intron_counts.tsv.gz"         
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
