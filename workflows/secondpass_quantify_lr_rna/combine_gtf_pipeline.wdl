version 1.0

workflow LongreadPipeline {
  input {
    Array[String] sample_ids          # List of sample identifiers
    Array[File] sample_gtfs           # List of GTF files from IsoQuant first pass
    String gffcompare_docker_image    # Docker image for gffcompare
    String sample_id                  # Sample identifier for SQANTI3 output
    String sqanti_docker_image        # SQANTI3 docker image name
    String sqanti_docker_image_hash_or_tag  # SQANTI3 docker image tag or hash
    File referenceGenome              # Reference genome FASTA (may be same as ref_fasta)
    File gencode_gtf                  # GENCODE annotation GTF
    Array[File] model_counts # counts values for transcript models
    String prefix #for gffcompare
    String process_tracking_gffread_docker_image 
    String process_gtf_docker_image 
  }
  
  Array[Pair[String, File]] zip_arrays = zip(sample_ids, sample_gtfs)

  scatter(pair in zip_arrays) {
    call FilterIsoQuant {
      input:
        raw_gtf = pair.right,
        sample_id = pair.left,
        docker_image = gffcompare_docker_image
    }

  }

  # Run gffcompare once on all filtered GTFs
  call RunGffcompare {
    input:
      ref_gtf = gencode_gtf,
      gtf_list = FilterIsoQuant.filtered_gtf,
      docker_image = gffcompare_docker_image,
      prefix = prefix 
  }

  call run_sqanti3 {
    input:
      sample_id = sample_id,
      isoquant_gtf = RunGffcompare.combined_gtf,  # Uses the combined GTF from step 1
      ref_annotation_gtf = gencode_gtf,
      ref_fasta = referenceGenome,
      docker_image = sqanti_docker_image,
      docker_image_hash_or_tag = sqanti_docker_image_hash_or_tag
  }

  call process_tracking_file {
    input:
      tracking_file = RunGffcompare.tracking_file,  # Uses tracking file from step 1
      sample_ids = sample_ids,
      model_counts = model_counts,
      docker_image = process_tracking_gffread_docker_image
  }

  call gffread {
    input:
      sample_id = sample_id,
      referenceGenome = referenceGenome,
      referenceAnnotation_full = RunGffcompare.combined_gtf,
      gencode_gtf = gencode_gtf,
      squanti_classification = run_sqanti3.sq_class,  # Uses SQANTI classification from step 2
      updated_tracking = process_tracking_file.updated_tracking,  # Uses processed tracking file
      docker_image = process_tracking_gffread_docker_image,
      prefix = prefix
  }

  call process_gtf {
    input:
      filtered_gtf = gffread.filtered_gtf,  # Uses filtered GTF from gffread
      sample_id = sample_id,
      referenceGenome = referenceGenome,
      referenceAnnotation_full = RunGffcompare.combined_gtf,
      gencode_gtf = gencode_gtf,
      docker_image = process_gtf_docker_image
  }

  output {
    File sq_class = run_sqanti3.sq_class
    File transcriptome_fasta = gffread.transcriptome_fasta
    File combined_sorted_gtf = process_gtf.sorted_gtf
    File updated_tracking_sq_filtered = gffread.updated_tracking_sq_filtered
  }
}

task FilterIsoQuant {
  input {
    File raw_gtf          # Raw GTF file from IsoQuant
    String sample_id      # Sample identifier
    String docker_image   # Docker image to use
  }

  command <<<
    echo "Filtering IsoQuant entries for sample: ~{sample_id}"

    # Decompress GTF file if needed
    if [[ "~{raw_gtf}" == *.gz ]]; then
      gunzip -c ~{raw_gtf} > raw.gtf
    else
      cp ~{raw_gtf} raw.gtf
    fi

    # Error handling for decompression issues
    if [[ ! -f raw.gtf ]]; then
      echo "ERROR: raw.gtf not found after decompression." >&2
      touch filtered.gtf
    else
      # Filter for lines where the source column (column 2) is "IsoQuant"
      awk '$2 == "IsoQuant"' raw.gtf > filtered.gtf
      if [[ ! -s filtered.gtf ]]; then
        echo "WARNING: No IsoQuant entries found for ~{sample_id}" >&2
        touch filtered.gtf
      fi
    fi
  >>>

  output {
    File filtered_gtf = "filtered.gtf"  # GTF file containing only IsoQuant entries
  }

  runtime {
    cpu: 1
    memory: "1G"
    docker: docker_image
  }
}


task RunGffcompare {
  input {
    File ref_gtf                # Reference GTF file
    Array[File] gtf_list        # List of filtered GTF files to compare
    String docker_image         # Docker image for gffcompare 
    String prefix
  }

  command <<<

    echo "Running gffcompare with the following GTFs:"
  
    gffcompare -r ~{ref_gtf} -X -p ~{prefix} -V -S -o gffcomp_out ~{sep=' ' gtf_list}

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


  >>>

  output {
    File combined_gtf = "renamed.gtf"  # Combined GTF from all samples
    File tracking_file = "gffcomp_out.newtracking"     # Tracking file showing transcript relationships
    File stats_file = "gffcomp_out.stats"           # Statistics about the comparison
  }

  runtime {
    cpu: 2
    memory: "64G"
    disks: "local-disk 500 SSD"
    docker: docker_image
  }
}

# Step 2 Task - SQANTI3 Quality Control and Classification
task run_sqanti3 {
  input {
    String sample_id              # Sample identifier
    File isoquant_gtf             # Combined GTF from gffcompare
    File ref_annotation_gtf       # Reference annotation GTF
    File ref_fasta                # Reference genome FASTA
    String docker_image           # SQANTI3 docker image name
    String docker_image_hash_or_tag # SQANTI3 docker image tag or hash
    Int cpu = 2                   # CPU cores to use
    Int mem_gb = 8                # Memory in GB
    Int preemptible = 2           # Number of preemptible attempts
    Int max_retries = 1           # Maximum number of retries
    Int additional_disk_gb = 0    # Additional disk space in GB
  }

  # Calculate required disk space based on input file sizes
  Int disk_space = (
    ceil(
      size(isoquant_gtf, "GiB") * 3 + size(ref_annotation_gtf, "GiB") + size(ref_fasta, "GiB")
    ) + 20 + additional_disk_gb
  )

  command <<<
    set -euo pipefail

    awk '{ if ($7 != ".") print }' ~{isoquant_gtf} > "~{sample_id}.unzipped.gtf"

    python /usr/local/src/SQANTI3-5.3.6/sqanti3_qc.py \
      --report both \
      --output ~{sample_id} \
      ~{sample_id}.unzipped.gtf \
      ~{ref_annotation_gtf} \
      ~{ref_fasta} 
  >>>

  output {
    File sq_class = "~{sample_id}_classification.txt"           #Transcript classification
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

# Step 3 Tasks
task gffread {
  input {
    String sample_id
    File referenceGenome
    File gencode_gtf
    File referenceAnnotation_full
    File squanti_classification
    File updated_tracking
    String prefix 
    Int cpu = 2
    Int memoryGB = 8
    Int diskSizeGB = 20
    String docker_image = "docker.io/hannharris/salmon-with-pandas@sha256:29adeea24ffb7d396e06c0bd7fd3c7ad258d65bd8fb9188f59e2c536023446f3"
  }

  command <<<
    set -euo pipefail

    # Step 1: Filter out features without a strand (column 7 = '.')
    if file ~{referenceAnnotation_full} | grep -q 'gzip compressed'; then
      zcat ~{referenceAnnotation_full} | awk '{ if ($7 != ".") print }' > annotation_filtered.gtf
    else
      awk '{ if ($7 != ".") print }' ~{referenceAnnotation_full} > annotation_filtered.gtf
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
    

    updated_tracking.to_csv("updated_tracking_sq_filtered.tsv", sep='\t', index=False)

    # Load GTF
    gtf = pd.read_csv("annotation_filtered.gtf", sep="\\t", comment='#', header=None)

    # Extract transcript_id from attributes column
    gtf['transcript_id'] = gtf[8].str.extract(r'transcript_id "([^"]+)"')

    # Filter GTF
    gtf_filtered = gtf[gtf['transcript_id'].isin(sq_filtered['isoform'])].copy()
    gtf_filtered.drop(columns=['transcript_id'], inplace=True)
    #gtf_filtered.drop(columns=['num_samples'], inplace=True)

    #gtf_filtered.drop(columns=['num_samples'], inplace=True)


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

    # Step 2: Concatenate and sort with gffread
    cat filtered1.gtf ~{gencode_gtf} > combined.gtf 

    cp combined.gtf sorted.gtf
    # Step 3: Restrict to contigs present in the reference genome
    grep '^>' ~{referenceGenome} | cut -d ' ' -f1 | sed 's/^>//' > contigs.txt
    awk 'NR==FNR {contigs[$1]; next} $1 in contigs' contigs.txt sorted.gtf > filtered.gtf

    gffread filtered.gtf -g ~{referenceGenome} -w ~{sample_id}_transcriptome.fa
  >>>

  output {
    File transcriptome_fasta = "~{sample_id}_transcriptome.fa"
    File filtered_gtf = "filtered.gtf"
    File updated_tracking_sq_filtered = "updated_tracking_sq_filtered.tsv"
  }

  runtime {
    cpu: "~{cpu}"
    memory: "~{memoryGB} GiB"
    disks: "local-disk ~{diskSizeGB} HDD"
    docker: "~{docker_image}"
    errorStrategy: "Continue"
  }
}

task process_gtf {
  input {
    String sample_id
    File referenceGenome
    File gencode_gtf
    File referenceAnnotation_full
    File filtered_gtf
    Int cpu = 2
    Int memoryGB = 8
    Int diskSizeGB = 20
    String docker_image = "quay.io/biocontainers/gawk:5.3.1"
  }

  command <<<
    set -euo pipefail
    grep -E 'HAVANA|ENSEMBL' ~{filtered_gtf} | \
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
    }
    ' ~{filtered_gtf} > updated.gtf
  >>>

  output {
    File sorted_gtf = "updated.gtf"
  }

  runtime {
    cpu: "~{cpu}"
    memory: "~{memoryGB} GiB"
    disks: "local-disk ~{diskSizeGB} HDD"
    docker: "~{docker_image}"
    errorStrategy: "Continue"
  }
}

task process_tracking_file {
  input {
    File tracking_file
    Array[String] sample_ids
    Array[File] model_counts
    Int cpu = 2
    Int memoryGB = 8
    Int diskSizeGB = 20
    String docker_image = "docker.io/hannharris/salmon-with-pandas@sha256:29adeea24ffb7d396e06c0bd7fd3c7ad258d65bd8fb9188f59e2c536023446f3"
    # String docker = "docker.io/hannharris/salmon-with-pandas@sha256:29adeea24ffb7d396e06c0bd7fd3c7ad258d65bd8fb9188f59e2c536023446f3"
  }

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
            # No TPM file for this sample, skip                                                                      
            print(f"Warning: No TPM file found for sample {sample}")                                                 
            continue                                                                                                 
                                                                                                                     
        model_counts_path = sample_to_tpm[sample]                                                                    
        print(f"Processing sample {sample} with TPM file {model_counts_path}")                                       
                                                                                                                     
        model_counts = pd.read_csv(model_counts_path, sep='\t', comment='#', header=None, compression='gzip')        
        model_counts = model_counts.rename(columns={0: 'id1', 1: 'tpm'})                                             
                                                                                                                     
        tracking_m_mask = tracking_m[tracking_m['sample'] == sample]                                                 
        tracking_m_mask = tracking_m_mask.merge(model_counts, on='id1', how='left')                                  
        tracking_m_mask = tracking_m_mask[tracking_m_mask['tpm'] >= 5]                                               
        tracking_m_mask['sample'] = sample                                                                           
        updated_tracking = pd.concat([updated_tracking, tracking_m_mask], ignore_index=True)                         
                                                                                                                     
    # Print summary of processed samples                                                                             
    sample_counts = updated_tracking['sample'].value_counts()                                                        
    print("Transcript counts per sample:")                                                                           
    for sample, count in sample_counts.items():                                                                      
        print(f"{sample}: {count} transcripts")                                                                       
                                                                                                                     
    updated_tracking.to_csv("updated_tracking.tsv", sep='\t', index=False) 
    
    EOF
    >>>

  output {
    File updated_tracking = "updated_tracking.tsv"
  }

  runtime {
    cpu: "~{cpu}"
    memory: "~{memoryGB} GiB"
    disks: "local-disk ~{diskSizeGB} HDD"
    docker: "~{docker_image}"
    errorStrategy: "Continue"
  }
}

