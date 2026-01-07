#!/bin/zsh

set -e

GENCODE_VERSION="38"

# download the compressed GTF
curl "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_$GENCODE_VERSION/gencode.v$GENCODE_VERSION.primary_assembly.annotation.gtf.gz" \
  -o "/tmp/gencode.v$GENCODE_VERSION.primary_assembly.annotation.gtf.gz"

# make the .db
uv run gffutils-cli \
  create \
  --disable-infer-genes \
  --disable-infer-transcripts \
  --output \
  "/tmp/gencode.v$GENCODE_VERSION.primary_assembly.annotation.db" \
  "/tmp/gencode.v$GENCODE_VERSION.primary_assembly.annotation.gtf.gz"

# make the .junc.bed
# see https://github.com/lh3/minimap2/blob/master/misc/README.md
k8-Darwin workflows/common/paftools.js gff2bed \
  "/tmp/gencode.v$GENCODE_VERSION.primary_assembly.annotation.gtf.gz" \
  > "/tmp/gencode.v$GENCODE_VERSION.primary_assembly.annotation.junc.bed"

# make the .gtf
gzip -kd "/tmp/gencode.v$GENCODE_VERSION.primary_assembly.annotation.gtf.gz"

gcloud storage cp \
  "/tmp/gencode.v$GENCODE_VERSION.primary_assembly.annotation.gtf" \
  "/tmp/gencode.v$GENCODE_VERSION.primary_assembly.annotation.db" \
  "/tmp/gencode.v$GENCODE_VERSION.primary_assembly.annotation.junc.bed" \
  "gs://ccleparams/long-read/gencode/v38/"
