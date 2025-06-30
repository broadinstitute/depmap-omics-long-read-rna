#!/bin/zsh

set -euo pipefail

uv pip compile pyproject.toml -o requirements.txt --emit-index-url > requirements.txt

gcloud functions deploy depmap-omics-long-read-rna \
  --gen2 \
  --runtime="python313" \
  --region="us-central1" \
  --source=. \
  --run-service-account="omics-pipeline-runner@depmap-omics.iam.gserviceaccount.com" \
  --entry-point="run" \
  --trigger-topic="run-depmap-omics-long-read-rna" \
  --timeout=540 \
  --memory="4GB" \
  --cpu=1 \
  --max-instances 1 \
  --concurrency 1
