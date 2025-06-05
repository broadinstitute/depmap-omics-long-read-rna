#!/bin/zsh

set -euo pipefail

uv pip compile pyproject.toml -o requirements.txt --emit-index-url > requirements.txt

gcloud functions deploy depmap-omics-long-read \
  --gen2 \
  --runtime=python312 \
  --region=us-central1 \
  --source=. \
  --run-service-account=dogspa-runner@depmap-omics.iam.gserviceaccount.com \
  --entry-point=run \
  --trigger-topic=run-depmap-omics-long-read \
  --set-secrets=/etc/secrets:/env=projects/201811582504/secrets/dogspa-secrets/versions/latest \
  --timeout=540 \
  --memory=4GB \
  --cpu=4 \
  --max-instances 1 \
  --concurrency 1
