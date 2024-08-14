DOGSPA: Depmap Omics Gumbo Syncing, Pre-processing, and Attestation (Long Reads)
---

This repository contains a single Python module to:

- copy sample data from a Terra workspace into Gumbo
- copy CRAM/CRAI/BAI files from the source GCS bucket to another one
- QC the metadata associated with the samples and their files in GCS
- upload new samples to the Gumbo sequencing table

It's deployed as a Google Cloud Function, which can be triggered by sending a message to the function's upstream pub/sub topic. Cloud Scheduler can then be used to trigger the topic on a schedule. 

# Architecture

The `entrypoint` function in `./dogspa/app.py` outlines all the quality control and data flow being performed:

1. Get `sample` table data from a Terra workspace.
2. Check for duplicate CRAM/CRAI (equivalently, BAI) GCS URIs across samples.
3. Get metadata (size, hash, update timestamp) about the CRAM files from GCS.
4. Check for missing CRAM/CRAI files.
5. Get the contents of the sequencing and profile tables from Gumbo.
6. Use CRAM file sizes as a key for checking which samples have already been uploaded to the Gumbo seq table. Existing samples are removed from the current batch (the corresponding rows in Gumbo **are not** replaced or updated).
7. Cross-reference sample SM-IDs to those in the profile table and assign a profile ID to each sample whenever the SM-ID uniquely identifies a profile.
8. Check samples for missing profile IDs.
9. Check whether CRAM file sizes are above the specified minimum.
10. Assign deterministic sequencing IDs (CDS-ZZZZZZ) by hashing relevant metadata.
11. Copy CRAM/CRAI files to our own GCS bucket, now with object keys determined by the CDS IDs.
12. Assign a version number to each sample based on how many samples already in the Gumbo seq table have the same profile ID.
13. Upload the samples to the Gumbo seq table.
14. Send a message to a Slack channel summarizing the execution results.

When a sample fails a quality control check or an operation on it fails (e.g. a BAM file was too small), the reason is recorded in the row's `issue` column and `blacklist` is set to `True`.

# Local development

## Installation

1. Install the required system dependencies:
   - [pyenv](https://github.com/pyenv/pyenv)
   - [Poetry](https://python-poetry.org/)
   - [gcloud CLI](https://cloud.google.com/sdk/gcloud)
   - [pre-commit](https://pre-commit.com/)
2. Install the required Python version (3.9.18):
	```shell
	pyenv install "$(cat .python-version)"
	```
3. Confirm that `python` maps to the correct version:
	```
	python --version
	```
4. Set the Poetry interpreter and install the Python dependencies:
	```shell
	poetry env use "$(pyenv which python)"
	poetry install
	```

## Credentials

Download a service account key for [dogspa-runner](https://console.cloud.google.com/iam-admin/serviceaccounts/details/112677449604656629898?project=depmap-omics) and set the `GOOGLE_APPLICATION_CREDENTIALS` environment variable to its location on your filesystem. This will simulate the permissions available inside the remote execution context. 

You can also set `SLACK_WEBHOOK_URL_ERRORS` and `SLACK_WEBHOOK_URL_STATS` environment variables if you're running with `dry_run=False` and want to send the errors and stats to Slack channels.

# Execution

The `entrypoint` function expects a [CloudEvent](https://cloud.google.com/eventarc/docs/workflows/cloudevents)-formatted message as its only argument.

## Local

The module's `__main__.py` file is configured for two example configurations:
```shell
poetry run python -m dogspa wgs
# or
poetry run python -m dogspa rna
```

## Remote

A GCP function expects a `main.py` function at the root of the source, so this exists as a thin wrapper around the `entrypoint` function.

# Deployment

The `deploy.sh` script will regenerate the `requirements.txt` file (GCP Functions don't know how to install dependencies with Poetry) and deploy the source as a Google Cloud function. 

To trigger remote execution of the function, send a message to the `run-dogspa` pub/sub topic:
```shell
gcloud pubsub topics publish run-dogspa --message='{}'
```

The `message` should be a JSON object conforming to the `DogspaConfig` Pydantic validation class:
```json5
{
    // data type of these samples ("rna" or "wgs"; only used to populate the `data_type`
    // column in Gumbo)
    "data_type": "rna",
  
    // a project name to use for the `source` column when uploading to Gumbo
    "source": "DEPMAP",
  
    // the Terra workspace we're checking for new samples
    "workspace": {
        "namespace": "terra-broad-cancer-prod",
        "name": "CCLE_DepMap_RNAseq"
    },
  
    // the column names in the Terra samples table representing BAI/CRAI URIs
    // (`bai_url`) and BAM URIs (`bam_url`)
    "terra_col_names": {
        "bai_url": "crai_or_bai_path",
        "bam_url": "cram_or_bam_path"
    },
  
    // the corresponding column names on the Gumbo side
    "gumbo_col_names": {
        "bai_url": "hg19_bam_filepath",
        "bam_url": "hg19_bai_filepath"
    },

    // (optional, default=false) upload new samples to Gumbo regardless of 
    // issue(s) identified instead of only the "bam too small" problem (for 
    // development purposes only; do not enable)
    "upload_all_failures": false,
      
    // minimum file size (in bytes) of the BAM file (if lower, upload to Gumbo as a
    // blacklisted sample)
    "min_file_size": 2000000000,
  
    // (optional, default=[]) manually ignore samples in Terra by `root_sample_id`
    "excluded_root_sample_ids": [],
  
    // the names of columns in the Gumbo profiles table that could contain SM-IDs (used
    // for mapping sequencings to profiles)
    "all_sm_id_col_names": [
        "smid_ordered",
        "smid_returned",
        "sm_id_matched"
    ],
  
    // GCP project this automation is running in 
    "gcp_project": "depmap-omics",
  
    // GCS bucket and path we're copying CRAM/CRAI/BAI objects to
    "gcs_destination": {
        "bucket": "cclebams",
        "prefix": "rna/"
    },
  
    // a UUIDv3 namespace to use for creating deterministic CDS-ZZZZZZ sequencing IDs
    // (arbitrary, but do not modify)
    "uuid_namespace": "00000000-0000-0000-0000-000000000000",
  
    // (optional, default=psutil.cpu_count()) number of CPUs to use in multiprocessing 
    "ncpus": 4,
  
    // (optional, default=true) don't actually do any copying of GCS objects, uploading
    // to Gumbo, or reporting of results to Slack
    "dry_run": true
}
```

Since invocation is asynchronous, you must use [Logs Explorer](https://console.cloud.google.com/logs/query;duration=PT1H;query=resource.type%3D%22cloud_run_revision%22%0Aresource.labels.service_name%3D%22dogspa%22?project=depmap-omics) to monitor progress.

## Development

Run `pre-commit run --all-files` to automatically format your code with [Ruff](https://docs.astral.sh/ruff/) and check static types with [Pyright](https://microsoft.github.io/pyright). There's also a hook to keep `requirements.txt` in sync with the current Poetry lockfile.

Whenever possible, function/method arguments and return values should be validated with Pydantic or Pandera (if a data frame).

### GraphQL code generation

This repo uses [ariadne-codegen](https://github.com/mirumee/ariadne-codegen) to generate the `gumbo_gql_client` module. It uses the folder of GraphQL queries (`./gql`) and the current GraphQL schema for a particular Gumbo environment to generate all of the Python classes, Pydantic models, and query/mutation methods for interacting with the [Gumbo GraphQL Service](https://github.com/broadinstitute/gumbo_client/tree/main/gumbo-gql-service). To regenerate the module using the current production schema:

```shell
HASURA_ADMIN_SECRET=... poetry run ariadne-codegen --config ariadne-prod.toml
```

The `HASURA_ADMIN_SECRET` variable should be the same as the secret/password for that endpoint's API (e.g. the [hasura-admin-secret-prod](https://console.cloud.google.com/security/secret-manager/secret/hasura-admin-secret-prod/versions?project=depmap-gumbo) secret). While developing or testing using a dev or staging Gumbo environment, change the `--config` option so that it uses one of the other available files.

### Testing

There's a minimal pytest suite available:
```shell
poetry run pytest
```
