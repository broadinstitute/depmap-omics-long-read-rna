DOGSPA: Depmap Omics Gumbo Syncing, Pre-processing, and Attestation (Long Reads)
---

This repo is a minor variation on [dogspa](https://github.com/broadinstitute/dogspa) for onboarding long-read RNA samples. Changes include:

1. Delivered BAM files are sourced from a Terra workspace's GCS bucket, but we're directly searching the bucket rather than using URLs found in a workspace data table.
2. Some "legacy" BAMs that were delivered to a slightly different location in that bucket. The `scratch/seed_legacy.py` script was usee to store those paths in this repo's `data` folder and the onboarding code appends these to the other BAMs found so that we aren't missing data during the initial run.
3. Lots of dogspa options are changed to hard-coded values.
4. In addition to uploading data to Gumbo, we update the `sample` data table in the Terra workspace with metadata and paths to the newly-copied BAMs in the `cclebams` bucket. 
