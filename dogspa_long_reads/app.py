import base64
import json
import logging
import os
from collections import OrderedDict

import functions_framework
import pandas as pd
from cloudevents.http import CloudEvent
from dotenv import load_dotenv

from dogspa_long_reads.utils.gcp import (
    check_file_sizes,
    copy_to_cclebams,
    list_bams,
    update_sample_file_uris,
)
from dogspa_long_reads.utils.metadata import (
    apply_col_map,
    assign_cds_ids,
    check_already_in_gumbo,
    explode_and_expand_models,
    id_bams,
    increment_sample_versions,
    join_metadata,
    join_short_read_metadata,
    upload_to_gumbo,
)
from dogspa_long_reads.utils.terra import TerraWorkspace
from dogspa_long_reads.utils.utils import (
    model_to_df,
    send_slack_message,
)
from dogspa_long_reads.utils.validators import (
    DogspaConfig,
    ModelsAndChildren,
    PanderaBaseSchema,
)
from gumbo_gql_client import GumboClient


@functions_framework.cloud_event
def entrypoint(cloud_event: CloudEvent) -> None:
    """
    Run the DOGSPA sample loader.

    :param cloud_event: the pub/sub CloudEvent payload
    """

    # decode payload into the module config
    ce_data = json.loads(base64.b64decode(cloud_event.data["message"]["data"]).decode())
    config = DogspaConfig(**ce_data)

    logging.debug("Full CloudEvent data:")
    logging.debug(json.dumps(ce_data, indent=4, sort_keys=True))

    # try to load secrets as ENV variables from attached Secrets Manager volume
    try:
        load_dotenv("/etc/secrets/env")
    except Exception as e:
        # we don't expect this to work if running locally
        logging.warning(f"Couldn't load attached secrets: {e}")

    # keep track of filtering and quality control checks performed
    stats = OrderedDict()
    report = OrderedDict()

    tw = TerraWorkspace(config.workspace.namespace, config.workspace.name)

    # get the sequencing and profile tables from Gumbo
    gumbo_client = GumboClient(
        url=os.environ["HASURA_URL"],
        headers={"X-Hasura-Admin-Secret": os.environ["HASURA_ADMIN_SECRET"]},
    )
    models = model_to_df(gumbo_client.get_models_and_children(), ModelsAndChildren)
    seq_table = explode_and_expand_models(models)

    # get source BAM files
    src_bams = list_bams(config.gcs_source.bucket, config.gcs_source.prefix)
    samples = id_bams(src_bams)

    # get existing destination BAM files
    # dest_bams = list_bams(config.gcs_destination.bucket, config.gcs_destination.prefix)
    # dest_bams = id_bams(dest_bams)

    # remove existing objects
    # samples = src_bams.loc[~src_bams["crc32c"].isin(dest_bams["crc32c"])]

    # compare file sizes to filter out samples that are in Gumbo already
    samples = check_already_in_gumbo(samples, seq_table, size_col_name="legacy_size")
    stats["n not yet in Gumbo"] = (~samples["already_in_gumbo"]).sum()
    samples = samples.loc[~samples["already_in_gumbo"]].drop(columns="already_in_gumbo")
    report["not yet in Gumbo"] = samples

    if len(samples) == 0:
        send_slack_message(
            os.getenv("SLACK_WEBHOOK_URL_ERRORS"),
            os.getenv("SLACK_WEBHOOK_URL_STATS"),
            stats,
            report,
            config,
        )
        return

    # join metadata to current samples
    samples = join_metadata(samples, seq_table)
    samples = join_short_read_metadata(samples, seq_table)

    # check that BAM file sizes are above minimum threshold
    samples, blacklisted = check_file_sizes(samples)
    stats["n with BAM file too small"] = blacklisted.sum()
    report["BAM file too small"] = samples.loc[blacklisted]

    # assign sequencing IDs
    samples = assign_cds_ids(samples, config)

    # copy files to our own bucket
    sample_files = copy_to_cclebams(samples, config)

    # replace URIs with ones for our bucket whenever the copy operation succeeded
    samples, blacklisted = update_sample_file_uris(samples, sample_files)
    stats["n with failed file copies"] = blacklisted.sum()
    report["failed copies"] = samples.loc[blacklisted]

    # rename some columns for Gumbo
    samples = apply_col_map(samples, config)

    # increment version numbers for samples with profile IDs already in seq table
    samples = increment_sample_versions(samples, seq_table, config)

    # finally upload the samples to the Gumbo sequencing table
    samples = upload_to_gumbo(gumbo_client, samples, config)
    stats["n successfully uploaded"] = len(samples)
    report["successfully uploaded"] = samples

    send_slack_message(
        os.getenv("SLACK_WEBHOOK_URL_ERRORS"),
        os.getenv("SLACK_WEBHOOK_URL_STATS"),
        stats,
        report,
        config,
    )
    logging.info("Done.")
