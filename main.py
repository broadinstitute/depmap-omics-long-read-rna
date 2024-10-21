import base64
import json
import logging

import functions_framework
import google.cloud.logging
import tomllib
from cloudevents.http import CloudEvent
from dotenv import load_dotenv
from nebelung.terra_workspace import TerraWorkspace

from depmap_omics_long_read_rna.types import GumboClient
from depmap_omics_long_read_rna.utils.bams import (
    do_delta_align_delivery_bams,
    do_upsert_delivery_bams,
)
from depmap_omics_long_read_rna.utils.onboarding import (
    do_join_short_read_data,
    do_onboard_samples,
)
from depmap_omics_long_read_rna.utils.utils import get_hasura_creds


@functions_framework.cloud_event
def run(cloud_event: CloudEvent) -> None:
    """
    Wrapper around the primary `entrypoint` function (needed for remote execution
    inside a GCP Function).

    :param cloud_event: the pub/sub CloudEvent payload
    """

    client = google.cloud.logging.Client()
    client.setup_logging(log_level=logging.INFO)

    ce_data = json.loads(base64.b64decode(cloud_event.data["message"]["data"]).decode())

    logging.debug("Full CloudEvent data:")
    logging.debug(json.dumps(ce_data, indent=4, sort_keys=True))

    # try to load secrets as ENV variables from attached Secrets Manager volume
    try:
        load_dotenv("/etc/secrets/env")
    except Exception as e:
        # we don't expect this to work if running locally
        logging.warning(f"Couldn't load attached secrets: {e}")

    # use same config loading as when calling the module CLI
    with open(ce_data["config_path"], "rb") as f:
        config = tomllib.load(f)

    # get URL and password for Gumbo GraphQL API
    hasura_creds = get_hasura_creds("prod")

    terra_workspace = TerraWorkspace(
        workspace_namespace=config["terra"]["workspace_namespace"],
        workspace_name=config["terra"]["workspace_name"],
    )

    gumbo_client = GumboClient(
        url=hasura_creds["url"],
        username="depmap-omics-long-read-rna",
        headers={"X-Hasura-Admin-Secret": hasura_creds["password"]},
    )

    if ce_data["cmd"] == "onboard-samples":
        do_onboard_samples(
            gcp_project_id=config["gcp_project_id"],
            unaligned_gcs_destination_bucket=config["onboarding"][
                "unaligned_gcs_destination"
            ]["bucket"],
            unaligned_gcs_destination_prefix=config["onboarding"][
                "unaligned_gcs_destination"
            ]["prefix"],
            aligned_gcs_destination_bucket=config["onboarding"][
                "aligned_gcs_destination"
            ]["bucket"],
            aligned_gcs_destination_prefix=config["onboarding"][
                "aligned_gcs_destination"
            ]["prefix"],
            terra_workspace=terra_workspace,
            gumbo_client=gumbo_client,
            dry_run=config["onboarding"]["dry_run"],
        )

        do_upsert_delivery_bams(
            uuid_namespace=config["uuid_namespace"],
            gcs_source_bucket=config["onboarding"]["gcs_source"]["bucket"],
            gcs_source_glob=config["onboarding"]["gcs_source"]["glob"],
            terra_workspace=config["terra_workspace"],
        )

        do_delta_align_delivery_bams(terra_workspace=terra_workspace)

        short_read_terra_workspace = TerraWorkspace(
            workspace_namespace=config["terra"]["short_read_workspace_namespace"],
            workspace_name=config["terra"]["short_read_workspace_name"],
        )

        do_join_short_read_data(
            terra_workspace=terra_workspace,
            short_read_terra_workspace=short_read_terra_workspace,
            gumbo_client=gumbo_client,
        )

    else:
        raise NotImplementedError(f"Invalid command: {ce_data['cmd']}")

    logging.info("Done.")
