import base64
import json
import logging
import os
from pathlib import Path

import functions_framework
import google.cloud.logging
from cloudevents.http import CloudEvent
from dotenv import load_dotenv
from nebelung.terra_workflow import TerraWorkflow
from nebelung.terra_workspace import TerraWorkspace

from dogspa_long_reads.types import GumboClient
from dogspa_long_reads.utils.bams import (
    do_delta_align_delivery_bams,
    do_upsert_delivery_bams,
)
from dogspa_long_reads.utils.onboarding import do_onboard_samples


@functions_framework.cloud_event
def run(cloud_event: CloudEvent) -> None:
    """
    Wrapper around the primary `entrypoint` function (needed for remote execution
    inside a GCP Function).

    :param cloud_event: the pub/sub CloudEvent payload
    """

    client = google.cloud.logging.Client()
    client.setup_logging(log_level=logging.INFO)

    config = json.loads(base64.b64decode(cloud_event.data["message"]["data"]).decode())

    logging.debug("Full CloudEvent data:")
    logging.debug(json.dumps(config, indent=4, sort_keys=True))

    # try to load secrets as ENV variables from attached Secrets Manager volume
    try:
        load_dotenv("/etc/secrets/env")
    except Exception as e:
        # we don't expect this to work if running locally
        logging.warning(f"Couldn't load attached secrets: {e}")

    terra_workspace = TerraWorkspace(
        workspace_namespace=config["terra"]["workspace_namespace"],
        workspace_name=config["terra"]["workspace_name"],
    )

    gumbo_client = GumboClient(
        url=os.environ["HASURA_URL"],
        username="snp_str_qc",
        headers={"X-Hasura-Admin-Secret": os.environ["HASURA_ADMIN_SECRET"]},
    )

    if config["cmd"] == "index-delivery-bams":
        do_upsert_delivery_bams(
            gcs_source_bucket=config["onboarding"]["gcs_source"]["bucket"],
            gcs_source_glob=config["onboarding"]["gcs_source"]["glob"],
            terra_workspace=config["terra_workspace"],
        )

        terra_workflow = TerraWorkflow(
            repo_namespace=config["terra"]["repo_namespace"],
            repo_method_name=config["terra"]["index_bam"]["repo_method_name"],
            method_config_name=config["terra"]["index_bam"]["method_config_name"],
            method_synopsis=config["terra"]["index_bam"]["method_synopsis"],
            workflow_wdl_path=Path(
                config["terra"]["index_bam"]["workflow_wdl_path"]
            ).resolve(),
            method_config_json_path=Path(
                config["terra"]["index_bam"]["method_config_json_path"]
            ).resolve(),
        )

        do_delta_align_delivery_bams(terra_workspace=terra_workspace)

    elif config["cmd"] == "onboard-samples":
        do_onboard_samples(
            gcp_project_id=config["gcp_project_id"],
            gcs_destination_bucket=config["gcs_destination_bucket"],
            gcs_destination_prefix=config["gcs_destination_prefix"],
            uuid_namespace=config["uuid_namespace"],
            terra_workspace=terra_workspace,
            gumbo_client=gumbo_client,
        )
    else:
        raise NotImplementedError(f"Invalid command: {config['cmd']}")

    logging.info("Done.")
