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
from depmap_omics_long_read_rna.utils.delivery_bams import do_upsert_delivery_bams
from depmap_omics_long_read_rna.utils.metadata import do_refresh_terra_samples
from depmap_omics_long_read_rna.utils.utils import (
    get_hasura_creds,
    make_workflow_from_config,
    submit_delta_job,
)


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

    terra_delivery_workspace = TerraWorkspace(
        workspace_namespace=config["terra"]["delivery_workspace_namespace"],
        workspace_name=config["terra"]["delivery_workspace_name"],
    )

    terra_workspace = TerraWorkspace(
        workspace_namespace=config["terra"]["workspace_namespace"],
        workspace_name=config["terra"]["workspace_name"],
    )

    gumbo_client = GumboClient(
        url=hasura_creds["url"],
        username="depmap-omics-long-read-rna",
        headers={"X-Hasura-Admin-Secret": hasura_creds["password"]},
    )

    if ce_data["cmd"] == "do-all":
        do_upsert_delivery_bams(
            gcs_source_bucket=config["alignment"]["gcs_source"]["bucket"],
            gcs_source_glob=config["alignment"]["gcs_source"]["glob"],
            uuid_namespace=config["uuid_namespace"],
            terra_workspace=terra_delivery_workspace,
            dry_run=config["onboarding"]["dry_run"],
        )

        submit_delta_job(
            terra_workspace=terra_delivery_workspace,
            terra_workflow=make_workflow_from_config(
                config, workflow_name="align_long_reads"
            ),
            entity_type="sample",
            entity_set_type="sample_set",
            entity_id_col="sample_id",
            expression="this.samples",
            dry_run=config["dry_run"],
        )

        do_refresh_terra_samples(
            terra_workspace=terra_workspace,
            short_read_terra_workspace=TerraWorkspace(
                workspace_namespace=config["terra"]["short_read_workspace_namespace"],
                workspace_name=config["terra"]["short_read_workspace_name"],
            ),
            gumbo_client=gumbo_client,
        )

        for workflow_name in ["quantify_long_reads", "call_fusions"]:
            submit_delta_job(
                terra_workspace=terra_workspace,
                terra_workflow=make_workflow_from_config(
                    config, workflow_name=workflow_name
                ),
                entity_type="sample",
                entity_set_type="sample_set",
                entity_id_col="sample_id",
                expression="this.samples",
                dry_run=config["dry_run"],
            )

    else:
        raise NotImplementedError(f"Invalid command: {ce_data['cmd']}")

    logging.info("Done.")
