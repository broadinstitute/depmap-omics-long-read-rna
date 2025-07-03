import base64
import json
import logging
import tomllib

import functions_framework
import google.cloud.logging
from cloudevents.http import CloudEvent
from dotenv import load_dotenv
from nebelung.terra_workspace import TerraWorkspace

from depmap_omics_long_read_rna.types import DeltaJob, GumboClient
from depmap_omics_long_read_rna.utils.aligned_bams import onboard_aligned_bams
from depmap_omics_long_read_rna.utils.delivery_bams import upsert_delivery_bams
from depmap_omics_long_read_rna.utils.metadata import refresh_terra_samples
from depmap_omics_long_read_rna.utils.utils import (
    get_hasura_creds,
    make_workflow_from_config,
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

    short_read_terra_workspace = TerraWorkspace(
        workspace_namespace=config["terra"]["short_read_workspace_namespace"],
        workspace_name=config["terra"]["short_read_workspace_name"],
    )

    gumbo_client = GumboClient(
        url=hasura_creds["url"],
        username="depmap-omics-long-read-rna",
        headers={"X-Hasura-Admin-Secret": hasura_creds["password"]},
    )

    if ce_data["cmd"] == "upsert-delivery-bams":
        upsert_delivery_bams(
            gcs_source_bucket=config["alignment"]["gcs_source"]["bucket"],
            gcs_source_glob=config["alignment"]["gcs_source"]["glob"],
            uuid_namespace=config["uuid_namespace"],
            terra_workspace=terra_delivery_workspace,
            dry_run=config["onboarding"]["dry_run"],
        )

    elif ce_data["cmd"] == "onboard-aligned-bams":
        onboard_aligned_bams(
            terra_workspace=terra_delivery_workspace,
            gumbo_client=gumbo_client,
            dry_run=config["onboarding"]["dry_run"],
        )

    elif ce_data["cmd"] == "refresh-terra-samples":
        refresh_terra_samples(
            terra_workspace=terra_workspace,
            short_read_terra_workspace=short_read_terra_workspace,
            gumbo_client=gumbo_client,
        )

    elif ce_data["cmd"] == "submit-delta-job":
        failed_submissions = []

        for x in ce_data["delta_jobs"]:
            # iterate over workflow names and their delta job submission attrs
            dj = DeltaJob.model_validate(x)

            try:
                terra_workspace.submit_delta_job(
                    terra_workflow=make_workflow_from_config(
                        config, workflow_name=dj.workflow_name
                    ),
                    entity_type=dj.entity_type,
                    entity_set_type=dj.entity_set_type,
                    entity_id_col=dj.entity_id_col,
                    expression=dj.expression,
                    input_cols=dj.input_cols,
                    output_cols=dj.output_cols,
                    resubmit_n_times=dj.resubmit_n_times,
                    force_retry=dj.force_retry,
                    use_callcache=dj.use_callcache,
                    use_reference_disks=dj.use_reference_disks,
                    delete_intermediate_output_files=dj.delete_intermediate_output_files,
                    memory_retry_multiplier=dj.memory_retry_multiplier,
                    per_workflow_cost_cap=dj.per_workflow_cost_cap,
                    workflow_failure_mode=dj.workflow_failure_mode,
                    user_comment=dj.user_comment,
                    max_n_entities=dj.max_n_entities,
                    dry_run=dj.dry_run,
                )
            except Exception as e:
                logging.error(f"Couldn't submit delta job: {e}")
                failed_submissions.append(dj.workflow_name)

        if len(failed_submissions) > 0:
            raise Exception(
                "Couldn't submit delta jobs for workflows: "
                + ", ".join(failed_submissions)
            )

    else:
        raise NotImplementedError(f"Invalid command: {ce_data['cmd']}")

    logging.info("Done.")
