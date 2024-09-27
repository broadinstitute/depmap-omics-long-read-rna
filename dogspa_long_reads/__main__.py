import json
import logging
import os
from pathlib import Path
from typing import Annotated, Any

import pandas as pd
import tomllib
import typer
from nebelung.terra_workflow import TerraWorkflow
from nebelung.terra_workspace import TerraWorkspace

from dogspa_long_reads.types import GumboClient
from dogspa_long_reads.utils.bams import (
    do_delta_index_delivery_bams,
    do_upsert_delivery_bams,
)
from dogspa_long_reads.utils.onboarding import do_onboard_samples

pd.set_option("display.max_columns", 30)
pd.set_option("display.max_colwidth", 50)
pd.set_option("display.max_info_columns", 30)
pd.set_option("display.max_info_rows", 20)
pd.set_option("display.max_rows", 20)
pd.set_option("display.max_seq_items", None)
pd.set_option("display.width", 200)
pd.set_option("expand_frame_repr", True)
pd.set_option("mode.chained_assignment", "warn")

app = typer.Typer()

config: dict[str, Any] = {}


# noinspection PyUnusedLocal
def done(*args, **kwargs):
    logging.info("Done.")


def make_workflow_from_config(
    repo_namespace: str, workflow_config: dict
) -> TerraWorkflow:
    """
    Make a TerraWorkflow object from a config entry.

    :param repo_namespace: the method repo namespace to store methods
    :param workflow_config: a dictionary with workflow/method config values
    :return: a TerraWorkflow instance
    """

    return TerraWorkflow(
        repo_namespace=repo_namespace,
        repo_method_name=workflow_config["repo_method_name"],
        method_config_name=workflow_config["method_config_name"],
        method_synopsis=workflow_config["method_synopsis"],
        workflow_wdl_path=Path(workflow_config["workflow_wdl_path"]).resolve(),
        method_config_json_path=Path(
            workflow_config["method_config_json_path"]
        ).resolve(),
    )


@app.callback(result_callback=done)
def main(
    ctx: typer.Context,
    config_path: Annotated[Path, typer.Option(exists=True)],
):
    logger = logging.getLogger()
    logger.addHandler(logging.StreamHandler())
    logger.setLevel(logging.INFO)

    with open(config_path, "rb") as f:
        config.update(tomllib.load(f))

    ctx.obj = {
        "terra_workspace": TerraWorkspace(
            workspace_namespace=config["terra"]["workspace_namespace"],
            workspace_name=config["terra"]["workspace_name"],
            owners=json.loads(os.environ["FIRECLOUD_OWNERS"]),
        ),
        "gumbo_client": GumboClient(
            url=os.environ["HASURA_URL"],
            username="snp_str_qc",
            headers={"X-Hasura-Admin-Secret": os.environ["HASURA_ADMIN_SECRET"]},
        ),
    }


@app.command()
def update_workflow(
    ctx: typer.Context, workflow_name: Annotated[str, typer.Option()]
) -> None:
    terra_workflow = make_workflow_from_config(
        repo_namespace=config["terra"]["repo_namespace"],
        workflow_config=config["terra"][workflow_name],
    )
    ctx.obj["terra_workspace"].update_workflow(terra_workflow=terra_workflow)


@app.command()
def upsert_delivery_bams(ctx: typer.Context) -> None:
    do_upsert_delivery_bams(
        gcs_source_bucket=config["onboarding"]["gcs_source"]["bucket"],
        gcs_source_glob=config["onboarding"]["gcs_source"]["glob"],
        terra_workspace=ctx.obj["terra_workspace"],
    )


@app.command()
def delta_index_delivery_bams(ctx: typer.Context) -> None:
    terra_workflow = make_workflow_from_config(
        repo_namespace=config["terra"]["repo_namespace"],
        workflow_config=config["terra"]["index_bam"],
    )

    do_delta_index_delivery_bams(
        terra_workspace=ctx.obj["terra_workspace"], terra_workflow=terra_workflow
    )


@app.command()
def onboard_samples(ctx: typer.Context) -> None:
    do_onboard_samples(
        gcp_project_id=config["gcp_project_id"],
        gcs_destination_bucket=config["onboarding"]["gcs_destination"]["bucket"],
        gcs_destination_prefix=config["onboarding"]["gcs_destination"]["prefix"],
        uuid_namespace=config["onboarding"]["uuid_namespace"],
        terra_workspace=ctx.obj["terra_workspace"],
        gumbo_client=ctx.obj["gumbo_client"],
        dry_run=config["onboarding"]["dry_run"],
    )


if __name__ == "__main__":
    app()
