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

from depmap_omics_long_read_rna.types import GumboClient
from depmap_omics_long_read_rna.utils.delivery_bams import (
    do_delta_align_delivery_bams,
    do_upsert_delivery_bams,
)
from depmap_omics_long_read_rna.utils.metadata import do_refresh_terra_samples
from depmap_omics_long_read_rna.utils.utils import get_hasura_creds, get_secret_from_sm

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

    if config["gumbo_env"] == "dev":
        gumbo_client = GumboClient(
            url="http://localhost:8080/v1/graphql",
            username="dogspa",
            headers={"X-Hasura-Admin-Secret": "secret"},
        )

    else:
        # get URL and password for Gumbo GraphQL API from secrets manager
        hasura_creds = get_hasura_creds(gumbo_env=config["gumbo_env"])

        gumbo_client = GumboClient(
            url=hasura_creds["url"],
            username="dogspa",
            headers={"X-Hasura-Admin-Secret": hasura_creds["password"]},
        )

    ctx.obj = {"gumbo_client": gumbo_client}


@app.command()
def update_workflow(
    ctx: typer.Context, workflow_name: Annotated[str, typer.Option()]
) -> None:
    # need a GitHub PAT for persisting WDL in gists
    github_pat = get_secret_from_sm(
        "projects/201811582504/secrets/github-pat-for-wdl-gists/versions/latest"
    )

    terra_workflow = TerraWorkflow(
        method_namespace=config["terra"][workflow_name]["method_namespace"],
        method_name=config["terra"][workflow_name]["method_name"],
        method_config_namespace=config["terra"][workflow_name][
            "method_config_namespace"
        ],
        method_config_name=config["terra"][workflow_name]["method_config_name"],
        method_synopsis=config["terra"][workflow_name]["method_synopsis"],
        workflow_wdl_path=Path(
            config["terra"][workflow_name]["workflow_wdl_path"]
        ).resolve(),
        method_config_json_path=Path(
            config["terra"][workflow_name]["method_config_json_path"]
        ).resolve(),
        github_pat=github_pat,
    )

    ctx.obj["terra_workspace"].update_workflow(terra_workflow=terra_workflow)


@app.command()
def upsert_delivery_bams(ctx: typer.Context) -> None:
    do_upsert_delivery_bams(
        gcs_source_bucket=config["alignment"]["gcs_source"]["bucket"],
        gcs_source_glob=config["alignment"]["gcs_source"]["glob"],
        uuid_namespace=config["uuid_namespace"],
        terra_workspace=TerraWorkspace(
            workspace_namespace=config["terra"]["delivery_workspace_namespace"],
            workspace_name=config["terra"]["delivery_workspace_name"],
            owners=json.loads(os.environ["FIRECLOUD_OWNERS"]),
        ),
        dry_run=config["dry_run"],
    )


@app.command()
def delta_align_delivery_bams(ctx: typer.Context) -> None:
    terra_workflow = TerraWorkflow(
        method_namespace=config["terra"]["align_long_reads"]["method_namespace"],
        method_name=config["terra"]["align_long_reads"]["method_name"],
        method_config_namespace=config["terra"]["align_long_reads"][
            "method_config_namespace"
        ],
        method_config_name=config["terra"]["align_long_reads"]["method_config_name"],
        method_synopsis=config["terra"]["align_long_reads"]["method_synopsis"],
        workflow_wdl_path=Path(
            config["terra"]["align_long_reads"]["workflow_wdl_path"]
        ).resolve(),
        method_config_json_path=Path(
            config["terra"]["align_long_reads"]["method_config_json_path"]
        ).resolve(),
    )

    do_delta_align_delivery_bams(
        terra_workspace=TerraWorkspace(
            workspace_namespace=config["terra"]["delivery_workspace_namespace"],
            workspace_name=config["terra"]["delivery_workspace_name"],
            owners=json.loads(os.environ["FIRECLOUD_OWNERS"]),
        ),
        terra_workflow=terra_workflow,
        dry_run=config["dry_run"],
    )


@app.command()
def refresh_terra_samples(ctx: typer.Context) -> None:
    terra_workspace = TerraWorkspace(
        workspace_namespace=config["terra"]["workspace_namespace"],
        workspace_name=config["terra"]["workspace_name"],
    )

    short_read_terra_workspace = TerraWorkspace(
        workspace_namespace=config["terra"]["short_read_workspace_namespace"],
        workspace_name=config["terra"]["short_read_workspace_name"],
    )

    do_refresh_terra_samples(
        terra_workspace=terra_workspace,
        short_read_terra_workspace=short_read_terra_workspace,
        gumbo_client=ctx.obj["gumbo_client"],
    )


if __name__ == "__main__":
    app()
