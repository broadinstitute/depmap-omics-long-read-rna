import json
import logging
import os
import tomllib
from pathlib import Path
from typing import Annotated, Any

import pandas as pd
import typer
from nebelung.terra_workspace import TerraWorkspace

import depmap_omics_long_read_rna.utils.aligned_bams as aligned_bams
import depmap_omics_long_read_rna.utils.delivery_bams as delivery_bams
import depmap_omics_long_read_rna.utils.metadata as metadata
from depmap_omics_long_read_rna.types import GumboClient
from depmap_omics_long_read_rna.utils.utils import (
    get_hasura_creds,
    get_secret_from_sm,
    make_workflow_from_config,
)

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

    def get_gumbo_client() -> GumboClient:
        if config["gumbo_env"] == "dev":
            return GumboClient(
                url="http://localhost:8080/v1/graphql",
                username="dogspa",
                headers={"X-Hasura-Admin-Secret": "secret"},
            )

        # get URL and password for Gumbo GraphQL API from secrets manager
        hasura_creds = get_hasura_creds(gumbo_env=config["gumbo_env"])

        return GumboClient(
            url=hasura_creds["url"],
            username="dogspa",
            headers={"X-Hasura-Admin-Secret": hasura_creds["password"]},
        )

    ctx.obj = {"get_gumbo_client": get_gumbo_client}


@app.command()
def update_workflow(workflow_name: Annotated[str, typer.Option()]) -> None:
    terra_workspace = TerraWorkspace(
        workspace_namespace=config["terra"]["workspace_namespace"],
        workspace_name=config["terra"]["workspace_name"],
        owners=json.loads(os.environ["FIRECLOUD_OWNERS"]),
    )

    # need a GitHub PAT for persisting WDL in gists
    github_pat = get_secret_from_sm(
        "projects/201811582504/secrets/github-pat-for-wdl-gists/versions/latest"
    )

    terra_workspace.update_workflow(
        terra_workflow=make_workflow_from_config(
            config, workflow_name, github_pat=github_pat
        )
    )


@app.command()
def upsert_delivery_bams() -> None:
    delivery_bams.upsert_delivery_bams(
        gcs_source_bucket=config["alignment"]["gcs_source"]["bucket"],
        gcs_source_glob=config["alignment"]["gcs_source"]["glob"],
        uuid_namespace=config["uuid_namespace"],
        terra_workspace=TerraWorkspace(
            workspace_namespace=config["terra"]["delivery_workspace_namespace"],
            workspace_name=config["terra"]["delivery_workspace_name"],
        ),
        dry_run=config["dry_run"],
    )


@app.command()
def onboard_aligned_bams(ctx: typer.Context) -> None:
    aligned_bams.onboard_aligned_bams(
        terra_workspace=TerraWorkspace(
            workspace_namespace=config["terra"]["workspace_namespace"],
            workspace_name=config["terra"]["workspace_name"],
        ),
        gumbo_client=ctx.obj["get_gumbo_client"](),
        dry_run=config["dry_run"],
    )


@app.command()
def refresh_terra_samples(ctx: typer.Context) -> None:
    metadata.refresh_terra_samples(
        terra_workspace=TerraWorkspace(
            workspace_namespace=config["terra"]["workspace_namespace"],
            workspace_name=config["terra"]["workspace_name"],
        ),
        short_read_terra_workspace=TerraWorkspace(
            workspace_namespace=config["terra"]["short_read_workspace_namespace"],
            workspace_name=config["terra"]["short_read_workspace_name"],
        ),
        gumbo_client=ctx.obj["get_gumbo_client"](),
    )


@app.command()
def delta_job(
    workflow_name: Annotated[str, typer.Option()],
    entity_type: Annotated[str, typer.Option()],
    entity_set_type: Annotated[str, typer.Option()],
    entity_id_col: Annotated[str, typer.Option()],
    expression: Annotated[str, typer.Option()],
    input_col: Annotated[list[str] | None, typer.Option()] = None,
    output_col: Annotated[list[str] | None, typer.Option()] = None,
) -> None:
    terra_workspace = TerraWorkspace(
        workspace_namespace=config["terra"]["workspace_namespace"],
        workspace_name=config["terra"]["workspace_name"],
    )

    input_cols_set = None
    output_cols_set = None

    if input_col is not None:
        input_cols_set = set(input_col)

    if output_col is not None:
        output_cols_set = set(output_col)

    terra_workspace.submit_delta_job(
        terra_workflow=make_workflow_from_config(config, workflow_name),
        entity_type=entity_type,
        entity_set_type=entity_set_type,
        entity_id_col=entity_id_col,
        expression=expression,
        dry_run=config["dry_run"],
        resubmit_n_times=10,
        input_cols=input_cols_set,
        output_cols=output_cols_set,
        max_n_entities=1,
    )


if __name__ == "__main__":
    app()
