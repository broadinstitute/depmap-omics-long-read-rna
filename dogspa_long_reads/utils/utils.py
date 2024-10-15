import json
import logging
import string
import uuid
from typing import List, Optional, OrderedDict, Type

import baseconv
import pandas as pd
import requests
from nebelung.terra_workspace import TerraWorkspace
from nebelung.types import PanderaBaseSchema
from nebelung.utils import type_data_frame
from pandera.typing import DataFrame as TypedDataFrame

from dogspa_long_reads.types import PydanticBaseModel


def model_to_df(
    model: PydanticBaseModel,
    pandera_schema: Type[PanderaBaseSchema],
    records_key: str = "records",
) -> TypedDataFrame[PanderaBaseSchema]:
    """
    Dump a Pydantic model and convert it to a data frame typed by a Pandera schema.

    :param model: a Pydandict model containing a list of objects keyed by `records_key`
    :param pandera_schema: the Pandera schema to cast the model to
    :param records_key: the key/method name in `model` containing the records
    """

    records = model.model_dump()[records_key]
    return type_data_frame(records, pandera_schema)


def df_to_model(
    df: pd.DataFrame, pydantic_schema: Type[PydanticBaseModel]
) -> List[PydanticBaseModel]:
    """
    Convert a Pandas data frame to a Pydantic model.

    :param df: a data frame
    :param pydantic_schema: the Pydantic schema to cast the data frame to
    :return: a Pydantic model
    """

    return [pydantic_schema(**x) for x in df.to_dict(orient="records")]


def uuid_to_base62(x: str) -> str:
    """
    Treat a UUID as a hexadecimal number and convert it to base62 (digits+uppercase and
    lowercase letters).

    :param x: a string representation of a UUID
    :return: the base62 representation of `x` as a string
    """

    converter = baseconv.BaseConverter("".join([string.digits, string.ascii_letters]))
    uuid_hex = x.replace("-", "")
    uuid4_as_int = int(uuid_hex, 16)
    return converter.encode(uuid4_as_int)


def send_slack_message(
    slack_webhook_url_errors: Optional[str],
    slack_webhook_url_stats: Optional[str],
    stats: OrderedDict[str, int],
    report: OrderedDict[str, pd.DataFrame],
    terra_workspace: TerraWorkspace,
    dry_run: bool,
) -> None:
    """
    Send a message to a Slack channel with content of `stats` and `report`.

    :param slack_webhook_url_errors: URL for a Slack Webhook
    :param slack_webhook_url_stats: URL for a Slack Webhook
    :param stats: ordered dictionary of metadata statistics
    :param report: ordered dictionary of relevant sample metadata
    :param terra_workspace: a TerraWorkspace instance
    :param dry_run: whether to skip updates to external data stores
    """

    if dry_run:
        logging.info("(skipping) Sending summary to Slack channel...")

        logging.info("Stats:")
        logging.info(json.dumps(dict(stats), default=int, indent=4))

        logging.info("Report:")
        logging.info(
            json.dumps(
                [
                    {
                        k: df.loc[:, df.columns.str.endswith("_id")].apply(
                            lambda x: " / ".join([s for s in x if s is not pd.NA]),
                            axis=1,
                        )
                    }
                    for k, df in report.items()
                    if len(df) > 0
                ],
                default=list,  # type: ignore
                indent=4,
            )
        )

        return

    if slack_webhook_url_errors is None or slack_webhook_url_stats is None:
        raise ValueError("Slack webhook URLs are required")

    logging.info("Sending stats to Slack channel...")

    terra_ws_name = "/".join(
        [terra_workspace.workspace_namespace, terra_workspace.workspace_name]
    )

    stats_blocks = [
        {
            "type": "header",
            "text": {
                "type": "plain_text",
                "text": f"Stats (Workspace: {terra_ws_name})",
                "emoji": False,
            },
        },
        {
            "type": "rich_text",
            "elements": [
                {
                    "type": "rich_text_list",
                    "style": "bullet",
                    "elements": [
                        {
                            "type": "rich_text_section",
                            "elements": [{"type": "text", "text": f"{k}: {v}"}],
                        }
                        for k, v in stats.items()
                    ],
                },
            ],
        },
    ]

    r = requests.post(
        slack_webhook_url_stats, data=json.dumps({"blocks": stats_blocks})
    )
    r.raise_for_status()

    if sum(len(df) for df in report.values()) == 0:
        return

    # there's at least one non-empty data frame in the report
    logging.info("Sending results to Slack channel...")

    results_blocks = [
        {
            "type": "header",
            "text": {
                "type": "plain_text",
                "text": f"Results (Workspace: {terra_ws_name})",
                "emoji": False,
            },
        }
    ]

    for k, v in report.items():
        if len(v) == 0:
            continue

        # print all available ID columns for each row in this report section
        texts = v.loc[:, v.columns.str.endswith("_id")].apply(
            lambda x: " / ".join([s for s in x if s is not pd.NA]),
            axis=1,
        )

        results_blocks.append(
            {
                "type": "rich_text",
                "elements": [
                    {
                        "type": "rich_text_section",
                        "elements": [{"type": "text", "text": f"{k} ({len(v)})"}],
                    },
                    {
                        "type": "rich_text_list",
                        "style": "bullet",
                        "indent": 1,
                        "elements": [
                            {
                                "type": "rich_text_section",
                                "elements": [{"type": "text", "text": x}],
                            }
                            for x in sorted(texts)
                        ],
                    },
                ],
            }
        )

    r = requests.post(
        slack_webhook_url_errors, data=json.dumps({"blocks": results_blocks})
    )
    r.raise_for_status()


def assign_hashed_uuids(
    df: pd.DataFrame, uuid_namespace: str, uuid_col_name: str, subset: list[str]
) -> pd.DataFrame:
    """
    Compute and add a consistent UUID-formatted ID column to a data frame.

    :param df: a data frame
    :param uuid_namespace: a namespace for generated UUIDv3s
    :param uuid_col_name: the name for the UUID column
    :param subset: the subset of columns to use for hashing
    :return: `df` with the new UUID column
    """

    df[uuid_col_name] = (
        df[sorted(subset)]
        .apply(lambda x: uuid.uuid3(uuid.UUID(uuid_namespace), x.to_json()), axis=1)
        .astype("string")
    )

    return df
