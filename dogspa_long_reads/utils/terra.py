import tempfile
from typing import Any, Callable, Type

import pandas as pd
import requests
from firecloud import api as firecloud_api
from pandera.typing import DataFrame as TypedDataFrame

from dogspa_long_reads.utils.utils import batch_evenly, maybe_retry
from dogspa_long_reads.utils.validators import PanderaBaseSchema


class TerraWorkspace:
    def __init__(self, workspace_namespace: str, workspace_name: str) -> None:
        self.workspace_namespace = workspace_namespace
        self.workspace_name = workspace_name

    def get_entities(
        self, entity_type: str, pandera_schema: Type[PanderaBaseSchema]
    ) -> TypedDataFrame[PanderaBaseSchema]:
        """
        Get a data frame of entities from a Terra data table.

        :param entity_type: the kind of entity (e.g. "sample")
        :param pandera_schema: a Pandera schema for the output data frame
        :return: a data frame of entities
        """

        print(f"Getting {entity_type} entities")
        j = call_firecloud_api(
            firecloud_api.get_entities,
            namespace=self.workspace_namespace,
            workspace=self.workspace_name,
            etype=entity_type,
        )

        records = [{f"{entity_type}_id": x["name"], **x["attributes"]} for x in j]

        return TypedDataFrame[pandera_schema](pd.DataFrame(records))

    def upload_entities(self, df: pd.DataFrame) -> None:
        """
        Upload a data frame of entities to a Terra data table.

        :param df: a data frame of entities
        """

        print(f"{len(df)} entities to upload to Terra")

        for batch in batch_evenly(df, max_batch_size=500):
            with tempfile.NamedTemporaryFile(suffix="tsv") as f:
                batch.to_csv(f, sep="\t", index=False)  # pyright: ignore
                f.flush()

                print(f"Upserting {len(batch)} entities to Terra")
                call_firecloud_api(
                    firecloud_api.upload_entities_tsv,
                    namespace=self.workspace_namespace,
                    workspace=self.workspace_name,
                    entities_tsv=f.name,
                    model="flexible",
                )


def call_firecloud_api(func: Callable, *args: Any, **kwargs: Any) -> Any:
    """
    Call a Firecloud API endpoint and check the response for a valid HTTP status code.

    :param func: a `firecloud.api` method
    :param args: arguments to `func`
    :param kwargs: keyword arguments to `func`
    :return: the API response, if any
    """

    res = maybe_retry(
        func,
        retryable_exceptions=(requests.ConnectionError, requests.ConnectTimeout),
        max_retries=4,
        *args,
        **kwargs,
    )

    if 200 <= res.status_code <= 299:
        try:
            return res.json()
        except requests.JSONDecodeError:
            return res.text

    try:
        raise requests.RequestException(f"HTTP {res.status_code} error: {res.json()}")
    except Exception as e:
        # it's returning HTML or we can't parse the JSON
        print(f"Error getting response as JSON: {e}")
        print(f"Response text: {res.text}")
        raise requests.RequestException(f"HTTP {res.status_code} error")
