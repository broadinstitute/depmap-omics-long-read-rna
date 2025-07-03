import string
import uuid
from pathlib import Path
from typing import Any, Callable, List, Type

import baseconv
import pandas as pd
from google.cloud import secretmanager_v1
from nebelung.terra_workflow import TerraWorkflow
from nebelung.types import PanderaBaseSchema
from nebelung.utils import type_data_frame
from pandera.typing import DataFrame as TypedDataFrame
from pydantic import BaseModel


def get_hasura_creds(gumbo_env: str) -> dict[str, str]:
    """
    Get URL and password for Hasura GraphQL API.

    :param gumbo_env: the Gumbo env to get credentials for ('staging' or 'prod')
    :return: a dictionary with the GraphQL API URL and password
    """

    return {
        "url": get_secret_from_sm(
            f"projects/814840278102/secrets/hasura-{gumbo_env}-api-url/versions/latest"
        ),
        "password": get_secret_from_sm(
            f"projects/814840278102/secrets/hasura-admin-secret-{gumbo_env}/versions/latest"
        ),
    }


def get_secret_from_sm(name: str) -> str:
    """
    Get the value of a secret from GCP Secret Manager.

    :param name: the fully-qualified name of a secret
    :return: the secret's decoded value
    """

    client = secretmanager_v1.SecretManagerServiceClient()
    request = secretmanager_v1.AccessSecretVersionRequest(mapping={"name": name})
    response = client.access_secret_version(request=request)
    return response.payload.data.decode()


def model_to_df(
    model: BaseModel,
    pandera_schema: Type[PanderaBaseSchema],
    records_key: str = "records",
    remove_unknown_cols: bool = False,
    mutator: Callable[[pd.DataFrame], pd.DataFrame] = lambda _: _,
) -> TypedDataFrame[PanderaBaseSchema]:
    """
    Dump a Pydantic model and convert it to a data frame typed by a Pandera schema.

    :param model: a Pydandict model containing a list of objects keyed by `records_key`
    :param pandera_schema: the Pandera schema to cast the model to
    :param records_key: the key/method name in `model` containing the records
    :param remove_unknown_cols: remove columns not specified in the schema
    :param mutator: an optional function to call on the data frame before typing (e.g.
    to rename columns to expected Pydantic field names)
    """

    records = model.model_dump()[records_key]

    df = pd.DataFrame(records)
    df = mutator(df)
    return type_data_frame(df, pandera_schema, remove_unknown_cols)


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
        .apply(
            lambda x: str(uuid.uuid3(uuid.UUID(uuid_namespace), x.to_json())), axis=1
        )
        .astype("string")
    )

    return df


def make_workflow_from_config(
    config: dict[str, Any], workflow_name: str, **kwargs: Any
) -> TerraWorkflow:
    """
    Make a TerraWorkflow object from a config entry.

    :param config: a config dictionary
    :param workflow_name: the name of the workflow referenced in the config
    :return: a TerraWorkflow instance
    """

    return TerraWorkflow(
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
        **kwargs,
    )
