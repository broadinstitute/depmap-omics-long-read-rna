import string
import uuid
from typing import List, Type

import baseconv
import pandas as pd
from google.cloud import secretmanager_v1
from nebelung.types import PanderaBaseSchema
from nebelung.utils import type_data_frame
from pandera.typing import DataFrame as TypedDataFrame

from depmap_omics_long_read_rna.types import PydanticBaseModel
from gumbo_gql_client import BaseModel


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
