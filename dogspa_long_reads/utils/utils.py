from __future__ import annotations

import json
import logging
import string
import uuid
from functools import partial
from math import ceil, sqrt
from time import sleep
from typing import (
    Callable,
    Generator,
    Iterable,
    List,
    Optional,
    OrderedDict,
    ParamSpec,
    Type,
    TypeVar,
)

import baseconv
import pandas as pd
import requests
from pandera.typing import DataFrame as TypedDataFrame

from dogspa_long_reads.utils.validators import (
    DogspaConfig,
    PanderaBaseSchema,
    PydanticBaseModel,
)

T = TypeVar("T")


def list_comprehender(x: list[T], i1: int, i2: int) -> list[T]:
    """
    Perform list comprehension to get a batch of items from a list.

    :param x: the list to extract a batch from
    :param i1: the lower index
    :param i2: the higher index
    :return: a batch from the list
    """

    return x[i1:i2]


def df_comprehender(x: pd.DataFrame, i1: int, i2: int) -> pd.DataFrame:
    """
    Perform list comprehension to get a batch of rows from a data frame.

    :param x: the data frame to extract a batch from
    :param i1: the lower index
    :param i2: the higher index
    :return: a batch from the data frame
    """

    return x.iloc[i1:i2, :]


def batch_evenly(
    items: Iterable[T] | pd.DataFrame, max_batch_size: int
) -> Generator[list[T] | pd.DataFrame, None, None]:
    """
    Yields evenly sized batches from an iterable or data frame such that each batch has
    at most `max_batch_size` items.

    :param items: the iterable or DataFrame to be batched
    :param max_batch_size: the maximum size of each batch
    :return: a generator yielding batches from the input items
    """

    if isinstance(items, pd.DataFrame):
        batchable_items = items
        comprehender = df_comprehender
    else:
        try:
            batchable_items = list(items)
            comprehender = list_comprehender
        except TypeError as e:
            raise TypeError(f"Cannot batch items of type {type(items)}: {e}")

    n_items = len(batchable_items)
    n_batches = 1 + n_items // max_batch_size
    batch_size = n_items / n_batches

    for i in range(n_batches):
        i1 = ceil(i * batch_size)
        i2 = ceil((i + 1) * batch_size)
        yield comprehender(batchable_items, i1, i2)  # pyright:ignore


P = ParamSpec("P")
R = TypeVar("R")


def generalized_fibonacci(n: int, *, f0: float = 1.0, f1: float = 1.0) -> float:
    """
    Calculate the nth number in a generalized Fibonacci sequence given two starting
    nonnegative real numbers. This generates a gradually increasing sequence that
    provides a good balance between linear and exponential functions for use as a
    backoff.

    :param n: the nth Fibonacci number to compute
    :param f0: the first starting value for the sequence
    :param f1: the second starting value for the sequence
    :return: the nth Fibonacci number
    """

    assert f0 >= 0, "f0 must be at least 0.0"
    assert f1 >= 0, "f1 must be at least 0.0"

    # compute constants for closed-form of Fibonacci sequence recurrence relation
    sqrt5 = sqrt(5)
    phi = (1 + sqrt5) / 2
    psi = 1 - phi
    a = (f1 - f0 * psi) / sqrt5
    b = (f0 * phi - f1) / sqrt5

    return max([0, a * phi**n + b * psi**n])


def maybe_retry(
    func: Callable[P, R],
    retryable_exceptions: tuple[Type[Exception], ...] = tuple([Exception]),
    max_retries: int = 0,
    waiter: Callable[..., float] = partial(generalized_fibonacci, f0=1.0, f1=1.0),
    *args: P.args,
    **kwargs: P.kwargs,
) -> R:
    """
    Call a function and optionally retry (at most `max_retries` times) if it raises
    certain exceptions.

    :param func: a function
    :param retryable_exceptions: a tuple of retryable exceptions
    :param max_retries: the maximum number of times to retry
    :param waiter: a function that returns the number of seconds to wait given how many
    tries have already happened
    :param kwargs: keyword arguments to `func`
    :return: the return value from `func`
    """

    if max_retries == 0:
        return func(*args, **kwargs)

    n_retries = 0

    while True:
        try:
            return func(*args, **kwargs)

        except retryable_exceptions as e:
            if n_retries == max_retries:
                raise e

            wait_seconds = round(waiter(n_retries + 1), 1)
            print(f"{e} (retrying in {wait_seconds}s)")
            sleep(wait_seconds)
            n_retries += 1


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
    return TypedDataFrame[pandera_schema](records)


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


def expand_dict_columns(
    df: pd.DataFrame,
    sep: str = "__",
    name_columns_with_parent: bool = True,
    parent_key: str = "",
    col_name_formatter: Callable[[str], str] = lambda _: _,
) -> pd.DataFrame:
    """
    Recursively expand columns in a data frame containing dictionaries into separate
    columns.

    :param df: a data frame
    :param sep: a separator character to use between `parent_key` and its column names
    :param name_columns_with_parent: whether to "namespace" nested column names using
    their parents' column names
    :param parent_key: the name of the parent column, applicable only if
    `name_columns_with_parent` is `True` (for recursion)
    :param col_name_formatter: an optional function to format resulting column names
    :return: a widened data frame
    """

    flattened_dict = {}

    for c, s in df.items():
        fvi = s.first_valid_index()

        if fvi is not None and isinstance(s.loc[fvi], dict):
            # if the column contains dictionaries, recursively flatten them
            nested_df = pd.json_normalize(s.tolist())
            nested_df.index = df.index

            if name_columns_with_parent:
                # e.g. if current column `c` is "foo" and the nested data contains a
                # field "bar", the resulting column name is "foo__bar"

                nested_df.columns = [
                    sep.join(
                        [
                            parent_key,
                            col_name_formatter(str(c)),
                            col_name_formatter(str(col)),
                        ]
                    )
                    if parent_key != ""
                    else sep.join(
                        [col_name_formatter(str(c)), col_name_formatter(str(col))]
                    )
                    for col in nested_df.columns
                ]

            # recurse on the nested data
            flattened_dict.update(
                expand_dict_columns(
                    nested_df,
                    sep=sep,
                    name_columns_with_parent=name_columns_with_parent,
                    parent_key=col_name_formatter(str(c)),
                )
            )

        else:
            # if not a dictionary, add the column as is
            flattened_dict[c] = s

    df = pd.DataFrame(flattened_dict)

    if parent_key == "":
        # make sure there are no duplicate column names after all expansion is done
        col_name_counts = df.columns.value_counts()

        if col_name_counts.gt(1).any():
            dup_names = set(
                col_name_counts[col_name_counts.gt(1)].index,  # pyright: ignore
            )
            raise NameError(
                f"Column names {dup_names} are duplicated. Try calling "
                "`expand_dict_columns` with `name_columns_with_parent=True`."
            )

    return df


def send_slack_message(
    slack_webhook_url_errors: Optional[str],
    slack_webhook_url_stats: Optional[str],
    stats: OrderedDict[str, int],
    report: OrderedDict[str, pd.DataFrame],
    config: DogspaConfig,
) -> None:
    """
    Send a message to a Slack channel with content of `stats` and `report`.

    :param slack_webhook_url_errors: URL for a Slack Webhook
    :param slack_webhook_url_stats: URL for a Slack Webhook
    :param stats: ordered dictionary of metadata statistics
    :param report: ordered dictionary of relevant sample metadata
    :param config: the dogspa configuration
    """

    if config.dry_run:
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

    terra_ws_name = "/".join([config.workspace.namespace, config.workspace.name])

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
