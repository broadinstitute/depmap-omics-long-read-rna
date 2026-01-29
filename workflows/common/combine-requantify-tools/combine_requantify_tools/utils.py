from typing import Type

import pandas as pd

from combine_requantify_tools.types import PanderaBaseSchema, TypedDataFrame


def type_data_frame(
    df: pd.DataFrame,
    pandera_schema: Type[PanderaBaseSchema],
    reorder_cols: bool = True,
    remove_unknown_cols: bool = False,
) -> TypedDataFrame[PanderaBaseSchema]:
    """
    Coerce a data frame into one specified by a Pandera schema and optionally reorder
    columns and remove unknown columns.

    :param df: a data frame
    :param pandera_schema: a Pandera schema
    :param reorder_cols: reorder columns as specified in the schema
    :param remove_unknown_cols: remove columns not specified in the schema
    :return: a data frame validated with the provided Pandera schema
    """

    if len(df) == 0:
        # make an empty data frame that conforms to the Pandera schema
        s = pandera_schema.to_schema()

        # `example` doesn't know how to instantiate columns with structured data
        dict_cols = []
        list_cols = []

        for c in s.columns:
            if s.columns[c].dtype.type is dict:
                dict_cols.append(c)
                s = s.remove_columns([c])
            elif s.columns[c].dtype.type is list:
                list_cols.append(c)
                s = s.remove_columns([c])

        df = pd.DataFrame(s.example(size=1))

        if len(dict_cols) > 0:
            for c in dict_cols:
                df[c] = [{}] * len(df)

        if len(list_cols) > 0:
            for c in list_cols:
                df[c] = [[]] * len(df)

        df = df.iloc[:0]
        return TypedDataFrame[pandera_schema](df)

    if not remove_unknown_cols and not reorder_cols:
        # can type and return
        return TypedDataFrame[pandera_schema](df)

    # we need to collect the current columns and schema columns (in original orders)
    current_cols = list(df.columns)
    schema_cols = list(pandera_schema.to_schema().columns.keys())

    if remove_unknown_cols:
        # drop excess columns (if any)
        excess_cols = list(set(current_cols) - set(schema_cols))

        if len(excess_cols) > 0:
            df = df.drop(columns=excess_cols)
            current_cols = list(df.columns)

    # `df` might contain extra columns, but we can still type it now
    df = TypedDataFrame[pandera_schema](df)

    if reorder_cols:
        # put columns in schema order, with extra columns in original order at the end
        all_cols = schema_cols.copy()
        all_cols.extend(current_cols)
        all_cols = list(dict.fromkeys(all_cols))
        df = TypedDataFrame[pandera_schema](df.loc[:, all_cols])

    return df
