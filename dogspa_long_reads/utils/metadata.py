from __future__ import annotations

import logging
import pathlib
import uuid

import pandas as pd
from pandera.typing import DataFrame as TypedDataFrame

from dogspa_long_reads.utils.terra import TerraWorkspace
from dogspa_long_reads.utils.utils import (
    df_to_model,
    expand_dict_columns,
    uuid_to_base62,
)
from dogspa_long_reads.utils.validators import (
    DogspaConfig,
    IdentifiedSrcBam,
    ModelsAndChildren,
    ObjectMetadata,
    SamplesForGumbo,
    SamplesMaybeInGumbo,
    SamplesWithCDSIDs,
    SamplesWithMetadata,
    SamplesWithRgUpdatedAt,
    SamplesWithShortReadMetadata,
    SeqTable,
    VersionedSamples,
)
from gumbo_gql_client import GumboClient, omics_sequencing_insert_input


def id_bams(
    bams: TypedDataFrame[ObjectMetadata],
) -> TypedDataFrame[IdentifiedSrcBam]:
    bams["model_id"] = (
        bams["url"]
        .apply(lambda x: pathlib.Path(x).name)
        .str.extract(r"^(ACH-[A-Z 0-9]+)")
    )

    if bams["model_id"].isna().any():
        raise ValueError("There are BAM files not named with model IDs (ACH-*)")

    bams = bams.rename(columns={"url": "bam_url"})

    # start tracking issues to store in Gumbo seq table
    bams["issue"] = pd.Series([set()] * len(bams))
    bams["blacklist"] = False

    return TypedDataFrame[IdentifiedSrcBam](bams)


def explode_and_expand_models(
    models: TypedDataFrame[ModelsAndChildren],
) -> TypedDataFrame[SeqTable]:
    for c in ["model_conditions", "omics_profiles", "omics_sequencings"]:
        models = models.explode(c).reset_index(drop=True)
        models = expand_dict_columns(models, name_columns_with_parent=False)

    models = models.dropna(subset=["model_condition_id", "profile_id"])

    models["is_main_sequencing_id"] = models["main_sequencing_id"].eq(
        models["sequencing_id"]
    )

    models[["blacklist_omics", "blacklist"]] = models[
        ["blacklist_omics", "blacklist"]
    ].fillna(False)

    return TypedDataFrame[SeqTable](models.drop(columns="main_sequencing_id"))


def join_metadata(
    samples: TypedDataFrame[SamplesWithShortReadMetadata],
    seq_table: TypedDataFrame[SeqTable],
) -> TypedDataFrame[SamplesWithMetadata]:
    seq_table = seq_table.loc[
        seq_table["datatype"].eq("long_read_rna")
        & ~seq_table["blacklist"]
        & ~seq_table["blacklist_omics"]
    ]

    samples = samples.merge(
        seq_table[["model_id", "model_condition_id", "profile_id"]],
        how="left",
        on="model_id",
    )

    return TypedDataFrame[SamplesWithMetadata](samples)


def join_short_read_metadata(
    samples: TypedDataFrame[SamplesWithMetadata], seq_table: TypedDataFrame[SeqTable]
) -> TypedDataFrame[SamplesWithShortReadMetadata]:
    # reproduce logic of `makeDefaultModelTable` in depmap_omics_upload
    source_priority = [
        "BROAD",
        "DEPMAP",
        "IBM",
        "CCLE2",
        "SANGER",
        "PRISM",
        "CCLF",
        "CHORDOMA",
    ]

    sp_df = pd.DataFrame(
        {
            "source": source_priority,
            "source_priority": list(range(len(source_priority))),
        }
    )

    seq_table = seq_table.loc[
        seq_table["model_id"].isin(samples["model_id"])
        & ~seq_table["blacklist"]
        & ~seq_table["blacklist_omics"]
    ].merge(sp_df, how="left", on="source")

    sr = seq_table.loc[seq_table["datatype"].eq("rna")]

    assert sr["source_priority"].notna().all()

    main_seq_ids = (
        sr.loc[sr["is_main_sequencing_id"], ["model_id", "sequencing_id"]]
        .drop_duplicates()
        .rename(columns={"sequencing_id": "main_sequencing_id"})
    )

    assert ~main_seq_ids["model_id"].duplicated().any()

    samples = samples.merge(main_seq_ids, how="left", on="model_id")

    sr_rna = sr.sort_values("source_priority").groupby("model_id").nth(0)

    sr_rna = sr_rna[["model_id", "profile_id", "bai_filepath", "bam_filepath"]].rename(
        columns={
            "profile_id": "sr_profile_id",
            "bai_filepath": "sr_bai_filepath",
            "bam_filepath": "sr_bam_filepath",
        }
    )

    samples = samples.merge(sr_rna, how="left", on="model_id")

    assert (
        samples[["sr_profile_id", "sr_bai_filepath", "sr_bam_filepath"]]
        .notna()
        .all(axis=None)
    )

    assert ~samples["sr_profile_id"].duplicated(keep=False).any()

    return TypedDataFrame[SamplesWithShortReadMetadata](samples)


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


def assign_cds_ids(
    samples: TypedDataFrame[IdentifiedSrcBam], config: DogspaConfig
) -> TypedDataFrame[SamplesWithCDSIDs]:
    """
    Assign a "CDS-" ID to each sample by hashing relevant columns.

    :param samples: the data frame of samples
    :param config: the dogspa configuration
    :return: the data frame with CDS IDs
    """

    logging.info("Assigning CDS IDs...")

    samples_w_ids = assign_hashed_uuids(
        samples,
        uuid_namespace=config.uuid_namespace,
        uuid_col_name="sequencing_id",
        subset=[
            "model_id",
            "crc32c",
            "size",
            "gcs_obj_updated_at",
        ],
    )

    # convert UUIDs to base62 and truncate to match existing ID format
    samples_w_ids["sequencing_id"] = "CDS-" + samples_w_ids["sequencing_id"].apply(
        uuid_to_base62
    ).str.slice(-6)

    return TypedDataFrame[SamplesWithCDSIDs](samples_w_ids)


def check_already_in_gumbo(
    samples: TypedDataFrame[IdentifiedSrcBam],
    seq_table: TypedDataFrame[SeqTable],
    size_col_name: str,
) -> TypedDataFrame[SamplesMaybeInGumbo]:
    """
    Mark samples that are already in Gumbo by comparing stored file sizes.

    :param samples: the data frame of samples
    :param seq_table: the data frame of Gumbo Omics sequencing data
    :param size_col_name: the Gumbo column containing file sizes
    :return: the samples data frame with `already_in_gumbo` column
    """

    logging.info("Checking Gumbo for existing samples...")

    samples["already_in_gumbo"] = (
        samples["size"]
        .isin(seq_table[size_col_name].dropna().astype("int64"))
        .astype("bool")
    )

    return TypedDataFrame[SamplesMaybeInGumbo](samples)


def apply_col_map(
    samples: TypedDataFrame[SamplesWithCDSIDs], config: DogspaConfig
) -> TypedDataFrame[SamplesForGumbo]:
    """
    Rename the columns in the samples data frame to their corresponding names in Gumbo.

    :param samples: the data frame of samples
    :param config: the dogspa configuration
    :return: the data frame with just the relevant renamed columns
    """

    logging.info("Renaming columns for Gumbo...")

    samples["source"] = "DEPMAP"
    samples["expected_type"] = "rna_long_read"
    samples["sequencing_date"] = samples["gcs_obj_updated_at"]
    samples["stranded"] = True

    # rename columns for Gumbo
    gumbo_samples = samples.rename(
        columns={
            "bam_url": "bam_filepath",
            "size": "legacy_size",
            "gcs_obj_updated_at": "update_time",
            "crc32c": "crc32c_hash",
        }
    ).drop(columns="sample_id")

    # concat issues for each sample into single string
    gumbo_samples["issue"] = gumbo_samples["issue"].apply(
        lambda x: "; ".join(x) if len(x) > 0 else pd.NA
    )

    return TypedDataFrame[SamplesForGumbo](gumbo_samples)


def upsert_terra_samples(tw: TerraWorkspace, samples: pd.DataFrame) -> None:
    """
    CellLineName
    LR_bam_filepath
    LongReadProfileID
    MainSequencingID
    ModelCondition
    ShortReadProfileID
    StrippedCellLineName
    participant
    """

    terra_samples = samples.rename(
        columns={
            "model_id": "entity:sample_id",
            "bam_url": "bam_filepath",
        }
    )


def increment_sample_versions(
    samples: TypedDataFrame[SamplesForGumbo],
    seq_table: TypedDataFrame[SeqTable],
    config: DogspaConfig,
) -> TypedDataFrame[VersionedSamples]:
    """
    For sample profile IDs that already exist in the Gumbo profiles table, increment
    the version numbers beyond the default value of 1.

    :param samples: the data frame of samples
    :param seq_table: the data frame of Gumbo sequence data
    :param config: the dogspa configuration
    :return: the samples data frame with optionally incremented version numbers
    """

    logging.info("Incrementing sample version numbers...")

    # count distinct profile IDs in existing seq table from among those observed in the
    # current batch of samples
    existing_profile_id_counts = (
        seq_table.loc[
            seq_table["expected_type"].eq(config.data_type)
            & seq_table["profile_id"].isin(samples["profile_id"]),
            "profile_id",
        ]
        .value_counts()
        .to_frame(name="version_n")
        .reset_index(names=["profile_id"])
    )

    # increment version numbers in samples if applicable
    versioned_samples = samples.merge(
        existing_profile_id_counts, how="left", on="profile_id"
    )
    versioned_samples["version_n"] = versioned_samples["version_n"].fillna(0)
    versioned_samples["version"] = 1 + versioned_samples["version_n"]
    versioned_samples = versioned_samples.drop(columns="version_n")

    return TypedDataFrame[VersionedSamples](versioned_samples)


def upload_to_gumbo(
    gumbo_client: GumboClient,
    samples: TypedDataFrame[VersionedSamples],
    config: DogspaConfig,
) -> TypedDataFrame[VersionedSamples]:
    """
    Upload the samples data frame to the Gumbo sequencing table.

    :param gumbo_client: an instance of the Gumbo GraphQL client
    :param samples: the data frame of samples
    :param config: the dogspa configuration
    :return: the samples data frame that was uploaded to Gumbo
    """

    # set custom hard-coded values if applicable
    for k, v in config.gumbo_custom_values.items():
        samples[k] = v

    # match historical behavior to exclude bad records except for too-small
    # BAM files (so that other types of bad records like missing GCS
    # objects that might get resolved later still get processed the next time
    # this automation is run)
    samples_to_upload = samples.loc[
        ~samples["blacklist"] | samples["issue"].eq("BAM file too small")
    ]

    if config.dry_run:
        logging.info(f"(skipping) Inserting {len(samples_to_upload)} samples")
        return samples_to_upload

    logging.info(f"Inserting {len(samples_to_upload)} samples to")

    if len(samples) > 0:
        objects = df_to_model(samples_to_upload, omics_sequencing_insert_input)
        res = gumbo_client.insert_omics_sequencings(username="dogspa", objects=objects)
        affected_rows = res.insert_omics_sequencing.affected_rows  # type: ignore
        logging.info(f"Inserted {affected_rows} samples to Gumbo")

    return samples_to_upload
