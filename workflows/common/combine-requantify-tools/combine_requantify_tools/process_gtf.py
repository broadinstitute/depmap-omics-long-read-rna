import logging
import re
from pathlib import Path

import pandas as pd

from combine_requantify_tools.types import Gtf, TypedDataFrame
from combine_requantify_tools.utils import read_gtf_from_file, type_data_frame

GENE_NAME_RE = re.compile(r'gene_name "([^"]+)"')
GENE_ID_RE = re.compile(r'gene_id "([^"]+)"')


def do_process_gtf(gtf_in: Path) -> TypedDataFrame[Gtf]:
    """
    Update GTF attributes by populating gene_id values.

    :param gtf_in: path to input GTF
    :return: processed GTF data frame
    """

    # load the GTF
    gtf = read_gtf_from_file(gtf_in)

    # parse attributes column to amke mapping of gene names to IDs
    logging.info(f"Mapping annotated gene IDs")
    annot_attributes = (
        gtf.loc[
            gtf["source"].str.contains("HAVANA|ENSEMBL", regex=True, na=False),
            "attribute",
        ]
        .dropna()
        .drop_duplicates()
    )

    gene_name_id_map = (
        pd.DataFrame(
            {
                "gene_name": annot_attributes.map(extract_gene_name),
                "gene_id": annot_attributes.map(extract_gene_id),
            }
        )
        .dropna()
        .drop_duplicates()
        .set_index("gene_name")["gene_id"]
        .to_dict()
    )

    # replace attributes using the mapping
    current_gname = None
    new_attrs = []

    for source, attr in zip(gtf["source"], gtf["attribute"]):
        gid = extract_gene_id(attr)
        gname = extract_gene_name(attr)

        # track last seen gene_name
        if gname is not None:
            current_gname = gname

        if source == "IsoQuant" and gid is not None:
            # we can potentially updated the gene ID using the mapping
            if gname is None and current_gname in gene_name_id_map:
                # no gene_name on this line; use most recent one
                attr = replace_gene_id(
                    attr, gene_id_old=gid, gene_id_new=gene_name_id_map[current_gname]
                )
            elif gname is not None and gname in gene_name_id_map:
                # gene_name present on this line
                attr = replace_gene_id(
                    attr, gene_id_old=gid, gene_id_new=gene_name_id_map[gname]
                )

        new_attrs.append(attr)

    gtf["attribute"] = new_attrs

    return type_data_frame(gtf, Gtf)


def extract_gene_name(attr: str) -> str | None:
    """
    Extract gene_name value from a GTF attribute string.

    :param attr: GTF attribute string containing key-value pairs
    :return: gene_name value if found, None otherwise
    """

    m = GENE_NAME_RE.search(attr)
    return m.group(1) if m else None


def extract_gene_id(attr: str) -> str | None:
    """
    Extract gene_id value from a GTF attribute string.

    :param attr: GTF attribute string containing key-value pairs
    :return: gene_id value if found, None otherwise
    """

    m = GENE_ID_RE.search(attr)
    return m.group(1) if m else None


def replace_gene_id(attr: str, gene_id_old: str, gene_id_new: str) -> str:
    """
    Replace gene_id value in a GTF attribute string.

    :param attr: GTF attribute string containing key-value pairs
    :param gene_id_old: existing gene_id value to replace
    :param gene_id_new: new gene_id value to use as replacement
    :return: updated attribute string with replaced gene_id
    """

    return re.sub(
        rf'gene_id "{re.escape(gene_id_old)}"',
        f'gene_id "{gene_id_new}"',
        attr,
        count=1,
    )
