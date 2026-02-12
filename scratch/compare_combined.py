import pandas as pd

pd.set_option("display.max_columns", 30)
pd.set_option("display.max_colwidth", 50)
pd.set_option("display.max_info_columns", 30)
pd.set_option("display.max_info_rows", 20)
pd.set_option("display.max_rows", 20)
pd.set_option("display.max_seq_items", None)
pd.set_option("display.width", 100)
pd.set_option("expand_frame_repr", True)
pd.set_option("mode.chained_assignment", "warn")


sq_old = pd.read_table("~/Desktop/combine/old/all-26q1_classification.tsv", sep="\t")
sq_new = pd.read_table("~/Desktop/combine/new/all-26q1_classification.tsv", sep="\t")

set(sq_old["isoform"]).symmetric_difference(set(sq_new["isoform"]))

sq_old["all_canonical"].value_counts()
sq_new["all_canonical"].value_counts()

sq_old["structural_category"].value_counts()
sq_new["structural_category"].value_counts()

sq_old["associated_gene"].value_counts()
sq_new["associated_gene"].value_counts()

tr_old = pd.read_table(
    "~/Desktop/combine/old/all-26q1_updated_tracking_sq_filtered.tsv", sep="\t"
)

tr_new = pd.read_parquet(
    "~/Desktop/combine/new/all-26q1_updated_tracking_sq_filtered.parquet"
)

set(tr_old["gene_id"]).symmetric_difference(set(tr_new["gene_id"]))

gtf_old = pd.read_table(
    "~/Desktop/combine/old/all-26q1_updated.gtf",
    sep="\t",
    comment="#",
    header=None,
    names=[
        "seqname",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attribute",
    ],
    dtype={
        "seqname": "string",
        "source": "string",
        "feature": "string",
        "start": "int64",
        "end": "int64",
        "score": "string",
        "strand": "string",
        "frame": "string",
        "attribute": "string",
    },
)

gtf_new = pd.read_table(
    "~/Desktop/combine/new/all-26q1_processed.gtf",
    sep="\t",
    comment="#",
    header=None,
    names=[
        "seqname",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attribute",
    ],
    dtype={
        "seqname": "string",
        "source": "string",
        "feature": "string",
        "start": "int64",
        "end": "int64",
        "score": "string",
        "strand": "string",
        "frame": "string",
        "attribute": "string",
    },
)

len(gtf_new)
len(gtf_old)

gtf_old["attribute"].values[0]
gtf_new["attribute"].values[0]
