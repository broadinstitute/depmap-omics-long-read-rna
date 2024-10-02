import pandas as pd

from dogspa_long_reads.utils.gcp import list_blobs

pd.set_option("display.max_columns", 30)
pd.set_option("display.max_colwidth", 100)
pd.set_option("display.max_info_columns", 30)
pd.set_option("display.max_info_rows", 20)
pd.set_option("display.max_rows", 20)
pd.set_option("display.max_seq_items", None)
pd.set_option("display.width", 200)
pd.set_option("expand_frame_repr", True)
pd.set_option("mode.chained_assignment", "warn")


bams = list_blobs(
    "fc-aaf4de93-c104-45c4-a01a-a036869119c6",
    glob="tag*/merge/*.bam",
)

bams = bams.rename(columns={"url": "bam_url"})
bams["model_id"] = bams["bam_url"].str.extract(r"(ACH-.+).merged.unaligned.bam$")
bams.to_csv("./data/tag.csv", index=False)
