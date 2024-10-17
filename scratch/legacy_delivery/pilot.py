import pandas as pd

from depmap_omics_long_read_rna.utils.gcp import list_blobs

pd.set_option("display.max_columns", 30)
pd.set_option("display.max_colwidth", 100)
pd.set_option("display.max_info_columns", 30)
pd.set_option("display.max_info_rows", 20)
pd.set_option("display.max_rows", 20)
pd.set_option("display.max_seq_items", None)
pd.set_option("display.width", 200)
pd.set_option("expand_frame_repr", True)
pd.set_option("mode.chained_assignment", "warn")


pilot = list_blobs(
    "fc-aaf4de93-c104-45c4-a01a-a036869119c6", glob="pilot/merge/*.sorted.bam"
)

pilot = pilot.rename(columns={"url": "bam_url"})
pilot["model_id"] = pilot["bam_url"].str.extract(r"(ACH-.+).sorted.bam$")
pilot.to_csv("./data/pilot.csv", index=False)
