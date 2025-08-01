import tomllib

import pandas as pd
from firecloud_api_cds import api as firecloud_api
from nebelung.terra_workspace import TerraWorkspace
from nebelung.types import Submissions
from nebelung.utils import call_firecloud_api, type_data_frame
from pd_flatten import pd_flatten
from tqdm import tqdm

pd.set_option("display.max_columns", 30)
pd.set_option("display.max_colwidth", 50)
pd.set_option("display.max_info_columns", 30)
pd.set_option("display.max_info_rows", 20)
pd.set_option("display.max_rows", 20)
pd.set_option("display.max_seq_items", None)
pd.set_option("display.width", 200)
pd.set_option("expand_frame_repr", True)
pd.set_option("mode.chained_assignment", "warn")

# use same config loading as when calling the module CLI
with open("config.toml", "rb") as f:
    config = tomllib.load(f)

workspace = TerraWorkspace(
    workspace_namespace=config["terra"]["workspace_namespace"],
    workspace_name=config["terra"]["workspace_name"],
)

submissions = type_data_frame(
    pd.DataFrame(
        call_firecloud_api(
            firecloud_api.list_submissions,
            namespace=workspace.workspace_namespace,
            workspace=workspace.workspace_name,
        )
    ),
    Submissions,
)

submissions = pd_flatten(submissions)

suc = submissions.loc[
    submissions["methodConfigurationName"].str.contains(r"quant|fusion")
    & submissions["workflowStatuses__Succeeded"].gt(0)
]

workflows = [
    call_firecloud_api(
        firecloud_api.get_submission,
        namespace=workspace.workspace_namespace,
        workspace=workspace.workspace_name,
        submission_id=x,
    )
    for x in suc["submissionId"]
]

wio = []

for s in tqdm(workflows):
    for w in s["workflows"]:
        if w["status"] != "Succeeded":
            continue

        r = call_firecloud_api(
            firecloud_api.get_workflow_metadata,
            namespace=workspace.workspace_namespace,
            workspace=workspace.workspace_name,
            submission_id=s["submissionId"],
            workflow_id=w["workflowId"],
            include_key=["inputs", "outputs", "workflowName"],
        )

        wio.append(
            {
                "submission_id": s["submissionId"],
                "workflow_id": w["workflowId"],
                "workflow_name": r["workflowName"],
                "ts": w["statusLastChangedDate"],
                "inputs": r["inputs"],
                "outputs": r["outputs"],
            }
        )

wio_df = pd.DataFrame(wio)
wio_df.to_parquet("./data/wio_df.parquet")

wio_df = pd_flatten(wio_df)

iodf = (
    wio_df.filter(regex=r"workflow_id|inputs__|outputs__")
    .melt(id_vars="workflow_id", var_name="iocol", value_name="url")
    .dropna()
)

iodf = iodf.loc[
    iodf["iocol"].isin(
        [
            "inputs__star_junctions",
            "inputs__sr_bam",
            "inputs__sr_cram_bam",
            "outputs__call_fusions.fusion_report",
            "outputs__call_lr_rna_fusions.fusion_report",
            "outputs__quantify_long_reads.gene_counts",
            "outputs__quantify_lr_rna.gene_counts",
        ]
    )
]

iodf = iodf.pivot(index="workflow_id", columns="iocol", values="url").reset_index()

iodf["fusion_sr_cram_bam_used"] = iodf["inputs__sr_bam"].fillna(
    iodf["inputs__sr_cram_bam"]
)

iodf["fusion_report"] = iodf["outputs__call_fusions.fusion_report"].fillna(
    iodf["outputs__call_lr_rna_fusions.fusion_report"]
)

iodf["gene_counts"] = iodf["outputs__quantify_long_reads.gene_counts"].fillna(
    iodf["outputs__quantify_lr_rna.gene_counts"]
)

iodf = iodf.rename(
    columns={"inputs__star_junctions": "quantify_sr_star_junctions_used"}
)

iodf = iodf.loc[
    :,
    [
        "workflow_id",
        "fusion_report",
        "gene_counts",
        "fusion_sr_cram_bam_used",
        "quantify_sr_star_junctions_used",
    ],
]

all_samples = workspace.get_entities("sample")

if "quantify_used_sr_evidence" not in all_samples.columns:
    all_samples["quantify_sr_star_junctions_used"] = pd.NA
    all_samples["quantify_used_sr_evidence"] = pd.NA

samples = all_samples.loc[
    :,
    [
        "sample_id",
        "fusion_report",
        "fusion_sr_cram_bam_used",
        "fusion_used_sr_evidence",
        "gene_counts",
        "quantify_sr_star_junctions_used",
        "quantify_used_sr_evidence",
    ],
].astype(
    {
        "sample_id": "string",
        "fusion_report": "string",
        "fusion_sr_cram_bam_used": "string",
        "fusion_used_sr_evidence": "boolean",
        "gene_counts": "string",
        "quantify_sr_star_junctions_used": "string",
        "quantify_used_sr_evidence": "boolean",
    }
)

# gs://fc-secure-d068ecea-bdcf-4e38-a449-edcd8933a71c/submissions/final-outputs/252275b2-1df9-4e56-b6ff-d0d59eb040db/call_lr_rna_fusions/993d269f-5bfe-45c7-a4c8-51b6d49e870b/call-ctat_lr_fusion/CDS-00nXiA.ctat-LR-fusion.fusion_predictions.tsv
samples["fusion_workflow_id"] = samples["fusion_report"].str.extract(
    r"\/([a-z 0-9 -]+)\/call-"
)
samples["quantify_workflow_id"] = samples["gene_counts"].str.extract(
    r"\/([a-z 0-9 -]+)\/call-"
)

iodf_fusion = iodf.loc[
    iodf["workflow_id"].isin(samples["fusion_workflow_id"]),
    ["workflow_id", "fusion_report", "fusion_sr_cram_bam_used"],
]

iodf_quantify = iodf.loc[
    iodf["workflow_id"].isin(samples["quantify_workflow_id"]),
    ["workflow_id", "gene_counts", "quantify_sr_star_junctions_used"],
]

samples = (
    samples.merge(
        iodf_fusion,
        how="left",
        left_on=["fusion_workflow_id", "fusion_report"],
        right_on=["workflow_id", "fusion_report"],
        suffixes=("", "_obs"),
    )
    .drop(columns="workflow_id")
    .merge(
        iodf_quantify,
        how="left",
        left_on=["quantify_workflow_id", "gene_counts"],
        right_on=["workflow_id", "gene_counts"],
        suffixes=("", "_obs"),
    )
    .drop(columns="workflow_id")
    .reset_index(drop=True)
)

fusion_updates_rows = (
    samples["fusion_used_sr_evidence"].isna()
    & samples["fusion_sr_cram_bam_used_obs"].notna()
)

samples.loc[fusion_updates_rows, "fusion_sr_cram_bam_used"] = samples.loc[
    fusion_updates_rows, "fusion_sr_cram_bam_used_obs"
]

samples.loc[fusion_updates_rows, "fusion_used_sr_evidence"] = samples.loc[
    fusion_updates_rows, "fusion_sr_cram_bam_used"
].notna()

quantify_updates_rows = (
    samples["quantify_used_sr_evidence"].isna()
    & samples["quantify_sr_star_junctions_used_obs"].notna()
)

samples.loc[quantify_updates_rows, "quantify_sr_star_junctions_used"] = samples.loc[
    quantify_updates_rows, "quantify_sr_star_junctions_used_obs"
]

samples.loc[quantify_updates_rows, "quantify_used_sr_evidence"] = samples.loc[
    quantify_updates_rows, "quantify_sr_star_junctions_used"
].notna()

samples = samples.loc[
    :,
    [
        "sample_id",
        "fusion_sr_cram_bam_used",
        "fusion_used_sr_evidence",
        "quantify_sr_star_junctions_used",
        "quantify_used_sr_evidence",
    ],
]

workspace.upload_entities(samples)
