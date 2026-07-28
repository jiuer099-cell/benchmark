from importlib.metadata import version
from pathlib import Path
import sys

import yaml
from snakemake.exceptions import WorkflowError
from snakemake.utils import min_version, validate


min_version("9.23.1")

outputflags: update

configfile: "config/config.example.yaml"
validate(config, "config/config.schema.yaml")

CONFIG_PATH = (
    str(workflow.overwrite_configfiles[0])
    if workflow.overwrite_configfiles
    else "config/config.example.yaml"
)
RUN_ID = config["project"]["run_id"]
SAMPLE_ID = config["sample"]["id"]
OFFICIAL_MODE = config["execution"]["official_score_mode"]
PANGENOME_ID = config["pangenome"]["id"]
SYNTHETIC_MODE = config.get("development", {}).get("synthetic_mode", False)
EVALUATION_MODE = "synthetic_smoke" if SYNTHETIC_MODE else "formal"
SNAKEMAKE_VERSION = version("snakemake")
PYTHON_EXECUTABLE = sys.executable
RULE_EXECUTOR = "workflow/scripts/pgbench_rule_exec.py"
PROVENANCE_LIBRARY = "workflow/scripts/pgbench_provenance.py"
METRICS_LIBRARY = "workflow/scripts/pgbench_metrics.py"
SCORING_LIBRARY = "workflow/scripts/pgbench_scoring.py"
CORE_ENV_SPEC = "workflow/envs/core.yaml"
VALIDATE_MANIFEST = "results/provenance/rules/validate_config/config.json"
CONTEXT_MANIFEST = (
    "results/provenance/rules/snapshot_run_context/context.json"
)
PANGENOME_RULE_MANIFEST = (
    f"results/provenance/rules/build_pangenome_manifest/"
    f"{PANGENOME_ID}.json"
)
GRAPH_ASSETS_ENABLED = bool(config["pangenome"]["build_graph_assets"])
GRAPH_ASSET_LOCK = (
    f"results/pangenome/{PANGENOME_ID}/graph-assets.lock.yaml"
)
GRAPH_ASSET_RULE_MANIFEST = (
    f"results/provenance/rules/lock_graph_assets/{PANGENOME_ID}.json"
)
CHALLENGE_RULE_MANIFEST = (
    f"results/provenance/rules/build_blinded_challenge_panel/"
    f"{SAMPLE_ID}.json"
)


def cli_repeated(flag, values):
    return [
        item
        for value in values
        for item in (flag, str(value))
    ]


def semantic_rule_manifest(rule_name, job_key):
    return f"results/provenance/rules/{rule_name}/{job_key}.json"

EXTERNAL_SETTINGS = []
for registration in config["external_plugins"]:
    tool_id = registration["id"]
    manifest_path = registration["manifest"]
    with open(manifest_path, encoding="utf-8") as handle:
        tool_manifest = yaml.safe_load(handle)
    if tool_manifest["id"] != tool_id:
        raise WorkflowError(
            f"external plugin ID {tool_id!r} does not match {manifest_path}"
        )
    output_dir = f"results/{SAMPLE_ID}/{OFFICIAL_MODE}/{tool_id}"
    output_relative = tool_manifest["outputs"]["vcf"]
    raw_output = f"{output_dir}/{output_relative}"
    settings = {
        "tool_id": tool_id,
        "semantic_rule_name": f"tool__{tool_id}__execute",
        "job_key": f"{SAMPLE_ID}.{OFFICIAL_MODE}",
        "tool_manifest": manifest_path,
        "executor": "workflow/scripts/pgbench_exec.py",
        "schema": "workflow/schemas/tool.schema.yaml",
        "rule_executor": RULE_EXECUTOR,
        "config_snapshot": CONFIG_PATH,
        "score_profile": config["catalogs"]["score_weights"],
        "run_context": "results/provenance/run-context.json",
        "pangenome_manifest": (
            f"results/pangenome/{PANGENOME_ID}/manifest.yaml"
        ),
        "graph_asset_manifest": None,
        "graph_asset_lock": None,
        "reference": config["reference"]["fasta"],
        "truth_profile": config["truth"]["primary"],
        "snakemake_version": SNAKEMAKE_VERSION,
        "execution_profile": "local",
        "random_seed": config["execution"]["random_seed"],
        "upstream_manifests": [
            PANGENOME_RULE_MANIFEST,
            CHALLENGE_RULE_MANIFEST,
        ],
        "mode": OFFICIAL_MODE,
        "run_id": RUN_ID,
        "sample_id": SAMPLE_ID,
        "output_dir": output_dir,
        "raw_output": raw_output,
        "output_vcf_relative": output_relative,
        "resolved_inputs": f"{output_dir}/meta/resolved_inputs.json",
        "attempt_record": f"{output_dir}/meta/attempt.json",
        "tool_log": (
            f"logs/tools/{tool_id}/{SAMPLE_ID}.{OFFICIAL_MODE}.log"
        ),
        "rule_manifest": semantic_rule_manifest(
            f"tool__{tool_id}__execute",
            f"{SAMPLE_ID}.{OFFICIAL_MODE}",
        ),
        "log": (
            f"logs/rules/tool__{tool_id}__execute/"
            f"{SAMPLE_ID}.{OFFICIAL_MODE}.log"
        ),
        "benchmark": (
            f"benchmarks/rules/tool__{tool_id}__execute/"
            f"{SAMPLE_ID}.{OFFICIAL_MODE}.jsonl"
        ),
        "threads": config["execution"].get("tool_threads", 1),
        "memory_mb": config["execution"].get("tool_memory_mb", 1024),
        "timeout_seconds": config["execution"].get(
            "tool_timeout_seconds", 120
        ),
        "execution_purpose": (
            "development_only" if SYNTHETIC_MODE else "formal"
        ),
        "alignment_kind": (
            config.get("caller_only", {}).get("alignment_kind")
            if OFFICIAL_MODE == "caller_only_shared_alignment"
            else None
        ),
        "inputs": {},
    }
    if OFFICIAL_MODE == "end_to_end_from_reads":
        canonical_fastq = (
            config["development"]["canonical_fastq"]
            if SYNTHETIC_MODE
            else config["sample"]["fastq"]
        )
        settings["inputs"] = {
            "canonical_fastq": canonical_fastq,
            "reference": config["reference"]["fasta"],
            "pangenome_manifest": (
                f"results/pangenome/{PANGENOME_ID}/manifest.yaml"
            ),
            "pangenome_panel": (
                f"results/pangenome/{PANGENOME_ID}/panel.vcf"
            ),
            "candidate_panel": (
                f"results/pangenome/{PANGENOME_ID}/challenge/"
                f"{SAMPLE_ID}.blinded.vcf"
            ),
        }
    else:
        shared_bam = config["caller_only"].get("shared_bam")
        settings["inputs"] = {
            "shared_alignment": shared_bam,
            "reference": config["reference"]["fasta"],
            "pangenome_manifest": (
                f"results/pangenome/{PANGENOME_ID}/manifest.yaml"
            ),
            "pangenome_panel": (
                f"results/pangenome/{PANGENOME_ID}/panel.vcf"
            ),
            "candidate_panel": (
                f"results/pangenome/{PANGENOME_ID}/challenge/"
                f"{SAMPLE_ID}.blinded.vcf"
            ),
        }
    mode_contract = tool_manifest["supported_modes"].get(OFFICIAL_MODE, {})
    graph_allowed = "graph_assets" in {
        *mode_contract.get("required_inputs", []),
        *mode_contract.get("optional_inputs", []),
    }
    if graph_allowed and GRAPH_ASSETS_ENABLED:
        graph_config = config["pangenome"]["graph_assets"]
        settings["inputs"]["graph_assets"] = str(
            Path(graph_config["gbz"]).parent
        )
        settings["graph_asset_manifest"] = graph_config["manifest"]
        settings["graph_asset_lock"] = GRAPH_ASSET_LOCK
        settings["upstream_manifests"].append(GRAPH_ASSET_RULE_MANIFEST)
    EXTERNAL_SETTINGS.append(settings)

RAW_OUTPUT_BY_TOOL = {
    settings["tool_id"]: settings["raw_output"]
    for settings in EXTERNAL_SETTINGS
}


def external_raw_vcf(wildcards):
    """Resolve the manifest-declared raw VCF for one registered tool."""

    try:
        return RAW_OUTPUT_BY_TOOL[wildcards.tool]
    except KeyError as error:
        raise WorkflowError(
            f"no external plugin is registered for tool {wildcards.tool!r}"
        ) from error

if SYNTHETIC_MODE and len(EXTERNAL_SETTINGS) != 1:
    raise WorkflowError(
        "Phase 1 synthetic workflow requires exactly one external plugin"
    )

REPORT_CARDS = [
    f"results/report/tool_cards/{settings['tool_id']}.html"
    for settings in EXTERNAL_SETTINGS
]
SCORE_JSONS = [
    f"results/summary/{settings['tool_id']}/score.json"
    for settings in EXTERNAL_SETTINGS
]
METRICS_JSONS = [
    f"results/summary/{settings['tool_id']}/metrics.json"
    for settings in EXTERNAL_SETTINGS
]
FINAL_SCORE_PACKAGES = [
    f"results/summary/{settings['tool_id']}/score-package.json"
    for settings in EXTERNAL_SETTINGS
]
TOOL_MANIFESTS = [
    settings["tool_manifest"]
    for settings in EXTERNAL_SETTINGS
]
if EXTERNAL_SETTINGS:
    FINAL_TARGETS = [
        "results/report/index.html",
        "results/summary/score.tsv",
        "results/summary/point_breakdown.tsv",
        "results/summary/metrics.long.tsv",
        "results/summary/metrics.json",
        *SCORE_JSONS,
        *FINAL_SCORE_PACKAGES,
        *METRICS_JSONS,
        *[
            f"results/{SAMPLE_ID}/{OFFICIAL_MODE}/"
            f"{settings['tool_id']}/canonical/linked.vcf"
            for settings in EXTERNAL_SETTINGS
        ],
        *[
            f"results/{SAMPLE_ID}/{OFFICIAL_MODE}/"
            f"{settings['tool_id']}/canonical/allele_links.tsv"
            for settings in EXTERNAL_SETTINGS
        ],
        *[
            f"results/provenance/{settings['tool_id']}/{name}"
            for settings in EXTERNAL_SETTINGS
            for name in (
                "pre-score-audit.json",
                "pre-score-lineage.json",
                "pre-score-lineage.tsv",
                "rule-lineage.json",
                "rule-lineage.tsv",
                "provenance-audit.json",
            )
        ],
    ]
else:
    FINAL_TARGETS = ["results/provenance/config.validated.json"]


rule all:
    input:
        FINAL_TARGETS


include: "workflow/rules/common.smk"
include: "workflow/rules/pangenome.smk"
include: "workflow/rules/panel.smk"
include: "workflow/rules/normalization.smk"
include: "workflow/rules/evaluation.smk"
include: "workflow/rules/provenance.smk"
include: "workflow/rules/scoring.smk"
include: "workflow/rules/report.smk"
include: "workflow/rules/obsidian.smk"


for settings in EXTERNAL_SETTINGS:
    module_name = f"external_{settings['tool_id']}"
    alias = f"{module_name}_"

    module:
        name: module_name
        snakefile: "workflow/modules/generic_external/Snakefile"
        config: {"generic_external": settings}

    use rule * from module_name as alias*
