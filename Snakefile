from importlib.metadata import version
import os
import sys
from pathlib import Path

import yaml
from snakemake.exceptions import WorkflowError
from snakemake.utils import min_version, validate


min_version("9.23.1")

# Conda activation emitted by Snakemake uses the Bash `source` builtin on
# Linux. Freeze the executor shell explicitly instead of inheriting `/bin/sh`
# (dash on Ubuntu), which otherwise fails before the first rule starts.
if os.name != "nt":
    shell.executable("/bin/bash")

outputflags: update

configfile: "config/config.example.yaml"
validate(config, "config/config.schema.yaml")

CONFIG_PATH = (
    str(workflow.overwrite_configfiles[0])
    if workflow.overwrite_configfiles
    else "config/config.example.yaml"
)
RUN_ID = config["project"]["run_id"]
RESULTS_ROOT = f"results/{RUN_ID}"
LOG_ROOT = f"logs/{RUN_ID}"
BENCHMARK_ROOT = f"benchmarks/{RUN_ID}"
SAMPLE_ID = config["sample"]["id"]
OFFICIAL_MODE = config["execution"]["official_score_mode"]
PANGENOME_ID = config["pangenome"]["id"]
SYNTHETIC_MODE = config.get("development", {}).get("synthetic_mode", False)
EVALUATION_MODE = "synthetic_smoke" if SYNTHETIC_MODE else "formal"
SNAKEMAKE_VERSION = version("snakemake")
# Resolve Python after Snakemake activates each rule's declared environment.
PYTHON_EXECUTABLE = sys.executable if os.name == "nt" else "python"
RULE_EXECUTOR = "workflow/scripts/pgbench_rule_exec.py"
PROVENANCE_LIBRARY = "workflow/scripts/pgbench_provenance.py"
METRICS_LIBRARY = "workflow/scripts/pgbench_metrics.py"
SCORING_LIBRARY = "workflow/scripts/pgbench_scoring.py"
CORE_ENV_SPEC = "workflow/envs/core.yaml"
VALIDATE_MANIFEST = RESULTS_ROOT + "/provenance/rules/validate_config/config.json"
CONTEXT_MANIFEST = (
    RESULTS_ROOT + "/provenance/rules/snapshot_run_context/context.json"
)
PANGENOME_RULE_MANIFEST = (
    f"{RESULTS_ROOT}/provenance/rules/build_pangenome_manifest/"
    f"{PANGENOME_ID}.json"
)
GRAPH_PROFILE = config["pangenome"]["graph_assets"]["profile"]


def selected_mode_uses_graph(registration):
    with open(registration["manifest"], encoding="utf-8") as handle:
        manifest = yaml.safe_load(handle)
    contract = manifest.get("supported_modes", {}).get(OFFICIAL_MODE, {})
    return "graph_assets" in {
        *contract.get("required_inputs", []),
        *contract.get("optional_inputs", []),
    }


GRAPH_ASSETS_ENABLED = bool(
    config["pangenome"]["build_graph_assets"]
    and any(
        selected_mode_uses_graph(registration)
        for registration in config["external_plugins"]
    )
)
GRAPH_ASSET_LOCK = (
    f"{RESULTS_ROOT}/pangenome/{PANGENOME_ID}/graph-assets.lock.yaml"
)
GRAPH_ASSET_RULE_MANIFEST = (
    f"{RESULTS_ROOT}/provenance/rules/lock_graph_assets/{PANGENOME_ID}.json"
)
CHALLENGE_RULE_MANIFEST = (
    f"{RESULTS_ROOT}/provenance/rules/build_blinded_challenge_panel/"
    f"{SAMPLE_ID}.json"
)
SEMANTIC_VALIDATION_JSON = RESULTS_ROOT + "/validation/evaluator-semantics.json"
SEMANTIC_VALIDATION_RULE_MANIFEST = (
    RESULTS_ROOT + "/provenance/rules/validate_evaluator_semantics/contract.json"
)
COVERAGE_MANIFEST = RESULTS_ROOT + "/coverage/coverage-manifest.json"
COVERAGE_RULE_MANIFEST = (
    RESULTS_ROOT + "/provenance/rules/build_coverage_datasets/coverage.json"
)
if SYNTHETIC_MODE:
    PRIMARY_TRUTH_VCF = config["development"]["truth_vcf"]
    PRIMARY_TRUTH_INDEX = None
    PRIMARY_TRUTH_AUDIT = None
    PRIMARY_TRUTH_RULE_MANIFEST = VALIDATE_MANIFEST
else:
    PRIMARY_TRUTH_VCF = (
        f"{RESULTS_ROOT}/truth/{config['truth']['primary']}/sv.truth.vcf.gz"
    )
    PRIMARY_TRUTH_INDEX = f"{PRIMARY_TRUTH_VCF}.tbi"
    PRIMARY_TRUTH_AUDIT = (
        f"{RESULTS_ROOT}/truth/{config['truth']['primary']}/audit.json"
    )
    PRIMARY_TRUTH_RULE_MANIFEST = (
        f"{RESULTS_ROOT}/provenance/rules/prepare_primary_truth/"
        f"{config['truth']['primary']}.json"
    )


def cli_repeated(flag, values):
    return [
        item
        for value in values
        for item in (flag, str(value))
    ]


def semantic_rule_manifest(rule_name, job_key):
    return f"{RESULTS_ROOT}/provenance/rules/{rule_name}/{job_key}.json"

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
    output_dir = f"{RESULTS_ROOT}/{SAMPLE_ID}/{OFFICIAL_MODE}/{tool_id}"
    output_relative = tool_manifest["outputs"]["vcf"]
    raw_output = f"{output_dir}/{output_relative}"
    settings = {
        "tool_id": tool_id,
        "semantic_rule_name": f"tool__{tool_id}__execute",
        "job_key": f"{SAMPLE_ID}.{tool_id}.{OFFICIAL_MODE}",
        "tool_manifest": manifest_path,
        "comparison_task": tool_manifest["comparison_task"],
        "executor": "workflow/scripts/pgbench_exec.py",
        "schema": "workflow/schemas/tool.schema.yaml",
        "rule_executor": RULE_EXECUTOR,
        "config_snapshot": CONFIG_PATH,
        "score_profile": config["catalogs"]["score_weights"],
        "run_context": RESULTS_ROOT + "/provenance/run-context.json",
        "pangenome_manifest": (
            f"{RESULTS_ROOT}/pangenome/{PANGENOME_ID}/manifest.yaml"
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
            f"{LOG_ROOT}/tools/{tool_id}/{SAMPLE_ID}.{OFFICIAL_MODE}.log"
        ),
        "rule_manifest": semantic_rule_manifest(
            f"tool__{tool_id}__execute",
            f"{SAMPLE_ID}.{tool_id}.{OFFICIAL_MODE}",
        ),
        "log": (
            f"{LOG_ROOT}/rules/tool__{tool_id}__execute/"
            f"{SAMPLE_ID}.{tool_id}.{OFFICIAL_MODE}.log"
        ),
        "benchmark": (
            f"{BENCHMARK_ROOT}/rules/tool__{tool_id}__execute/"
            f"{SAMPLE_ID}.{tool_id}.{OFFICIAL_MODE}.jsonl"
        ),
        "threads": config["execution"].get("tool_threads", 1),
        "memory_mb": config["execution"].get("tool_memory_mb", 1024),
        "timeout_seconds": config["execution"].get(
            "tool_timeout_seconds", 120
        ),
        "benchmark_repeats": config["execution"]["benchmark_repeats"],
        "cache_policy": config["execution"]["cache_policy"],
        "execution_purpose": (
            "development_only" if SYNTHETIC_MODE else "formal"
        ),
        "inputs": {},
    }
    mode_contract = tool_manifest["supported_modes"].get(OFFICIAL_MODE, {})
    allowed_inputs = {
        *mode_contract.get("required_inputs", []),
        *mode_contract.get("optional_inputs", []),
    }
    input_candidates = {
        "reference": config["reference"]["fasta"],
        "reference_index": config["reference"].get("fai"),
        "pangenome_manifest": (
            f"{RESULTS_ROOT}/pangenome/{PANGENOME_ID}/manifest.yaml"
        ),
        "pangenome_panel": (
            f"{RESULTS_ROOT}/pangenome/{PANGENOME_ID}/panel.vcf"
        ),
        "candidate_panel": (
            f"{RESULTS_ROOT}/pangenome/{PANGENOME_ID}/challenge/"
            f"{SAMPLE_ID}.blinded.vcf"
        ),
        "tool_index": registration.get("tool_index"),
        "short_fastq_r1": config["sample"]["fastq_r1"],
        "short_fastq_r2": config["sample"]["fastq_r2"],
    }
    settings["inputs"] = {
        name: path
        for name, path in input_candidates.items()
        if name in allowed_inputs and path
    }
    graph_allowed = "graph_assets" in {
        *mode_contract.get("required_inputs", []),
        *mode_contract.get("optional_inputs", []),
    }
    if graph_allowed and GRAPH_ASSETS_ENABLED:
        graph_config = config["pangenome"]["graph_assets"]
        graph_path = next(
            (
                graph_config.get(name)
                for name in ("gbz", "gfa", "xg", "min", "dist")
                if graph_config.get(name)
            ),
            None,
        )
        if graph_path is None:
            raise WorkflowError("enabled graph profile has no graph asset path")
        settings["inputs"]["graph_assets"] = str(
            Path(graph_path).parent
        )
        settings["graph_asset_manifest"] = graph_config["manifest"]
        settings["graph_asset_lock"] = GRAPH_ASSET_LOCK
        settings["graph_profile"] = GRAPH_PROFILE
        settings["upstream_manifests"].append(GRAPH_ASSET_RULE_MANIFEST)
    EXTERNAL_SETTINGS.append(settings)

TOOL_SETTINGS_BY_ID = {
    settings["tool_id"]: settings for settings in EXTERNAL_SETTINGS
}


def tool_comparison_task(wildcards):
    try:
        return TOOL_SETTINGS_BY_ID[wildcards.tool]["comparison_task"]
    except KeyError as error:
        raise WorkflowError(
            f"no comparison task for tool {wildcards.tool!r}"
        ) from error


def evaluation_truth_vcf(wildcards):
    return PRIMARY_TRUTH_VCF


def evaluation_truth_index(wildcards):
    return PRIMARY_TRUTH_INDEX


def evaluation_truth_rule_manifest(wildcards):
    return PRIMARY_TRUTH_RULE_MANIFEST


def evaluation_truth_profile(wildcards):
    return config["truth"]["primary"]

RAW_OUTPUT_BY_TOOL = {
    settings["tool_id"]: settings["raw_output"]
    for settings in EXTERNAL_SETTINGS
}
TOOL_BENCHMARK_BY_TOOL = {
    settings["tool_id"]: settings["benchmark"]
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
    f"{RESULTS_ROOT}/report/tool_cards/{settings['tool_id']}.html"
    for settings in EXTERNAL_SETTINGS
]
SCORE_JSONS = [
    f"{RESULTS_ROOT}/summary/{settings['tool_id']}/score.json"
    for settings in EXTERNAL_SETTINGS
]
METRICS_JSONS = [
    f"{RESULTS_ROOT}/summary/{settings['tool_id']}/metrics.json"
    for settings in EXTERNAL_SETTINGS
]
FINAL_SCORE_PACKAGES = [
    f"{RESULTS_ROOT}/summary/{settings['tool_id']}/score-package.json"
    for settings in EXTERNAL_SETTINGS
]
TOOL_MANIFESTS = [
    settings["tool_manifest"]
    for settings in EXTERNAL_SETTINGS
]
if EXTERNAL_SETTINGS:
    FINAL_TARGETS = [
        RESULTS_ROOT + "/report/index.html",
        RESULTS_ROOT + "/summary/score.tsv",
        RESULTS_ROOT + "/summary/point_breakdown.tsv",
        RESULTS_ROOT + "/summary/metrics.long.tsv",
        RESULTS_ROOT + "/summary/metrics.json",
        SEMANTIC_VALIDATION_JSON,
        COVERAGE_MANIFEST,
        *SCORE_JSONS,
        *FINAL_SCORE_PACKAGES,
        *METRICS_JSONS,
        *[
            f"{RESULTS_ROOT}/provenance/{settings['tool_id']}/contracts/{name}"
            for settings in EXTERNAL_SETTINGS
            for name in (
                "allowed_inputs.json",
                "parameter_manifest.json",
                "external_resource_hashes.json",
                "training_or_tuning_status.json",
            )
        ],
        *[
            f"{RESULTS_ROOT}/{SAMPLE_ID}/{OFFICIAL_MODE}/"
            f"{settings['tool_id']}/canonical/all-sites.vcf"
            for settings in EXTERNAL_SETTINGS
        ],
        *[
            f"{RESULTS_ROOT}/{SAMPLE_ID}/{OFFICIAL_MODE}/"
            f"{settings['tool_id']}/canonical/allele_links.tsv"
            for settings in EXTERNAL_SETTINGS
        ],
        *[
            f"{RESULTS_ROOT}/provenance/{settings['tool_id']}/{name}"
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
    FINAL_TARGETS = [RESULTS_ROOT + "/provenance/config.validated.json"]


rule all:
    input:
        FINAL_TARGETS


include: "workflow/rules/common.smk"
include: "workflow/rules/pangenome.smk"
include: "workflow/rules/panel.smk"
include: "workflow/rules/contracts.smk"
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
