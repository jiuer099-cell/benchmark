"""The pre-score audit and its sealing replay must audit the same job set.

``audit_score_inputs`` builds the stored audit from
``expected_pre_score_jobs()``: that list is the denominator of
``manifest_completeness`` and decides which manifests are actually validated.
``finalize_score_provenance`` later replays the stored audit with the audit's
own lineage nodes (i.e. the manifests ``pre_score_manifest_paths()`` supplied)
as the expected-job list, and rejects the seal when
``audited_job_count`` / ``valid_manifest_log_benchmark_count`` do not reproduce.

The two lists therefore have to cover exactly the same (rule, job) pairs.  They
did not: ``pre_score_manifest_paths`` gained ``validate_evaluator_semantics``,
``freeze_information_contract`` and ``materialize_evaluation_query`` while
``expected_pre_score_jobs`` did not, so the audit counted 13 jobs and the replay
16 and the seal could not be written.  These tests pin that invariant down.

The rule functions are plain Python defined in a Snakemake file, so they are
extracted from source and evaluated against a small stub namespace.  That keeps
the check structural - it compares what the two functions actually return for
the same wildcards - without needing Snakemake itself.
"""

from __future__ import annotations

import re
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

RULES = Path(__file__).resolve().parents[2] / "workflow" / "rules"
PROVENANCE_SMK = RULES / "provenance.smk"

MANIFEST_PATH_RE = re.compile(r"provenance/rules/(?P<rule>[^/]+)/(?P<job>[^/]+)\.json$")

RUN_ID = "pangenie_sr_production_v2_20260915_cachefixed"
RESULTS_ROOT = f"results/{RUN_ID}"
SAMPLE_ID = "HG002"
OFFICIAL_MODE = "end_to_end_from_reads"
PANGENOME_ID = "HG002_LOO_HPRC_GRCh38_SV_v1"
TRUTH_PROFILE = "giab_hg002_grch38_v5_0q"
TOOL = "pangenie"
CORE_JOB = f"{SAMPLE_ID}.{TOOL}.{OFFICIAL_MODE}"
FORMAL_EVALUATORS = ("truvari", "aardvark", "vcfdist")


def _extract_function(source: str, name: str) -> str:
    """Slice one top-level ``def`` out of a Snakemake file."""

    match = re.search(rf"(?m)^def {re.escape(name)}\(", source)
    if match is None:
        raise AssertionError(f"{name} is not defined in the source")
    lines = source[match.start() :].splitlines(keepends=True)
    body = [lines[0]]
    for line in lines[1:]:
        if line[:1] not in ("", " ", "\t") and (
            line.startswith("def ") or line.startswith("rule ") or line.startswith("@")
        ):
            break
        body.append(line)
    return "".join(body)


def _namespace(*, synthetic: bool, graph_assets: bool) -> dict[str, Any]:
    def semantic_rule_manifest(rule_name: str, job_key: str) -> str:
        return f"{RESULTS_ROOT}/provenance/rules/{rule_name}/{job_key}.json"

    validate_manifest = semantic_rule_manifest("validate_config", "config")
    primary_truth_manifest = semantic_rule_manifest(
        "prepare_primary_truth", TRUTH_PROFILE
    )
    return {
        "RESULTS_ROOT": RESULTS_ROOT,
        "SAMPLE_ID": SAMPLE_ID,
        "OFFICIAL_MODE": OFFICIAL_MODE,
        "PANGENOME_ID": PANGENOME_ID,
        "SYNTHETIC_MODE": synthetic,
        "GRAPH_ASSETS_ENABLED": graph_assets,
        "FORMAL_EVALUATORS": list(FORMAL_EVALUATORS),
        "config": {"truth": {"primary": TRUTH_PROFILE}},
        "VALIDATE_MANIFEST": validate_manifest,
        "CONTEXT_MANIFEST": semantic_rule_manifest("snapshot_run_context", "context"),
        "GRAPH_ASSET_RULE_MANIFEST": semantic_rule_manifest(
            "lock_graph_assets", PANGENOME_ID
        ),
        "PANGENOME_RULE_MANIFEST": semantic_rule_manifest(
            "build_pangenome_manifest", PANGENOME_ID
        ),
        "CHALLENGE_RULE_MANIFEST": semantic_rule_manifest(
            "build_blinded_challenge_panel", SAMPLE_ID
        ),
        "SEMANTIC_VALIDATION_RULE_MANIFEST": semantic_rule_manifest(
            "validate_evaluator_semantics", "contract"
        ),
        "PRIMARY_TRUTH_RULE_MANIFEST": (
            validate_manifest if synthetic else primary_truth_manifest
        ),
        "semantic_rule_manifest": semantic_rule_manifest,
        "evaluation_truth_rule_manifest": lambda wildcards: (
            validate_manifest if synthetic else primary_truth_manifest
        ),
    }


def _load_rules(*, synthetic: bool, graph_assets: bool) -> tuple[Any, Any]:
    source = PROVENANCE_SMK.read_text(encoding="utf-8")
    namespace = _namespace(synthetic=synthetic, graph_assets=graph_assets)
    for name in ("pre_score_manifest_paths", "expected_pre_score_jobs"):
        exec(_extract_function(source, name), namespace)  # noqa: S102 - test harness
    return namespace["pre_score_manifest_paths"], namespace["expected_pre_score_jobs"]


def _dynamic_metric_upstreams(wildcards: SimpleNamespace) -> list[str]:
    """Evaluate scoring.smk's fuse_evaluator_metrics upstream manifest list."""

    source = (RULES / "scoring.smk").read_text(encoding="utf-8")
    namespace = dict(_namespace(synthetic=False, graph_assets=False))
    namespace["formal_evaluator_manifests"] = lambda wc: [
        f"{RESULTS_ROOT}/provenance/rules/evaluate_{evaluator}/"
        f"{SAMPLE_ID}.{wc.tool}.{OFFICIAL_MODE}.json"
        for evaluator in FORMAL_EVALUATORS
    ]
    exec(  # noqa: S102 - test harness
        _extract_function(source, "dynamic_metric_upstreams"), namespace
    )
    return namespace["dynamic_metric_upstreams"](wildcards)


def _job_pairs(paths: list[str]) -> set[tuple[str, str]]:
    pairs: set[tuple[str, str]] = set()
    for path in paths:
        match = MANIFEST_PATH_RE.search(path)
        if match is None:
            raise AssertionError(f"manifest path is not rule-scoped: {path}")
        pairs.add((match.group("rule"), match.group("job")))
    return pairs


def _expected_pairs(values: list[str]) -> set[tuple[str, str]]:
    pairs: set[tuple[str, str]] = set()
    for value in values:
        rule_name, separator, job_key = value.partition("=")
        if not separator:
            raise AssertionError(f"expected job must use RULE=JOB_KEY: {value}")
        pairs.add((rule_name, job_key))
    return pairs


@pytest.mark.parametrize("graph_assets", [False, True])
def test_formal_manifest_and_expected_job_lists_cover_the_same_jobs(
    graph_assets: bool,
) -> None:
    manifests, expected = _load_rules(synthetic=False, graph_assets=graph_assets)
    wildcards = SimpleNamespace(tool=TOOL)

    manifest_pairs = _job_pairs(manifests(wildcards))
    expected_pairs = _expected_pairs(expected(wildcards))

    assert manifest_pairs == expected_pairs


@pytest.mark.parametrize("graph_assets", [False, True])
def test_synthetic_manifest_and_expected_job_lists_cover_the_same_jobs(
    graph_assets: bool,
) -> None:
    manifests, expected = _load_rules(synthetic=True, graph_assets=graph_assets)
    wildcards = SimpleNamespace(tool=TOOL)

    assert _job_pairs(manifests(wildcards)) == _expected_pairs(expected(wildcards))


def test_formal_ancestry_covers_every_manifest_fuse_consumes() -> None:
    """Every manifest fuse_evaluator_metrics consumes must also be audited.

    scoring.smk hands ``dynamic_metric_upstreams`` to fuse_evaluator_metrics, so
    those manifest IDs appear in its lineage closure.  If the seal does not load
    them the closure reports them as missing nodes and the seal is refused -
    the defect that made the first freezing replay abort with an empty
    "manifest path validation failed:" message.  This test reads the real
    scoring.smk function rather than a copy, so the two files cannot drift.
    """

    manifests, _expected = _load_rules(synthetic=False, graph_assets=False)
    wildcards = SimpleNamespace(tool=TOOL)
    audited = _job_pairs(manifests(wildcards))

    scored = _job_pairs(_dynamic_metric_upstreams(wildcards))

    assert scored <= audited
    assert {
        ("validate_evaluator_semantics", "contract"),
        ("freeze_information_contract", TOOL),
        ("materialize_evaluation_query", CORE_JOB),
    } <= scored


def test_extracted_lists_are_non_trivial() -> None:
    """Guard the source extractor itself against silently returning nothing."""

    manifests, expected = _load_rules(synthetic=False, graph_assets=False)
    wildcards = SimpleNamespace(tool=TOOL)

    assert len(_job_pairs(manifests(wildcards))) >= 16
    assert len(expected(wildcards)) >= 16
