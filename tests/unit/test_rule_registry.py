from __future__ import annotations

import re
from pathlib import Path

import yaml
from jsonschema import Draft202012Validator

ROOT = Path(__file__).resolve().parents[2]
REGISTRY_PATH = ROOT / "workflow" / "rule-registry.yaml"
SCHEMA_PATH = ROOT / "workflow" / "schemas" / "rule-registry.schema.yaml"


EXPECTED_RULE_IDS = {
    "validate_config",
    "snapshot_run_context",
    "inspect_input_bam",
    "prepare_reference",
    "prepare_primary_truth",
    "prepare_legacy_truth",
    "prepare_stratifications",
    "extract_canonical_reads",
    "split_haplotype_reads",
    "prepare_population_sv_source",
    "normalize_population_panel",
    "assign_pangenome_allele_ids",
    "lock_graph_assets",
    "build_pangenome_manifest",
    "audit_truth_leakage",
    "build_canonical_graph_recipe",
    "build_shared_alignment__{kind}",
    "match_hidden_panel_truth",
    "build_blinded_challenge_panel",
    "audit_challenge_panel",
    "tool__{tool}__preflight",
    "tool__{tool}__execute",
    "tool__{tool}__prepare_assets",
    "tool__{tool}__index",
    "tool__{tool}__map",
    "tool__{tool}__discover",
    "tool__{tool}__assemble",
    "tool__{tool}__integrate",
    "tool__{tool}__call",
    "tool__{tool}__genotype",
    "tool__{tool}__postprocess",
    "tool__{tool}__validate_raw",
    "canonicalize_vcf",
    "link_pangenome_alleles",
    "make_truvari_view",
    "make_aardvark_view",
    "make_vcfdist_view",
    "evaluate_truvari",
    "evaluate_aardvark",
    "evaluate_vcfdist",
    "parse_truvari_metrics",
    "parse_aardvark_metrics",
    "parse_vcfdist_metrics",
    "stratify_evaluator_metrics",
    "fuse_evaluator_metrics",
    "collect_tool_resources",
    "audit_score_inputs",
    "compute_pgbench_score",
    "finalize_score_provenance",
    "build_rule_lineage",
    "render_report",
    "render_report_index",
    "check_obsidian_sync",
    "sync_obsidian_design",
}


def _load_yaml(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        value = yaml.safe_load(handle)
    assert isinstance(value, dict)
    return value


def test_registry_matches_schema_and_confirmed_rule_catalog() -> None:
    registry = _load_yaml(REGISTRY_PATH)
    schema = _load_yaml(SCHEMA_PATH)
    Draft202012Validator(schema).validate(registry)
    ids = [entry["id"] for entry in registry["rules"]]
    assert len(ids) == len(set(ids)) == 54
    assert set(ids) == EXPECTED_RULE_IDS


def test_all_upstream_rules_resolve_and_all_jobs_are_auditable() -> None:
    registry = _load_yaml(REGISTRY_PATH)
    rules = {entry["id"]: entry for entry in registry["rules"]}
    for rule in rules.values():
        assert rule["manifest_required"] is True
        assert rule["log_required"] is True
        assert rule["benchmark_required"] is True
        for upstream in rule["upstream"]:
            assert upstream in rules, (rule["id"], upstream)


def test_partial_or_verified_sources_exist_when_declared() -> None:
    registry = _load_yaml(REGISTRY_PATH)
    for rule in registry["rules"]:
        source = rule.get("source_path")
        if source is None:
            continue
        if rule["implementation_status"] in {"partial", "implemented", "verified"}:
            assert (ROOT / source).exists(), (rule["id"], source)


def test_every_executable_core_rule_is_registered_and_not_planned() -> None:
    registry = _load_yaml(REGISTRY_PATH)
    rules = {entry["id"]: entry for entry in registry["rules"]}
    source_paths = [ROOT / "Snakefile", *(ROOT / "workflow" / "rules").glob("*.smk")]
    executable = {
        match.group(1)
        for path in source_paths
        for match in re.finditer(
            r"(?m)^rule\s+([A-Za-z_][A-Za-z0-9_]*)\s*:",
            path.read_text(encoding="utf-8"),
        )
    }
    executable.discard("all")
    # One Snakemake wildcard rule implements three semantically registered
    # evaluator rules; the wrapped executor records evaluate_<evaluator>.
    executable.discard("run_formal_evaluator")
    assert executable
    assert executable <= set(rules)
    for rule_id in executable:
        assert rules[rule_id]["implementation_status"] != "planned", rule_id

    external_execute = rules["tool__{tool}__execute"]
    assert external_execute["implementation_status"] in {
        "partial",
        "implemented",
        "verified",
    }
    for evaluator in ("truvari", "aardvark", "vcfdist"):
        assert rules[f"evaluate_{evaluator}"]["implementation_status"] in {
            "implemented",
            "verified",
        }
