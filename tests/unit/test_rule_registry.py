from __future__ import annotations

from pathlib import Path

import jsonschema
import yaml

ROOT = Path(__file__).resolve().parents[2]


def test_registry_is_valid_and_only_contains_existing_implementations() -> None:
    registry = yaml.safe_load((ROOT / "workflow/rule-registry.yaml").read_text(encoding="utf-8"))
    schema = yaml.safe_load((ROOT / "workflow/schemas/rule-registry.schema.yaml").read_text(encoding="utf-8"))
    jsonschema.Draft202012Validator(schema).validate(registry)
    ids = [entry["id"] for entry in registry["rules"]]
    assert len(ids) == len(set(ids))
    assert {"validate_evaluator_semantics", "freeze_information_contract", "build_coverage_datasets", "materialize_all_sites", "materialize_canonical_units", "normalize_pgf1_evidence", "compute_pg_f1", "validate_pgf1_leaderboard_admission", "finalize_pgf1_release"}.issubset(ids)
    for entry in registry["rules"]:
        source = entry.get("source_path")
        if source:
            assert (ROOT / source).is_file(), source


def test_registry_exposes_one_primary_score_path() -> None:
    registry = yaml.safe_load((ROOT / "workflow/rule-registry.yaml").read_text(encoding="utf-8"))
    by_id = {entry["id"]: entry for entry in registry["rules"]}
    assert by_id["compute_pg_f1"]["upstream"] == ["normalize_pgf1_evidence"]
    assert by_id["validate_pgf1_leaderboard_admission"]["upstream"] == ["compute_pg_f1"]
    assert by_id["finalize_pgf1_release"]["upstream"] == ["validate_pgf1_leaderboard_admission"]


def test_pgf1_release_path_targets_only_track_local_artifacts() -> None:
    pgf1_rules = (ROOT / "workflow/rules/pgf1.smk").read_text(encoding="utf-8")
    snakefile = (ROOT / "Snakefile").read_text(encoding="utf-8")
    assert "finalize_pgf1_release" in pgf1_rules
    assert '"/summary/{tool}/score.json"' not in pgf1_rules
    assert "release/score-package.json" in snakefile
