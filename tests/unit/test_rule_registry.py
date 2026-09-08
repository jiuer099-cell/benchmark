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
    assert {"validate_evaluator_semantics", "freeze_information_contract", "build_coverage_datasets", "materialize_all_sites", "compute_me_f1"}.issubset(ids)
    for entry in registry["rules"]:
        source = entry.get("source_path")
        if source:
            assert (ROOT / source).is_file(), source


def test_registry_exposes_one_primary_score_path() -> None:
    registry = yaml.safe_load((ROOT / "workflow/rule-registry.yaml").read_text(encoding="utf-8"))
    by_id = {entry["id"]: entry for entry in registry["rules"]}
    assert by_id["compute_me_f1"]["upstream"] == [
        "evaluate_truvari",
        "evaluate_aardvark",
        "evaluate_vcfdist",
        "freeze_information_contract",
    ]
