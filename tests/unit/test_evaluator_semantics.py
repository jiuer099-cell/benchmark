from __future__ import annotations

import csv
import sys
from pathlib import Path

import yaml

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from validate_evaluator_semantics import EVALUATORS, validate  # noqa: E402


def test_all_three_evaluators_must_match_the_same_gt_contract(tmp_path: Path) -> None:
    cases = Path(__file__).resolve().parents[1] / "fixtures/evaluator_semantics/cases.yaml"
    suite = yaml.safe_load(cases.read_text(encoding="utf-8"))
    ledgers: dict[str, Path] = {}
    for evaluator in EVALUATORS:
        path = tmp_path / f"{evaluator}.tsv"
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["case_id", "tp", "fp", "fn"],
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            for case in suite["cases"]:
                writer.writerow({"case_id": case["id"], **case["expected"]})
        ledgers[evaluator] = path
    result = validate(cases, ledgers)
    assert result["status"] == "valid"
    assert set(result["evaluators"]) == set(EVALUATORS)
