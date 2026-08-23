from __future__ import annotations

import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from validate_config import (  # noqa: E402
    ConfigValidationError,
    _mode_semantics,
    validate_configuration,
)


def _example_config() -> dict:
    return yaml.safe_load(
        (ROOT / "config" / "config.example.yaml").read_text(encoding="utf-8")
    )


def _named_config(name: str) -> dict:
    return yaml.safe_load(
        (ROOT / "config" / name).read_text(encoding="utf-8")
    )


def test_synthetic_config_may_explicitly_allow_missing_bam(tmp_path: Path) -> None:
    config = _example_config()
    config["development"] = {
        "synthetic_mode": True,
        "canonical_fastq": "tests/fixtures/synthetic/reads.fastq",
        "population_vcf": "tests/fixtures/synthetic/population.vcf",
        "truth_vcf": "tests/fixtures/synthetic/truth.vcf",
        "benchmark_bed": "tests/fixtures/synthetic/benchmark.bed",
        "evaluator_fixture_dir": "tests/fixtures/synthetic/evaluators",
        "allow_missing_bam": True,
    }
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
    validated, plugins = validate_configuration(
        config_path,
        config_schema_path=ROOT / "config" / "config.schema.yaml",
        tool_schema_path=ROOT / "workflow" / "schemas" / "tool.schema.yaml",
        repo_root=ROOT,
    )
    assert validated["development"]["synthetic_mode"] is True
    assert set(plugins) == {"kanpig"}


def test_formal_config_rejects_missing_bam(tmp_path: Path) -> None:
    config = _example_config()
    config["sample"]["bam"] = str(tmp_path / "missing.bam")
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
    with pytest.raises(ConfigValidationError, match="does not exist"):
        validate_configuration(
            config_path,
            config_schema_path=ROOT / "config" / "config.schema.yaml",
            tool_schema_path=ROOT / "workflow" / "schemas" / "tool.schema.yaml",
            repo_root=ROOT,
        )


def test_mode_semantics_reject_overlapping_input_sets() -> None:
    manifest = {
        "id": "bad",
        "tasks": ["genotyping"],
        "inputs": {"read_sequences_auxiliary": False},
        "supported_modes": {
            "caller_only_shared_alignment": {
                "required_inputs": ["shared_alignment"],
                "optional_inputs": ["reference"],
                "forbidden_inputs": [
                    "reference",
                    "canonical_fastq",
                    "original_input_bam",
                ],
                "billable_stages": ["genotype"],
            }
        },
    }
    with pytest.raises(ConfigValidationError, match="overlapping"):
        _mode_semantics(manifest)


def test_mode_semantics_rejects_unbilled_declared_task() -> None:
    manifest = {
        "id": "bad",
        "tasks": ["genotyping"],
        "inputs": {"read_sequences_auxiliary": False},
        "supported_modes": {
            "caller_only_shared_alignment": {
                "required_inputs": ["shared_alignment"],
                "optional_inputs": [],
                "forbidden_inputs": [
                    "canonical_fastq",
                    "original_input_bam",
                ],
                "billable_stages": ["map"],
            }
        },
    }
    with pytest.raises(ConfigValidationError, match="undeclared stages"):
        _mode_semantics(manifest)


def test_pangenie_requires_complete_paired_short_reads(tmp_path: Path) -> None:
    config = _named_config("config.pangenie.example.yaml")
    config["sample"]["fastq_r2"] = None
    config_path = tmp_path / "pangenie.yaml"
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
    with pytest.raises(ConfigValidationError, match=r"sample\.fastq_r2"):
        validate_configuration(
            config_path,
            config_schema_path=ROOT / "config" / "config.schema.yaml",
            tool_schema_path=ROOT / "workflow" / "schemas" / "tool.schema.yaml",
            repo_root=ROOT,
        )


def test_plugin_technology_must_match_sample(tmp_path: Path) -> None:
    config = _named_config("config.pangenie.example.yaml")
    config["sample"]["technology"] = "pacbio_clr"
    config_path = tmp_path / "technology.yaml"
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
    with pytest.raises(ConfigValidationError, match="does not support sample technology"):
        validate_configuration(
            config_path,
            config_schema_path=ROOT / "config" / "config.schema.yaml",
            tool_schema_path=ROOT / "workflow" / "schemas" / "tool.schema.yaml",
            repo_root=ROOT,
        )


def test_vg_requires_accepted_graph_profile(tmp_path: Path) -> None:
    config = _named_config("config.vg.example.yaml")
    config["pangenome"]["build_graph_assets"] = False
    config["pangenome"]["graph_assets"]["profile"] = "none"
    for name in ("manifest", "gbz", "xg", "min", "dist", "sample_list"):
        config["pangenome"]["graph_assets"][name] = None
    config_path = tmp_path / "vg.yaml"
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
    with pytest.raises(ConfigValidationError, match="requires graph assets"):
        validate_configuration(
            config_path,
            config_schema_path=ROOT / "config" / "config.schema.yaml",
            tool_schema_path=ROOT / "workflow" / "schemas" / "tool.schema.yaml",
            repo_root=ROOT,
        )


def test_novel_truth_profile_must_exist(tmp_path: Path) -> None:
    config = _example_config()
    config["catalogs"]["novel_truth_profile"] = "config/missing-novel.yaml"
    config["development"] = {
        "synthetic_mode": True,
        "canonical_fastq": "tests/fixtures/synthetic/reads.fastq",
        "population_vcf": "tests/fixtures/synthetic/population.vcf",
        "truth_vcf": "tests/fixtures/synthetic/truth.vcf",
        "benchmark_bed": "tests/fixtures/synthetic/benchmark.bed",
        "evaluator_fixture_dir": "tests/fixtures/synthetic/evaluators",
        "allow_missing_bam": True,
    }
    config_path = tmp_path / "missing-novel.yaml"
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
    with pytest.raises(ConfigValidationError, match="novel_truth_profile"):
        validate_configuration(
            config_path,
            config_schema_path=ROOT / "config" / "config.schema.yaml",
            tool_schema_path=ROOT / "workflow" / "schemas" / "tool.schema.yaml",
            repo_root=ROOT,
        )
