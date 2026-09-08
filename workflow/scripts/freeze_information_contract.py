#!/usr/bin/env python3
"""Freeze the allowed-information and parameter contract for one adapter run."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import yaml  # type: ignore[import-untyped]


class InformationContractError(ValueError):
    """Raised when a run exposes undeclared or target-specific information."""


ALLOWED_INPUT_NAMES = {
    "short_fastq_r1", "short_fastq_r2", "reference", "reference_index",
    "pangenome_manifest", "pangenome_panel", "candidate_panel", "allele_fasta",
    "graph_assets", "tool_index",
}
REQUIRED_INFORMATION_FIELDS = {
    "allowed_inputs_manifested", "target_truth_used_for_calling",
    "target_assembly_used", "target_family_genotypes_used",
    "target_specific_external_callset", "tuning_mode", "parameter_source",
}


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise InformationContractError(f"{path} must contain an object")
    return value


def _load_yaml(path: Path) -> dict[str, Any]:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise InformationContractError(f"{path} must contain a mapping")
    return value


def _sha256_json(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _reject_plugin_scores(value: Any, location: str = "tool") -> None:
    if isinstance(value, Mapping):
        for key, item in value.items():
            normalized = str(key).casefold().replace("-", "_")
            if normalized in {"me_f1", "benchmark_score", "primary_score"}:
                raise InformationContractError(
                    f"plugins cannot provide benchmark scores: {location}.{key}"
                )
            _reject_plugin_scores(item, f"{location}.{key}")
    elif isinstance(value, Sequence) and not isinstance(value, str):
        for index, item in enumerate(value):
            _reject_plugin_scores(item, f"{location}[{index}]")


def freeze(
    tool_manifest: Path,
    resolved_inputs: Path,
    allowed_information_policy: Path,
    tuning_policy_path: Path,
) -> dict[str, dict[str, Any]]:
    tool = _load_yaml(tool_manifest)
    _reject_plugin_scores(tool)
    resolved = _load_json(resolved_inputs)
    policy = _load_yaml(allowed_information_policy)
    tuning_policy = _load_yaml(tuning_policy_path)
    if policy.get("contract") != "pgbench_allowed_information_policy_v1":
        raise InformationContractError("unsupported allowed-information policy")
    if tuning_policy.get("contract") != "pgbench_tuning_policy_v1":
        raise InformationContractError("unsupported tuning policy")
    information = tool.get("information_contract")
    if not isinstance(information, Mapping) or set(information) != REQUIRED_INFORMATION_FIELDS:
        raise InformationContractError("tool information_contract is incomplete")
    forbidden_flags = {
        "target_truth_used_for_calling", "target_assembly_used",
        "target_family_genotypes_used", "target_specific_external_callset",
    }
    if any(information.get(field) is not False for field in forbidden_flags):
        raise InformationContractError("target-specific information is forbidden")
    declared_information = tool.get("information_inputs")
    if (
        not isinstance(declared_information, Sequence)
        or isinstance(declared_information, str)
        or not declared_information
        or len(set(declared_information)) != len(declared_information)
    ):
        raise InformationContractError("tool information_inputs must be unique and non-empty")
    allowed_policy = policy.get("allowed")
    forbidden_policy = policy.get("forbidden")
    resolved_mapping = policy.get("resolved_input_information")
    if not all(
        isinstance(value, Mapping)
        for value in (allowed_policy, forbidden_policy, resolved_mapping)
    ):
        raise InformationContractError("allowed-information policy is incomplete")
    unknown_declared = sorted(
        str(name)
        for name in declared_information
        if allowed_policy.get(name) is not True
    )
    if unknown_declared:
        raise InformationContractError(
            "tool declares information outside the allowlist: "
            + ", ".join(unknown_declared)
        )
    inputs = resolved.get("inputs")
    if not isinstance(inputs, list):
        raise InformationContractError("resolved input manifest has no inputs")
    by_name: dict[str, dict[str, Any]] = {}
    for item in inputs:
        if not isinstance(item, dict) or item.get("name") not in ALLOWED_INPUT_NAMES:
            raise InformationContractError("resolved input is outside the allowlist")
        name = str(item["name"])
        if name in by_name:
            raise InformationContractError(f"duplicate resolved input: {name}")
        digest = item.get("sha256")
        if not isinstance(digest, str) or len(digest) != 64:
            raise InformationContractError(f"input {name} has no frozen SHA-256")
        by_name[name] = {
            "sha256": digest,
            "size_bytes": item.get("size_bytes"),
            "path_type": item.get("path_type"),
        }
    if not {"short_fastq_r1", "short_fastq_r2", "candidate_panel"}.issubset(by_name):
        raise InformationContractError("paired FASTQ and candidate panel are mandatory")
    required_information = {
        str(resolved_mapping[name])
        for name in by_name
        if name in resolved_mapping
    }
    missing_declarations = sorted(required_information - set(declared_information))
    if missing_declarations:
        raise InformationContractError(
            "resolved inputs use undeclared information: "
            + ", ".join(missing_declarations)
        )
    unmapped_inputs = sorted(set(by_name) - set(resolved_mapping))
    if unmapped_inputs:
        raise InformationContractError(
            "resolved inputs have no information-policy mapping: "
            + ", ".join(unmapped_inputs)
        )
    parameter_contract = tool.get("parameter_contract")
    if not isinstance(parameter_contract, Mapping):
        parameter_contract = {
            "contract": "pgbench_tool_parameters_v1",
            "source_release": {"tag": str(tool.get("version"))},
            "explicit": {},
            "implicit": {},
        }
    common = {
        "schema_version": 1,
        "tool_id": tool.get("id"),
        "tool_version": tool.get("version"),
        "tool_manifest_sha256": hashlib.sha256(tool_manifest.read_bytes()).hexdigest(),
    }
    allowed = {
        **common,
        "contract": "pgbench_allowed_information_v1",
        "allowed_inputs": sorted(by_name),
        "declared_information_inputs": sorted(map(str, declared_information)),
        "resolved_information_inputs": sorted(required_information),
        "input_contract": dict(information),
        "policy_sha256": hashlib.sha256(
            allowed_information_policy.read_bytes()
        ).hexdigest(),
        "undeclared_information": [],
        "forbidden_information": [],
        "status": "valid",
    }
    parameters = {
        **common,
        "contract": "pgbench_parameter_freeze_v1",
        "parameters": dict(parameter_contract),
        "parameter_sha256": _sha256_json(parameter_contract),
        "status": "frozen",
    }
    resources = {
        **common,
        "contract": "pgbench_external_resources_v1",
        "resources": by_name,
        "resources_sha256": _sha256_json(by_name),
        "status": "frozen",
    }
    tuning_mode = information["tuning_mode"]
    allowed_modes = tuning_policy.get("allowed_modes")
    if not isinstance(allowed_modes, Mapping) or tuning_mode not in allowed_modes:
        raise InformationContractError("tool tuning_mode is not allowed")
    mode_policy = allowed_modes[tuning_mode]
    if not isinstance(mode_policy, Mapping):
        raise InformationContractError("tuning mode policy must be a mapping")
    parameter_source = information["parameter_source"]
    if parameter_source != mode_policy.get("parameter_source_required"):
        raise InformationContractError(
            "tool parameter_source does not match its tuning mode"
        )
    tuning = {
        **common,
        "contract": "pgbench_tuning_status_v1",
        "tuning_mode": tuning_mode,
        "parameter_source": parameter_source,
        "parameter_hash": parameters["parameter_sha256"],
        "tuning_policy_sha256": hashlib.sha256(
            tuning_policy_path.read_bytes()
        ).hexdigest(),
        "target_truth_inspected": False,
        "final_me_f1_inspected_before_freeze": False,
        "status": "frozen_before_test",
    }
    return {
        "allowed_inputs": allowed,
        "parameter_manifest": parameters,
        "external_resource_hashes": resources,
        "training_or_tuning_status": tuning,
    }


def _write(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tool-manifest", required=True, type=Path)
    parser.add_argument("--resolved-inputs", required=True, type=Path)
    parser.add_argument(
        "--allowed-information-policy", required=True, type=Path
    )
    parser.add_argument("--tuning-policy", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        documents = freeze(
            args.tool_manifest,
            args.resolved_inputs,
            args.allowed_information_policy,
            args.tuning_policy,
        )
        names = {
            "allowed_inputs": "allowed_inputs.json",
            "parameter_manifest": "parameter_manifest.json",
            "external_resource_hashes": "external_resource_hashes.json",
            "training_or_tuning_status": "training_or_tuning_status.json",
        }
        for key, filename in names.items():
            _write(args.output_dir / filename, documents[key])
    except (OSError, json.JSONDecodeError, yaml.YAMLError, InformationContractError) as exc:
        print(f"freeze_information_contract: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
