#!/usr/bin/env python3
"""Validate PGBench YAML configuration and registered external tools."""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import yaml  # type: ignore[import-untyped]
from jsonschema import Draft202012Validator, FormatChecker


class ConfigValidationError(ValueError):
    """Raised when schema-valid data still violates semantic constraints."""


TASK_TO_STAGE = {
    "prepare_assets": "prepare_assets",
    "index": "index",
    "mapping": "map",
    "discovery": "discover",
    "assembly": "assemble",
    "integration": "integrate",
    "calling": "call",
    "genotyping": "genotype",
    "postprocess": "postprocess",
}
READ_INPUTS = {"canonical_fastq", "short_fastq_r1", "short_fastq_r2"}
GRAPH_PROFILE_ASSETS = {
    "none": set(),
    "vg_gbz_min_dist": {"manifest", "gbz", "min", "dist", "sample_list"},
    "vg_giraffe_shortread": {
        "manifest", "gbz", "min", "zipcodes", "dist", "sample_list"
    },
    "vg_legacy_xg": {"manifest", "gbz", "xg", "min", "dist", "sample_list"},
}


def load_yaml(path: Path) -> dict[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as handle:
            loaded = yaml.safe_load(handle)
    except (OSError, yaml.YAMLError) as exc:
        raise ConfigValidationError(f"cannot load YAML {path}: {exc}") from exc
    if not isinstance(loaded, dict):
        raise ConfigValidationError(f"{path} must contain a YAML mapping")
    return loaded


def validate_json_schema(instance: Mapping[str, Any], schema_path: Path) -> None:
    schema = load_yaml(schema_path)
    validator = Draft202012Validator(schema, format_checker=FormatChecker())
    errors = sorted(
        validator.iter_errors(instance),
        key=lambda error: [str(item) for item in error.absolute_path],
    )
    if errors:
        formatted = []
        for error in errors:
            location = ".".join(str(item) for item in error.absolute_path) or "$"
            formatted.append(f"{location}: {error.message}")
        raise ConfigValidationError("; ".join(formatted))


def _mode_semantics(tool: Mapping[str, Any]) -> None:
    inputs = tool["inputs"]
    for mode, contract in tool["supported_modes"].items():
        required = set(contract["required_inputs"])
        optional = set(contract["optional_inputs"])
        forbidden = set(contract["forbidden_inputs"])
        overlaps = {
            "required/optional": required & optional,
            "required/forbidden": required & forbidden,
            "optional/forbidden": optional & forbidden,
        }
        conflicts = {name: values for name, values in overlaps.items() if values}
        if conflicts:
            raise ConfigValidationError(
                f"tool {tool['id']} mode {mode} has overlapping input sets: {conflicts}"
            )
        allowed = required | optional
        if mode == "caller_only_shared_alignment":
            if "shared_alignment" not in required:
                raise ConfigValidationError(
                    f"tool {tool['id']} caller-only mode requires shared_alignment"
                )
            if "original_input_bam" not in forbidden:
                raise ConfigValidationError(
                    f"tool {tool['id']} caller-only mode must forbid original_input_bam"
                )
            declares_paired_contract = bool(
                {"short_fastq_r1", "short_fastq_r2"}
                & (required | optional | forbidden)
            )
            if (
                not inputs.get("read_sequences_auxiliary", False)
                and declares_paired_contract
                and not READ_INPUTS.issubset(forbidden)
            ):
                raise ConfigValidationError(
                    f"tool {tool['id']} caller-only mode must forbid all FASTQ inputs"
                )
        if mode == "end_to_end_from_reads":
            has_single = "canonical_fastq" in required
            has_pair = {"short_fastq_r1", "short_fastq_r2"}.issubset(required)
            has_partial_pair = bool(
                {"short_fastq_r1", "short_fastq_r2"} & required
            ) and not has_pair
            if not has_single and not has_pair:
                raise ConfigValidationError(
                    f"tool {tool['id']} end-to-end mode requires canonical_fastq "
                    "or the complete short_fastq_r1/short_fastq_r2 pair"
                )
            if has_partial_pair:
                raise ConfigValidationError(
                    f"tool {tool['id']} must require both paired FASTQ inputs"
                )
            for forbidden_input in ("original_input_bam", "shared_alignment"):
                if forbidden_input not in forbidden:
                    raise ConfigValidationError(
                        f"tool {tool['id']} end-to-end mode must forbid "
                        f"{forbidden_input}"
                    )
        if allowed & forbidden:
            raise ConfigValidationError(
                f"tool {tool['id']} mode {mode} exposes forbidden inputs"
            )

        declared_stages = set(contract["billable_stages"])
        permitted_stages = {TASK_TO_STAGE[task] for task in tool["tasks"]}
        if not declared_stages.issubset(permitted_stages):
            raise ConfigValidationError(
                f"tool {tool['id']} mode {mode} bills undeclared stages: "
                f"{sorted(declared_stages - permitted_stages)}"
            )
        if not permitted_stages.intersection(declared_stages):
            raise ConfigValidationError(
                f"tool {tool['id']} mode {mode} has no billable declared task"
            )


def validate_external_plugins(
    config: Mapping[str, Any],
    *,
    repo_root: Path,
    tool_schema_path: Path,
) -> dict[str, dict[str, Any]]:
    enabled_tools = set(config["tools"]["enabled"])
    planned_tools = set(config["tools"]["planned_tools"])
    if enabled_tools & planned_tools:
        raise ConfigValidationError(
            "tools.enabled and tools.planned_tools must be disjoint"
        )

    loaded: dict[str, dict[str, Any]] = {}
    for registration in config["external_plugins"]:
        plugin_id = registration["id"]
        if plugin_id in loaded:
            raise ConfigValidationError(f"duplicate external plugin ID {plugin_id}")
        if plugin_id in enabled_tools or plugin_id in planned_tools:
            raise ConfigValidationError(
                f"external plugin ID conflicts with built-in tool ID {plugin_id}"
            )
        manifest_path = (repo_root / registration["manifest"]).resolve()
        manifest = load_yaml(manifest_path)
        validate_json_schema(manifest, tool_schema_path)
        if manifest["id"] != plugin_id:
            raise ConfigValidationError(
                f"config plugin ID {plugin_id!r} does not match manifest "
                f"ID {manifest['id']!r}"
            )
        _mode_semantics(manifest)
        loaded[plugin_id] = manifest
    return loaded


def _validate_graph_profile(config: Mapping[str, Any]) -> None:
    pangenome = config["pangenome"]
    graph = pangenome["graph_assets"]
    profile = graph["profile"]
    if profile not in GRAPH_PROFILE_ASSETS:
        raise ConfigValidationError(f"unsupported graph asset profile {profile!r}")
    enabled = bool(pangenome["build_graph_assets"])
    if enabled and profile == "none":
        raise ConfigValidationError(
            "build_graph_assets=true requires a non-none graph profile"
        )
    if not enabled and profile != "none":
        raise ConfigValidationError(
            "build_graph_assets=false requires graph_assets.profile=none"
        )
    missing = sorted(
        name
        for name in GRAPH_PROFILE_ASSETS[profile]
        if not isinstance(graph.get(name), str) or not graph[name]
    )
    if missing:
        raise ConfigValidationError(
            f"graph profile {profile} is missing required assets: {missing}"
        )


def _validate_plugin_compatibility(
    config: Mapping[str, Any],
    plugins: Mapping[str, Mapping[str, Any]],
) -> None:
    technology = config["sample"]["technology"]
    mode = config["execution"]["official_score_mode"]
    synthetic = bool(config.get("development", {}).get("synthetic_mode", False))
    graph_profile = config["pangenome"]["graph_assets"]["profile"]
    for plugin_id, manifest in plugins.items():
        supported_technologies = set(manifest["capabilities"]["technology"])
        if technology not in supported_technologies:
            raise ConfigValidationError(
                f"tool {plugin_id} does not support sample technology {technology}; "
                f"supported={sorted(supported_technologies)}"
            )
        contract = manifest["supported_modes"].get(mode)
        if not isinstance(contract, Mapping):
            raise ConfigValidationError(
                f"tool {plugin_id} does not support official mode {mode}"
            )
        required = set(contract["required_inputs"])
        if mode == "end_to_end_from_reads":
            canonical = (
                config["development"].get("canonical_fastq")
                if synthetic
                else config["sample"].get("fastq")
            )
            if "canonical_fastq" in required and not canonical:
                raise ConfigValidationError(
                    f"tool {plugin_id} requires sample.fastq"
                )
            for field, contract_name in (
                ("fastq_r1", "short_fastq_r1"),
                ("fastq_r2", "short_fastq_r2"),
            ):
                if contract_name in required and not config["sample"].get(field):
                    raise ConfigValidationError(
                        f"tool {plugin_id} requires sample.{field}"
                    )
        if "graph_assets" in required:
            if graph_profile == "none":
                raise ConfigValidationError(
                    f"tool {plugin_id} requires graph assets"
                )
            accepted = set(
                manifest["inputs"]["graph_assets"].get("accepted_profiles", [])
            )
            if accepted and graph_profile not in accepted:
                raise ConfigValidationError(
                    f"tool {plugin_id} rejects graph profile {graph_profile}; "
                    f"accepted={sorted(accepted)}"
                )


def validate_configuration(
    config_path: Path,
    *,
    config_schema_path: Path,
    tool_schema_path: Path,
    repo_root: Path,
) -> tuple[dict[str, Any], dict[str, dict[str, Any]]]:
    config = load_yaml(config_path)
    validate_json_schema(config, config_schema_path)
    plugins = validate_external_plugins(
        config,
        repo_root=repo_root,
        tool_schema_path=tool_schema_path,
    )
    _validate_graph_profile(config)
    _validate_plugin_compatibility(config, plugins)
    evaluator_profile = (repo_root / config["catalogs"]["evaluator_profile"]).resolve()
    if not evaluator_profile.is_file():
        raise ConfigValidationError(
            f"evaluator profile does not exist: {evaluator_profile}"
        )
    development = config.get("development", {})
    needs_bam = (
        config["execution"]["official_score_mode"]
        == "caller_only_shared_alignment"
    )
    bam_value = config["sample"].get("bam")
    bam_path = Path(bam_value) if isinstance(bam_value, str) else None
    if (
        needs_bam
        and (
            not development.get("synthetic_mode", False)
            or not development.get("allow_missing_bam", False)
        )
    ) and (bam_path is None or not bam_path.exists()):
        raise ConfigValidationError(
            f"sample BAM does not exist: {bam_path}; use synthetic_mode only "
            "for local fixtures"
        )
    return config, plugins


def _atomic_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True, ensure_ascii=False)
        handle.write("\n")
    temporary.replace(path)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate PGBench configuration.")
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument(
        "--config-schema",
        type=Path,
        default=Path("config/config.schema.yaml"),
    )
    parser.add_argument(
        "--tool-schema",
        type=Path,
        default=Path("workflow/schemas/tool.schema.yaml"),
    )
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--output", type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        config, plugins = validate_configuration(
            args.config,
            config_schema_path=args.config_schema,
            tool_schema_path=args.tool_schema,
            repo_root=args.repo_root.resolve(),
        )
        if args.output is not None:
            _atomic_json(
                args.output,
                {
                    "config": config,
                    "external_plugins": plugins,
                    "status": "valid",
                },
            )
    except ConfigValidationError as exc:
        print(f"validate_config: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
