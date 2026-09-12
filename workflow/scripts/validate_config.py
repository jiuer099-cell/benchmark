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
GRAPH_METADATA_FIELDS = {"profile", "reference_path"}


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
    for mode, contract in tool["supported_modes"].items():
        required = set(contract["required_inputs"])
        optional = set(contract["optional_inputs"])
        if required & optional:
            raise ConfigValidationError(
                f"tool {tool['id']} mode {mode} has overlapping required/optional inputs"
            )
        if mode != "end_to_end_from_reads":
            raise ConfigValidationError(f"tool {tool['id']} declares unsupported mode {mode}")
        capabilities = tool.get("capabilities", {})
        read_class = (
            capabilities.get("read_class", "short")
            if isinstance(capabilities, Mapping)
            else "short"
        )
        evidence = (
            {"short_fastq_r1", "short_fastq_r2"}
            if read_class == "short"
            else {"long_reads_fastq"}
        )
        if not evidence.issubset(required):
            raise ConfigValidationError(
                f"tool {tool['id']} must require {sorted(evidence)} for its "
                f"{read_class}-read channel"
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
    loaded: dict[str, dict[str, Any]] = {}
    for registration in config["external_plugins"]:
        plugin_id = registration["id"]
        if plugin_id in loaded:
            raise ConfigValidationError(f"duplicate external plugin ID {plugin_id}")
        manifest_path = (repo_root / registration["manifest"]).resolve()
        manifest = load_yaml(manifest_path)
        validate_json_schema(manifest, tool_schema_path)
        if manifest["id"] != plugin_id:
            raise ConfigValidationError(
                f"config plugin ID {plugin_id!r} does not match manifest "
                f"ID {manifest['id']!r}"
            )
        _mode_semantics(manifest)
        if manifest["source"] != "external":
            raise ConfigValidationError(f"tool {plugin_id} must be an external adapter")
        if manifest["outputs"].get("candidate_output_contract") != "all_sites":
            raise ConfigValidationError(
                f"tool {plugin_id} must produce canonical all-sites output"
            )
        output_path = Path(str(manifest["outputs"].get("vcf", "")))
        reserved_names = {
            "all-sites.vcf", "all-sites.vcf.gz",
            "evaluation-query.vcf", "evaluation-query.vcf.gz",
            "candidate-status.tsv", "addressability.audit.json",
        }
        if "canonical" in {part.casefold() for part in output_path.parts} or (
            output_path.name.casefold() in reserved_names
        ):
            raise ConfigValidationError(
                f"tool {plugin_id} output uses a benchmark-core reserved path"
            )
        information = manifest.get("information_contract")
        if not isinstance(information, Mapping):
            raise ConfigValidationError(
                f"tool {plugin_id} is missing its allowed-information contract"
            )
        loaded[plugin_id] = manifest
    return loaded


def _validate_graph_profile(config: Mapping[str, Any]) -> None:
    pangenome = config["pangenome"]
    graph = pangenome["graph_assets"]
    profile = graph["profile"]
    enabled = bool(pangenome["build_graph_assets"])
    if enabled and profile == "none":
        raise ConfigValidationError(
            "build_graph_assets=true requires a non-none graph profile"
        )
    if not enabled and profile != "none":
        raise ConfigValidationError(
            "build_graph_assets=false requires graph_assets.profile=none"
        )
    configured_assets = {
        name: value
        for name, value in graph.items()
        if name not in GRAPH_METADATA_FIELDS and isinstance(value, str) and value
    }
    if enabled and "manifest" not in configured_assets:
        raise ConfigValidationError(
            f"graph profile {profile} requires a source manifest"
        )
    if enabled and set(configured_assets) == {"manifest"}:
        raise ConfigValidationError(
            f"graph profile {profile} declares no graph assets"
        )


def _validate_plugin_compatibility(
    config: Mapping[str, Any],
    plugins: Mapping[str, Mapping[str, Any]],
) -> None:
    technology = config["sample"]["technology"]
    benchmark_track = config["benchmark_contract"]["track"]
    track_class = config["benchmark_contract"]["read_class"]
    mode = config["execution"]["official_score_mode"]
    graph_profile = config["pangenome"]["graph_assets"]["profile"]
    native_contract = config["pangenome"].get("native_input_contract")
    scoring_contract = config["pangenome"]["canonical_scoring_contract"]
    if (
        scoring_contract.get("canonical_candidate_count") != 18164
        or scoring_contract.get("tool_specific_denominator") != "forbidden"
        or scoring_contract.get("deterministic_native_projection_required") is not True
    ):
        raise ConfigValidationError(
            "canonical scoring contract must freeze the 18,164-candidate universe "
            "and deterministic native projection"
        )
    for plugin_id, manifest in plugins.items():
        supported_technologies = set(manifest["capabilities"]["technology"])
        if technology not in supported_technologies:
            raise ConfigValidationError(
                f"tool {plugin_id} does not support sample technology {technology}; "
                f"supported={sorted(supported_technologies)}"
            )
        if manifest["capabilities"]["read_class"] != track_class:
            # Registered adapters for the other channel are explicitly
            # ineligible rather than counted as failed/zero-score tools.
            continue
        contract = manifest["supported_modes"].get(mode)
        if not isinstance(contract, Mapping):
            raise ConfigValidationError(
                f"tool {plugin_id} does not support official mode {mode}"
            )
        required = set(contract["required_inputs"])
        for field, contract_name in (
            ("fastq_r1", "short_fastq_r1"),
            ("fastq_r2", "short_fastq_r2"),
            ("fastq", "long_reads_fastq"),
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
        if {"shared_shortread_alignment", "shared_shortread_alignment_index"}.intersection(required):
            alignment = (
                native_contract.get("shared_shortread_alignment")
                if isinstance(native_contract, Mapping)
                else None
            )
            if alignment is not None and not isinstance(alignment, Mapping):
                raise ConfigValidationError(
                    f"tool {plugin_id} shared alignment contract must be a mapping"
                )
            if not {"shared_shortread_alignment", "shared_shortread_alignment_index"}.issubset(required):
                raise ConfigValidationError(
                    f"tool {plugin_id} must require the benchmark shared BAM and BAI"
                )
        native_assets = manifest["native_assets"]
        bundle = config["pangenome"]["frozen_haplotype_source_bundle"]
        if native_assets["uses_population_or_pangenome_information"]:
            if native_assets["frozen_haplotype_source_bundle"] != bundle["id"]:
                raise ConfigValidationError(
                    f"tool {plugin_id} must derive its native assets from "
                    "the configured Frozen Haplotype Source Bundle"
                )
            if bundle["source_cohort_id"] != config["pangenome"]["population_panel"]:
                raise ConfigValidationError("Frozen Haplotype Source Bundle cohort mismatch")
            if bundle["target_family_excluded"] is not True:
                raise ConfigValidationError("Frozen Haplotype Source Bundle requires family exclusion")
        elif native_assets["frozen_haplotype_source_bundle"] is not None:
            raise ConfigValidationError(
                f"tool {plugin_id} declares a frozen bundle without consuming population information"
            )
        accepted = set(manifest["inputs"]["graph_assets"].get("accepted_profiles", []))
        if accepted and graph_profile not in accepted:
            raise ConfigValidationError(
                f"tool {plugin_id} rejects graph profile {graph_profile}; "
                f"accepted={sorted(accepted)}"
            )
        registered_assets = config["external_plugins"][list(plugins).index(plugin_id)].get("adapter_assets", {})
        adapter_assets = {name for name in required if name.startswith("adapter_asset.")}
        missing_assets = sorted(adapter_assets - set(registered_assets))
        if missing_assets:
            raise ConfigValidationError(
                f"tool {plugin_id} is missing registered adapter assets: {missing_assets}"
            )

    if benchmark_track not in {"short_read_fixed_panel_genotyping", "long_read_fixed_panel_genotyping"}:
        raise ConfigValidationError(f"unsupported benchmark track {benchmark_track}")


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
    catalogues: dict[str, dict[str, Any]] = {}
    for catalog_name in (
        "truthsets",
        "stratifications",
        "score_weights",
        "evaluator_profile",
        "metric_dictionary",
        "allowed_information",
        "tuning_policy",
        "panel_provenance",
    ):
        catalog_path = (repo_root / config["catalogs"][catalog_name]).resolve()
        if not catalog_path.is_file():
            raise ConfigValidationError(
                f"{catalog_name} does not exist: {catalog_path}"
            )
        catalogues[catalog_name] = load_yaml(catalog_path)
    if (
        catalogues["allowed_information"].get("contract")
        != "pgbench_allowed_information_policy_v1"
    ):
        raise ConfigValidationError("unsupported allowed-information policy")
    if catalogues["tuning_policy"].get("contract") != "pgbench_tuning_policy_v1":
        raise ConfigValidationError("unsupported tuning policy")
    panel_status = catalogues["panel_provenance"].get("status")
    if panel_status not in {"verified", "unverified_information"}:
        raise ConfigValidationError("invalid panel provenance status")
    if (
        config["benchmark_contract"]["score_status"] == "valid"
        and panel_status != "verified"
    ):
        raise ConfigValidationError(
            "score_status=valid requires verified panel source provenance"
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
