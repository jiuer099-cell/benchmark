#!/usr/bin/env python3
"""Execute a validated external tool under the isolated PGBench contract."""

from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import signal
import subprocess
import sys
import tempfile
import threading
import time
import uuid
from collections.abc import Iterator, Mapping, Sequence
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path, PurePosixPath
from itertools import zip_longest
from typing import Any, TextIO, cast

import jsonschema
import yaml  # type: ignore[import-untyped]

try:
    from pgbench_provenance import (
        atomic_write_json,
        fingerprint_path,
        sha256_file,
    )
    from validate_tool_output import (
        ToolOutputValidation,
        ToolOutputValidationError,
        validate_tool_output,
    )
except ModuleNotFoundError:  # pragma: no cover - package-style invocation
    from .pgbench_provenance import (
        atomic_write_json,
        fingerprint_path,
        sha256_file,
    )
    from .validate_tool_output import (
        ToolOutputValidation,
        ToolOutputValidationError,
        validate_tool_output,
    )


SUPPORTED_MODES = {"end_to_end_from_reads"}
MODE_INPUTS = {
    "short_fastq_r1",
    "short_fastq_r2",
    "reference",
    "reference_index",
    "pangenome_manifest",
    "pangenome_panel",
    "pangenie_private_phased_panel",
    "pangenie_private_biallelic_panel",
    "pangenie_biallelic_converter",
    "canonical_allele_projection",
    "candidate_panel",
    "allele_fasta",
    "graph_assets",
    "tool_index",
}
TRANSPORT_INPUTS = MODE_INPUTS
TRANSPORT_TO_CONTRACT = {name: name for name in MODE_INPUTS}


def _fastq_records(path: Path) -> Iterator[tuple[str, int]]:
    """Yield normalized read names and base counts after strict FASTQ checks."""

    opener = gzip.open if path.name.casefold().endswith(".gz") else open
    try:
        with opener(path, "rt", encoding="utf-8", errors="strict", newline="") as handle:
            record_number = 0
            while True:
                header = handle.readline()
                if header == "":
                    return
                sequence = handle.readline()
                plus = handle.readline()
                quality = handle.readline()
                record_number += 1
                if not sequence or not plus or not quality:
                    raise ToolContractError(
                        f"truncated FASTQ record {record_number}: {path}"
                    )
                header = header.rstrip("\r\n")
                sequence = sequence.rstrip("\r\n")
                plus = plus.rstrip("\r\n")
                quality = quality.rstrip("\r\n")
                if not header.startswith("@"):
                    raise ToolContractError(
                        f"FASTQ record {record_number} lacks @ header: {path}"
                    )
                if not plus.startswith("+"):
                    raise ToolContractError(
                        f"FASTQ record {record_number} lacks + separator: {path}"
                    )
                if not sequence or len(sequence) != len(quality):
                    raise ToolContractError(
                        f"FASTQ sequence/quality length mismatch at record "
                        f"{record_number}: {path}"
                    )
                name = header[1:].split(maxsplit=1)[0]
                if name.endswith("/1") or name.endswith("/2"):
                    name = name[:-2]
                if not name:
                    raise ToolContractError(
                        f"FASTQ record {record_number} has an empty read name: {path}"
                    )
                yield name, len(sequence)
    except (OSError, EOFError, UnicodeError) as exc:
        raise ToolContractError(f"cannot fully decode FASTQ {path}: {exc}") from exc


def validate_read_evidence(supplied_inputs: Mapping[str, Path]) -> dict[str, Any]:
    """Validate complete FASTQ streams and mate correspondence before execution."""

    r1 = supplied_inputs.get("short_fastq_r1")
    r2 = supplied_inputs.get("short_fastq_r2")
    if (r1 is None) != (r2 is None):
        raise ToolContractError(
            "short-read evidence must supply both short_fastq_r1 and short_fastq_r2"
        )
    result: dict[str, Any] = {"contract": "pgbench_fastq_validation_v1"}
    if r1 is not None and r2 is not None:
        read_count = 0
        read_bases = 0
        for pair_number, pair in enumerate(
            zip_longest(_fastq_records(r1), _fastq_records(r2)), start=1
        ):
            left, right = pair
            if left is None or right is None:
                raise ToolContractError(
                    "paired FASTQ files contain different record counts"
                )
            if left[0] != right[0]:
                raise ToolContractError(
                    f"paired FASTQ read-name mismatch at pair {pair_number}: "
                    f"{left[0]!r} != {right[0]!r}"
                )
            read_count += 1
            read_bases += left[1] + right[1]
        if read_count == 0:
            raise ToolContractError("paired FASTQ evidence is empty")
        result["paired_fastq"] = {
            "read_pairs": read_count,
            "read_count": read_count * 2,
            "read_bases": read_bases,
            "mate_names_match": True,
        }
    return result
INPUT_ENVIRONMENT = {
    "short_fastq_r1": "PGBENCH_INPUT_FASTQ_R1",
    "short_fastq_r2": "PGBENCH_INPUT_FASTQ_R2",
    "reference": "PGBENCH_REFERENCE_FASTA",
    "reference_index": "PGBENCH_REFERENCE_INDEX",
    "candidate_panel": "PGBENCH_CANDIDATE_VCF",
    "pangenome_manifest": "PGBENCH_PANGENOME_MANIFEST",
    "pangenome_panel": "PGBENCH_PANEL_VCF",
    "pangenie_private_phased_panel": "PGBENCH_PANGENIE_PRIVATE_PHASED_PANEL",
    "pangenie_private_biallelic_panel": "PGBENCH_PANGENIE_PRIVATE_BIALLELIC_PANEL",
    "pangenie_biallelic_converter": "PGBENCH_PANGENIE_BIALLELIC_CONVERTER",
    "canonical_allele_projection": "PGBENCH_CANONICAL_ALLELE_PROJECTION",
    "allele_fasta": "PGBENCH_ALLELES_FASTA",
    "graph_assets": "PGBENCH_GRAPH_DIR",
    "tool_index": "PGBENCH_INDEX_DIR",
}
BASE_PGBENCH_ENVIRONMENT = {
    "PGBENCH_RUN_ID",
    "PGBENCH_SAMPLE_ID",
    "PGBENCH_RESOLVED_INPUTS",
    "PGBENCH_OUTPUT_DIR",
    "PGBENCH_OUTPUT_VCF",
    "PGBENCH_THREADS",
    "PGBENCH_MEMORY_MB",
    "PGBENCH_CACHE_POLICY",
    "PGBENCH_ALLELE_NAMESPACE",
    "PGBENCH_RANDOM_SEED",
    "PGBENCH_GRAPH_REFERENCE_PATH",
}
OUTPUT_INDEX_SUFFIXES = (".tbi", ".csi", ".idx")
GRAPH_ASSET_LOCK_NAME = "graph-assets.lock.yaml"


class ToolContractError(ValueError):
    """Raised when a tool manifest or invocation violates the interface."""


class ToolExecutionError(RuntimeError):
    """Raised after an external tool attempt has been recorded as failed."""


class ToolTimeoutError(ToolExecutionError):
    """Raised after the complete external-tool process group times out."""

    def __init__(self, message: str, *, returncode: int | None = None) -> None:
        self.returncode = returncode
        super().__init__(message)


class ToolInterruptedError(ToolExecutionError):
    """Raised when SIGTERM requests interruption of the current tool attempt."""

    def __init__(self, message: str, *, returncode: int | None = None) -> None:
        self.returncode = returncode
        super().__init__(message)


@dataclass(frozen=True)
class ResolvedInput:
    """One mode-authorized input and its immutable content identity."""

    name: str
    contract_name: str
    path: str
    sha256: str
    size_bytes: int
    mtime_ns: int
    path_type: str
    read_only: bool
    environment_variable: str | None

    def to_dict(self) -> dict[str, str | int | bool | None]:
        return asdict(self)


@dataclass(frozen=True)
class ToolExecutionResult:
    """Successful external-tool execution metadata."""

    attempt_id: str
    command: tuple[str, ...]
    output: ToolOutputValidation
    resolved_inputs_path: str
    attempt_record_path: str
    attempt_archive_path: str
    log_path: str
    log_archive_path: str


@dataclass(frozen=True)
class PreviousOutputArchive:
    """Recoverable location of one provenance-backed previous output."""

    path: str
    source_attempt_id: str
    files: tuple[str, ...]
    output_sha256: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _utc_now() -> str:
    return datetime.now(UTC).isoformat().replace("+00:00", "Z")


def _load_yaml_mapping(path: Path, label: str) -> dict[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as handle:
            loaded = yaml.safe_load(handle)
    except (OSError, yaml.YAMLError) as exc:
        raise ToolContractError(f"cannot load {label} {path}: {exc}") from exc
    if not isinstance(loaded, dict):
        raise ToolContractError(f"{label} must contain a YAML mapping: {path}")
    return cast(dict[str, Any], loaded)


def _json_path(error: jsonschema.ValidationError) -> str:
    return ".".join(str(part) for part in error.absolute_path) or "<root>"


def _semantic_manifest_errors(manifest: Mapping[str, Any]) -> list[str]:
    errors: list[str] = []
    if not isinstance(manifest.get("version"), str) or not manifest.get("version"):
        errors.append("tools must declare a non-empty version")

    execution = manifest.get("execution")
    if isinstance(execution, Mapping) and "module" in execution:
        errors.append(
            "generic external execution accepts runner-only plugins; "
            "custom modules require a separate trusted-reviewed path"
        )

    outputs = manifest.get("outputs")
    if isinstance(outputs, Mapping):
        output_value = outputs.get("vcf")
        if isinstance(output_value, str):
            pure_output = PurePosixPath(output_value)
            if pure_output.is_absolute() or ".." in pure_output.parts:
                errors.append("outputs.vcf must be a relative path without '..'")
            if not (output_value.endswith(".vcf") or output_value.endswith(".vcf.gz")):
                errors.append("outputs.vcf must end in .vcf or .vcf.gz")

    supported_modes = manifest.get("supported_modes")
    inputs = manifest.get("inputs")
    paradigm = manifest.get("paradigm")
    if not isinstance(supported_modes, Mapping) or not isinstance(inputs, Mapping):
        return errors

    globally_required: set[str] = set()
    for key in ("reference", "pangenome_manifest", "pangenome_panel"):
        if inputs.get(key) == "required":
            globally_required.add(key)
    if inputs.get("candidates") == "required":
        globally_required.add("candidate_panel")
    graph = inputs.get("graph_assets")
    if isinstance(graph, Mapping) and graph.get("required") is True:
        globally_required.add("graph_assets")

    for mode, raw_contract in supported_modes.items():
        if mode not in SUPPORTED_MODES or not isinstance(raw_contract, Mapping):
            continue
        required = set(cast(Sequence[str], raw_contract.get("required_inputs", [])))
        optional = set(cast(Sequence[str], raw_contract.get("optional_inputs", [])))
        overlap = sorted(required & optional)
        if overlap:
            errors.append(
                f"supported_modes.{mode} has required/optional overlap: "
                + ", ".join(overlap)
            )
        for asset in sorted(globally_required):
            if asset not in required:
                errors.append(
                    f"supported_modes.{mode}.required_inputs must include {asset}"
                )

        if mode != "end_to_end_from_reads":
            errors.append(f"unsupported execution mode: {mode}")
        if not {"short_fastq_r1", "short_fastq_r2"}.issubset(required):
            errors.append("end_to_end_from_reads must require paired FASTQ inputs")

        if paradigm == "genotyping_only" and "candidate_panel" not in required:
            errors.append(
                f"supported_modes.{mode}.required_inputs must include "
                "candidate_panel for genotyping_only"
            )
    return errors


def load_tool_manifest(
    tool_manifest_path: Path,
    schema_path: Path,
) -> dict[str, Any]:
    """Load and validate a tool manifest against schema and mode semantics."""

    manifest = _load_yaml_mapping(tool_manifest_path, "tool manifest")
    schema = _load_yaml_mapping(schema_path, "tool schema")
    try:
        jsonschema.Draft202012Validator.check_schema(schema)
    except jsonschema.SchemaError as exc:
        raise ToolContractError(f"invalid tool schema {schema_path}: {exc}") from exc
    validation_errors = sorted(
        jsonschema.Draft202012Validator(schema).iter_errors(manifest),
        key=lambda error: list(error.absolute_path),
    )
    errors = [f"{_json_path(error)}: {error.message}" for error in validation_errors]
    errors.extend(_semantic_manifest_errors(manifest))
    errors.extend(_execution_path_errors(tool_manifest_path, manifest))
    if errors:
        raise ToolContractError("invalid tool.yaml: " + "; ".join(errors))
    return manifest


def _execution_path_errors(
    tool_manifest_path: Path,
    manifest: Mapping[str, Any],
) -> list[str]:
    errors: list[str] = []
    execution = manifest.get("execution")
    if not isinstance(execution, Mapping):
        return errors
    try:
        plugin_root = tool_manifest_path.parent.resolve(strict=True)
    except OSError as exc:
        return [f"plugin root cannot be resolved: {exc}"]
    for key in ("runner", "environment", "rule_registry"):
        value = execution.get(key)
        if value is None:
            continue
        if not isinstance(value, str) or not value:
            continue
        raw = Path(value).expanduser()
        candidate = raw if raw.is_absolute() else plugin_root / raw
        try:
            resolved = candidate.resolve(strict=True)
        except OSError:
            errors.append(f"execution.{key} does not exist: {value}")
            continue
        try:
            resolved.relative_to(plugin_root)
        except ValueError:
            errors.append(f"execution.{key} must stay within plugin root: {value}")
            continue
        if not resolved.is_file():
            errors.append(f"execution.{key} is not a regular file: {value}")
    return errors


def _parse_inputs(values: Sequence[str]) -> dict[str, Path]:
    parsed: dict[str, Path] = {}
    for value in values:
        name, separator, path_text = value.partition("=")
        if not separator or not name or not path_text:
            raise ToolContractError(
                f"--input must use NAME=PATH syntax, found {value!r}"
            )
        if name not in TRANSPORT_INPUTS:
            raise ToolContractError(f"unknown or non-exposable input name: {name}")
        if name in parsed:
            raise ToolContractError(f"duplicate input name: {name}")
        parsed[name] = Path(path_text)
    return parsed


def _contract_input_is_present(name: str, supplied: Mapping[str, Path]) -> bool:
    return name in supplied


def validate_mode_inputs(
    manifest: Mapping[str, Any],
    mode: str,
    supplied_inputs: Mapping[str, Path],
) -> None:
    """Enforce required and declared short-read inputs for one mode."""

    supported_modes = cast(Mapping[str, Mapping[str, Any]], manifest["supported_modes"])
    if mode not in supported_modes:
        raise ToolContractError(
            f"tool {manifest['id']} does not support execution mode {mode}"
        )
    contract = supported_modes[mode]
    required = set(cast(Sequence[str], contract["required_inputs"]))
    optional = set(cast(Sequence[str], contract.get("optional_inputs", [])))

    missing = sorted(
        name
        for name in required
        if not _contract_input_is_present(name, supplied_inputs)
    )
    if missing:
        raise ToolContractError(
            f"missing required input(s) for {mode}: {', '.join(missing)}"
        )

    for transport_name, path in supplied_inputs.items():
        contract_name = TRANSPORT_TO_CONTRACT[transport_name]
        if not path.exists():
            raise ToolContractError(
                f"declared input does not exist: {transport_name}={path}"
            )
        if contract_name not in required | optional:
            raise ToolContractError(
                f"input {contract_name} is not whitelisted by supported_modes.{mode}"
            )

def resolve_inputs(
    *,
    manifest: Mapping[str, Any],
    tool_manifest_path: Path,
    mode: str,
    supplied_inputs: Mapping[str, Path],
    resolved_inputs_path: Path,
    run_id: str,
) -> tuple[ResolvedInput, ...]:
    """Resolve and atomically freeze authorized inputs with content hashes."""

    mode_contract = cast(
        Mapping[str, Any],
        cast(Mapping[str, Any], manifest["supported_modes"])[mode],
    )
    resolved: list[ResolvedInput] = []
    for name in sorted(supplied_inputs):
        candidate = supplied_inputs[name].expanduser().resolve(strict=True)
        fingerprint = fingerprint_path(candidate)
        resolved.append(
            ResolvedInput(
                name=name,
                contract_name=TRANSPORT_TO_CONTRACT[name],
                path=str(candidate),
                sha256=fingerprint.sha256,
                size_bytes=fingerprint.size,
                mtime_ns=fingerprint.mtime_ns,
                path_type=fingerprint.path_type,
                read_only=True,
                environment_variable=INPUT_ENVIRONMENT.get(name),
            )
        )

    payload = {
        "schema_version": 1,
        "run_id": run_id,
        "tool_id": manifest["id"],
        "tool_manifest_path": str(tool_manifest_path.resolve(strict=True)),
        "tool_manifest_sha256": sha256_file(tool_manifest_path.resolve(strict=True)),
        "mode": mode,
        "billable_stages": list(cast(Sequence[str], mode_contract["billable_stages"])),
        "required_inputs": list(cast(Sequence[str], mode_contract["required_inputs"])),
        "optional_inputs": list(
            cast(Sequence[str], mode_contract.get("optional_inputs", []))
        ),
        "inputs": [item.to_dict() for item in resolved],
        "read_validation": validate_read_evidence(supplied_inputs),
    }
    atomic_write_json(resolved_inputs_path, payload)
    return tuple(resolved)


def _ensure_safe_directory_chain(root: Path, relative: Path) -> Path:
    """Create a directory chain below ``root`` without traversing symlinks."""

    current = root
    for part in relative.parts:
        if part in {"", "."}:
            continue
        current = current / part
        if current.is_symlink():
            raise ToolContractError(
                f"assigned output path traverses a symlink: {current}"
            )
        if current.exists():
            if not current.is_dir():
                raise ToolContractError(
                    f"assigned output parent is not a directory: {current}"
                )
            continue
        try:
            # Microsoft Store Python can translate POSIX-style restrictive
            # modes into unusable Windows ACLs.  Keep owner-only permissions
            # on POSIX and inherit the normal parent ACL on Windows.
            if os.name == "nt":
                current.mkdir()
            else:
                current.mkdir(mode=0o700)
        except FileExistsError as exc:
            if current.is_symlink() or not current.is_dir():
                raise ToolContractError(
                    f"assigned output parent is not a safe directory: {current}"
                ) from exc
    return current


def _safe_output_path(output_dir: Path, relative_output: str) -> Path:
    relative = PurePosixPath(relative_output)
    if relative.is_absolute() or ".." in relative.parts:
        raise ToolContractError(
            "tool outputs.vcf must stay below the assigned output directory"
        )
    requested_root = output_dir.expanduser()
    if requested_root.is_symlink():
        raise ToolContractError("assigned output directory must not be a symlink")
    requested_root.mkdir(parents=True, exist_ok=True)
    root = requested_root.resolve(strict=True)
    relative_path = Path(*relative.parts)
    _ensure_safe_directory_chain(root, relative_path.parent)
    output = root / relative_path
    try:
        output.relative_to(root)
    except ValueError as exc:
        raise ToolContractError(
            "tool outputs.vcf escapes the assigned output directory"
        ) from exc
    if output.is_symlink():
        raise ToolContractError("assigned output VCF must not be a symlink")
    return output


def _resolve_runner(tool_manifest_path: Path, runner_value: str) -> Path:
    plugin_root = tool_manifest_path.parent.resolve(strict=True)
    raw = Path(runner_value).expanduser()
    candidate = raw if raw.is_absolute() else plugin_root / raw
    try:
        resolved = candidate.resolve(strict=True)
    except OSError as exc:
        raise ToolContractError(
            f"external runner does not exist: {runner_value}"
        ) from exc
    try:
        resolved.relative_to(plugin_root)
    except ValueError as exc:
        raise ToolContractError(
            "external runner must stay within the registered plugin root"
        ) from exc
    if not resolved.is_file():
        raise ToolContractError(f"external runner is not a file: {runner_value}")
    return resolved


def _runner_command(
    runner: Path,
    *,
    sandbox_backend: str,
) -> tuple[str, ...]:
    if runner.suffix == ".py":
        interpreter = "python" if sandbox_backend == "apptainer" else sys.executable
        return (interpreter, str(runner))
    if not os.access(runner, os.X_OK):
        raise ToolContractError(f"external runner is not executable: {runner}")
    return (str(runner),)


def _locked_graph_reference_path(
    resolved_inputs: Sequence[ResolvedInput],
) -> str | None:
    """Read the reference selector from the resolved, content-locked graph bundle."""

    graph_inputs = [item for item in resolved_inputs if item.name == "graph_assets"]
    if not graph_inputs:
        return None
    if len(graph_inputs) != 1:  # defensive assertion against interface drift
        raise ToolContractError("exactly one graph_assets input is supported")

    graph_input = graph_inputs[0]
    graph_dir = Path(graph_input.path)
    if graph_input.path_type != "directory" or not graph_dir.is_dir():
        raise ToolContractError("graph_assets must resolve to a directory")

    lock_path = graph_dir / GRAPH_ASSET_LOCK_NAME
    if lock_path.is_symlink():
        raise ToolContractError(
            f"graph asset lock must not be a symlink: {lock_path}"
        )
    if not lock_path.is_file():
        raise ToolContractError(
            "graph_assets is missing its frozen reference selector: "
            f"{lock_path}"
        )
    lock = _load_yaml_mapping(lock_path, "graph asset lock")
    reference_path = lock.get("reference_path")
    if (
        not isinstance(reference_path, str)
        or not reference_path.strip()
        or "\n" in reference_path
        or "\r" in reference_path
    ):
        raise ToolContractError(
            "graph asset lock reference_path must be a non-empty single-line string"
        )
    return reference_path


def build_tool_environment(
    *,
    base_environment: Mapping[str, str],
    manifest: Mapping[str, Any],
    run_id: str,
    sample_id: str,
    resolved_inputs: Sequence[ResolvedInput],
    resolved_inputs_path: Path,
    output_dir: Path,
    output_vcf: Path,
    threads: int,
    memory_mb: int,
    alignment_kind: str | None = None,
    cache_policy: str = "isolated_empty_tool_cache",
    attempt_work_dir: Path,
    random_seed: int = 0,
) -> dict[str, str]:
    """Build a minimal deterministic environment with PGBENCH values scrubbed."""

    del alignment_kind  # mapping is tool-owned in the only supported execution mode

    inherited_allowlist = {
        "PATH",
        "CONDA_PREFIX",
        "VIRTUAL_ENV",
        "LD_LIBRARY_PATH",
        "DYLD_LIBRARY_PATH",
        "SYSTEMROOT",
    }
    environment = {
        key: value
        for key, value in base_environment.items()
        if key in inherited_allowlist
    }
    home_dir = attempt_work_dir / "home"
    temporary_dir = attempt_work_dir / "tmp"
    cache_dir = attempt_work_dir / "cache"
    for directory in (home_dir, temporary_dir, cache_dir):
        directory.mkdir(parents=True, exist_ok=True)
    environment.update(
        {
            "HOME": str(home_dir.resolve()),
            "PYTHONDONTWRITEBYTECODE": "1",
            "TMPDIR": str(temporary_dir.resolve()),
            "XDG_CACHE_HOME": str(cache_dir.resolve()),
            "TZ": "UTC",
            "LC_ALL": "C",
            "LANG": "C",
            "PYTHONHASHSEED": "0",
            "PYTHONDONTWRITEBYTECODE": "1",
        }
    )
    namespaces = cast(
        Sequence[str],
        cast(Mapping[str, Any], manifest["capabilities"])["allele_namespace"],
    )
    pgbench_values = {
        "PGBENCH_RUN_ID": run_id,
        "PGBENCH_SAMPLE_ID": sample_id,
        "PGBENCH_RESOLVED_INPUTS": str(resolved_inputs_path.resolve()),
        "PGBENCH_OUTPUT_DIR": str(output_dir.resolve()),
        "PGBENCH_OUTPUT_VCF": str(output_vcf.resolve()),
        "PGBENCH_THREADS": str(threads),
        "PGBENCH_MEMORY_MB": str(memory_mb),
        "PGBENCH_CACHE_POLICY": cache_policy,
        "PGBENCH_ALLELE_NAMESPACE": namespaces[0],
        "PGBENCH_RANDOM_SEED": str(random_seed),
    }
    for item in resolved_inputs:
        if item.environment_variable is not None:
            pgbench_values[item.environment_variable] = item.path
    graph_reference_path = _locked_graph_reference_path(resolved_inputs)
    if graph_reference_path is not None:
        pgbench_values["PGBENCH_GRAPH_REFERENCE_PATH"] = graph_reference_path

    allowed = BASE_PGBENCH_ENVIRONMENT | {
        item.environment_variable
        for item in resolved_inputs
        if item.environment_variable is not None
    }
    unexpected = sorted(set(pgbench_values) - allowed)
    if unexpected:  # defensive assertion against interface drift
        raise ToolContractError(
            "runner attempted to expose unapproved PGBENCH variables: "
            + ", ".join(unexpected)
        )
    environment.update(pgbench_values)
    return environment


def _prepare_log(log_path: Path) -> tuple[Path, TextIO]:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=log_path.parent,
        prefix=f".{log_path.name}.",
        suffix=".attempt.tmp",
        text=True,
    )
    return Path(temporary_name), os.fdopen(
        descriptor, "w", encoding="utf-8", newline=""
    )


def _fsync_directory(path: Path) -> None:
    try:
        descriptor = os.open(path, os.O_RDONLY)
    except OSError:
        return
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _recorded_relative_path(
    value: Any,
    recorded_root: Any,
    label: str,
) -> Path:
    """Return a recorded artifact path relative to its recorded output root.

    Attempt records are portable provenance: an output tree may be copied from
    macOS or Linux to Windows before a rerun.  Validate the lexical relationship
    recorded at creation time, then compare that relative location with the
    current, already-resolved output tree.  File identity is still protected by
    the recorded SHA-256 and size checks below.
    """

    if not isinstance(value, str) or not value:
        raise ToolContractError(
            f"previous attempt does not contain a valid {label} path"
        )
    if not isinstance(recorded_root, str) or not recorded_root:
        raise ToolContractError(
            "previous attempt does not contain a valid output_dir path"
        )
    recorded_path = PurePosixPath(value.replace("\\", "/"))
    root_path = PurePosixPath(recorded_root.replace("\\", "/"))
    try:
        relative = recorded_path.relative_to(root_path)
    except ValueError as exc:
        raise ToolContractError(
            f"previous attempt {label} escapes its recorded output_dir"
        ) from exc
    if relative.is_absolute() or ".." in relative.parts:
        raise ToolContractError(
            f"previous attempt {label} escapes its recorded output_dir"
        )
    return Path(*relative.parts)


def _verified_previous_output(
    *,
    attempt_record_path: Path,
    final_output_vcf: Path,
    output_root: Path,
    tool_id: str,
    sample_id: str,
    mode: str,
) -> tuple[dict[str, Any], str]:
    """Validate that an existing VCF is the last successful owned output."""

    if final_output_vcf.is_symlink():
        raise ToolContractError("assigned output VCF must not be a symlink")
    if not final_output_vcf.exists():
        raise ToolContractError("assigned output VCF disappeared before archival")
    if not final_output_vcf.is_file():
        raise ToolContractError("assigned output VCF is not a regular file")
    try:
        final_output_vcf.resolve(strict=True).relative_to(output_root)
    except (OSError, ValueError) as exc:
        raise ToolContractError(
            "assigned output VCF escapes the assigned output directory"
        ) from exc

    if attempt_record_path.is_symlink():
        raise ToolContractError("previous attempt record must not be a symlink")
    if not attempt_record_path.is_file():
        raise ToolContractError(
            "assigned output VCF already exists but has no previous attempt record"
        )
    try:
        with attempt_record_path.open("r", encoding="utf-8") as handle:
            previous = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise ToolContractError(
            "assigned output VCF already exists but its previous attempt record "
            "cannot be loaded"
        ) from exc
    if not isinstance(previous, dict):
        raise ToolContractError("previous attempt record must be a JSON object")

    expected_identity = {
        "status": "success",
        "tool_id": tool_id,
        "sample_id": sample_id,
        "mode": mode,
    }
    mismatched = [
        field
        for field, expected in expected_identity.items()
        if previous.get(field) != expected
    ]
    if mismatched:
        raise ToolContractError(
            "assigned output VCF already exists but previous attempt identity "
            "does not match: " + ", ".join(mismatched)
        )

    expected_path = final_output_vcf.resolve(strict=True)
    expected_relative = expected_path.relative_to(output_root)
    recorded_output_dir = previous.get("output_dir")
    if _recorded_relative_path(
        previous.get("expected_output_vcf"),
        recorded_output_dir,
        "expected_output_vcf",
    ) != expected_relative:
        raise ToolContractError(
            "assigned output VCF already exists but previous expected path differs"
        )

    output = previous.get("output")
    if not isinstance(output, Mapping):
        raise ToolContractError(
            "assigned output VCF already exists but previous output metadata is absent"
        )
    if _recorded_relative_path(
        output.get("path"),
        recorded_output_dir,
        "output.path",
    ) != expected_relative:
        raise ToolContractError(
            "assigned output VCF already exists but previous output path differs"
        )
    current = fingerprint_path(final_output_vcf)
    metadata_matches = (
        current.path_type == "file"
        and output.get("sha256") == current.sha256
        and output.get("size_bytes") == current.size
    )
    if not metadata_matches:
        raise ToolContractError(
            "assigned output VCF already exists but no longer matches the "
            "previous successful attempt"
        )
    source_attempt_id = previous.get("attempt_id")
    if not isinstance(source_attempt_id, str) or not source_attempt_id:
        raise ToolContractError("previous attempt record has no attempt_id")
    return cast(dict[str, Any], previous), source_attempt_id


def _archive_previous_output(
    *,
    attempt_record_path: Path,
    final_output_vcf: Path,
    output_root: Path,
    attempt_id: str,
    tool_id: str,
    sample_id: str,
    mode: str,
) -> PreviousOutputArchive | None:
    """Atomically publish a recoverable archive of a trusted prior VCF."""

    if not final_output_vcf.exists() and not final_output_vcf.is_symlink():
        return None
    previous, source_attempt_id = _verified_previous_output(
        attempt_record_path=attempt_record_path,
        final_output_vcf=final_output_vcf,
        output_root=output_root,
        tool_id=tool_id,
        sample_id=sample_id,
        mode=mode,
    )

    sources = [final_output_vcf]
    sources.extend(
        Path(f"{final_output_vcf}{suffix}") for suffix in OUTPUT_INDEX_SUFFIXES
    )
    existing_sources: list[Path] = []
    for source in sources:
        if not source.exists() and not source.is_symlink():
            continue
        if source.is_symlink():
            raise ToolContractError(
                f"previous output artifact must not be a symlink: {source}"
            )
        if not source.is_file():
            raise ToolContractError(
                f"previous output artifact is not a regular file: {source}"
            )
        try:
            source.resolve(strict=True).relative_to(output_root)
        except (OSError, ValueError) as exc:
            raise ToolContractError(
                f"previous output artifact escapes output_dir: {source}"
            ) from exc
        existing_sources.append(source)

    history_root = _ensure_safe_directory_chain(
        output_root, Path("meta") / "previous_outputs"
    )
    destination = history_root / attempt_id
    if destination.exists() or destination.is_symlink():
        raise ToolContractError(
            f"previous-output archive already exists; refusing overwrite: {destination}"
        )
    staging = history_root / f".{attempt_id}.{uuid.uuid4().hex}.tmp"
    if staging.exists() or staging.is_symlink():
        raise ToolContractError(
            f"previous-output staging path already exists: {staging}"
        )
    if os.name == "nt":
        staging.mkdir()
    else:
        staging.mkdir(mode=0o700)

    moved: list[tuple[Path, Path]] = []
    relative_files: list[str] = []
    try:
        for source in existing_sources:
            relative = source.relative_to(output_root)
            staged = staging / relative
            staged.parent.mkdir(parents=True, exist_ok=True)
            os.replace(source, staged)
            moved.append((source, staged))
            relative_files.append(relative.as_posix())
        if destination.exists() or destination.is_symlink():
            raise ToolContractError(
                "previous-output archive appeared during archival; refusing overwrite"
            )
        os.replace(staging, destination)
    except BaseException as exc:
        rollback_errors: list[str] = []
        for source, staged in reversed(moved):
            if not staged.exists() and not staged.is_symlink():
                continue
            if source.exists() or source.is_symlink():
                rollback_errors.append(str(source))
                continue
            try:
                os.replace(staged, source)
            except OSError:
                rollback_errors.append(str(source))
        if staging.exists() and not staging.is_symlink():
            shutil.rmtree(staging)
        if rollback_errors:
            raise ToolExecutionError(
                "failed to roll back previous-output archival for: "
                + ", ".join(rollback_errors)
            ) from exc
        if isinstance(exc, ToolContractError):
            raise
        raise ToolContractError(
            f"cannot archive previous output safely: {exc}"
        ) from exc

    _fsync_directory(history_root)
    previous_output = cast(Mapping[str, Any], previous["output"])
    return PreviousOutputArchive(
        path=str(destination),
        source_attempt_id=source_attempt_id,
        files=tuple(relative_files),
        output_sha256=cast(str, previous_output["sha256"]),
    )


def _commit_log(
    temporary_path: Path,
    handle: TextIO,
    log_path: Path,
    archive_path: Path,
) -> None:
    handle.flush()
    os.fsync(handle.fileno())
    handle.close()
    os.replace(temporary_path, archive_path)
    link_path = log_path.with_name(f".{log_path.name}.{uuid.uuid4().hex}.tmp")
    try:
        os.link(archive_path, link_path)
        os.replace(link_path, log_path)
        _fsync_directory(log_path.parent)
    except BaseException:
        link_path.unlink(missing_ok=True)
        raise


def _archive_path(path: Path, attempt_id: str) -> Path:
    return path.with_name(f"{path.stem}.{attempt_id}{path.suffix}")


def _write_attempt_records(
    *,
    current_path: Path,
    archive_path: Path,
    payload: Mapping[str, Any],
) -> None:
    atomic_write_json(archive_path, payload)
    atomic_write_json(current_path, payload)


def _mount_parent_arguments(paths: Sequence[Path]) -> list[str]:
    parents: set[Path] = set()
    for path in paths:
        for parent in path.parents:
            if parent == Path("/"):
                continue
            parents.add(parent)
    arguments: list[str] = []
    for parent in sorted(parents, key=lambda value: len(value.parts)):
        arguments.extend(("--dir", str(parent)))
    return arguments


def _sandbox_command(
    *,
    tool_command: Sequence[str],
    backend: str,
    manifest: Mapping[str, Any],
    plugin_root: Path,
    resolved_inputs: Sequence[ResolvedInput],
    attempt_work_dir: Path,
    additional_read_only_paths: Sequence[Path],
) -> tuple[str, ...]:
    if backend == "none":
        return tuple(tool_command)
    executable = shutil.which(backend)
    if executable is None:
        raise ToolContractError(
            f"sandbox backend {backend!r} is required but not installed"
        )

    input_paths = [
        *[Path(item.path) for item in resolved_inputs],
        *additional_read_only_paths,
    ]
    if backend == "bwrap":
        runtime_roots = [Path(sys.prefix)]
        for variable in ("CONDA_PREFIX", "VIRTUAL_ENV"):
            value = os.environ.get(variable)
            if value:
                runtime_roots.append(Path(value))
        system_roots: list[Path] = []
        seen_roots: set[Path] = set()
        for path in (
            Path("/usr"),
            Path("/bin"),
            Path("/lib"),
            Path("/lib64"),
            Path("/etc/alternatives"),
            *runtime_roots,
        ):
            resolved = path.resolve(strict=False)
            if path.exists() and resolved not in seen_roots:
                seen_roots.add(resolved)
                system_roots.append(path)
        mount_targets = [
            *system_roots,
            plugin_root,
            *input_paths,
            attempt_work_dir,
        ]
        command: list[str] = [
            executable,
            "--die-with-parent",
            "--new-session",
            "--unshare-all",
            "--unshare-net",
            "--proc",
            "/proc",
            "--dev",
            "/dev",
            *_mount_parent_arguments(mount_targets),
        ]
        for root in system_roots:
            command.extend(("--ro-bind", str(root), str(root)))
        command.extend(("--ro-bind", str(plugin_root), str(plugin_root)))
        for path in input_paths:
            command.extend(("--ro-bind", str(path), str(path)))
        command.extend(
            (
                "--bind",
                str(attempt_work_dir),
                str(attempt_work_dir),
                "--chdir",
                str(attempt_work_dir),
                "--",
                *tool_command,
            )
        )
        return tuple(command)

    if backend == "apptainer":
        execution = cast(Mapping[str, Any], manifest["execution"])
        container = execution.get("container")
        if not isinstance(container, str) or not container:
            raise ToolContractError("apptainer sandbox requires execution.container")
        command = [
            executable,
            "exec",
            "--containall",
            "--cleanenv",
            "--no-home",
            "--net",
            "--network",
            "none",
            "--no-mount",
            "home",
            "--pwd",
            str(attempt_work_dir),
            "--bind",
            f"{plugin_root}:{plugin_root}:ro",
        ]
        for path in input_paths:
            command.extend(("--bind", f"{path}:{path}:ro"))
        command.extend(
            (
                "--bind",
                f"{attempt_work_dir}:{attempt_work_dir}:rw",
                container,
                *tool_command,
            )
        )
        return tuple(command)

    raise ToolContractError(f"unsupported sandbox backend: {backend}")


def _terminate_process_group(process: subprocess.Popen[str]) -> int | None:
    if os.name == "nt":
        process.poll()
        if process.returncode is not None:
            return process.returncode
        try:
            subprocess.run(
                ["taskkill", "/PID", str(process.pid), "/T", "/F"],
                stdin=subprocess.DEVNULL,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=False,
                timeout=10,
            )
        except (OSError, subprocess.TimeoutExpired):
            process.kill()
        try:
            return process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            process.kill()
            return process.wait()

    def group_exists() -> bool:
        try:
            os.killpg(process.pid, 0)
        except ProcessLookupError:
            return False
        except PermissionError:
            return True
        return True

    process.poll()
    if not group_exists():
        return process.returncode
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        return process.poll()
    deadline = time.monotonic() + 5
    while group_exists() and time.monotonic() < deadline:
        process.poll()
        time.sleep(0.05)
    if group_exists():
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
    if process.poll() is None:
        process.wait()
    return process.returncode


@contextmanager
def _termination_as_exception() -> Iterator[None]:
    if threading.current_thread() is not threading.main_thread():
        yield
        return
    previous = signal.getsignal(signal.SIGTERM)

    def handle_sigterm(signum: int, frame: Any) -> None:
        del signum, frame
        raise ToolInterruptedError("external runner interrupted by SIGTERM")

    signal.signal(signal.SIGTERM, handle_sigterm)
    try:
        yield
    finally:
        signal.signal(signal.SIGTERM, previous)


def _run_command(
    *,
    command: Sequence[str],
    cwd: Path,
    environment: Mapping[str, str],
    log_handle: TextIO,
    timeout_seconds: int,
) -> int:
    popen_options: dict[str, Any] = {
        "cwd": cwd,
        "env": dict(environment),
        "stdin": subprocess.DEVNULL,
        "stdout": log_handle,
        "stderr": subprocess.STDOUT,
        "shell": False,
        "text": True,
        "start_new_session": True,
        "close_fds": True,
    }
    # POSIX umask is part of the sandbox boundary.  Passing it on Windows is
    # not portable and can create directories with an unusable DACL under
    # Microsoft Store Python, preventing even the desktop user from auditing
    # or cleaning a completed run.
    if os.name != "nt":
        popen_options["umask"] = 0o077
    try:
        process = subprocess.Popen(
            list(command),
            **popen_options,
        )
    except OSError as exc:
        raise ToolExecutionError(f"cannot start external runner: {exc}") from exc
    try:
        return process.wait(timeout=timeout_seconds)
    except subprocess.TimeoutExpired as exc:
        returncode = _terminate_process_group(process)
        raise ToolTimeoutError(
            f"external runner exceeded timeout of {timeout_seconds} seconds",
            returncode=returncode,
        ) from exc
    except BaseException as exc:
        returncode = _terminate_process_group(process)
        if isinstance(exc, ToolInterruptedError):
            exc.returncode = returncode
        raise


def _verify_immutable_inputs(resolved_inputs: Sequence[ResolvedInput]) -> None:
    changed: list[str] = []
    for item in resolved_inputs:
        current = fingerprint_path(item.path)
        if (
            current.sha256 != item.sha256
            or current.size != item.size_bytes
            or current.mtime_ns != item.mtime_ns
        ):
            changed.append(item.name)
    if changed:
        raise ToolExecutionError(
            "plugin modified read-only input(s): " + ", ".join(sorted(changed))
        )


def execute_tool(
    *,
    tool_manifest_path: Path,
    schema_path: Path,
    mode: str,
    run_id: str,
    sample_id: str,
    supplied_inputs: Mapping[str, Path],
    output_dir: Path,
    resolved_inputs_path: Path,
    attempt_record_path: Path,
    log_path: Path,
    threads: int,
    memory_mb: int,
    cache_policy: str = "isolated_empty_tool_cache",
    execution_purpose: str = "development_only",
    timeout_seconds: int = 3600,
    random_seed: int = 0,
) -> ToolExecutionResult:
    """Execute plugin code now and return only after validating its new VCF."""

    if threads < 1:
        raise ToolContractError("threads must be at least 1")
    if memory_mb < 1:
        raise ToolContractError("memory_mb must be at least 1")
    if cache_policy != "isolated_empty_tool_cache":
        raise ToolContractError(
            "only cache_policy=isolated_empty_tool_cache is supported"
        )
    if timeout_seconds < 1:
        raise ToolContractError("timeout_seconds must be at least 1")
    if execution_purpose not in {"development_only", "formal"}:
        raise ToolContractError("execution_purpose must be development_only or formal")
    manifest = load_tool_manifest(tool_manifest_path, schema_path)
    validate_mode_inputs(
        manifest,
        mode,
        supplied_inputs,
    )
    execution_contract = cast(Mapping[str, Any], manifest["execution"])
    sandbox_backend = cast(str, execution_contract.get("sandbox_backend", "none"))
    trust_level = cast(str, execution_contract.get("trust_level", "untrusted"))
    if sandbox_backend not in {"none", "bwrap", "apptainer"}:
        raise ToolContractError(f"unsupported sandbox backend: {sandbox_backend}")
    if execution_purpose == "formal" and sandbox_backend == "none":
        raise ToolContractError(
            "formal execution refuses sandbox_backend=none; bwrap or "
            "apptainer is required for every plugin, including locally "
            "reviewed plugins"
        )
    isolation_status = (
        "development_only_contract_enforcement"
        if sandbox_backend == "none"
        else f"sandboxed_{sandbox_backend}"
    )
    resolved_inputs = resolve_inputs(
        manifest=manifest,
        tool_manifest_path=tool_manifest_path,
        mode=mode,
        supplied_inputs=supplied_inputs,
        resolved_inputs_path=resolved_inputs_path,
        run_id=run_id,
    )

    output_contract = cast(Mapping[str, Any], manifest["outputs"])
    relative_output = cast(str, output_contract["vcf"])
    final_output_vcf = _safe_output_path(output_dir, relative_output)
    output_root = output_dir.resolve(strict=True)
    for item in resolved_inputs:
        input_path = Path(item.path)
        if input_path == final_output_vcf:
            raise ToolContractError("assigned output VCF aliases a declared input")
        if input_path.is_dir():
            try:
                output_root.relative_to(input_path)
            except ValueError:
                pass
            else:
                raise ToolContractError(
                    f"assigned output directory is nested inside input {item.name}"
                )
    runner_value = cast(str, execution_contract["runner"])
    runner = _resolve_runner(tool_manifest_path, runner_value)
    tool_command = _runner_command(
        runner,
        sandbox_backend=sandbox_backend,
    )
    plugin_root = tool_manifest_path.parent.resolve(strict=True)
    runner_sha256 = sha256_file(runner)
    manifest_sha256 = sha256_file(tool_manifest_path.resolve(strict=True))
    plugin_root_sha256 = fingerprint_path(plugin_root).sha256
    resolved_inputs_sha256 = sha256_file(resolved_inputs_path.resolve(strict=True))

    attempt_id = uuid.uuid4().hex
    attempt_archive_path = _archive_path(attempt_record_path, attempt_id)
    log_archive_path = _archive_path(log_path, attempt_id)
    attempt_work_dir = output_dir.resolve() / ".pgbench_attempts" / attempt_id
    attempt_output_dir = attempt_work_dir / "output"
    attempt_output_dir.mkdir(parents=True, exist_ok=False)
    attempt_output_vcf = _safe_output_path(attempt_output_dir, relative_output)
    environment = build_tool_environment(
        base_environment=os.environ,
        manifest=manifest,
        run_id=run_id,
        sample_id=sample_id,
        resolved_inputs=resolved_inputs,
        resolved_inputs_path=resolved_inputs_path,
        output_dir=attempt_output_dir,
        output_vcf=attempt_output_vcf,
        threads=threads,
        memory_mb=memory_mb,
        cache_policy=cache_policy,
        attempt_work_dir=attempt_work_dir,
        random_seed=random_seed,
    )
    command = _sandbox_command(
        tool_command=tool_command,
        backend=sandbox_backend,
        manifest=manifest,
        plugin_root=plugin_root,
        resolved_inputs=resolved_inputs,
        attempt_work_dir=attempt_work_dir,
        additional_read_only_paths=(
            resolved_inputs_path.resolve(strict=True),
            # Permit adapters to inspect only their own immutable history.
            # The current attempt is rebound read/write by ``_sandbox_command``;
            # its siblings remain read-only.  This supports safe reuse of a
            # completed, content-verified checkpoint after an adapter-only fix
            # without exposing other tools' results or the host filesystem.
            attempt_work_dir.parent,
        ),
    )
    if sandbox_backend == "apptainer":
        environment.update(
            {
                f"APPTAINERENV_{key}": value
                for key, value in environment.items()
                if key.startswith("PGBENCH_")
                or key
                in {
                    "HOME",
                    "TMPDIR",
                    "XDG_CACHE_HOME",
                    "TZ",
                    "LC_ALL",
                    "LANG",
                    "PYTHONHASHSEED",
                    "PYTHONDONTWRITEBYTECODE",
                }
            }
        )

    previous_output_archive = _archive_previous_output(
        attempt_record_path=attempt_record_path,
        final_output_vcf=final_output_vcf,
        output_root=output_root,
        attempt_id=attempt_id,
        tool_id=cast(str, manifest["id"]),
        sample_id=sample_id,
        mode=mode,
    )

    started_at = _utc_now()
    started_at_ns = time.time_ns()
    attempt: dict[str, Any] = {
        "schema_version": 1,
        "attempt_id": attempt_id,
        "run_id": run_id,
        "sample_id": sample_id,
        "tool_id": manifest["id"],
        "tool_version": manifest.get("version"),
        "mode": mode,
        "execution_purpose": execution_purpose,
        "trust_level": trust_level,
        "sandbox_backend": sandbox_backend,
        "isolation_status": isolation_status,
        "formal_score_eligible": execution_purpose == "formal",
        "command": list(command),
        "tool_command": list(tool_command),
        "runner_path": str(runner),
        "runner_sha256": runner_sha256,
        "tool_manifest_path": str(tool_manifest_path.resolve(strict=True)),
        "tool_manifest_sha256": manifest_sha256,
        "plugin_root_sha256": plugin_root_sha256,
        "resolved_inputs_path": str(resolved_inputs_path.resolve(strict=True)),
        "resolved_inputs_sha256": resolved_inputs_sha256,
        "output_dir": str(output_dir.resolve()),
        "attempt_work_dir": str(attempt_work_dir),
        "attempt_output_vcf": str(attempt_output_vcf),
        "expected_output_vcf": str(final_output_vcf),
        "previous_output_archive": (
            previous_output_archive.to_dict()
            if previous_output_archive is not None
            else None
        ),
        "timeout_seconds": timeout_seconds,
        "cache_policy": cache_policy,
        "random_seed": random_seed,
        "started_at": started_at,
        "started_at_ns": started_at_ns,
        "finished_at": None,
        "exit_code": None,
        "status": "running",
        "error": None,
        "output": None,
        "log_path": str(log_path.resolve()),
        "attempt_archive_path": str(attempt_archive_path.resolve()),
        "log_archive_path": str(log_archive_path.resolve()),
    }
    _write_attempt_records(
        current_path=attempt_record_path,
        archive_path=attempt_archive_path,
        payload=attempt,
    )

    temporary_log, log_handle = _prepare_log(log_path)
    exit_code: int | None = None
    validation: ToolOutputValidation | None = None
    failure: BaseException | None = None
    try:
        with _termination_as_exception():
            exit_code = _run_command(
                command=command,
                cwd=attempt_work_dir,
                environment=environment,
                log_handle=log_handle,
                timeout_seconds=timeout_seconds,
            )
        _verify_immutable_inputs(resolved_inputs)
        if sha256_file(runner) != runner_sha256:
            raise ToolExecutionError("plugin runner changed during execution")
        if sha256_file(tool_manifest_path.resolve(strict=True)) != manifest_sha256:
            raise ToolExecutionError("tool manifest changed during execution")
        if fingerprint_path(plugin_root).sha256 != plugin_root_sha256:
            raise ToolExecutionError("plugin package changed during execution")
        if (
            sha256_file(resolved_inputs_path.resolve(strict=True))
            != resolved_inputs_sha256
        ):
            raise ToolExecutionError("resolved_inputs.json changed during execution")
        if exit_code != 0:
            raise ToolExecutionError(f"external runner exited with status {exit_code}")
        candidate_vcf = supplied_inputs.get("candidate_panel")
        validation = validate_tool_output(
            output_vcf=attempt_output_vcf,
            output_dir=attempt_output_dir,
            started_at_ns=started_at_ns,
            sample_id=sample_id,
            candidate_output_contract=cast(
                str, output_contract["candidate_output_contract"]
            ),
            candidate_vcf=(
                candidate_vcf.resolve(strict=True)
                if candidate_vcf is not None
                else None
            ),
        )
        if (
            execution_purpose == "formal"
            and validation.full_compression_validation_status
            == "infrastructure_required"
        ):
            raise ToolExecutionError(
                "formal .vcf.gz acceptance requires external BGZF/index "
                "validation infrastructure"
            )
        if final_output_vcf.exists() or final_output_vcf.is_symlink():
            raise ToolExecutionError(
                "assigned output VCF appeared during execution; refusing overwrite"
            )
        os.link(attempt_output_vcf, final_output_vcf)
        attempt_output_vcf.unlink()
        validation = validate_tool_output(
            output_vcf=final_output_vcf,
            output_dir=output_dir,
            started_at_ns=started_at_ns,
            sample_id=sample_id,
            candidate_output_contract=cast(
                str, output_contract["candidate_output_contract"]
            ),
            candidate_vcf=(
                candidate_vcf.resolve(strict=True)
                if candidate_vcf is not None
                else None
            ),
        )
    except BaseException as exc:
        if isinstance(exc, ToolTimeoutError | ToolInterruptedError):
            exit_code = exc.returncode
        failure = exc
    finally:
        try:
            _commit_log(
                temporary_log,
                log_handle,
                log_path,
                log_archive_path,
            )
        except BaseException as log_exc:
            if failure is None:
                failure = log_exc
        if failure is not None and (
            final_output_vcf.is_file() or final_output_vcf.is_symlink()
        ):
            final_output_vcf.unlink()
        attempt.update(
            {
                "finished_at": _utc_now(),
                "exit_code": exit_code,
                "termination_signal": (
                    -exit_code if isinstance(exit_code, int) and exit_code < 0 else None
                ),
                "timed_out": isinstance(failure, ToolTimeoutError),
                "status": "success" if failure is None else "failed",
                "error": (
                    None if failure is None else f"{type(failure).__name__}: {failure}"
                ),
                "output": validation.to_dict() if validation is not None else None,
                "formal_score_eligible": (
                    execution_purpose == "formal" and failure is None
                ),
            }
        )
        _write_attempt_records(
            current_path=attempt_record_path,
            archive_path=attempt_archive_path,
            payload=attempt,
        )

    if failure is not None:
        if isinstance(failure, ToolContractError | ToolExecutionError):
            raise failure
        if isinstance(failure, ToolOutputValidationError):
            raise ToolExecutionError(
                f"tool output validation failed: {failure}"
            ) from failure
        if isinstance(failure, KeyboardInterrupt | SystemExit):
            raise failure
        raise ToolExecutionError(
            f"external tool execution failed: {failure}"
        ) from failure
    assert validation is not None
    return ToolExecutionResult(
        attempt_id=attempt_id,
        command=command,
        output=validation,
        resolved_inputs_path=str(resolved_inputs_path.resolve(strict=True)),
        attempt_record_path=str(attempt_record_path.resolve(strict=True)),
        attempt_archive_path=str(attempt_archive_path.resolve(strict=True)),
        log_path=str(log_path.resolve(strict=True)),
        log_archive_path=str(log_archive_path.resolve(strict=True)),
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tool-manifest", required=True, type=Path)
    parser.add_argument("--schema", required=True, type=Path)
    parser.add_argument("--mode", required=True, choices=sorted(SUPPORTED_MODES))
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--input", action="append", default=[], metavar="NAME=PATH")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--resolved-inputs", required=True, type=Path)
    parser.add_argument("--attempt-record", required=True, type=Path)
    parser.add_argument("--log", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--memory-mb", type=int, default=1024)
    parser.add_argument(
        "--cache-policy",
        choices=("isolated_empty_tool_cache",),
        default="isolated_empty_tool_cache",
    )
    parser.add_argument(
        "--execution-purpose",
        choices=("development_only", "formal"),
        default="development_only",
    )
    parser.add_argument("--timeout-seconds", type=int, default=3600)
    parser.add_argument("--random-seed", type=int, default=0)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        result = execute_tool(
            tool_manifest_path=args.tool_manifest,
            schema_path=args.schema,
            mode=args.mode,
            run_id=args.run_id,
            sample_id=args.sample_id,
            supplied_inputs=_parse_inputs(args.input),
            output_dir=args.output_dir,
            resolved_inputs_path=args.resolved_inputs,
            attempt_record_path=args.attempt_record,
            log_path=args.log,
            threads=args.threads,
            memory_mb=args.memory_mb,
            cache_policy=args.cache_policy,
            execution_purpose=args.execution_purpose,
            timeout_seconds=args.timeout_seconds,
            random_seed=args.random_seed,
        )
    except (ToolContractError, ToolExecutionError) as exc:
        raise SystemExit(f"pgbench_exec: {exc}") from exc
    print(json.dumps(asdict(result), ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
