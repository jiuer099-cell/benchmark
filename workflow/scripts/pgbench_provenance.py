"""Deterministic provenance primitives for the PGBench workflow.

This module is intentionally independent of Snakemake and of command
execution.  It provides the stable data model used by the workflow runner,
manifest writer, lineage builder, and provenance auditor.
"""

from __future__ import annotations

import csv
import fnmatch
import hashlib
import io
import json
import os
import re
import tempfile
from collections import defaultdict, deque
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

MANIFEST_SCHEMA_VERSION = "pgbench.rule_manifest.v1"
LINEAGE_SCHEMA_VERSION = "pgbench.rule_lineage.v1"
AUDIT_SCHEMA_VERSION = "pgbench.provenance_audit.v1"

# Before the run-context snapshot became the authoritative owner of evaluator
# and scoring profiles, validate_config also recorded these files as inputs.
# validate_config only verifies that the evaluator profile path exists; it does
# not consume either profile's contents.  Old manifests can therefore report
# harmless drift after a profile is frozen correctly by snapshot_run_context.
_LEGACY_VALIDATE_PROFILE_SUFFIXES = (
    "/config/consensus_scoring.yaml",
    "/config/evaluator_profile.yaml",
)

MANIFEST_ID_FIELDS = (
    "manifest_schema_version",
    "run_id",
    "rule_name",
    "job_key",
    "attempt_id",
    "rule_source_sha256",
)

MANIFEST_REQUIRED_FIELDS = (
    "manifest_schema_version",
    "manifest_id",
    "attempt_id",
    "run_id",
    "rule_name",
    "job_key",
    "wildcards",
    "module_or_tool_id",
    "snakefile_path",
    "rule_source_sha256",
    "script_or_wrapper_path",
    "script_or_wrapper_sha256",
    "command",
    "params",
    "threads",
    "requested_resources",
    "input_paths",
    "input_sha256",
    "input_size",
    "input_mtime",
    "output_paths",
    "output_sha256",
    "config_snapshot_sha256",
    "score_profile_sha256",
    "pangenome_manifest_sha256",
    "reference_sha256",
    "truth_profile",
    "conda_lock_sha256",
    "container_uri",
    "container_digest",
    "git_head",
    "git_dirty",
    "git_diff_sha256",
    "snakemake_version",
    "execution_profile",
    "hardware_fingerprint_sha256",
    "random_seed",
    "upstream_manifest_ids",
    "started_at",
    "finished_at",
    "exit_code",
    "status",
)

DEFAULT_CORE_RULE_PATTERNS = (
    "tool__*__*",
    "canonicalize_vcf",
    "link_pangenome_alleles",
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
)

_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_SAFE_KEY_RE = re.compile(r"^[A-Za-z0-9_.-]+$")
_FINAL_STATUSES = {"success", "failed", "skipped", "not_applicable"}
_ALLOWED_STATUSES = _FINAL_STATUSES | {"running"}
_PATH_HASH_FIELDS = ("input_sha256", "output_sha256")
_PATH_METADATA_FIELDS = ("input_size", "input_mtime")
_TOP_LEVEL_SHA256_FIELDS = (
    "rule_source_sha256",
    "script_or_wrapper_sha256",
    "config_snapshot_sha256",
    "score_profile_sha256",
    "pangenome_manifest_sha256",
    "reference_sha256",
    "conda_lock_sha256",
    "git_diff_sha256",
    "hardware_fingerprint_sha256",
)


class ProvenanceError(ValueError):
    """Base class for deterministic provenance failures."""


class PathHashError(ProvenanceError):
    """Raised when an input cannot be represented by the hashing contract."""


class ManifestValidationError(ProvenanceError):
    """Raised when a rule manifest violates the versioned contract."""

    def __init__(self, errors: Sequence[str]) -> None:
        self.errors = tuple(errors)
        super().__init__("; ".join(self.errors))


class LineageError(ProvenanceError):
    """Raised when lineage inputs are ambiguous or malformed."""


@dataclass(frozen=True)
class PathFingerprint:
    """Stable content and lstat metadata for one declared path."""

    path_type: str
    sha256: str
    size: int
    mtime_ns: int
    link_target: str | None = None
    target_sha256: str | None = None

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable record."""

        return asdict(self)


@dataclass(frozen=True)
class ExpectedJob:
    """One executed job expected to have a formal provenance package."""

    rule_name: str
    job_key: str
    manifest_id: str | None = None
    core: bool = False


def canonical_json_bytes(value: Any) -> bytes:
    """Serialize JSON deterministically for hashing and identity generation."""

    try:
        text = json.dumps(
            value,
            ensure_ascii=False,
            allow_nan=False,
            separators=(",", ":"),
            sort_keys=True,
        )
    except (TypeError, ValueError) as exc:
        raise ProvenanceError(f"value is not canonical JSON: {exc}") from exc
    return text.encode("utf-8")


def sha256_bytes(data: bytes) -> str:
    """Return a lowercase hexadecimal SHA-256 digest."""

    return hashlib.sha256(data).hexdigest()


def sha256_file(path: str | os.PathLike[str], *, chunk_size: int = 1024 * 1024) -> str:
    """Hash a regular file without loading it into memory."""

    candidate = Path(path)
    if candidate.is_symlink() or not candidate.is_file():
        raise PathHashError(f"not a regular file: {candidate}")
    digest = hashlib.sha256()
    with candidate.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def _sort_path_key(relative_path: str) -> bytes:
    return os.fsencode(relative_path)


def _path_identity(path: Path) -> tuple[int, int]:
    stat_result = path.stat()
    return stat_result.st_dev, stat_result.st_ino


def _symlink_record(path: Path, active: frozenset[tuple[int, int]]) -> dict[str, Any]:
    target_text = os.readlink(path)
    target = Path(target_text)
    if not target.is_absolute():
        target = path.parent / target
    if not target.exists():
        raise PathHashError(f"broken symlink cannot record target content: {path}")
    target_sha256 = _content_sha256(target, active)
    return {
        "type": "symlink",
        "target": target_text,
        "target_sha256": target_sha256,
    }


def _directory_entries(
    root: Path, active: frozenset[tuple[int, int]]
) -> tuple[list[dict[str, Any]], int]:
    entries: list[dict[str, Any]] = []
    total_size = 0
    stack = [root]
    while stack:
        directory = stack.pop()
        try:
            children = list(os.scandir(directory))
        except OSError as exc:
            raise PathHashError(f"cannot scan directory {directory}: {exc}") from exc
        children.sort(
            key=lambda item: _sort_path_key(
                Path(item.path).relative_to(root).as_posix()
            ),
            reverse=True,
        )
        for child in children:
            child_path = Path(child.path)
            relative = child_path.relative_to(root).as_posix()
            try:
                stat_result = child_path.lstat()
            except OSError as exc:
                raise PathHashError(f"cannot stat {child_path}: {exc}") from exc
            if child.is_symlink():
                record = _symlink_record(child_path, active)
                record["path"] = relative
                record["size"] = stat_result.st_size
                entries.append(record)
                total_size += stat_result.st_size
            elif child.is_file(follow_symlinks=False):
                entries.append(
                    {
                        "path": relative,
                        "type": "file",
                        "sha256": sha256_file(child_path),
                        "size": stat_result.st_size,
                    }
                )
                total_size += stat_result.st_size
            elif child.is_dir(follow_symlinks=False):
                entries.append({"path": relative, "type": "directory"})
                stack.append(child_path)
            else:
                raise PathHashError(f"unsupported filesystem entry: {child_path}")
    entries.sort(key=lambda item: _sort_path_key(str(item["path"])))
    return entries, total_size


def _content_sha256(path: Path, active: frozenset[tuple[int, int]]) -> str:
    if path.is_symlink():
        return sha256_bytes(canonical_json_bytes(_symlink_record(path, active)))
    if not path.exists():
        raise PathHashError(f"path does not exist: {path}")
    if path.is_file():
        return sha256_file(path)
    if path.is_dir():
        identity = _path_identity(path)
        if identity in active:
            raise PathHashError(f"symlink cycle encountered while hashing: {path}")
        entries, _ = _directory_entries(path, active | {identity})
        tree = {"tree_schema": "pgbench.merkle_tree.v1", "entries": entries}
        return sha256_bytes(canonical_json_bytes(tree))
    raise PathHashError(f"unsupported path type: {path}")


def fingerprint_path(path: str | os.PathLike[str]) -> PathFingerprint:
    """Fingerprint a file, directory, or symlink under the fixed hash contract."""

    candidate = Path(path)
    try:
        stat_result = candidate.lstat()
    except OSError as exc:
        raise PathHashError(f"cannot stat {candidate}: {exc}") from exc

    if candidate.is_symlink():
        record = _symlink_record(candidate, frozenset())
        digest = sha256_bytes(canonical_json_bytes(record))
        return PathFingerprint(
            path_type="symlink",
            sha256=digest,
            size=stat_result.st_size,
            mtime_ns=stat_result.st_mtime_ns,
            link_target=str(record["target"]),
            target_sha256=str(record["target_sha256"]),
        )
    if candidate.is_file():
        return PathFingerprint(
            path_type="file",
            sha256=sha256_file(candidate),
            size=stat_result.st_size,
            mtime_ns=stat_result.st_mtime_ns,
        )
    if candidate.is_dir():
        identity = _path_identity(candidate)
        entries, total_size = _directory_entries(candidate, frozenset({identity}))
        tree = {"tree_schema": "pgbench.merkle_tree.v1", "entries": entries}
        return PathFingerprint(
            path_type="directory",
            sha256=sha256_bytes(canonical_json_bytes(tree)),
            size=total_size,
            mtime_ns=stat_result.st_mtime_ns,
        )
    raise PathHashError(f"unsupported path type: {candidate}")


def sha256_path(path: str | os.PathLike[str]) -> str:
    """Return the deterministic content hash for a declared path."""

    return fingerprint_path(path).sha256


def sha256_directory(path: str | os.PathLike[str]) -> str:
    """Hash a directory and reject non-directory inputs."""

    candidate = Path(path)
    if candidate.is_symlink() or not candidate.is_dir():
        raise PathHashError(f"not a directory: {candidate}")
    return fingerprint_path(candidate).sha256


def fingerprint_paths(
    paths: Sequence[str | os.PathLike[str]],
    *,
    base_dir: str | os.PathLike[str] | None = None,
) -> tuple[dict[str, str], dict[str, int], dict[str, int]]:
    """Fingerprint declared paths while preserving their manifest path keys."""

    root = Path(base_dir) if base_dir is not None else Path.cwd()
    hashes: dict[str, str] = {}
    sizes: dict[str, int] = {}
    mtimes: dict[str, int] = {}
    for path_value in paths:
        key = os.fspath(path_value)
        if key in hashes:
            raise PathHashError(f"duplicate declared path: {key}")
        candidate = Path(key)
        if not candidate.is_absolute():
            candidate = root / candidate
        fingerprint = fingerprint_path(candidate)
        hashes[key] = fingerprint.sha256
        sizes[key] = fingerprint.size
        mtimes[key] = fingerprint.mtime_ns
    return hashes, sizes, mtimes


def atomic_write_text(
    path: str | os.PathLike[str], text: str, *, encoding: str = "utf-8"
) -> None:
    """Atomically replace a text file and fsync both file and parent directory."""

    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=destination.parent,
        prefix=f".{destination.name}.",
        suffix=".tmp",
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding=encoding, newline="") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
        try:
            directory_descriptor = os.open(destination.parent, os.O_RDONLY)
        except OSError:
            return
        try:
            os.fsync(directory_descriptor)
        finally:
            os.close(directory_descriptor)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def atomic_write_json(
    path: str | os.PathLike[str], payload: Any, *, indent: int = 2
) -> None:
    """Write JSON through a same-directory temporary file and ``os.replace``."""

    try:
        text = json.dumps(
            payload,
            ensure_ascii=False,
            allow_nan=False,
            indent=indent,
            sort_keys=True,
        )
    except (TypeError, ValueError) as exc:
        raise ProvenanceError(f"payload is not valid JSON: {exc}") from exc
    atomic_write_text(path, f"{text}\n")


def load_json(path: str | os.PathLike[str]) -> Any:
    """Load one UTF-8 JSON document with a provenance-specific error."""

    source = Path(path)
    try:
        with source.open("r", encoding="utf-8") as handle:
            return json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise ProvenanceError(f"cannot load JSON {source}: {exc}") from exc


def build_manifest_id(
    manifest: Mapping[str, Any] | None = None,
    *,
    manifest_schema_version: str | None = None,
    run_id: str | None = None,
    rule_name: str | None = None,
    job_key: str | None = None,
    attempt_id: str | None = None,
    rule_source_sha256: str | None = None,
) -> str:
    """Derive the stable manifest ID from the six normative identity fields."""

    explicit = {
        "manifest_schema_version": manifest_schema_version,
        "run_id": run_id,
        "rule_name": rule_name,
        "job_key": job_key,
        "attempt_id": attempt_id,
        "rule_source_sha256": rule_source_sha256,
    }
    identity: dict[str, Any] = {}
    for field in MANIFEST_ID_FIELDS:
        value = manifest.get(field) if manifest is not None else explicit[field]
        if not isinstance(value, str) or not value:
            raise ManifestValidationError(
                [f"{field} must be a non-empty string to build manifest_id"]
            )
        identity[field] = value
    return sha256_bytes(canonical_json_bytes(identity))


def _parse_timestamp(value: Any, field: str, errors: list[str]) -> datetime | None:
    if value is None:
        return None
    if not isinstance(value, str) or not value:
        errors.append(f"{field} must be an RFC 3339 timestamp or null")
        return None
    normalized = value[:-1] + "+00:00" if value.endswith("Z") else value
    try:
        timestamp = datetime.fromisoformat(normalized)
    except ValueError:
        errors.append(f"{field} is not a valid RFC 3339 timestamp")
        return None
    if timestamp.tzinfo is None:
        errors.append(f"{field} must include a timezone")
        return None
    return timestamp


def _null_reason(manifest: Mapping[str, Any], field: str) -> str | None:
    reasons = manifest.get("not_applicable_reason")
    if isinstance(reasons, str) and reasons:
        return reasons
    if isinstance(reasons, Mapping):
        reason = reasons.get(field)
        if isinstance(reason, str) and reason:
            return reason
    return None


def manifest_validation_errors(
    manifest: Mapping[str, Any],
    *,
    verify_paths: bool = False,
    base_dir: str | os.PathLike[str] | None = None,
    require_success: bool = False,
) -> list[str]:
    """Return all deterministic contract violations for one rule manifest."""

    errors: list[str] = []
    missing = [field for field in MANIFEST_REQUIRED_FIELDS if field not in manifest]
    errors.extend(f"missing required field: {field}" for field in missing)
    if missing:
        return errors

    for field in MANIFEST_ID_FIELDS:
        value = manifest.get(field)
        if not isinstance(value, str) or not value:
            errors.append(f"{field} must be a non-empty string")
    if manifest.get("manifest_schema_version") != MANIFEST_SCHEMA_VERSION:
        errors.append(f"manifest_schema_version must be {MANIFEST_SCHEMA_VERSION!r}")
    for field in ("rule_name", "job_key"):
        value = manifest.get(field)
        if isinstance(value, str) and not _SAFE_KEY_RE.fullmatch(value):
            errors.append(f"{field} contains unsafe path characters")
    if not errors:
        expected_id = build_manifest_id(manifest)
        if manifest.get("manifest_id") != expected_id:
            errors.append("manifest_id does not match the normative identity fields")
    manifest_id = manifest.get("manifest_id")
    if not isinstance(manifest_id, str) or not _SHA256_RE.fullmatch(manifest_id):
        errors.append("manifest_id must be a lowercase 64-character SHA-256")

    mapping_fields = (
        "wildcards",
        "params",
        "requested_resources",
        "input_sha256",
        "input_size",
        "input_mtime",
        "output_sha256",
    )
    for field in mapping_fields:
        if not isinstance(manifest.get(field), Mapping):
            errors.append(f"{field} must be a mapping")
    for field in ("input_paths", "output_paths", "upstream_manifest_ids"):
        value = manifest.get(field)
        if not isinstance(value, list) or not all(
            isinstance(item, str) and item for item in value
        ):
            errors.append(f"{field} must be a list of non-empty strings")
        elif len(value) != len(set(value)):
            errors.append(f"{field} must not contain duplicates")

    threads = manifest.get("threads")
    if isinstance(threads, bool) or not isinstance(threads, int) or threads < 1:
        errors.append("threads must be an integer greater than zero")
    git_dirty = manifest.get("git_dirty")
    if not isinstance(git_dirty, bool):
        errors.append("git_dirty must be boolean")
    random_seed = manifest.get("random_seed")
    if random_seed is not None and (
        isinstance(random_seed, bool) or not isinstance(random_seed, int)
    ):
        errors.append("random_seed must be an integer or null")
    exit_code = manifest.get("exit_code")
    if exit_code is not None and (
        isinstance(exit_code, bool) or not isinstance(exit_code, int)
    ):
        errors.append("exit_code must be an integer or null")

    status = manifest.get("status")
    if status not in _ALLOWED_STATUSES:
        errors.append(f"status must be one of {sorted(_ALLOWED_STATUSES)}")
    if require_success and status != "success":
        errors.append("formal rule manifest must have status='success'")
    if status == "success" and exit_code != 0:
        errors.append("successful manifest must have exit_code=0")
    if status == "failed" and (not isinstance(exit_code, int) or exit_code == 0):
        errors.append("failed manifest must have a non-zero exit_code")

    started = _parse_timestamp(manifest.get("started_at"), "started_at", errors)
    finished = _parse_timestamp(manifest.get("finished_at"), "finished_at", errors)
    if status in _FINAL_STATUSES and finished is None:
        errors.append("finished_at is required for a final status")
    if started is not None and finished is not None and finished < started:
        errors.append("finished_at must not precede started_at")

    for field in _TOP_LEVEL_SHA256_FIELDS:
        value = manifest.get(field)
        if value is not None and (
            not isinstance(value, str) or not _SHA256_RE.fullmatch(value)
        ):
            errors.append(f"{field} must be a lowercase SHA-256 or null")
    score_profile_sha256 = manifest.get("score_profile_sha256")
    if not isinstance(score_profile_sha256, str) or not _SHA256_RE.fullmatch(
        score_profile_sha256
    ):
        errors.append("score_profile_sha256 must be a lowercase SHA-256")
    container_digest = manifest.get("container_digest")
    if container_digest is not None and (
        not isinstance(container_digest, str)
        or not re.fullmatch(r"sha256:[0-9a-f]{64}", container_digest)
    ):
        errors.append("container_digest must have the form sha256:<64 lowercase hex>")

    for field in _PATH_HASH_FIELDS:
        value = manifest.get(field)
        if isinstance(value, Mapping):
            for path_key, digest in value.items():
                if not isinstance(path_key, str) or not path_key:
                    errors.append(f"{field} keys must be non-empty strings")
                if not isinstance(digest, str) or not _SHA256_RE.fullmatch(digest):
                    errors.append(f"{field}[{path_key!r}] must be a lowercase SHA-256")
    for field in _PATH_METADATA_FIELDS:
        value = manifest.get(field)
        if isinstance(value, Mapping):
            for path_key, number in value.items():
                if not isinstance(path_key, str) or not path_key:
                    errors.append(f"{field} keys must be non-empty strings")
                if (
                    isinstance(number, bool)
                    or not isinstance(number, int)
                    or number < 0
                ):
                    errors.append(
                        f"{field}[{path_key!r}] must be a non-negative integer"
                    )

    input_paths = manifest.get("input_paths")
    output_paths = manifest.get("output_paths")
    if isinstance(input_paths, list):
        expected_keys = set(input_paths)
        for field in ("input_sha256", "input_size", "input_mtime"):
            value = manifest.get(field)
            if isinstance(value, Mapping) and set(value) != expected_keys:
                errors.append(f"{field} keys must exactly match input_paths")
    if isinstance(output_paths, list):
        expected_keys = set(output_paths)
        output_hashes = manifest.get("output_sha256")
        if (
            status == "success"
            and isinstance(output_hashes, Mapping)
            and set(output_hashes) != expected_keys
        ):
            errors.append(
                "output_sha256 keys must exactly match output_paths on success"
            )

    upstream = manifest.get("upstream_manifest_ids")
    if isinstance(upstream, list):
        if manifest_id in upstream:
            errors.append("manifest cannot name itself as an upstream manifest")
        for upstream_id in upstream:
            if not _SHA256_RE.fullmatch(upstream_id):
                errors.append(
                    "upstream_manifest_ids values must be lowercase SHA-256 IDs"
                )

    nullable_exemptions = {"finished_at", "exit_code"}
    if status != "running":
        nullable_exemptions.clear()
    for field in MANIFEST_REQUIRED_FIELDS:
        if (
            manifest.get(field) is None
            and field not in nullable_exemptions
            and _null_reason(manifest, field) is None
        ):
            errors.append(
                f"{field} is null but not_applicable_reason does not explain it"
            )

    if verify_paths:
        root = Path(base_dir) if base_dir is not None else Path.cwd()
        for prefix in ("input", "output"):
            if prefix == "output" and status != "success":
                continue
            declared_paths = manifest.get(f"{prefix}_paths")
            declared_hashes = manifest.get(f"{prefix}_sha256")
            if not isinstance(declared_paths, list) or not isinstance(
                declared_hashes, Mapping
            ):
                continue
            for declared_path in declared_paths:
                candidate = Path(declared_path)
                if not candidate.is_absolute():
                    candidate = root / candidate
                try:
                    actual = sha256_path(candidate)
                except PathHashError as exc:
                    errors.append(f"cannot verify {prefix} path {declared_path}: {exc}")
                    continue
                if declared_hashes.get(declared_path) != actual:
                    errors.append(f"{prefix} hash mismatch: {declared_path}")
    return errors


def validate_manifest(
    manifest: Mapping[str, Any],
    *,
    verify_paths: bool = False,
    base_dir: str | os.PathLike[str] | None = None,
    require_success: bool = False,
) -> None:
    """Raise ``ManifestValidationError`` when a manifest is not valid."""

    errors = manifest_validation_errors(
        manifest,
        verify_paths=verify_paths,
        base_dir=base_dir,
        require_success=require_success,
    )
    if errors:
        raise ManifestValidationError(errors)


def is_valid_manifest(
    manifest: Mapping[str, Any],
    *,
    verify_paths: bool = False,
    base_dir: str | os.PathLike[str] | None = None,
    require_success: bool = False,
) -> bool:
    """Return whether a manifest satisfies the requested validation level."""

    return not manifest_validation_errors(
        manifest,
        verify_paths=verify_paths,
        base_dir=base_dir,
        require_success=require_success,
    )


def prepare_rule_manifest(
    payload: Mapping[str, Any],
    *,
    base_dir: str | os.PathLike[str] | None = None,
    calculate_path_metadata: bool = True,
    manifest_output_path: str | os.PathLike[str] | None = None,
) -> dict[str, Any]:
    """Populate deterministic IDs and declared path fingerprints."""

    try:
        manifest = json.loads(canonical_json_bytes(payload))
    except json.JSONDecodeError as exc:  # pragma: no cover - canonical JSON is valid
        raise ProvenanceError(str(exc)) from exc
    if not isinstance(manifest, dict):
        raise ProvenanceError("manifest payload must be a JSON object")
    manifest.setdefault("manifest_schema_version", MANIFEST_SCHEMA_VERSION)

    for field in ("input_paths", "output_paths"):
        paths = manifest.get(field)
        if not isinstance(paths, list):
            raise ManifestValidationError([f"{field} must be a list"])
        if not all(isinstance(path, str) and path for path in paths):
            raise ManifestValidationError(
                [f"{field} must contain only non-empty strings"]
            )

    if manifest_output_path is not None:
        destination = Path(manifest_output_path).absolute()
        root = Path(base_dir).absolute() if base_dir is not None else Path.cwd()
        for output_path in manifest["output_paths"]:
            candidate = Path(output_path)
            if not candidate.is_absolute():
                candidate = root / candidate
            if candidate.absolute() == destination:
                raise ManifestValidationError(
                    ["formal manifest must not include itself in output_paths"]
                )

    if calculate_path_metadata:
        input_hashes, input_sizes, input_mtimes = fingerprint_paths(
            manifest["input_paths"], base_dir=base_dir
        )
        output_hashes, _, _ = fingerprint_paths(
            manifest["output_paths"], base_dir=base_dir
        )
        manifest["input_sha256"] = input_hashes
        manifest["input_size"] = input_sizes
        manifest["input_mtime"] = input_mtimes
        manifest["output_sha256"] = output_hashes
    manifest["manifest_id"] = build_manifest_id(manifest)
    return manifest


def write_rule_manifest(
    payload: Mapping[str, Any],
    output_path: str | os.PathLike[str],
    *,
    base_dir: str | os.PathLike[str] | None = None,
    calculate_path_metadata: bool = True,
    require_success: bool = True,
) -> dict[str, Any]:
    """Prepare, validate, and atomically write one formal rule manifest."""

    manifest = prepare_rule_manifest(
        payload,
        base_dir=base_dir,
        calculate_path_metadata=calculate_path_metadata,
        manifest_output_path=output_path,
    )
    validate_manifest(
        manifest,
        verify_paths=calculate_path_metadata,
        base_dir=base_dir,
        require_success=require_success,
    )
    atomic_write_json(output_path, manifest)
    return manifest


def load_manifest(path: str | os.PathLike[str]) -> dict[str, Any]:
    """Load one manifest JSON object."""

    payload = load_json(path)
    if not isinstance(payload, dict):
        raise ManifestValidationError([f"manifest must be a JSON object: {path}"])
    return payload


def load_manifests(paths: Iterable[str | os.PathLike[str]]) -> list[dict[str, Any]]:
    """Load manifests in stable path order."""

    ordered_paths = sorted(paths, key=lambda item: os.fsencode(os.fspath(item)))
    return [load_manifest(path) for path in ordered_paths]


def discover_manifest_paths(
    manifest_paths: Iterable[str | os.PathLike[str]] = (),
    manifest_directories: Iterable[str | os.PathLike[str]] = (),
) -> list[Path]:
    """Resolve explicit files and recursively discovered JSON manifests."""

    discovered = {Path(path) for path in manifest_paths}
    for directory in manifest_directories:
        root = Path(directory)
        if not root.is_dir():
            raise ProvenanceError(f"manifest directory does not exist: {root}")
        discovered.update(path for path in root.rglob("*.json") if path.is_file())
    return sorted(discovered, key=lambda path: os.fsencode(os.fspath(path)))


def _index_manifests(
    manifests: Sequence[Mapping[str, Any]],
) -> dict[str, Mapping[str, Any]]:
    index: dict[str, Mapping[str, Any]] = {}
    for position, manifest in enumerate(manifests):
        manifest_id = manifest.get("manifest_id")
        if not isinstance(manifest_id, str) or not manifest_id:
            raise LineageError(f"manifest at index {position} has no manifest_id")
        if manifest_id in index:
            raise LineageError(f"duplicate manifest_id: {manifest_id}")
        index[manifest_id] = manifest
    return index


def _find_cycles(
    nodes: set[str], upstream_by_node: Mapping[str, Sequence[str]]
) -> list[list[str]]:
    visiting: set[str] = set()
    visited: set[str] = set()
    path: list[str] = []
    cycles: set[tuple[str, ...]] = set()

    def canonical_cycle(cycle: list[str]) -> tuple[str, ...]:
        body = cycle[:-1]
        rotations = [tuple(body[index:] + body[:index]) for index in range(len(body))]
        return min(rotations)

    def visit(node: str) -> None:
        if node in visited:
            return
        if node in visiting:
            start = path.index(node)
            cycles.add(canonical_cycle(path[start:] + [node]))
            return
        visiting.add(node)
        path.append(node)
        for upstream_id in sorted(upstream_by_node.get(node, ())):
            if upstream_id in nodes:
                visit(upstream_id)
        path.pop()
        visiting.remove(node)
        visited.add(node)

    for node in sorted(nodes):
        visit(node)
    return [list(cycle) + [cycle[0]] for cycle in sorted(cycles)]


def _shared_artifacts(
    upstream: Mapping[str, Any], downstream: Mapping[str, Any]
) -> list[dict[str, str]]:
    output_hashes = upstream.get("output_sha256")
    input_hashes = downstream.get("input_sha256")
    if not isinstance(output_hashes, Mapping) or not isinstance(input_hashes, Mapping):
        return []
    outputs_by_hash: dict[str, list[str]] = defaultdict(list)
    for path, digest in output_hashes.items():
        if isinstance(path, str) and isinstance(digest, str):
            outputs_by_hash[digest].append(path)
    shared: list[dict[str, str]] = []
    for input_path, digest in input_hashes.items():
        if not isinstance(input_path, str) or not isinstance(digest, str):
            continue
        for output_path in outputs_by_hash.get(digest, ()):
            shared.append(
                {
                    "sha256": digest,
                    "upstream_output_path": output_path,
                    "downstream_input_path": input_path,
                }
            )
    return sorted(
        shared,
        key=lambda item: (
            item["sha256"],
            item["upstream_output_path"],
            item["downstream_input_path"],
        ),
    )


def _manifest_attestations(
    upstream: Mapping[str, Any], downstream: Mapping[str, Any]
) -> list[dict[str, str]]:
    """Find cryptographic edges created by consuming an upstream manifest.

    Formal orchestration rules such as lineage builders and auditors consume
    the manifest document itself rather than one of the rule's biological
    outputs.  That is a valid hash link when the downstream input fingerprint
    equals the deterministic on-disk serialization of the supplied upstream
    manifest.  Requiring the exact digest prevents a bare
    ``upstream_manifest_ids`` declaration from being treated as evidence.
    """

    input_hashes = downstream.get("input_sha256")
    manifest_id = upstream.get("manifest_id")
    if not isinstance(input_hashes, Mapping) or not isinstance(manifest_id, str):
        return []
    try:
        serialized = json.dumps(
            upstream,
            ensure_ascii=False,
            allow_nan=False,
            indent=2,
            sort_keys=True,
        )
    except (TypeError, ValueError):
        return []
    digest = sha256_bytes(f"{serialized}\n".encode("utf-8"))
    return sorted(
        [
            {
                "sha256": digest,
                "upstream_output_path": f"manifest:{manifest_id}",
                "downstream_input_path": input_path,
            }
            for input_path, input_digest in input_hashes.items()
            if isinstance(input_path, str) and input_digest == digest
        ],
        key=lambda item: item["downstream_input_path"],
    )


def build_lineage(
    manifests: Sequence[Mapping[str, Any]],
    *,
    target_manifest_ids: Sequence[str] | None = None,
    validate_manifests: bool = True,
) -> dict[str, Any]:
    """Build a deterministic reverse lineage rooted at final target manifests."""

    if validate_manifests:
        for position, manifest in enumerate(manifests):
            errors = manifest_validation_errors(manifest, require_success=True)
            if errors:
                raise LineageError(
                    f"manifest at index {position} is invalid: {'; '.join(errors)}"
                )
    index = _index_manifests(manifests)
    referenced = {
        upstream_id
        for manifest in manifests
        for upstream_id in manifest.get("upstream_manifest_ids", ())
        if isinstance(upstream_id, str)
    }
    if target_manifest_ids is None:
        targets = sorted(set(index) - referenced)
        if not targets and index:
            # A non-empty graph without a sink is cyclic. Selecting every node
            # ensures that an omitted target argument cannot hide that cycle.
            targets = sorted(index)
    else:
        targets = sorted(set(target_manifest_ids))

    missing_targets = sorted(set(targets) - set(index))
    selected: set[str] = set()
    missing_upstream: set[str] = set()
    queue = deque(target for target in targets if target in index)
    while queue:
        manifest_id = queue.popleft()
        if manifest_id in selected:
            continue
        selected.add(manifest_id)
        upstream_ids = index[manifest_id].get("upstream_manifest_ids", ())
        if not isinstance(upstream_ids, list):
            continue
        for upstream_id in upstream_ids:
            if upstream_id in index:
                queue.append(upstream_id)
            else:
                missing_upstream.add(upstream_id)

    upstream_by_node = {
        manifest_id: [
            upstream_id
            for upstream_id in index[manifest_id].get("upstream_manifest_ids", ())
            if upstream_id in selected
        ]
        for manifest_id in selected
    }
    cycles = _find_cycles(selected, upstream_by_node)

    children: dict[str, set[str]] = defaultdict(set)
    indegree = {manifest_id: 0 for manifest_id in selected}
    edges: list[dict[str, Any]] = []
    for downstream_id in sorted(selected):
        for upstream_id in sorted(upstream_by_node[downstream_id]):
            children[upstream_id].add(downstream_id)
            indegree[downstream_id] += 1
            shared = _shared_artifacts(index[upstream_id], index[downstream_id])
            # Prefer the biological/intermediate artifact link already used by
            # historical audits.  A manifest-document attestation is a strict
            # fallback for orchestration-only edges, not an additional edge
            # representation that would rewrite otherwise identical lineage.
            if not shared:
                shared = _manifest_attestations(
                    index[upstream_id], index[downstream_id]
                )
            shared = sorted(
                shared,
                key=lambda item: (
                    item["sha256"],
                    item["upstream_output_path"],
                    item["downstream_input_path"],
                ),
            )
            edges.append(
                {
                    "upstream_manifest_id": upstream_id,
                    "downstream_manifest_id": downstream_id,
                    "hash_linked": bool(shared),
                    "shared_artifacts": shared,
                }
            )

    ready = deque(sorted(node for node, degree in indegree.items() if degree == 0))
    topological_order: list[str] = []
    while ready:
        node = ready.popleft()
        topological_order.append(node)
        for child in sorted(children.get(node, ())):
            indegree[child] -= 1
            if indegree[child] == 0:
                ready.append(child)
        ready = deque(sorted(ready))
    remaining = sorted(selected - set(topological_order))
    node_order = topological_order + remaining

    nodes = []
    for manifest_id in node_order:
        manifest = index[manifest_id]
        nodes.append(
            {
                "manifest_id": manifest_id,
                "run_id": manifest.get("run_id"),
                "rule_name": manifest.get("rule_name"),
                "job_key": manifest.get("job_key"),
                "attempt_id": manifest.get("attempt_id"),
                "status": manifest.get("status"),
                "input_sha256": manifest.get("input_sha256"),
                "output_sha256": manifest.get("output_sha256"),
                "upstream_manifest_ids": manifest.get("upstream_manifest_ids"),
            }
        )

    missing_ids = sorted(set(missing_targets) | missing_upstream)
    # ``all([])`` is true in Python, but a multi-node provenance graph with no
    # edges is not a lineage.  Requiring at least one edge for more than one
    # selected manifest prevents a collection of unrelated manifests from
    # receiving a vacuous ``hash_lineage_complete`` result.
    has_required_edges = len(selected) <= 1 or bool(edges)
    hash_lineage_complete = (
        not missing_ids
        and not cycles
        and has_required_edges
        and all(edge["hash_linked"] for edge in edges)
    )
    return {
        "lineage_schema_version": LINEAGE_SCHEMA_VERSION,
        "status": "complete" if hash_lineage_complete else "incomplete",
        "target_manifest_ids": targets,
        "topological_order": topological_order,
        "nodes": nodes,
        "edges": edges,
        "missing_manifest_ids": missing_ids,
        "cycles": cycles,
        "hash_lineage_complete": hash_lineage_complete,
    }


def lineage_tsv(lineage: Mapping[str, Any]) -> str:
    """Render deterministic lineage edges, including root-only nodes, as TSV."""

    nodes = {
        node["manifest_id"]: node
        for node in lineage.get("nodes", ())
        if isinstance(node, Mapping) and isinstance(node.get("manifest_id"), str)
    }
    rows: list[dict[str, Any]] = []
    downstream_ids = {
        edge.get("downstream_manifest_id")
        for edge in lineage.get("edges", ())
        if isinstance(edge, Mapping)
    }
    for manifest_id, node in sorted(nodes.items()):
        if manifest_id not in downstream_ids:
            rows.append(
                {
                    "upstream_manifest_id": "",
                    "downstream_manifest_id": manifest_id,
                    "upstream_rule_name": "",
                    "downstream_rule_name": node.get("rule_name", ""),
                    "hash_linked": "",
                    "shared_artifact_count": 0,
                }
            )
    for edge in lineage.get("edges", ()):
        if not isinstance(edge, Mapping):
            continue
        upstream_id = str(edge.get("upstream_manifest_id", ""))
        downstream_id = str(edge.get("downstream_manifest_id", ""))
        shared = edge.get("shared_artifacts", ())
        rows.append(
            {
                "upstream_manifest_id": upstream_id,
                "downstream_manifest_id": downstream_id,
                "upstream_rule_name": nodes.get(upstream_id, {}).get("rule_name", ""),
                "downstream_rule_name": nodes.get(downstream_id, {}).get(
                    "rule_name", ""
                ),
                "hash_linked": str(bool(edge.get("hash_linked"))).lower(),
                "shared_artifact_count": len(shared) if isinstance(shared, list) else 0,
            }
        )
    fieldnames = (
        "upstream_manifest_id",
        "downstream_manifest_id",
        "upstream_rule_name",
        "downstream_rule_name",
        "hash_linked",
        "shared_artifact_count",
    )
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(
        stream, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(
        sorted(
            rows,
            key=lambda row: (
                str(row["downstream_manifest_id"]),
                str(row["upstream_manifest_id"]),
            ),
        )
    )
    return stream.getvalue()


def _expected_jobs(
    manifests: Sequence[Mapping[str, Any]],
    expected_jobs: Sequence[ExpectedJob | Mapping[str, Any]] | None,
) -> list[ExpectedJob]:
    if expected_jobs is None:
        return [
            ExpectedJob(
                rule_name=str(manifest.get("rule_name", "")),
                job_key=str(manifest.get("job_key", "")),
                manifest_id=(
                    str(manifest["manifest_id"])
                    if isinstance(manifest.get("manifest_id"), str)
                    else None
                ),
            )
            for manifest in manifests
        ]
    normalized: list[ExpectedJob] = []
    for item in expected_jobs:
        if isinstance(item, ExpectedJob):
            normalized.append(item)
            continue
        if not isinstance(item, Mapping):
            raise ProvenanceError("expected jobs must be objects")
        rule_name = item.get("rule_name")
        job_key = item.get("job_key")
        manifest_id = item.get("manifest_id")
        core = item.get("core", False)
        if not isinstance(rule_name, str) or not rule_name:
            raise ProvenanceError("expected job rule_name must be non-empty")
        if not isinstance(job_key, str) or not job_key:
            raise ProvenanceError("expected job job_key must be non-empty")
        if not _SAFE_KEY_RE.fullmatch(rule_name):
            raise ProvenanceError("expected job rule_name contains unsafe characters")
        if not _SAFE_KEY_RE.fullmatch(job_key):
            raise ProvenanceError("expected job job_key contains unsafe characters")
        if manifest_id is not None and not isinstance(manifest_id, str):
            raise ProvenanceError("expected job manifest_id must be a string or null")
        if not isinstance(core, bool):
            raise ProvenanceError("expected job core must be boolean")
        normalized.append(ExpectedJob(rule_name, job_key, manifest_id, core))
    keys = [(job.rule_name, job.job_key) for job in normalized]
    if len(keys) != len(set(keys)):
        raise ProvenanceError("expected jobs contain duplicate rule_name/job_key pairs")
    return sorted(normalized, key=lambda job: (job.rule_name, job.job_key))


def _matches_core(rule_name: str, patterns: Sequence[str]) -> bool:
    return any(fnmatch.fnmatchcase(rule_name, pattern) for pattern in patterns)


def _companion_path_candidates(
    root: Path,
    job: ExpectedJob,
    *,
    run_id: Any,
) -> tuple[tuple[Path, ...], tuple[Path, ...]]:
    """Return preferred run-scoped and legacy companion locations.

    Current workflows isolate companions below ``logs/<run_id>`` and
    ``benchmarks/<run_id>``.  The unscoped paths remain a read-only fallback so
    manifests produced by releases before run isolation can still be audited.
    Unsafe run identifiers are never interpolated into a filesystem path.
    """

    legacy = (
        root / "logs" / "rules" / job.rule_name / f"{job.job_key}.log",
        root / "benchmarks" / "rules" / job.rule_name / f"{job.job_key}.jsonl",
    )
    if not isinstance(run_id, str) or not _SAFE_KEY_RE.fullmatch(run_id):
        return ((legacy[0],), (legacy[1],))
    scoped = (
        root / "logs" / run_id / "rules" / job.rule_name / f"{job.job_key}.log",
        root
        / "benchmarks"
        / run_id
        / "rules"
        / job.rule_name
        / f"{job.job_key}.jsonl",
    )
    return ((scoped[0], legacy[0]), (scoped[1], legacy[1]))


def audit_manifests(
    manifests: Sequence[Mapping[str, Any]],
    *,
    expected_jobs: Sequence[ExpectedJob | Mapping[str, Any]] | None = None,
    target_manifest_ids: Sequence[str] | None = None,
    core_rule_patterns: Sequence[str] | None = None,
    workspace_root: str | os.PathLike[str] | None = None,
    require_companions: bool = False,
    verify_paths: bool = False,
) -> dict[str, Any]:
    """Audit formal manifests and derive the four scoring traceability fields."""

    effective_core_patterns = (
        DEFAULT_CORE_RULE_PATTERNS
        if core_rule_patterns is None
        else tuple(core_rule_patterns)
    )
    jobs = _expected_jobs(manifests, expected_jobs)
    manifest_by_id: dict[str, Mapping[str, Any]] = {}
    manifests_by_key: dict[tuple[str, str], list[Mapping[str, Any]]] = defaultdict(list)
    issues: list[dict[str, Any]] = []
    duplicate_manifest_id = False
    for manifest in manifests:
        manifest_id = manifest.get("manifest_id")
        if isinstance(manifest_id, str):
            if manifest_id in manifest_by_id:
                duplicate_manifest_id = True
                issues.append(
                    {
                        "code": "duplicate_manifest_id",
                        "severity": "error",
                        "manifest_id": manifest_id,
                        "message": "manifest_id appears more than once",
                    }
                )
            else:
                manifest_by_id[manifest_id] = manifest
        key = (str(manifest.get("rule_name", "")), str(manifest.get("job_key", "")))
        manifests_by_key[key].append(manifest)

    valid_packages = 0
    invalid_core = duplicate_manifest_id
    valid_manifests: list[Mapping[str, Any]] = []
    selected_by_job: dict[tuple[str, str], Mapping[str, Any]] = {}
    root = Path(workspace_root) if workspace_root is not None else Path.cwd()

    # A legacy validate_config manifest may contain the scoring/evaluator
    # profiles even though validate_config never reads their contents.  Treat
    # drift as superseded only when snapshot_run_context independently froze
    # the *current* content of that exact profile path.  This keeps the
    # compatibility exception narrow and still rejects unfrozen profile drift.
    frozen_legacy_profiles: dict[str, str] = {}
    for manifest in manifests:
        if manifest.get("rule_name") != "snapshot_run_context":
            continue
        input_hashes = manifest.get("input_sha256")
        if not isinstance(input_hashes, Mapping):
            continue
        for raw_path, digest in input_hashes.items():
            normalized = "/" + str(raw_path).replace("\\", "/").lstrip("/")
            if not normalized.endswith(_LEGACY_VALIDATE_PROFILE_SUFFIXES):
                continue
            candidate = Path(str(raw_path))
            if not candidate.is_absolute():
                candidate = root / candidate
            if (
                isinstance(digest, str)
                and candidate.is_file()
                and sha256_file(candidate) == digest
            ):
                frozen_legacy_profiles[str(raw_path)] = digest

    for job in jobs:
        core = job.core or _matches_core(job.rule_name, effective_core_patterns)
        matches = manifests_by_key.get((job.rule_name, job.job_key), [])
        selected_manifest = (
            manifest_by_id.get(job.manifest_id)
            if job.manifest_id is not None
            else (matches[0] if len(matches) == 1 else None)
        )
        if selected_manifest is not None and (
            selected_manifest.get("rule_name") != job.rule_name
            or selected_manifest.get("job_key") != job.job_key
        ):
            issues.append(
                {
                    "code": "manifest_job_mismatch",
                    "severity": "error" if core else "warning",
                    "manifest_id": selected_manifest.get("manifest_id"),
                    "rule_name": job.rule_name,
                    "job_key": job.job_key,
                    "message": "manifest_id resolves to a different rule_name/job_key",
                }
            )
            selected_manifest = None
        if selected_manifest is None:
            issues.append(
                {
                    "code": "missing_manifest",
                    "severity": "error" if core else "warning",
                    "rule_name": job.rule_name,
                    "job_key": job.job_key,
                    "message": "executed job has no unambiguous manifest",
                }
            )
            invalid_core = invalid_core or core
            continue
        selected_by_job[(job.rule_name, job.job_key)] = selected_manifest
        validation_errors = manifest_validation_errors(
            selected_manifest,
            verify_paths=verify_paths,
            base_dir=root,
            require_success=True,
        )
        retained_validation_errors: list[str] = []
        for message in validation_errors:
            prefix = "input hash mismatch: "
            drift_path = message[len(prefix) :] if message.startswith(prefix) else None
            if (
                job.rule_name == "validate_config"
                and drift_path in frozen_legacy_profiles
                and (
                    "/" + str(drift_path).replace("\\", "/").lstrip("/")
                ).endswith(_LEGACY_VALIDATE_PROFILE_SUFFIXES)
            ):
                issues.append(
                    {
                        "code": "superseded_validation_profile",
                        "severity": "warning",
                        "manifest_id": selected_manifest.get("manifest_id"),
                        "rule_name": job.rule_name,
                        "job_key": job.job_key,
                        "path": drift_path,
                        "message": (
                            "legacy non-semantic validate_config profile drift "
                            "is superseded by the frozen run-context profile"
                        ),
                    }
                )
                continue
            retained_validation_errors.append(message)
        validation_errors = retained_validation_errors
        package_valid = not validation_errors
        for message in validation_errors:
            issues.append(
                {
                    "code": "invalid_manifest",
                    "severity": "error" if core else "warning",
                    "manifest_id": selected_manifest.get("manifest_id"),
                    "rule_name": job.rule_name,
                    "job_key": job.job_key,
                    "message": message,
                }
            )
        if not validation_errors:
            valid_manifests.append(selected_manifest)
        if require_companions:
            log_paths, benchmark_paths = _companion_path_candidates(
                root,
                job,
                run_id=selected_manifest.get("run_id"),
            )
            for kind, candidates in (
                ("log", log_paths),
                ("benchmark", benchmark_paths),
            ):
                if not any(candidate.is_file() for candidate in candidates):
                    package_valid = False
                    issues.append(
                        {
                            "code": f"missing_{kind}",
                            "severity": "error" if core else "warning",
                            "manifest_id": selected_manifest.get("manifest_id"),
                            "rule_name": job.rule_name,
                            "job_key": job.job_key,
                            "path": str(candidates[0]),
                            "searched_paths": [
                                str(candidate) for candidate in candidates
                            ],
                            "message": (
                                f"required {kind} record is missing from "
                                "run-scoped and legacy locations"
                            ),
                        }
                    )
        if package_valid:
            valid_packages += 1
        elif core:
            invalid_core = True

    try:
        lineage = build_lineage(
            list(manifest_by_id.values()),
            target_manifest_ids=target_manifest_ids,
            validate_manifests=False,
        )
    except LineageError as exc:
        lineage = {
            "lineage_schema_version": LINEAGE_SCHEMA_VERSION,
            "status": "incomplete",
            "target_manifest_ids": sorted(set(target_manifest_ids or ())),
            "topological_order": [],
            "nodes": [],
            "edges": [],
            "missing_manifest_ids": [],
            "cycles": [],
            "hash_lineage_complete": False,
        }
        issues.append(
            {
                "code": "lineage_error",
                "severity": "error",
                "message": str(exc),
            }
        )
        invalid_core = True

    selected_core_ids = {
        str(manifest.get("manifest_id"))
        for job in jobs
        if job.core or _matches_core(job.rule_name, effective_core_patterns)
        for manifest in [selected_by_job.get((job.rule_name, job.job_key))]
        if manifest is not None
    }
    for edge in lineage.get("edges", ()):
        if (
            isinstance(edge, Mapping)
            and not edge.get("hash_linked")
            and (
                edge.get("upstream_manifest_id") in selected_core_ids
                or edge.get("downstream_manifest_id") in selected_core_ids
            )
        ):
            invalid_core = True
    if lineage.get("missing_manifest_ids") or lineage.get("cycles"):
        invalid_core = True

    denominator = len(jobs)
    manifest_completeness = valid_packages / denominator if denominator else 0.0
    hash_lineage_complete = bool(lineage.get("hash_lineage_complete"))
    environment_complete = bool(valid_manifests) and all(
        manifest.get("conda_lock_sha256") is not None
        or manifest.get("container_digest") is not None
        for manifest in valid_manifests
    )
    run_context_complete = bool(valid_manifests) and all(
        isinstance(manifest.get("config_snapshot_sha256"), str)
        and bool(_SHA256_RE.fullmatch(manifest["config_snapshot_sha256"]))
        and isinstance(manifest.get("score_profile_sha256"), str)
        and bool(_SHA256_RE.fullmatch(manifest["score_profile_sha256"]))
        and isinstance(manifest.get("run_id"), str)
        and bool(manifest.get("run_id"))
        and isinstance(manifest.get("truth_profile"), str)
        and bool(manifest.get("truth_profile"))
        and isinstance(manifest.get("git_head"), str)
        and bool(re.fullmatch(r"[0-9a-f]{40}", manifest["git_head"]))
        and isinstance(manifest.get("git_dirty"), bool)
        and isinstance(manifest.get("git_diff_sha256"), str)
        and bool(_SHA256_RE.fullmatch(manifest["git_diff_sha256"]))
        and isinstance(manifest.get("snakemake_version"), str)
        and bool(manifest.get("snakemake_version"))
        and isinstance(manifest.get("hardware_fingerprint_sha256"), str)
        and bool(_SHA256_RE.fullmatch(manifest["hardware_fingerprint_sha256"]))
        and isinstance(manifest.get("execution_profile"), str)
        and bool(manifest.get("execution_profile"))
        and isinstance(manifest.get("random_seed"), int)
        and not isinstance(manifest.get("random_seed"), bool)
        for manifest in valid_manifests
    )
    if run_context_complete:
        # Every audited job must describe the same frozen run context.  A
        # complete-looking manifest cannot silently substitute a different
        # config, Git state, hardware identity, profile, or random seed.
        context_fields = (
            "run_id",
            "config_snapshot_sha256",
            "score_profile_sha256",
            "truth_profile",
            "git_head",
            "git_dirty",
            "git_diff_sha256",
            "snakemake_version",
            "hardware_fingerprint_sha256",
            "execution_profile",
            "random_seed",
        )
        mismatched_context_fields = [
            field
            for field in context_fields
            if len({manifest.get(field) for manifest in valid_manifests}) != 1
        ]
        if mismatched_context_fields:
            run_context_complete = False
            invalid_core = True
            for field in mismatched_context_fields:
                issues.append(
                    {
                        "code": "run_context_mismatch",
                        "severity": "error",
                        "field": field,
                        "message": (
                            "audited manifests disagree on frozen run-context "
                            f"field {field}"
                        ),
                    }
                )
    scoring_manifests = [
        manifest
        for manifest in valid_manifests
        if manifest.get("rule_name") == "compute_pgbench_score"
    ]
    if scoring_manifests:
        scoring_profile_declared = all(
            isinstance(manifest.get("params"), Mapping)
            and isinstance(manifest["params"].get("score_profile"), str)
            and bool(manifest["params"].get("score_profile"))
            for manifest in scoring_manifests
        )
        run_context_complete = run_context_complete and scoring_profile_declared
        if not scoring_profile_declared:
            invalid_core = True
            issues.append(
                {
                    "code": "missing_scoring_profile",
                    "severity": "error",
                    "message": (
                        "compute_pgbench_score must declare params.score_profile"
                    ),
                }
            )

    core_provenance_valid = not invalid_core
    if not core_provenance_valid:
        status = "invalid"
    elif (
        manifest_completeness == 1.0
        and hash_lineage_complete
        and environment_complete
        and run_context_complete
    ):
        status = "valid"
    else:
        status = "provisional"

    return {
        "audit_schema_version": AUDIT_SCHEMA_VERSION,
        "status": status,
        "audited_job_count": denominator,
        "valid_manifest_log_benchmark_count": valid_packages,
        "manifest_completeness": manifest_completeness,
        "manifest_points": 2.0 * manifest_completeness,
        "hash_lineage_complete": hash_lineage_complete,
        "environment_complete": environment_complete,
        "run_context_complete": run_context_complete,
        "core_provenance_valid": core_provenance_valid,
        "lineage": lineage,
        "issues": sorted(
            issues,
            key=lambda issue: (
                str(issue.get("severity", "")),
                str(issue.get("code", "")),
                str(issue.get("rule_name", "")),
                str(issue.get("job_key", "")),
                str(issue.get("message", "")),
            ),
        ),
    }
