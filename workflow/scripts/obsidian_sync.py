"""Deterministic Obsidian design-note synchronization for PGBench."""

from __future__ import annotations

import os
from collections.abc import Mapping
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import yaml  # type: ignore[import-untyped]

try:
    from pgbench_provenance import (
        atomic_write_json,
        atomic_write_text,
        canonical_json_bytes,
        sha256_bytes,
        sha256_file,
    )
    from snapshot_run_context import capture_run_context
except ModuleNotFoundError:  # pragma: no cover - package-style invocation
    from .pgbench_provenance import (
        atomic_write_json,
        atomic_write_text,
        canonical_json_bytes,
        sha256_bytes,
        sha256_file,
    )
    from .snapshot_run_context import capture_run_context


BEGIN_MARKER = "<!-- PGBENCH:BEGIN GENERATED DESIGN -->"
END_MARKER = "<!-- PGBENCH:END GENERATED DESIGN -->"
SYNC_SCHEMA_VERSION = "pgbench.obsidian_sync.v1"
CODE_STATES = {"design-only", "partial", "implemented", "verified"}


class ObsidianSyncError(RuntimeError):
    """Raised when the binding, note, or generated region is invalid."""


@dataclass(frozen=True)
class SyncPaths:
    repo_root: Path
    design: Path
    config: Path
    rule_registry: Path
    tool_root: Path
    binding: Path
    vault_root: Path
    canonical_note: Path
    score_profile: Path
    pangenome_schema: Path


def _utc_now() -> str:
    return datetime.now(UTC).isoformat().replace("+00:00", "Z")


def _load_yaml_mapping(path: Path, label: str) -> dict[str, Any]:
    try:
        value = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise ObsidianSyncError(f"cannot load {label} {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ObsidianSyncError(f"{label} must be a YAML mapping: {path}")
    return value


def _split_frontmatter(text: str) -> tuple[dict[str, Any], str]:
    if not text.startswith("---\n"):
        raise ObsidianSyncError("canonical note must start with YAML frontmatter")
    end = text.find("\n---\n", 4)
    if end < 0:
        raise ObsidianSyncError("canonical note frontmatter is not closed")
    raw = text[4:end]
    try:
        metadata = yaml.safe_load(raw)
    except yaml.YAMLError as exc:
        raise ObsidianSyncError(f"invalid canonical note frontmatter: {exc}") from exc
    if not isinstance(metadata, dict):
        raise ObsidianSyncError("canonical note frontmatter must be a mapping")
    return metadata, text[end + 5 :]


def _render_note(metadata: Mapping[str, Any], body: str) -> str:
    frontmatter = yaml.safe_dump(
        dict(metadata),
        allow_unicode=True,
        default_flow_style=False,
        sort_keys=False,
    ).rstrip()
    return f"---\n{frontmatter}\n---\n\n{body.lstrip()}"


def _replace_generated_design(body: str, design_text: str) -> str:
    begin_count = body.count(BEGIN_MARKER)
    end_count = body.count(END_MARKER)
    if begin_count == 0 and end_count == 0:
        separator = "" if not body.strip() else "\n\n"
        return (
            body.rstrip()
            + separator
            + BEGIN_MARKER
            + "\n"
            + design_text.rstrip()
            + "\n"
            + END_MARKER
            + "\n"
        )
    if begin_count != end_count:
        raise ObsidianSyncError(
            "canonical note has an unbalanced generated-design marker set"
        )
    # The generated design describes this contract and therefore contains the
    # marker literals inside a fenced example. The outer managed region is the
    # first BEGIN marker through the last END marker.
    begin = body.index(BEGIN_MARKER)
    end = body.rindex(END_MARKER)
    if end < begin:
        raise ObsidianSyncError("generated-design end marker precedes begin marker")
    prefix = body[: begin + len(BEGIN_MARKER)]
    suffix = body[end:]
    return f"{prefix}\n{design_text.rstrip()}\n{suffix}"


def _generated_design(body: str) -> str | None:
    begin_count = body.count(BEGIN_MARKER)
    end_count = body.count(END_MARKER)
    if begin_count < 1 or begin_count != end_count:
        return None
    begin = body.index(BEGIN_MARKER) + len(BEGIN_MARKER)
    end = body.rindex(END_MARKER)
    if end < begin:
        return None
    return body[begin:end].strip()


def _binding_frontmatter(path: Path) -> dict[str, Any]:
    metadata, _ = _split_frontmatter(path.read_text(encoding="utf-8"))
    return metadata


def resolve_paths(
    *,
    repo_root: Path,
    design: Path,
    config: Path,
    rule_registry: Path,
    tool_root: Path,
    binding: Path,
    vault_root: Path | None = None,
    canonical_note: Path | None = None,
    pangenome_schema: Path = Path("workflow/schemas/pangenome-manifest.schema.yaml"),
) -> SyncPaths:
    root = repo_root.resolve(strict=True)

    def repo_path(value: Path) -> Path:
        candidate = value if value.is_absolute() else root / value
        return candidate.resolve(strict=True)

    binding_path = repo_path(binding)
    binding_data = _binding_frontmatter(binding_path)
    config_path = repo_path(config)
    config_data = _load_yaml_mapping(config_path, "benchmark config")
    catalogs = config_data.get("catalogs")
    pangenome = config_data.get("pangenome")
    obsidian = config_data.get("obsidian", {})
    if not isinstance(catalogs, Mapping) or not isinstance(pangenome, Mapping):
        raise ObsidianSyncError("config requires catalogs and pangenome mappings")
    if not isinstance(obsidian, Mapping):
        raise ObsidianSyncError("config.obsidian must be a mapping")

    selected_vault = vault_root
    if selected_vault is None:
        env_name = obsidian.get("vault_root_env")
        if isinstance(env_name, str) and os.environ.get(env_name):
            selected_vault = Path(os.environ[env_name])
    if selected_vault is None and isinstance(obsidian.get("vault_root"), str):
        selected_vault = Path(str(obsidian["vault_root"]))
    if selected_vault is None and isinstance(binding_data.get("vault_root"), str):
        selected_vault = Path(str(binding_data["vault_root"]))
    if selected_vault is None:
        raise ObsidianSyncError(
            "cannot resolve vault root from CLI, environment, config, or binding"
        )
    selected_vault = selected_vault.expanduser().resolve(strict=True)

    selected_note = canonical_note
    if selected_note is None and isinstance(obsidian.get("canonical_note"), str):
        selected_note = Path(str(obsidian["canonical_note"]))
    if selected_note is None:
        binding_note = binding_data.get("canonical_note")
        if isinstance(binding_note, str):
            binding_path_value = Path(binding_note)
            project_prefix = Path("Research") / str(
                binding_data.get("project_slug", "")
            )
            try:
                selected_note = binding_path_value.relative_to(project_prefix)
            except ValueError:
                selected_note = binding_path_value
    if selected_note is None:
        raise ObsidianSyncError("cannot resolve canonical note path")
    if selected_note.is_absolute():
        note_path = selected_note
        project_root = selected_note.parent
    else:
        direct_note = selected_vault / selected_note
        project_relpath = obsidian.get("project_relpath")
        nested_note = (
            selected_vault / str(project_relpath) / selected_note
            if isinstance(project_relpath, str) and project_relpath
            else None
        )
        if direct_note.is_file():
            note_path = direct_note
            project_root = selected_vault
        elif nested_note is not None and nested_note.is_file():
            note_path = nested_note
            project_root = selected_vault / str(project_relpath)
        else:
            raise ObsidianSyncError(
                "canonical note does not exist under either the selected "
                "project root or config.obsidian.project_relpath"
            )

    score_value = catalogs.get("score_weights")
    if not isinstance(score_value, str):
        raise ObsidianSyncError("config.catalogs.score_weights must be a path")
    return SyncPaths(
        repo_root=root,
        design=repo_path(design),
        config=config_path,
        rule_registry=repo_path(rule_registry),
        tool_root=repo_path(tool_root),
        binding=binding_path,
        vault_root=project_root.resolve(strict=True),
        canonical_note=note_path.resolve(strict=True),
        score_profile=repo_path(Path(score_value)),
        pangenome_schema=repo_path(pangenome_schema),
    )


def tool_registry_sha256(tool_root: Path, repo_root: Path) -> str | None:
    paths = sorted(
        {
            *tool_root.rglob("tool.yaml"),
            *tool_root.rglob("rule-registry.yaml"),
        },
        key=lambda path: os.fsencode(path.relative_to(repo_root).as_posix()),
    )
    if not paths:
        return None
    inventory = [
        {
            "path": path.relative_to(repo_root).as_posix(),
            "sha256": sha256_file(path),
        }
        for path in paths
    ]
    return sha256_bytes(canonical_json_bytes(inventory))


def expected_metadata(
    paths: SyncPaths,
    *,
    code_state: str,
    snakemake_version: str,
) -> dict[str, Any]:
    if code_state not in CODE_STATES:
        raise ObsidianSyncError(f"code_state must be one of {sorted(CODE_STATES)}")
    config = _load_yaml_mapping(paths.config, "benchmark config")
    context = capture_run_context(
        repo_root=paths.repo_root,
        score_profile=paths.score_profile,
        random_seed=int(config["execution"]["random_seed"]),
        snakemake_version=snakemake_version,
        execution_profile="obsidian-sync",
    )
    pangenome = config["pangenome"]
    score = config["score"]
    design_relative = paths.design.relative_to(paths.repo_root).as_posix()
    return {
        "repo_root": str(paths.repo_root),
        "repo_design_path": design_relative,
        "repo_commit": context["git_head"],
        "git_head": context["git_head"],
        "git_dirty": context["git_dirty"],
        "git_diff_sha256": context["git_diff_sha256"],
        "git_status_sha256": context["git_status_sha256"],
        "untracked_files_sha256": context["untracked_files_sha256"],
        "repo_state_sha256": context["repo_state_sha256"],
        "tracked_tree_sha256": context["tracked_tree_sha256"],
        "design_sha256": sha256_file(paths.design),
        "code_commit": context["git_head"],
        "code_state": code_state,
        "score_profile": score["profile"],
        "score_profile_sha256": sha256_file(paths.score_profile),
        "pangenome_id": pangenome["id"],
        "pangenome_manifest_schema_sha256": sha256_file(paths.pangenome_schema),
        "rule_registry_sha256": sha256_file(paths.rule_registry),
        "tool_registry_sha256": tool_registry_sha256(paths.tool_root, paths.repo_root),
        "sync_status": "stale" if context["git_dirty"] else "current",
    }


def inspect_sync(
    paths: SyncPaths,
    *,
    code_state: str,
    snakemake_version: str,
) -> dict[str, Any]:
    expected = expected_metadata(
        paths,
        code_state=code_state,
        snakemake_version=snakemake_version,
    )
    note_text = paths.canonical_note.read_text(encoding="utf-8")
    metadata, body = _split_frontmatter(note_text)
    mismatches: list[dict[str, Any]] = []
    for key, expected_value in expected.items():
        actual = metadata.get(key)
        if actual != expected_value:
            mismatches.append(
                {
                    "field": key,
                    "expected": expected_value,
                    "actual": actual,
                }
            )
    design_text = paths.design.read_text(encoding="utf-8").strip()
    generated = _generated_design(body)
    if generated != design_text:
        mismatches.append(
            {
                "field": "generated_design",
                "expected_sha256": sha256_bytes(design_text.encode("utf-8")),
                "actual_sha256": (
                    sha256_bytes(generated.encode("utf-8"))
                    if generated is not None
                    else None
                ),
            }
        )
    status = (
        "current"
        if not mismatches and expected["sync_status"] == "current"
        else "stale"
    )
    return {
        "schema_version": SYNC_SCHEMA_VERSION,
        "checked_at": _utc_now(),
        "status": status,
        "repo_root": str(paths.repo_root),
        "vault_root": str(paths.vault_root),
        "canonical_note": str(paths.canonical_note),
        "canonical_note_sha256": sha256_file(paths.canonical_note),
        "expected": expected,
        "mismatches": mismatches,
    }


def synchronize_design(
    paths: SyncPaths,
    *,
    code_state: str,
    snakemake_version: str,
) -> dict[str, Any]:
    note_text = paths.canonical_note.read_text(encoding="utf-8")
    previous_note_sha256 = sha256_bytes(note_text.encode("utf-8"))
    metadata, body = _split_frontmatter(note_text)
    expected = expected_metadata(
        paths,
        code_state=code_state,
        snakemake_version=snakemake_version,
    )
    metadata.update(expected)
    metadata["updated"] = datetime.now(UTC).date().isoformat()
    metadata["updated_at"] = _utc_now()
    design_text = paths.design.read_text(encoding="utf-8")
    body = _replace_generated_design(body, design_text)
    rendered = _render_note(metadata, body)

    temporary = paths.canonical_note.with_name(
        f".{paths.canonical_note.name}.pgbench-sync.tmp"
    )
    atomic_write_text(temporary, rendered)
    verification_text = temporary.read_text(encoding="utf-8")
    verification_metadata, verification_body = _split_frontmatter(verification_text)
    if _generated_design(verification_body) != design_text.strip():
        temporary.unlink(missing_ok=True)
        raise ObsidianSyncError("temporary note failed generated-design check")
    for key, value in expected.items():
        if verification_metadata.get(key) != value:
            temporary.unlink(missing_ok=True)
            raise ObsidianSyncError(f"temporary note failed metadata check for {key}")
    os.replace(temporary, paths.canonical_note)
    result = inspect_sync(
        paths,
        code_state=code_state,
        snakemake_version=snakemake_version,
    )
    result["previous_canonical_note_sha256"] = previous_note_sha256
    return result


def write_check(path: Path, payload: Mapping[str, Any]) -> None:
    atomic_write_json(path, dict(payload))
