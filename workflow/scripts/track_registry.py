"""Declarative PGBench technology-track registry.

The registry is deliberately data-driven: adding a technology requires a
registry entry and an adapter declaration, never a Core code branch.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml  # type: ignore[import-untyped]


class TrackRegistryError(ValueError):
    """Raised for an invalid or internally inconsistent track registry."""


@dataclass(frozen=True)
class Track:
    id: str
    release_id: str
    read_class: str
    technology: str
    legacy_ids: tuple[str, ...]
    frozen_input: Mapping[str, Any] | None


def default_registry_path(repo_root: Path) -> Path:
    return repo_root / "config" / "track_registry.yaml"


def load_track_registry(path: Path) -> dict[str, Track]:
    try:
        loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise TrackRegistryError(f"cannot load track registry {path}: {exc}") from exc
    if not isinstance(loaded, Mapping) or loaded.get("schema_version") != 1:
        raise TrackRegistryError("track registry must declare schema_version: 1")
    policy = loaded.get("policy")
    if not isinstance(policy, Mapping) or any(
        policy.get(name) is not False
        for name in (
            "cross_track_ranking",
            "cross_technology_ranking",
            "cross_track_runtime_ranking",
        )
    ):
        raise TrackRegistryError("track registry must forbid all cross-track ranking")
    raw_tracks = loaded.get("tracks")
    if not isinstance(raw_tracks, Mapping) or not raw_tracks:
        raise TrackRegistryError("track registry must contain at least one track")
    result: dict[str, Track] = {}
    identifiers: set[str] = set()
    technologies: set[str] = set()
    for track_id, raw in raw_tracks.items():
        if not isinstance(track_id, str) or not isinstance(raw, Mapping):
            raise TrackRegistryError("track registry entries must be mappings")
        fields = ("release_id", "read_class", "technology", "legacy_ids")
        if any(not raw.get(field) for field in fields[:3]) or not isinstance(raw.get("legacy_ids"), list):
            raise TrackRegistryError(f"track {track_id} is missing required fields")
        read_class = raw["read_class"]
        technology = raw["technology"]
        legacy_ids = raw["legacy_ids"]
        if read_class not in {"short", "long"}:
            raise TrackRegistryError(f"track {track_id} has invalid read_class")
        if not isinstance(technology, str) or not isinstance(raw["release_id"], str):
            raise TrackRegistryError(f"track {track_id} has invalid identity fields")
        if not all(isinstance(item, str) and item for item in legacy_ids):
            raise TrackRegistryError(f"track {track_id} has invalid legacy_ids")
        frozen_input = raw.get("frozen_input")
        if frozen_input is not None:
            if not isinstance(frozen_input, Mapping):
                raise TrackRegistryError(f"track {track_id} frozen_input must be a mapping")
            required_input_fields = ("dataset_id", "object_name", "read_selection_policy")
            if any(
                not isinstance(frozen_input.get(field), str) or not frozen_input[field]
                for field in required_input_fields
            ):
                raise TrackRegistryError(
                    f"track {track_id} frozen_input is missing an immutable input identity"
                )
            prohibited = frozen_input.get("prohibited_read_selection", [])
            if not isinstance(prohibited, list) or not all(
                isinstance(item, str) and item for item in prohibited
            ):
                raise TrackRegistryError(
                    f"track {track_id} has invalid prohibited_read_selection"
                )
        if track_id in identifiers or track_id in legacy_ids or any(item in identifiers for item in legacy_ids):
            raise TrackRegistryError(f"track {track_id} has a duplicate track identifier")
        if technology in technologies:
            raise TrackRegistryError(f"technology {technology} is assigned to multiple tracks")
        identifiers.update((track_id, *legacy_ids))
        technologies.add(technology)
        result[track_id] = Track(
            track_id,
            raw["release_id"],
            read_class,
            technology,
            tuple(legacy_ids),
            dict(frozen_input) if frozen_input is not None else None,
        )
    return result


def resolve_track(registry: Mapping[str, Track], track_id: str) -> Track:
    if track_id in registry:
        return registry[track_id]
    for track in registry.values():
        if track_id in track.legacy_ids:
            return track
    raise TrackRegistryError(f"unsupported benchmark track {track_id}")


def technologies(registry: Mapping[str, Track]) -> set[str]:
    return {track.technology for track in registry.values()}
