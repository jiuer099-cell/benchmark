"""Regression gates for the tool-neutral external-adapter boundary."""

from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CORE_SOURCE_ROOTS = (
    ROOT / "workflow" / "rules",
    ROOT / "workflow" / "modules",
    ROOT / "workflow" / "scripts",
    ROOT / "workflow" / "schemas",
)


def _core_sources() -> tuple[Path, ...]:
    """Return every executable/declarative Core source, not a hand-picked subset."""

    sources = {
        ROOT / "Snakefile",
        ROOT / "config" / "config.schema.yaml",
        ROOT / "workflow" / "rule-registry.yaml",
    }
    for source_root in CORE_SOURCE_ROOTS:
        sources.update(
            path
            for path in source_root.rglob("*")
            if path.is_file()
            and (path.suffix in {".py", ".smk", ".yaml", ".yml"} or path.name == "Snakefile")
        )
    return tuple(sorted(sources))


def _adapter_manifests() -> dict[str, dict[str, object]]:
    """Discover adapters rather than teaching this CI gate each tool name."""

    import yaml

    manifests: dict[str, dict[str, object]] = {}
    for path in sorted((ROOT / "plugins").glob("*/tool.yaml")):
        manifest = yaml.safe_load(path.read_text(encoding="utf-8"))
        assert isinstance(manifest, dict), f"invalid adapter manifest: {path}"
        tool_id = manifest.get("id")
        assert isinstance(tool_id, str) and tool_id, f"missing adapter id: {path}"
        manifests[tool_id] = manifest
    return manifests


def test_core_boundary_has_no_current_adapter_identity() -> None:
    """Adding or changing an adapter must not reintroduce Core branching."""

    adapter_ids = set(_adapter_manifests())
    for source in _core_sources():
        text = source.read_text(encoding="utf-8").casefold()
        hits = sorted(
            tool for tool in adapter_ids
            if re.search(rf"(?<![a-z0-9]){re.escape(tool)}(?![a-z0-9])", text)
        )
        assert not hits, f"tool-specific Core logic in {source}: {hits}"


def test_core_source_names_have_no_current_adapter_identity() -> None:
    """Tool-owned helpers must live below their adapter, never in Core."""

    adapter_ids = set(_adapter_manifests())
    for source in _core_sources():
        normalized_name = source.name.casefold().replace("-", "_")
        hits = sorted(
            tool
            for tool in adapter_ids
            if re.search(
                rf"(?:^|_){re.escape(tool.replace('-', '_'))}(?:_|\.|$)",
                normalized_name,
            )
        )
        assert not hits, f"tool-specific Core filename {source}: {hits}"


def test_every_discovered_adapter_declares_external_contract() -> None:
    """A new tool is admitted by an adapter manifest alone, never a Core branch."""

    manifests = _adapter_manifests()
    assert manifests
    for tool_id, manifest in manifests.items():
        assert manifest["source"] == "external"
        assert manifest["capabilities"]["read_class"] in {"short", "long"}
        assert "native_assets" in manifest
        assert manifest["interface_version"] == 2
