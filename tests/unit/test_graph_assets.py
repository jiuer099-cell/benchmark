from __future__ import annotations

import sys
from pathlib import Path

import pytest
import yaml

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from build_pangenome_manifest import (  # noqa: E402
    PangenomeManifestError,
    load_graph_assets_lock,
)
from validate_graph_assets import (  # noqa: E402
    GraphAssetError,
    lock_graph_assets,
    sha256_file,
)


def _write_bundle(tmp_path: Path) -> dict[str, Path]:
    root = tmp_path / "graph"
    root.mkdir()
    paths = {
        "source_manifest": root / "graph-assets.lock.yaml",
        "gbz": root / "graph.gbz",
        "xg": root / "graph.xg",
        "min_index": root / "graph.min",
        "dist": root / "graph.dist",
        "sample_list": root / "samples.txt",
    }
    paths["source_manifest"].write_text(
        "producer: vg autoindex\nreference_path: GRCh38\n",
        encoding="utf-8",
    )
    for key in ("gbz", "xg", "min_index", "dist"):
        paths[key].write_bytes(f"synthetic-{key}".encode())
    paths["sample_list"].write_text("HG001\nHG003\n", encoding="utf-8")
    return paths


def test_graph_bundle_is_content_locked(tmp_path: Path) -> None:
    paths = _write_bundle(tmp_path)
    lock = lock_graph_assets(**paths, reference_path="GRCh38")

    assert lock["schema_version"] == 1
    assert lock["asset_root"] == str((tmp_path / "graph").resolve())
    assert lock["reference_path"] == "GRCh38"
    assert lock["sample_count"] == 2
    assert set(lock["assets"]) == {"gbz", "xg", "min", "dist", "sample_list"}
    assert lock["assets"]["gbz"]["sha256"] == sha256_file(paths["gbz"])
    assert lock["source_manifest"]["path"] == str(paths["source_manifest"])


def test_giraffe_shortread_indexes_are_content_locked(tmp_path: Path) -> None:
    paths = _write_bundle(tmp_path)
    paths["min_index"] = tmp_path / "graph" / "graph.shortread.withzip.min"
    paths["min_index"].write_bytes(b"synthetic-minimizer-with-zipcodes")
    paths["zipcodes"] = tmp_path / "graph" / "graph.shortread.zipcodes"
    paths["zipcodes"].write_bytes(b"synthetic-zipcodes")

    lock = lock_graph_assets(
        **paths,
        reference_path="GRCh38",
        profile="vg_giraffe_shortread",
    )

    assert set(lock["assets"]) == {
        "gbz",
        "min",
        "zipcodes",
        "dist",
        "sample_list",
    }


def test_source_lock_checksum_mismatch_is_rejected(tmp_path: Path) -> None:
    paths = _write_bundle(tmp_path)
    paths["source_manifest"].write_text(
        "reference_path: GRCh38\n"
        "assets:\n"
        "  gbz:\n"
        f"    sha256: {'0' * 64}\n",
        encoding="utf-8",
    )
    with pytest.raises(GraphAssetError, match="checksum mismatch"):
        lock_graph_assets(**paths, reference_path="GRCh38")


def test_assets_must_use_one_conventional_directory(tmp_path: Path) -> None:
    paths = _write_bundle(tmp_path)
    outside = tmp_path / "graph.xg"
    outside.write_bytes(b"xg")
    paths["xg"] = outside
    with pytest.raises(GraphAssetError, match="one directory"):
        lock_graph_assets(**paths, reference_path="GRCh38")


def test_duplicate_graph_samples_are_rejected(tmp_path: Path) -> None:
    paths = _write_bundle(tmp_path)
    paths["sample_list"].write_text("HG001\nHG001\n", encoding="utf-8")
    with pytest.raises(GraphAssetError, match="duplicate"):
        lock_graph_assets(**paths, reference_path="GRCh38")


def test_hg002_alias_in_graph_is_rejected(tmp_path: Path) -> None:
    paths = _write_bundle(tmp_path)
    paths["sample_list"].write_text(
        "HG001\nHG002#1#chr1\n",
        encoding="utf-8",
    )
    with pytest.raises(GraphAssetError, match="excluded HG002"):
        lock_graph_assets(**paths, reference_path="GRCh38")


def test_pangenome_embedding_rechecks_locked_content(tmp_path: Path) -> None:
    paths = _write_bundle(tmp_path)
    payload = lock_graph_assets(**paths, reference_path="GRCh38")
    lock_path = tmp_path / "results" / "graph-assets.lock.yaml"
    lock_path.parent.mkdir()
    lock_path.write_text(
        yaml.safe_dump(payload, sort_keys=False),
        encoding="utf-8",
    )

    embedded = load_graph_assets_lock(lock_path)
    assert embedded["manifest"]["sha256"] == sha256_file(paths["source_manifest"])
    assert embedded["lock_manifest"]["sha256"] == sha256_file(lock_path)
    assert embedded["asset_root"] == str((tmp_path / "graph").resolve())
    assert embedded["sample_count"] == 2
    assert embedded["excluded_samples"] == ["HG002", "NA24385"]

    paths["gbz"].write_bytes(b"changed-after-lock")
    with pytest.raises(PangenomeManifestError, match="no longer matches"):
        load_graph_assets_lock(lock_path)
