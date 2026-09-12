from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))

from freeze_haplotype_source_bundle import freeze, parse_software  # noqa: E402


def test_bundle_lock_seals_required_assets_and_build_software(tmp_path: Path) -> None:
    names = (
        "gfa_or_gbz",
        "population_vcf",
        "sample_roster",
        "haplotype_roster",
        "family_exclusion_manifest",
        "reference",
    )
    assets = []
    for name in names:
        path = tmp_path / name
        path.write_text(name, encoding="utf-8")
        assets.append((name, path))
    tool_sha = "a" * 64

    lock = freeze(
        bundle_id="bundle_v1",
        release="release_v1",
        reference="grch38",
        excluded_samples=["HG002", "HG003", "HG004"],
        assets=assets,
        software=[parse_software(f"vg=1.55.0={tool_sha}")],
    )

    assert lock["family_exclusion"] == ["HG002", "HG003", "HG004"]
    assert lock["assets"]["population_vcf"]["sha256"] == hashlib.sha256(
        b"population_vcf"
    ).hexdigest()
    assert lock["software"] == [
        {"name": "vg", "version": "1.55.0", "sha256": tool_sha}
    ]


def test_bundle_lock_rejects_missing_asset_or_unaddressed_software(tmp_path: Path) -> None:
    asset = tmp_path / "graph.gbz"
    asset.write_text("graph", encoding="utf-8")
    with pytest.raises(ValueError, match="missing required assets"):
        freeze(
            bundle_id="bundle_v1",
            release="release_v1",
            reference="grch38",
            excluded_samples=[],
            assets=[("gfa_or_gbz", asset)],
            software=[],
        )
    with pytest.raises(ValueError, match="NAME=VERSION=SHA256"):
        parse_software("vg=1.55.0")
