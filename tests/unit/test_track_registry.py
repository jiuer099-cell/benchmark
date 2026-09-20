from __future__ import annotations

import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))

from track_registry import (  # noqa: E402
    TrackRegistryError,
    load_track_registry,
    resolve_track,
)


def test_registry_declares_four_isolated_tracks_and_ont_pass_policy() -> None:
    registry = load_track_registry(ROOT / "config" / "track_registry.yaml")
    assert set(registry) == {"sr_illumina", "lr_hifi", "lr_clr", "lr_ont"}
    assert resolve_track(registry, "short_read_fixed_panel_genotyping").id == "sr_illumina"
    ont = registry["lr_ont"]
    assert ont.frozen_input == {
        "dataset_id": "HG002_ONT_R9.4.1_Guppy5.0.6_SUP_pass_v1",
        "object_name": "basecalls.fastq.gz",
        "read_selection_policy": "official_qscore_pass_only",
        "prohibited_read_selection": [
            "pass_plus_fail", "fail_only", "adapter_selected_subset"
        ],
    }


def test_registry_rejects_an_incomplete_frozen_input(tmp_path: Path) -> None:
    registry = (ROOT / "config" / "track_registry.yaml").read_text(encoding="utf-8")
    tmp = tmp_path / "tracks.yaml"
    tmp.write_text(
        registry.replace("      object_name: basecalls.fastq.gz\n", ""),
        encoding="utf-8",
    )
    with pytest.raises(TrackRegistryError, match="immutable input identity"):
        load_track_registry(tmp)
