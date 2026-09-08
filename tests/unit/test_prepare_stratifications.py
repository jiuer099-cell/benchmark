from __future__ import annotations

import gzip
import hashlib
import json
import sys
from pathlib import Path

import pytest
import yaml

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from prepare_stratifications import StratificationError, prepare  # noqa: E402


def test_existing_context_bed_is_hashed_and_frozen(tmp_path: Path) -> None:
    bed = tmp_path / "Union" / "context.bed.gz"
    bed.parent.mkdir()
    with gzip.open(bed, "wt", encoding="utf-8") as handle:
        handle.write("chr1\t0\t100\n")
    catalogue = tmp_path / "strata.yaml"
    catalogue.write_text(
        yaml.safe_dump(
            {
                "schema_version": 2,
                "profile_id": "test",
                "reference_id": "grch38",
                "release": "v1",
                "source_root": "https://example.invalid",
                "local_root": str(tmp_path),
                "strata": {
                    "difficult": {
                        "relative_path": "Union/context.bed.gz",
                        "md5": hashlib.md5(bed.read_bytes()).hexdigest(),  # noqa: S324
                        "category": "difficulty",
                    }
                },
            }
        ),
        encoding="utf-8",
    )
    output = tmp_path / "lock.json"

    payload = prepare(catalogue, output, download=False)

    assert payload["assets"]["difficult"]["interval_count"] == 1
    assert len(payload["assets"]["difficult"]["sha256"]) == 64
    assert json.loads(output.read_text())["contract"] == (
        "pgbench_context_stratifications_v1"
    )


def test_stratification_path_traversal_is_rejected(tmp_path: Path) -> None:
    catalogue = tmp_path / "bad.yaml"
    catalogue.write_text(
        yaml.safe_dump(
            {
                "schema_version": 2,
                "profile_id": "test",
                "reference_id": "grch38",
                "release": "v1",
                "source_root": "https://example.invalid",
                "local_root": str(tmp_path),
                "strata": {
                    "bad": {
                        "relative_path": "../escape.bed.gz",
                        "md5": "0" * 32,
                    }
                },
            }
        ),
        encoding="utf-8",
    )
    with pytest.raises(StratificationError, match="unsafe"):
        prepare(catalogue, tmp_path / "lock.json", download=False)
