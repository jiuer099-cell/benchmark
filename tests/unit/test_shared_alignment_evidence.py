from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from materialize_formal_metrics import evidence_profile  # noqa: E402


def test_evidence_profile_records_shared_bam_as_primary_short_read_input(
    tmp_path: Path,
) -> None:
    resolved = tmp_path / "resolved-inputs.json"
    resolved.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "mode": "end_to_end_from_reads",
                "alignment_kind": "bam",
                "inputs": [
                    {
                        "name": "shared_shortread_alignment",
                        "sha256": "a" * 64,
                        "size_bytes": 9,
                        "path_type": "file",
                    },
                    {
                        "name": "shared_shortread_alignment_index",
                        "sha256": "b" * 64,
                        "size_bytes": 5,
                        "path_type": "file",
                    },
                ],
            }
        ),
        encoding="utf-8",
    )

    profile = evidence_profile(
        resolved,
        technology="illumina_pe",
        library_id="HG002_Illumina_canonical",
        source_evidence_id="HG002_Illumina_paired_fastq",
        coverage_x=None,
        read_count=None,
        read_bases=None,
        downsampling_seed=None,
    )

    assert profile["evidence_kind"] == "shared_alignment"
    assert set(profile["input_assets"]) == {
        "shared_shortread_alignment",
        "shared_shortread_alignment_index",
    }
