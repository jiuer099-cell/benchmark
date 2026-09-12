from __future__ import annotations

import argparse
import hashlib
import importlib.util
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "workflow" / "scripts" / "validate_shared_alignment.py"
SPEC = importlib.util.spec_from_file_location("shared_alignment", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_shared_alignment_lock_binds_one_bam_to_the_frozen_fastq_pair(tmp_path: Path) -> None:
    r1, r2, reference, bam, bai = [tmp_path / name for name in ("r1.fq", "r2.fq", "ref.fa", "shared.bam", "shared.bam.bai")]
    for path, content in zip((r1, r2, reference, bam, bai), (b"r1", b"r2", b"ref", b"bam", b"bai"), strict=True):
        path.write_bytes(content)
    args = argparse.Namespace(
        fastq_r1=r1, fastq_r2=r2, reference=reference, bam=bam, bai=bai,
        source_fastq_sha256=MODULE.paired_fastq_sha256(r1, r2),
        reference_sha256=_sha(reference), bam_sha256=_sha(bam), bai_sha256=_sha(bai),
        aligner="bwa-mem2", aligner_version="2.2.1", command_sha256="a" * 64,
    )
    lock = MODULE.validate(args)
    assert lock["status"] == "verified"
    assert lock["bam_sha256"] == _sha(bam)
    args.bam_sha256 = "0" * 64
    with pytest.raises(MODULE.SharedAlignmentError, match="bam_sha256"):
        MODULE.validate(args)
