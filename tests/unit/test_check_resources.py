from __future__ import annotations

import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from check_resources import (  # noqa: E402
    CALLER_ONLY,
    END_TO_END,
    build_inventory,
    main,
)


def _touch(root: Path, relative: str, content: bytes = b"x") -> str:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)
    return relative


def _config(root: Path) -> dict:
    fasta = _touch(root, "reference.fa", b">chr1\nA\n")
    fai = _touch(root, "reference.fa.fai")
    sequence_dict = _touch(root, "reference.dict")
    bam = _touch(root, "reads/HG002.bam", b"bam")
    _touch(root, "reads/HG002.bam.bai", b"bai")
    fastq = _touch(root, "reads/HG002.fastq.gz", b"fastq")
    truth = _touch(root, "truth/truth.vcf.gz", b"vcf")
    _touch(root, "truth/truth.vcf.gz.tbi", b"tbi")
    bed = _touch(root, "truth/benchmark.bed", b"chr1\t0\t1\n")
    population = _touch(root, "pangenome/population.vcf.gz", b"vcf")
    _touch(root, "pangenome/population.vcf.gz.tbi", b"tbi")
    graph_assets = {
        name: _touch(root, f"pangenome/graph.{name}", name.encode())
        for name in ("manifest", "gbz", "xg", "min", "dist", "sample_list")
    }
    return {
        "execution": {"tracks": [CALLER_ONLY, END_TO_END]},
        "sample": {"bam": bam, "fastq": fastq},
        "caller_only": {"shared_bam": bam},
        "reference": {"fasta": fasta, "fai": fai, "dict": sequence_dict},
        "evaluation": {"truth_vcf": truth, "benchmark_bed": bed},
        "pangenome": {
            "population_vcf": population,
            "build_graph_assets": True,
            "graph_assets": graph_assets,
        },
        "development": {"synthetic_mode": False},
    }


def _by_id(report: dict) -> dict[str, dict]:
    return {asset["id"]: asset for asset in report["assets"]}


def test_all_tracks_are_ready_and_report_exact_sizes(tmp_path: Path) -> None:
    report = build_inventory(_config(tmp_path), repo_root=tmp_path)
    assets = _by_id(report)

    assert report["status"] == "ready"
    assert report["tracks"] == [CALLER_ONLY, END_TO_END]
    assert assets["sample.bam"]["required"] is True
    assert assets["sample.bam"]["size_bytes"] == 3
    assert assets["sample.bai"]["status"] == "ok"
    assert assets["sample.fastq"]["required"] is True
    assert assets["evaluation.truth_vcf_tbi"]["status"] == "ok"
    assert assets["pangenome.population_vcf_tbi"]["status"] == "ok"


def test_track_override_changes_bam_and_fastq_requiredness(tmp_path: Path) -> None:
    config = _config(tmp_path)
    (tmp_path / config["sample"]["fastq"]).unlink()

    caller_report = build_inventory(
        config,
        repo_root=tmp_path,
        track_selection="caller-only",
    )
    end_to_end_report = build_inventory(
        config,
        repo_root=tmp_path,
        track_selection="end-to-end",
    )

    caller_assets = _by_id(caller_report)
    end_to_end_assets = _by_id(end_to_end_report)
    assert caller_report["status"] == "ready"
    assert caller_assets["sample.fastq"]["required"] is False
    assert end_to_end_report["status"] == "incomplete"
    assert end_to_end_assets["sample.bam"]["required"] is False
    assert end_to_end_assets["sample.fastq"]["required"] is True


def test_empty_required_file_fails_preflight(tmp_path: Path) -> None:
    config = _config(tmp_path)
    (tmp_path / config["reference"]["fai"]).write_bytes(b"")

    report = build_inventory(config, repo_root=tmp_path)

    assert report["status"] == "incomplete"
    assert _by_id(report)["reference.fai"]["status"] == "empty"
    assert "reference.fai" in report["missing_required"]


def test_synthetic_smoke_does_not_require_reference_indexes(
    tmp_path: Path,
) -> None:
    config = _config(tmp_path)
    config["development"]["synthetic_mode"] = True
    config["reference"].pop("fai")
    config["reference"].pop("dict")
    (tmp_path / "reference.fa.fai").unlink()
    (tmp_path / "reference.dict").unlink()

    report = build_inventory(config, repo_root=tmp_path)
    assets = _by_id(report)

    assert report["status"] == "ready"
    assert assets["reference.fai"]["required"] is False
    assert assets["reference.dict"]["required"] is False


def test_cli_is_read_only_and_returns_one_for_missing_assets(
    tmp_path: Path,
    capsys,
) -> None:
    config = _config(tmp_path)
    missing_bed = tmp_path / config["evaluation"]["benchmark_bed"]
    missing_bed.unlink()
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
    before = sorted(
        (path.relative_to(tmp_path), path.stat().st_size)
        for path in tmp_path.rglob("*")
        if path.is_file()
    )

    exit_code = main(
        [
            "--config",
            str(config_path),
            "--repo-root",
            str(tmp_path),
            "--json",
        ]
    )
    after = sorted(
        (path.relative_to(tmp_path), path.stat().st_size)
        for path in tmp_path.rglob("*")
        if path.is_file()
    )

    assert exit_code == 1
    assert before == after
    assert '"status": "incomplete"' in capsys.readouterr().out
