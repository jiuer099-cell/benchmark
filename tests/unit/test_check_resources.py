from __future__ import annotations

import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))

from check_resources import TRACK, ResourceCheckError, build_inventory, main  # noqa: E402


def _touch(root: Path, relative: str, content: bytes = b"x") -> str:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)
    return relative


def _config(root: Path) -> dict:
    fasta = _touch(root, "reference.fa", b">chr1\nA\n")
    fai = _touch(root, "reference.fa.fai")
    sequence_dict = _touch(root, "reference.dict")
    r1 = _touch(root, "reads/HG002.R1.fastq.gz", b"fastq1")
    r2 = _touch(root, "reads/HG002.R2.fastq.gz", b"fastq2")
    truth = _touch(root, "truth/truth.vcf.gz", b"vcf")
    _touch(root, "truth/truth.vcf.gz.tbi", b"tbi")
    bed = _touch(root, "truth/benchmark.bed", b"chr1\t0\t1\n")
    population = _touch(root, "pangenome/population.vcf.gz", b"vcf")
    _touch(root, "pangenome/population.vcf.gz.tbi", b"tbi")
    return {
        "execution": {"tracks": [TRACK]},
        "sample": {"fastq_r1": r1, "fastq_r2": r2},
        "reference": {"fasta": fasta, "fai": fai, "dict": sequence_dict},
        "evaluation": {"truth_vcf": truth, "benchmark_bed": bed},
        "pangenome": {
            "population_vcf": population,
            "build_graph_assets": False,
            "graph_assets": {"profile": "none"},
        },
        "development": {"synthetic_mode": False},
    }


def _by_id(report: dict) -> dict[str, dict]:
    return {asset["id"]: asset for asset in report["assets"]}


def test_short_read_inventory_is_ready(tmp_path: Path) -> None:
    report = build_inventory(_config(tmp_path), repo_root=tmp_path)
    assets = _by_id(report)
    assert report["status"] == "ready"
    assert report["tracks"] == [TRACK]
    assert assets["sample.fastq_r1"]["status"] == "ok"
    assert assets["sample.fastq_r2"]["status"] == "ok"
    assert assets["evaluation.truth_vcf_tbi"]["status"] == "ok"


def test_missing_mate_fails_preflight(tmp_path: Path) -> None:
    config = _config(tmp_path)
    (tmp_path / config["sample"]["fastq_r2"]).unlink()
    report = build_inventory(config, repo_root=tmp_path)
    assert report["status"] == "not_ready"
    assert "sample.fastq_r2" in report["missing_required"]


def test_only_fixed_short_read_track_is_accepted(tmp_path: Path) -> None:
    config = _config(tmp_path)
    config["execution"]["tracks"] = ["unsupported_track"]
    try:
        build_inventory(config, repo_root=tmp_path)
    except ResourceCheckError as exc:
        assert "must be exactly" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("unsupported track was accepted")


def test_synthetic_smoke_does_not_require_reference_indexes(tmp_path: Path) -> None:
    config = _config(tmp_path)
    config["development"]["synthetic_mode"] = True
    config["reference"].pop("fai")
    config["reference"].pop("dict")
    report = build_inventory(config, repo_root=tmp_path)
    assert report["status"] == "ready"
    assert _by_id(report)["reference.fai"]["required"] is False


def test_cli_is_read_only_and_reports_missing_assets(tmp_path: Path, capsys) -> None:
    config = _config(tmp_path)
    (tmp_path / config["evaluation"]["benchmark_bed"]).unlink()
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
    before = sorted(
        (path.relative_to(tmp_path), path.stat().st_size)
        for path in tmp_path.rglob("*")
        if path.is_file()
    )
    exit_code = main(
        ["--config", str(config_path), "--repo-root", str(tmp_path), "--json"]
    )
    after = sorted(
        (path.relative_to(tmp_path), path.stat().st_size)
        for path in tmp_path.rglob("*")
        if path.is_file()
    )
    assert exit_code == 1
    assert before == after
    assert '"status": "not_ready"' in capsys.readouterr().out
