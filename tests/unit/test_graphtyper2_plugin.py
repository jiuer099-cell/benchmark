from __future__ import annotations

import gzip
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RUNNER = ROOT / "plugins" / "graphtyper2" / "run.py"
SPEC = importlib.util.spec_from_file_location("graphtyper2_runner", RUNNER)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _candidate(path: Path) -> None:
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t10\tCAND_1\tA\t<DEL>\t.\tPASS\tEND=20\n"
        "chr1\t10000001\tCAND_2\tC\t<INS>\t.\tPASS\tEND=10000001\n",
        encoding="utf-8",
    )


def test_regions_cover_candidate_chunks(tmp_path: Path) -> None:
    candidate = tmp_path / "candidate.vcf"
    regions = tmp_path / "regions.txt"
    _candidate(candidate)
    MODULE.write_regions(candidate, regions)
    assert regions.read_text(encoding="utf-8").splitlines() == [
        "chr1:1-10000000",
        "chr1:10000001-20000000",
    ]


def test_projection_preserves_all_candidates_and_no_calls(tmp_path: Path) -> None:
    candidate = tmp_path / "candidate.vcf"
    _candidate(candidate)
    generated = tmp_path / "generated.vcf.gz"
    with gzip.open(generated, "wt", encoding="utf-8") as handle:
        handle.write(
            "##fileformat=VCFv4.2\n"
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
            "chr1\t10\told\tA\t<DEL>\t99\tPASS\tEND=20\tGT\t0/1\n"
        )
    output = tmp_path / "calls.vcf"
    MODULE.project_calls(candidate, [generated], output, "HG002")
    records = [
        line.split("\t")
        for line in output.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]
    assert [record[2] for record in records] == ["CAND_1", "CAND_2"]
    assert records[0][9] == "0/1"
    assert records[1][9] == "./."
