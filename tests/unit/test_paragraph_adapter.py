from __future__ import annotations

import sys
from pathlib import Path

import pytest

PLUGIN = Path(__file__).resolve().parents[2] / "plugins" / "paragraph"
sys.path.insert(0, str(PLUGIN))

import run as paragraph_adapter  # noqa: E402

META = [
    "##fileformat=VCFv4.2\n",
    '##INFO=<ID=PANGENOME_ALLELE_ID,Number=1,Type=String,Description="x">\n',
    '##INFO=<ID=AF,Number=A,Type=Float,Description="af">\n',
]
HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
FULL_HEADER = HEADER + "\tFORMAT\tSAMPLE"

IDS = {
    ("chr1", 100): "PGSV_a1",
    ("chr1", 200): "PGSV_a2",
    ("chr1", 300): "PGSV_a3",
    ("chr2", 100): "PGSV_a4",
    ("chr3", 100): "PGSV_a5",
}


def _record(contig: str, pos: int, allele_id: str) -> str:
    return (
        f"{contig}\t{pos}\t{IDS[(contig, pos)]}\tA\t<TYPE>\t.\t.\t"
        f"SVTYPE=DEL;END={pos + 50};AF=0.1;PANGENOME_ALLELE_ID={allele_id}\n"
    )


def _native_record(contig: str, pos: int, allele_id: str, gt: str) -> str:
    return (
        f"{contig}\t{pos}\t{IDS[(contig, pos)]}\tA\t<TYPE>\t.\t.\t"
        f"SVTYPE=DEL;END={pos + 50};AF=0.1;PANGENOME_ALLELE_ID={allele_id}\tGT\t{gt}\n"
    )


def _panel_path(tmp_path: Path, lines: list[str]) -> Path:
    path = tmp_path / "candidates.vcf"
    path.write_text("".join(META) + FULL_HEADER + "\n" + "".join(lines), encoding="utf-8")
    return path


def _native_path(tmp_path: Path, name: str, records: list[str]) -> Path:
    path = tmp_path / name
    path.write_text("".join(META) + FULL_HEADER + "\n" + "".join(records), encoding="utf-8")
    return path


PANEL_LINES = [
    _record("chr1", 100, "allele-1"),
    _record("chr1", 200, "allele-2"),
    _record("chr1", 300, "allele-3"),
    _record("chr2", 100, "allele-4"),
    _record("chr3", 100, "allele-5"),
]


def test_load_candidate_panel_groups_contigs_in_first_appearance_order(tmp_path):
    panel = paragraph_adapter.load_candidate_panel(_panel_path(tmp_path, PANEL_LINES))
    assert list(panel.contigs) == ["chr1", "chr2", "chr3"]
    assert panel.contig_records["chr1"] == ("PGSV_a1", "PGSV_a2", "PGSV_a3")
    assert panel.contig_records["chr2"] == ("PGSV_a4",)
    assert panel.by_allele["allele-5"] == "PGSV_a5"
    assert panel.by_key[("chr3", "100", "A", "<TYPE>")] == "PGSV_a5"


def test_write_chunks_respects_per_chunk_budget(tmp_path):
    panel = paragraph_adapter.load_candidate_panel(_panel_path(tmp_path, PANEL_LINES))
    chunks = paragraph_adapter.write_chunks(panel, tmp_path / "chunks", max_per_chunk=2)
    # chr1: 3 records -> 2 chunks (2+1); chr2: 1 chunk; chr3: 1 chunk -> 4 total
    assert len(chunks) == 4
    body = [line for line in chunks[0].read_text(encoding="utf-8").splitlines() if not line.startswith("#")]
    assert len(body) == 2
    for chunk in chunks:
        text = chunk.read_text(encoding="utf-8")
        assert text.count("##fileformat") == 1
        assert FULL_HEADER in text


def test_write_chunks_contig_filter_selects_only_requested_contigs(tmp_path):
    panel = paragraph_adapter.load_candidate_panel(_panel_path(tmp_path, PANEL_LINES))
    chunks = paragraph_adapter.write_chunks(
        panel, tmp_path / "chunks", max_per_chunk=10, contig_filter=frozenset({"chr2"})
    )
    assert len(chunks) == 1
    body = [line for line in chunks[0].read_text(encoding="utf-8").splitlines() if not line.startswith("#")]
    assert len(body) == 1
    assert body[0].split("\t")[0] == "chr2"


def test_write_chunks_rejects_unknown_pilot_contig(tmp_path):
    panel = paragraph_adapter.load_candidate_panel(_panel_path(tmp_path, PANEL_LINES))
    with pytest.raises(RuntimeError, match="absent from candidate panel"):
        paragraph_adapter.write_chunks(
            panel, tmp_path / "chunks", max_per_chunk=10, contig_filter=frozenset({"chrY"})
        )


def test_paragraph_private_input_left_anchors_only_incompatible_records(tmp_path):
    reference = tmp_path / "reference.fa"
    # Write LF bytes directly so the fixture's FAI offsets are portable to
    # Windows, where text-mode writes otherwise translate newlines.
    reference.write_bytes(b">chr1\nACGTACGT\n")
    # ``chr1`` bases begin at byte 6; lines contain eight bases plus one LF.
    (tmp_path / "reference.fa.fai").write_text("chr1\t8\t6\t8\t9\n", encoding="utf-8")
    panel_path = _panel_path(
        tmp_path,
        [
            "chr1\t2\tPGSV_a1\tC\tTTT\t.\t.\tPANGENOME_ALLELE_ID=allele-1\n",
            "chr1\t4\tPGSV_a2\tT\tTGG\t.\t.\tPANGENOME_ALLELE_ID=allele-2\n",
        ],
    )
    panel = paragraph_adapter.load_candidate_panel(panel_path)

    prepared, changed = paragraph_adapter.paragraph_compatible_records(panel, reference)

    assert changed == 1
    assert prepared["PGSV_a1"][1:5] == ["1", "PGSV_a1", "AC", "ATTT"]
    assert prepared["PGSV_a2"][1:5] == ["4", "PGSV_a2", "T", "TGG"]
    assert panel.records["PGSV_a1"][1:5] == ["2", "PGSV_a1", "C", "TTT"]


def test_project_all_sites_merges_disjoint_chunks_and_fills_no_calls(tmp_path):
    panel = paragraph_adapter.load_candidate_panel(_panel_path(tmp_path, PANEL_LINES))
    native_a = _native_path(tmp_path, "a.vcf", [_native_record("chr1", 100, "allele-1", "0/1")])
    native_b = _native_path(
        tmp_path,
        "b.vcf",
        [
            _native_record("chr2", 100, "allele-4", "1/1"),
            _native_record("chr3", 100, "allele-5", "0/0"),
        ],
    )
    out = tmp_path / "out" / "all_sites.vcf"
    matched, no_call = paragraph_adapter.project_all_sites(
        [native_a, native_b], panel, out, "SAMPLE", pilot_contigs=["chr1"]
    )
    assert (matched, no_call) == (3, 2)
    text = out.read_text(encoding="utf-8")
    assert "##PGBENCH_Paragraph_PilotContigs=chr1\n" in text
    body = {
        line.split("\t")[0] + ":" + line.split("\t")[1]: line.split("\t")[9]
        for line in text.splitlines()
        if line and not line.startswith("#")
    }
    assert body["chr1:100"] == "0/1"
    assert body["chr1:200"] == "./."
    assert body["chr2:100"] == "1/1"
    assert body["chr3:100"] == "0/0"


def test_project_all_sites_rejects_duplicate_candidate_across_chunks(tmp_path):
    panel = paragraph_adapter.load_candidate_panel(_panel_path(tmp_path, PANEL_LINES))
    duplicate = _native_record("chr1", 100, "allele-1", "0/1")
    native = _native_path(tmp_path, "a.vcf", [duplicate, duplicate])
    with pytest.raises(RuntimeError, match="duplicates"):
        paragraph_adapter.project_all_sites([native], panel, tmp_path / "out.vcf", "SAMPLE")


def test_read_length_from_bam_streams_and_stops_early(tmp_path, monkeypatch):
    calls = {"killed": False}

    class FakeProcess:
        def __init__(self):
            self.stdout = iter(
                [
                    "q1\t99\tchr1\t100\t60\t10M\t=\t200\t300\t" + "A" * 150 + "\t*\n",
                ]
            )
            self.returncode = None

        def kill(self):
            calls["killed"] = True

        def wait(self):
            return 0

    monkeypatch.setattr(paragraph_adapter.subprocess, "Popen", lambda *a, **k: FakeProcess())
    length = paragraph_adapter.read_length_from_bam(tmp_path / "in.bam")
    assert length == 150
    assert calls["killed"] is True
