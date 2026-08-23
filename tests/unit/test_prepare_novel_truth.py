from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from prepare_novel_truth import (  # noqa: E402
    IndexedReference,
    _walk_parts,
    graph_alleles,
    load_graph_bubbles,
    load_graph_segments,
    materialize_plain_novel_truth,
)
from sv_matching import load_evaluator_profile  # noqa: E402


def _reference(tmp_path: Path) -> tuple[Path, Path]:
    fasta = tmp_path / "reference.fa"
    fasta.write_bytes(b">chr1\n" + b"A" * 400 + b"\n")
    fai = tmp_path / "reference.fa.fai"
    fai.write_text("chr1\t400\t6\t400\t401\n", encoding="utf-8")
    return fasta, fai


def test_graph_exclusion_builds_a_nonleaking_novel_truth_fixture(
    tmp_path: Path,
) -> None:
    fasta, fai = _reference(tmp_path)
    insertion = "ACGT" * 15
    reverse_source = "A" * 20 + "C" * 20 + "G" * 20
    reverse_insertion = "C" * 20 + "G" * 20 + "T" * 20
    gfa = tmp_path / "graph.gfa"
    gfa.write_text(
        "H\tVN:Z:1.0\n"
        f"S\tins\t{insertion}\n"
        f"S\trev\t{reverse_source}\n",
        encoding="utf-8",
    )
    calls = tmp_path / "graph.variation.calls.bed"
    calls.write_text(
        "#CHROM\tSTART\tEND\tINFO\tFORMAT\n"
        "chr1\t20\t20\tNS=2;NA=3;ALEN=0,60,60;"
        "AWALK=*,>ins,<rev\tGT:CSTRAND\t0:+\n"
        "chr1\t100\t160\tNS=2;NA=2;ALEN=60,0;"
        "AWALK=>unused,*\tGT:CSTRAND\t0:+\n",
        encoding="utf-8",
    )
    truth = tmp_path / "truth.vcf"
    truth.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        f"chr1\t20\tKNOWN_INS\tA\tA{insertion}\t.\tPASS\t"
        "SVTYPE=INS;END=20;SVLEN=60\tGT\t0/1\n"
        f"chr1\t20\tKNOWN_REV\tA\tA{reverse_insertion}\t.\tPASS\t"
        "SVTYPE=INS;END=20;SVLEN=60\tGT\t0/1\n"
        f"chr1\t100\tKNOWN_DEL\t{'A' * 61}\tA\t.\tPASS\t"
        "SVTYPE=DEL;END=160;SVLEN=-60\tGT\t0/1\n"
        f"chr1\t250\tNOVEL_INS\tA\tA{'G' * 60}\t.\tPASS\t"
        "SVTYPE=INS;END=250;SVLEN=60\tGT\t0/1\n",
        encoding="utf-8",
    )
    evaluator = load_evaluator_profile(ROOT / "config" / "evaluator_profile.yaml")
    bubbles = load_graph_bubbles(calls)
    required = {
        name
        for bubble in bubbles
        for allele_index, walk in enumerate(bubble.allele_walks)
        if allele_index != bubble.reference_allele
        for _, name in _walk_parts(walk)
    }
    segments = load_graph_segments(gfa, required)
    records = graph_alleles(
        bubbles=bubbles,
        segments=segments,
        reference=IndexedReference(fasta, fai),
        evaluator_profile=evaluator,
    )
    output = tmp_path / "novel.vcf"
    ledger = tmp_path / "exclusion.tsv"

    counts = materialize_plain_novel_truth(
        truth_vcf=truth,
        graph_records=records,
        evaluator_profile=evaluator,
        output_vcf=output,
        exclusion_ledger=ledger,
    )

    assert counts == {
        "source_truth_records": 4,
        "graph_eligible_alleles": 3,
        "excluded_graph_represented": 3,
        "novel_truth_records": 1,
        "ambiguous_truth_matches": 0,
    }
    text = output.read_text(encoding="utf-8")
    assert "NOVEL_INS" in text
    assert "KNOWN_INS" not in text
    assert "KNOWN_REV" not in text
    assert "KNOWN_DEL" not in text
    ledger_text = ledger.read_text(encoding="utf-8")
    assert ledger_text.count("graph_represented") == 3
    assert "NOVEL_INS\tnovel" in ledger_text
