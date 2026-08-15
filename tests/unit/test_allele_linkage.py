from __future__ import annotations

import csv
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from link_pangenome_alleles import (  # noqa: E402
    Allele,
    Call,
    link_call,
    link_vcf,
    load_ledger,
)


def _allele(
    allele_id: str,
    *,
    pos: int,
    end: int,
    alt: str = "<DEL>",
    svtype: str = "DEL",
    svlen: int = -10,
) -> Allele:
    return Allele(
        allele_id=allele_id,
        chrom="1",
        pos=pos,
        end=end,
        svtype=svtype,
        svlen=svlen,
        ref="A",
        alt=alt,
        af="0.01",
        graph_class="simple_biallelic",
    )


def _call(
    *,
    pos: int,
    end: int,
    alt: str = "<DEL>",
    svtype: str = "DEL",
    svlen: int = -10,
    claimed_id: str | None = None,
) -> Call:
    return Call(
        canon_id="CANON_1",
        chrom="1",
        pos=pos,
        end=end,
        svtype=svtype,
        svlen=svlen,
        ref="A",
        alt=alt,
        claimed_id=claimed_id,
    )


def test_link_statuses_cover_exact_equivalent_and_out_of_panel() -> None:
    allele = _allele("PGSV_A", pos=100, end=110)
    by_id = {allele.allele_id: allele}
    assert link_call(_call(pos=100, end=110), by_id, [allele]).status == (
        "in_panel_exact"
    )
    assert link_call(_call(pos=105, end=115), by_id, [allele]).status == (
        "in_panel_equivalent"
    )
    assert link_call(_call(pos=1000, end=1010), by_id, [allele]).status == (
        "out_of_panel"
    )


def test_link_statuses_cover_invalid_wrong_ambiguous_and_unresolved() -> None:
    first = _allele("PGSV_A", pos=100, end=110)
    second = _allele("PGSV_B", pos=105, end=115)
    by_id = {first.allele_id: first, second.allele_id: second}
    assert (
        link_call(
            _call(pos=100, end=110, claimed_id="MISSING"), by_id, [first, second]
        ).status
        == "invalid_allele_id"
    )
    assert (
        link_call(
            _call(pos=105, end=115, claimed_id="PGSV_A"), by_id, [first, second]
        ).status
        == "wrong_link"
    )
    assert link_call(_call(pos=103, end=113), by_id, [first, second]).status == (
        "ambiguous"
    )
    bnd = _call(
        pos=200,
        end=200,
        alt="<BND>",
        svtype="BND",
        svlen=0,
    )
    assert link_call(bnd, by_id, [first, second]).status == "unresolved"


def test_link_vcf_writes_auditable_ledger(tmp_path: Path) -> None:
    ledger = tmp_path / "alleles.tsv"
    ledger.write_text(
        "PANGENOME_ALLELE_ID\tCHROM\tPOS\tEND\tSVTYPE\tSVLEN\tREF\tALT\tAF\t"
        "GRAPH_COMPLEXITY_CLASS\tSOURCE_RECORD_ID\n"
        "PGSV_A\t1\t100\t110\tDEL\t-10\tA\t<DEL>\t0.01\t"
        "simple_biallelic\ta\n",
        encoding="utf-8",
    )
    canonical = tmp_path / "canonical.vcf"
    canonical.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        "1\t100\tCANON_1\tA\t<DEL>\t.\tPASS\t"
        "CANON_ID=CANON_1;SVTYPE=DEL;END=110;SVLEN=-10\tGT\t0/1\n",
        encoding="utf-8",
    )
    linked = tmp_path / "linked.vcf"
    links = tmp_path / "links.tsv"
    assert link_vcf(canonical, ledger, linked, links) == 1
    assert "PANGENOME_LINK_STATUS=in_panel_exact" in linked.read_text()
    with links.open(encoding="utf-8") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["PANGENOME_ALLELE_ID"] == "PGSV_A"
    assert row["LINK_STATUS"] == "in_panel_exact"


def test_link_vcf_declares_retained_info_fields(tmp_path: Path) -> None:
    ledger = tmp_path / "alleles.tsv"
    ledger.write_text(
        "PANGENOME_ALLELE_ID\tCHROM\tPOS\tEND\tSVTYPE\tSVLEN\tREF\tALT\tAF\t"
        "GRAPH_COMPLEXITY_CLASS\tSOURCE_RECORD_ID\n"
        "PGSV_A\t1\t100\t110\tDEL\t-10\tA\t<DEL>\t0.01\t"
        "simple_biallelic\ta\n",
        encoding="utf-8",
    )
    canonical = tmp_path / "canonical.vcf"
    canonical.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        "1\t100\tCANON_1\tA\t<DEL>\t.\tPASS\t"
        "PANGENOME_ALLELE_ID=PGSV_A;CONFLICT=source;UNDECLARED=value;"
        "SVTYPE=DEL;END=110;SVLEN=-10\tGT\t0/1\n",
        encoding="utf-8",
    )
    linked = tmp_path / "linked.vcf"
    link_vcf(canonical, ledger, linked, tmp_path / "links.tsv")
    text = linked.read_text(encoding="utf-8")
    assert "##INFO=<ID=PANGENOME_ALLELE_ID,Number=1,Type=String" in text
    assert "##INFO=<ID=CONFLICT,Number=.,Type=String" in text
    assert "##INFO=<ID=UNDECLARED,Number=.,Type=String" in text
    assert text.index("##INFO=<ID=UNDECLARED") < text.index("#CHROM")


def test_relevant_ledger_scan_accepts_large_irrelevant_alleles(
    tmp_path: Path,
) -> None:
    ledger = tmp_path / "alleles.tsv"
    ledger.write_text(
        "PANGENOME_ALLELE_ID\tCHROM\tPOS\tEND\tSVTYPE\tSVLEN\tREF\tALT\tAF\t"
        "GRAPH_COMPLEXITY_CLASS\tSOURCE_RECORD_ID\n"
        "PGSV_HUGE\t2\t1\t200002\tINS\t200001\tA\t"
        + ("T" * 200_000)
        + "\t0.01\tcomplex\thuge\n"
        "PGSV_A\t1\t100\t110\tDEL\t-10\tA\t<DEL>\t0.01\t"
        "simple_biallelic\ta\n",
        encoding="utf-8",
    )
    call = _call(pos=100, end=110, claimed_id="PGSV_A")
    by_id, alleles = load_ledger(ledger, relevant_calls=[call])
    assert set(by_id) == {"PGSV_A"}
    assert [allele.allele_id for allele in alleles] == ["PGSV_A"]
