from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from materialize_all_sites import materialize  # noqa: E402


def test_projects_every_candidate_and_keeps_missing_distinct(tmp_path: Path) -> None:
    panel = tmp_path / "panel.vcf"
    panel.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        "chr1\t100\tPGSV_A\tA\t<DEL>\t.\tPASS\t"
        "SVTYPE=DEL;END=200;SVLEN=-100;PANGENOME_ALLELE_ID=A\tGT\t./.\n"
        "chr1\t300\tPGSV_B\tA\tAT\t.\tPASS\t"
        "SVTYPE=INS;END=300;SVLEN=80;PANGENOME_ALLELE_ID=B\tGT\t./.\n",
        encoding="utf-8",
    )
    query = tmp_path / "query.vcf"
    query.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        "chr1\t100\tnative1\tA\t<DEL>\t.\tPASS\t"
        "SVTYPE=DEL;END=200;SVLEN=-100;PANGENOME_LINKED_ID=A\tGT\t0/1\n",
        encoding="utf-8",
    )
    output = tmp_path / "all-sites.vcf"
    status = tmp_path / "addressability.tsv"
    audit = tmp_path / "audit.json"
    result = materialize(
        canonical_panel=panel,
        linked_query=query,
        output_vcf=output,
        addressability_tsv=status,
        audit_json=audit,
    )
    records = [
        line for line in output.read_text().splitlines() if not line.startswith("#")
    ]
    assert len(records) == 2
    assert records[0].endswith("\tGT\t0/1")
    assert "PGBENCH_OUTPUT_STATUS=missing_output" in records[1]
    assert records[1].endswith("\tGT\t./.")
    with status.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["output_status"] for row in rows] == ["called", "missing_output"]
    assert result["silent_candidate_deletion_count"] == 0
    assert json.loads(audit.read_text())["canonical_candidate_count"] == 2
