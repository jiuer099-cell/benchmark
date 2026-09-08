from __future__ import annotations

import csv
import gzip
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from materialize_evaluation_query import materialize_evaluation_query  # noqa: E402


def test_evaluation_query_is_a_deterministic_non_reference_projection(
    tmp_path: Path,
) -> None:
    all_sites = tmp_path / "all-sites.vcf.gz"
    with gzip.open(all_sites, "wt", encoding="utf-8") as handle:
        handle.write(
            "##fileformat=VCFv4.2\n"
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
            "chr1\t10\tA\tA\tAT\t0\tPASS\tPGBENCH_OUTPUT_STATUS=called\tGT\t0/1\n"
            "chr1\t20\tB\tA\tAT\t999\tPASS\tPGBENCH_OUTPUT_STATUS=called\tGT\t0/0\n"
            "chr1\t30\tC\tA\tAT\t.\tPASS\tPGBENCH_OUTPUT_STATUS=explicit_no_call\tGT\t./.\n"
            "chr1\t40\tD\tA\tAT\t.\tPASS\tPGBENCH_OUTPUT_STATUS=unsupported_representation\tGT\t./.\n"
        )
    status = tmp_path / "candidate-status.tsv"
    fieldnames = [
        "candidate_id", "tool_status", "link_status",
        "addressability_status", "raw_gt", "canonical_gt", "final_status",
    ]
    with status.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(
            [
                {"candidate_id": "A", "tool_status": "called", "link_status": "linked", "addressability_status": "addressable", "raw_gt": "0/1", "canonical_gt": "0/1", "final_status": "addressable_called"},
                {"candidate_id": "B", "tool_status": "called", "link_status": "linked", "addressability_status": "addressable", "raw_gt": "0/0", "canonical_gt": "0/0", "final_status": "addressable_called"},
                {"candidate_id": "C", "tool_status": "explicit_no_call", "link_status": "linked", "addressability_status": "addressable", "raw_gt": "./.", "canonical_gt": "./.", "final_status": "explicit_no_call"},
                {"candidate_id": "D", "tool_status": "unsupported_representation", "link_status": "not_applicable", "addressability_status": "unsupported_representation", "raw_gt": ".", "canonical_gt": "./.", "final_status": "unsupported_representation"},
            ]
        )
    output = tmp_path / "evaluation-query.vcf.gz"
    audit = materialize_evaluation_query(
        canonical_panel=all_sites,
        all_sites_vcf=all_sites,
        candidate_status_tsv=status,
        output_vcf=output,
        audit_json=tmp_path / "audit.json",
    )
    with gzip.open(output, "rt", encoding="utf-8") as handle:
        records = [line for line in handle if not line.startswith("#")]
    assert len(records) == 1
    assert records[0].split("\t")[2] == "A"
    assert records[0].split("\t")[5] == "0"  # QUAL is not filtered.
    assert audit["canonical_candidate_count"] == 4
    assert audit["non_reference_query_count"] == 1
    assert audit["reference_genotype_excluded_count"] == 1
    assert audit["unique_candidate_mapping"] is True
    assert audit["canonical_alleles_verified"] is True
    assert audit["extra_filtering_applied"] is False

    second_output = tmp_path / "second-name.vcf.gz"
    materialize_evaluation_query(
        canonical_panel=all_sites,
        all_sites_vcf=all_sites,
        candidate_status_tsv=status,
        output_vcf=second_output,
        audit_json=tmp_path / "second-audit.json",
    )
    assert output.read_bytes() == second_output.read_bytes()
