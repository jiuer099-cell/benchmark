from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))

from prepare_pgf1_small_replay import main


def test_small_replay_is_stratified_and_never_production_eligible(tmp_path: Path):
    header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
    panel = tmp_path / "panel.vcf"
    linked = tmp_path / "linked.vcf"
    all_sites = tmp_path / "all-sites.vcf.gz"
    truth = tmp_path / "truth.vcf.gz"
    hidden = tmp_path / "hidden.tsv"
    panel_rows = []
    linked_rows = []
    all_rows = []
    hidden_rows = []
    for index, gt in enumerate(("0/1", "1/1", "0/0", "."), start=1):
        candidate = f"C{index}"
        pos = 100 * index
        panel_rows.append(f"chr1\t{pos}\t{candidate}\tN\t<DEL>\t.\tPASS\tEND={pos + 50};SVTYPE=DEL\n")
        linked_rows.append(f"chr1\t{pos}\tN{index}\tN\t<DEL>\t.\tPASS\tEND={pos + 50};PANGENOME_LINKED_ID={candidate};PANGENOME_CANDIDATE_COUNT=1\tGT\t{gt}\n")
        all_rows.append(f"chr1\t{pos}\t{candidate}\tN\t<DEL>\t.\tPASS\tEND={pos + 50}\tGT\t{gt}\n")
        hidden_rows.append(f"{candidate}\t{gt if gt != '.' else '0/1'}\t1\n")
    panel.write_text(header.replace("\tFORMAT\tHG002", "") + "".join(panel_rows), encoding="utf-8")
    linked.write_text(header + "".join(linked_rows), encoding="utf-8")
    with gzip.open(all_sites, "wt", encoding="utf-8") as handle:
        handle.write(header + "".join(all_rows))
    with gzip.open(truth, "wt", encoding="utf-8") as handle:
        handle.write(header + "chr1\t100\tT1\tN\t<DEL>\t.\tPASS\tEND=450\tGT\t0/1\n")
    hidden.write_text("candidate_id\ttruth_gt\ttruth_scorable\n" + "".join(hidden_rows), encoding="utf-8")
    output = tmp_path / "out"
    args = [
        "--canonical-panel", str(panel), "--hidden-truth-ledger", str(hidden),
        "--linked-vcf", str(linked), "--all-sites-vcf", str(all_sites),
        "--truth-vcf", str(truth), "--output-dir", str(output),
        "--positive-units", "2", "--reference-units", "1", "--no-call-units", "1",
    ]
    assert main(args) == 0
    manifest = json.loads((output / "selection-manifest.json").read_text(encoding="utf-8"))
    assert manifest["selected_genotype_counts"] == {"no_call": 1, "positive": 2, "reference": 1}
    assert manifest["production_eligible"] is False
    assert manifest["leaderboard_admissible"] is False
    assert manifest["selected_unit_count"] == 4
    assert (output / "truth.vcf").read_text(encoding="utf-8").startswith("#CHROM")
