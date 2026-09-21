#!/usr/bin/env python3
"""Materialize fixed release invariants used by leaderboard admission."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import yaml
from pgbench_pgf1_contract import PG_F1_SCORE_CONTRACT_SHA256
from pgf1 import PGF1Error,json_dump,sha256_file
def main(argv:list[str]|None=None)->int:
    p=argparse.ArgumentParser();p.add_argument("--config",required=True,type=Path);p.add_argument("--panel",required=True,type=Path);p.add_argument("--units",required=True,type=Path);p.add_argument("--evaluator-profile",required=True,type=Path);p.add_argument("--output",required=True,type=Path);a=p.parse_args(argv);c=yaml.safe_load(a.config.read_text(encoding="utf-8"));
    release_id=c["benchmark_contract"].get("release_id")
    if not isinstance(release_id,str) or not release_id:raise PGF1Error("LEADERBOARD_ADMISSION_FAILED: PG-F1 release_id must be frozen")
    paths=[c["evaluation"]["truth_vcf"],c["evaluation"]["benchmark_bed"],c["reference"]["fasta"]]
    normalizer_sha=sha256_file(Path(__file__).with_name("normalize_evaluator_evidence.py"))
    data={"track_registry_id":str(c["benchmark_contract"]["track"]),"release_id":release_id,"canonical_panel_sha256":"", "canonical_unit_set_sha256":sha256_file(a.units),"truth_vcf_sha256":sha256_file(Path(paths[0])),"benchmark_bed_sha256":sha256_file(Path(paths[1])),"reference_sha256":sha256_file(Path(paths[2])),"score_contract_sha256":PG_F1_SCORE_CONTRACT_SHA256,"evaluator_bundle_sha256":sha256_file(a.evaluator_profile),"evaluator_query_contract_sha256":sha256_file(Path(__file__).with_name("materialize_pgf1_evaluator_query.py")),"normalizer_contract_sha256":normalizer_sha,"normalizer_sha256s":{"truvari":normalizer_sha,"aardvark":normalizer_sha,"vcfdist":normalizer_sha}}
    data["canonical_panel_sha256"]=sha256_file(a.panel)
    json_dump(a.output,data);return 0
if __name__=="__main__":raise SystemExit(main())
