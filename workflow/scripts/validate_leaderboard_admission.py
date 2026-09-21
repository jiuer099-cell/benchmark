#!/usr/bin/env python3
"""Fail-closed admission record for a PG-F1 leaderboard release."""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
from pgbench_pgf1_contract import PG_F1_SCORE_CONTRACT_SHA256
from pgf1 import PGF1Error,json_dump,sha256_file
REQUIRED=("track_registry_id","release_id","canonical_panel_sha256","canonical_unit_set_sha256","truth_vcf_sha256","benchmark_bed_sha256","reference_sha256","score_contract_sha256","evaluator_bundle_sha256","evaluator_query_contract_sha256","normalizer_contract_sha256","normalizer_sha256s")
STRING_INVARIANTS=tuple(key for key in REQUIRED if key!="normalizer_sha256s")
def load(p:Path)->dict:
    x=json.loads(p.read_text(encoding="utf-8"))
    if not isinstance(x,dict):raise PGF1Error("LEADERBOARD_ADMISSION_FAILED")
    return x
def validate(score:Path,invariants:Path,out:Path)->dict:
    s,i=load(score),load(invariants); missing=[k for k in STRING_INVARIANTS if not isinstance(i.get(k),str) or not i[k]]
    if missing or not isinstance(i.get("normalizer_sha256s"),dict) or set(i["normalizer_sha256s"])!={"truvari","aardvark","vcfdist"} or any(not isinstance(v,str) or len(v)!=64 for v in i["normalizer_sha256s"].values()):raise PGF1Error("LEADERBOARD_ADMISSION_FAILED: missing/invalid invariant")
    if s.get("score_status")!="valid" or s.get("score_name")!="PG-F1" or s.get("score_contract_sha256")!=PG_F1_SCORE_CONTRACT_SHA256 or i["score_contract_sha256"]!=PG_F1_SCORE_CONTRACT_SHA256 or s.get("canonical_unit_set_sha256")!=i["canonical_unit_set_sha256"]:raise PGF1Error("LEADERBOARD_ADMISSION_FAILED: invariant mismatch")
    r={"status":"PASS",**{k:i[k] for k in REQUIRED},"score_sha256":sha256_file(score)};json_dump(out,r);return r
def main(argv:list[str]|None=None)->int:
    p=argparse.ArgumentParser();p.add_argument("--score",required=True,type=Path);p.add_argument("--invariants",required=True,type=Path);p.add_argument("--output",required=True,type=Path);a=p.parse_args(argv)
    try:validate(a.score,a.invariants,a.output)
    except (OSError,json.JSONDecodeError,PGF1Error) as e:print(f"validate_leaderboard_admission: {e}",file=sys.stderr);return 2
    return 0
if __name__=="__main__":raise SystemExit(main())
