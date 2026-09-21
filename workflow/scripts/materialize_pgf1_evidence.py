#!/usr/bin/env python3
"""Materialize the auditable PG-F1 ledger and sole formal score."""
from __future__ import annotations
import argparse
from pathlib import Path
from pgbench_pgf1_contract import PG_F1_SCORE_CONTRACT,PG_F1_SCORE_CONTRACT_SHA256
from pgf1 import EVALUATORS,PGF1Error,PGF1_LEDGER_FIELDS,build_ledger,json_dump,load_evidence,load_units,open_text,sha256_file,write_tsv
def gt(f:list[str])->str:
    if len(f)<10:return "./."
    k,v=f[8].split(":"),f[9].split(":");return v[k.index("GT")] if "GT" in k and k.index("GT")<len(v) and v[k.index("GT")] else "./."
def materialize(*,units:Path,all_sites:Path,evidence:dict[str,Path],ledger:Path,audit:Path,score:Path,run_id:str,sample:str,tool:str,track:str)->dict[str,object]:
    u=load_units(units);calls={}
    with open_text(all_sites,"rt") as h:
        for line in h:
            if line.strip() and not line.startswith("#"):
                f=line.rstrip("\r\n").split("\t")
                if len(f)<10 or f[2] in {"","."} or f[2] in calls:raise PGF1Error("CANONICAL_MAPPING_FAILURE")
                calls[f[2]]=gt(f)
    if set(calls)!={x["candidate_id"] for x in u.values()}:raise PGF1Error("UNIT_SET_MISMATCH")
    e={x:load_evidence(evidence[x],x,set(u)) for x in EVALUATORS};rows,metrics=build_ledger(units=u,query_gt={i:calls[x["candidate_id"]] for i,x in u.items()},evidence=e);write_tsv(ledger,PGF1_LEDGER_FIELDS,rows)
    au={"contract":"pgbench_pgf1_evidence_v1","status":"valid","canonical_unit_set_sha256":sha256_file(units),"all_sites_sha256":sha256_file(all_sites),"evidence_sha256":{x:sha256_file(evidence[x]) for x in EVALUATORS},"pgf1_evidence_sha256":sha256_file(ledger),"score_contract_sha256":PG_F1_SCORE_CONTRACT_SHA256,"normalizer_contract":"pgf1-evaluator-normalizer-v1","scored_units":metrics["scored_units"]};json_dump(audit,au)
    result={"score_name":"PG-F1","score_profile":"pgbench_pgf1_v1","score_status":"valid","benchmark_score":metrics["pg_f1"],"pg_f1":metrics["pg_f1"],"precision":metrics["precision"],"recall":metrics["recall"],"tp":metrics["tp"],"fp":metrics["fp"],"fn":metrics["fn"],"tuple":{"run_id":run_id,"sample":sample,"tool":tool,"track":track},"score_contract":PG_F1_SCORE_CONTRACT,"score_contract_sha256":PG_F1_SCORE_CONTRACT_SHA256,"canonical_unit_set_sha256":au["canonical_unit_set_sha256"],"pgf1_evidence_sha256":au["pgf1_evidence_sha256"],"raw_evaluator_metrics_diagnostic_only":True};json_dump(score,result);return result
def main(argv:list[str]|None=None)->int:
    p=argparse.ArgumentParser();p.add_argument("--units",required=True,type=Path);p.add_argument("--all-sites",required=True,type=Path);[p.add_argument(f"--{x}-evidence",required=True,type=Path) for x in EVALUATORS];p.add_argument("--ledger-output",required=True,type=Path);p.add_argument("--audit-output",required=True,type=Path);p.add_argument("--score-output",required=True,type=Path);p.add_argument("--run-id",required=True);p.add_argument("--sample-id",required=True);p.add_argument("--tool-id",required=True);p.add_argument("--track",required=True);a=p.parse_args(argv)
    try:materialize(units=a.units,all_sites=a.all_sites,evidence={x:getattr(a,x+"_evidence") for x in EVALUATORS},ledger=a.ledger_output,audit=a.audit_output,score=a.score_output,run_id=a.run_id,sample=a.sample_id,tool=a.tool_id,track=a.track)
    except (OSError,PGF1Error,ValueError) as e:print(f"materialize_pgf1_evidence: {e}",file=__import__('sys').stderr);return 2
    return 0
if __name__=="__main__":raise SystemExit(main())
