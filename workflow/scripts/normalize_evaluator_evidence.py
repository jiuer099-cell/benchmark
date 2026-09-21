#!/usr/bin/env python3
"""Normalize direct evaluator records to PG-F1 per-unit binary evidence."""
from __future__ import annotations
import argparse,hashlib,sys
from pathlib import Path
from pgf1 import EVIDENCE_FIELDS,PGF1Error,canonical_alt_dosage,load_units,open_text,read_tsv,sha256_file,write_tsv
VERSION="pgf1-evaluator-normalizer-v1"
def all_site_gts(path:Path)->dict[str,str]:
    calls={}
    with open_text(path,"rt") as h:
        for line in h:
            if not line.strip() or line.startswith("#"):continue
            f=line.rstrip("\r\n").split("\t")
            if len(f)<10 or f[2] in {"","."} or f[2] in calls:raise PGF1Error("CANONICAL_MAPPING_FAILURE")
            keys,values=f[8].split(":"),f[9].split(":");calls[f[2]]=values[keys.index("GT")] if "GT" in keys and keys.index("GT")<len(values) else "./."
    return calls
def normalize(*,units:Path,all_sites:Path,raw:Path,evaluator:str,release_id:str,run_id:str,sample_id:str,tool_id:str,out:Path)->dict[str,int]:
    unit=load_units(units);required={"result_id","evaluator","evaluator_accept","evaluator_resolved","mapping_status","truth_event_id"};by={}
    for x in read_tsv(raw,required):
        c=x["result_id"]
        if x["evaluator"]!=evaluator or not c or c in by:raise PGF1Error("NORMALIZER_RECONCILIATION_FAILURE")
        if x["mapping_status"] in {"partial_complex","ambiguous","many_to_many","one_to_many","many_to_one"}:raise PGF1Error("AMBIGUOUS_COMPLEX_MAPPING")
        by[c]=x
    calls=all_site_gts(all_sites); candidates={x["candidate_id"] for x in unit.values()}
    if set(calls)!=candidates or not set(by).issubset(candidates):raise PGF1Error("UNIT_SET_MISMATCH")
    raw_sha=sha256_file(raw); script_sha=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(); rows=[]
    for uid,x in sorted(unit.items()):
        r=by.get(x["candidate_id"]); query_positive=canonical_alt_dosage(calls[x["candidate_id"]]) not in {None,0}
        # A reference/no-call candidate has no submitted evaluator event. It is
        # still an explicit evidence row, not a missing row. Conversely a
        # positive query without a raw direct decision is a formal failure.
        if r is None:
            if query_positive:raise PGF1Error("CONSENSUS_EVIDENCE_INCOMPLETE")
            r={"truth_event_id":"","result_id":"","evaluator_accept":"0","mapping_status":"core_absence"};event=False;resolved=False
        else:
            event=r.get("event_eligible","0")=="1"; resolved=r["evaluator_resolved"]=="1"
            if query_positive and (not event or not resolved):raise PGF1Error("CONSENSUS_EVIDENCE_INCOMPLETE")
        rel="AM" if resolved and r["evaluator_accept"]=="1" else "NONE"
        rows.append({"release_id":release_id,"run_id":run_id,"sample_id":sample_id,"tool_id":tool_id,"unit_id":uid,"evaluator":evaluator,"truth_record_id":r.get("truth_event_id",""),"query_record_id":r.get("result_id",""),"relation":rel,"binary_vote":int(rel=="AM"),"mapping_cardinality":"one_to_one" if event else "not_applicable","mapping_status":"resolved" if resolved else "not_applicable","match_group_id":r.get("truth_event_id",""),"superlocus_id":"","native_decision":r.get("evaluator_accept","0"),"native_reason":r.get("mapping_status",""),"raw_evidence_sha256":raw_sha,"normalizer_version":VERSION,"normalizer_sha256":script_sha})
    write_tsv(out,EVIDENCE_FIELDS,rows);return {"units":len(rows),"raw_rows":len(by)}
def main(argv:list[str]|None=None)->int:
    p=argparse.ArgumentParser();p.add_argument("--units",required=True,type=Path);p.add_argument("--all-sites",required=True,type=Path);p.add_argument("--raw-ledger",required=True,type=Path);p.add_argument("--evaluator",required=True,choices=("truvari","aardvark","vcfdist"));p.add_argument("--release-id",required=True);p.add_argument("--run-id",required=True);p.add_argument("--sample-id",required=True);p.add_argument("--tool-id",required=True);p.add_argument("--output",required=True,type=Path);a=p.parse_args(argv)
    try:normalize(units=a.units,all_sites=a.all_sites,raw=a.raw_ledger,evaluator=a.evaluator,release_id=a.release_id,run_id=a.run_id,sample_id=a.sample_id,tool_id=a.tool_id,out=a.output)
    except (OSError,PGF1Error) as e:print(f"normalize_evaluator_evidence: {e}",file=sys.stderr);return 2
    return 0
if __name__=="__main__":raise SystemExit(main())
