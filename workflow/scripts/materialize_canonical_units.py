#!/usr/bin/env python3
"""Freeze ALT-specific canonical PG-F1 units without changing the panel."""
from __future__ import annotations
import argparse,csv,json,sys
from pathlib import Path
from pgf1 import PGF1Error,UNIT_FIELDS,open_text,sha256_file,write_tsv
def info(s:str,k:str,d:str="")->str:
    for x in s.split(";"):
        if x.startswith(k+"="):return x.split("=",1)[1]
    return d
def materialize(panel:Path,hidden:Path,out:Path,sha:Path,audit:Path)->dict[str,object]:
    truth={}
    with hidden.open(encoding="utf-8",newline="") as h:
        r=csv.DictReader(h,delimiter="\t")
        if r.fieldnames is None or not {"candidate_id","truth_gt"}.issubset(r.fieldnames):raise PGF1Error("hidden truth ledger lacks required fields")
        for x in r:
            if not x["candidate_id"] or x["candidate_id"] in truth:raise PGF1Error("UNIT_SET_MISMATCH")
            truth[x["candidate_id"]]=x
    seen=set(); rows=[]
    with open_text(panel,"rt") as h:
        for n,line in enumerate(h,1):
            if not line.strip() or line.startswith("#"):continue
            f=line.rstrip("\r\n").split("\t")
            if len(f)<8 or f[2] in {"","."} or f[2] in seen or f[2] not in truth:raise PGF1Error(f"UNIT_SET_MISMATCH at panel line {n}")
            seen.add(f[2]);alts=f[4].split(",")
            if not alts or any(not a or a=="." for a in alts):raise PGF1Error("invalid ALT")
            status="SCORABLE" if truth[f[2]].get("truth_scorable","true").lower() in {"true","1","yes"} else "UNSCORABLE"
            for i,alt in enumerate(alts,1):rows.append({"unit_id":f"{f[2]}:A{i}","candidate_id":f[2],"allele_id":f"{info(f[7],'PANGENOME_ALLELE_ID',f[2])}:A{i}","chrom":f[0],"pos":f[1],"end":info(f[7],"END",f[1]),"svtype":info(f[7],"SVTYPE","UNKNOWN"),"canonical_ref":f[3],"canonical_alt":alt,"truth_gt":truth[f[2]]["truth_gt"],"scorable_status":status})
    if set(truth)!=seen:raise PGF1Error("UNIT_SET_MISMATCH")
    write_tsv(out,UNIT_FIELDS,rows); digest=sha256_file(out);sha.parent.mkdir(parents=True,exist_ok=True);sha.write_text(f"{digest}  {out.name}\n",encoding="utf-8")
    result={"contract":"pgbench_canonical_units_v1","status":"valid","panel_sha256":sha256_file(panel),"hidden_truth_ledger_sha256":sha256_file(hidden),"canonical_unit_set_sha256":digest,"candidate_count":len(seen),"unit_count":len(rows),"multiallelic_candidates":sum(1 for r in rows if r["unit_id"].endswith(":A2"))};audit.parent.mkdir(parents=True,exist_ok=True);audit.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8");return result
def main(argv:list[str]|None=None)->int:
    p=argparse.ArgumentParser();p.add_argument("--canonical-panel",required=True,type=Path);p.add_argument("--hidden-truth-ledger",required=True,type=Path);p.add_argument("--output-units",required=True,type=Path);p.add_argument("--output-sha256",required=True,type=Path);p.add_argument("--audit-json",required=True,type=Path);a=p.parse_args(argv)
    try:materialize(a.canonical_panel,a.hidden_truth_ledger,a.output_units,a.output_sha256,a.audit_json)
    except (OSError,PGF1Error) as e:print(f"materialize_canonical_units: {e}",file=sys.stderr);return 2
    return 0
if __name__=="__main__":raise SystemExit(main())
