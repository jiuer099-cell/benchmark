#!/usr/bin/env python3
"""Fail closed before production; a failed adapter never becomes PG-F1=0."""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
from pgf1 import PGF1Error,json_dump
GATES=("schema","environment","fixed_mini_panel","canonical_mapping","pgf1_evidence","real_input")
def validate(checks:Path,out:Path)->dict:
    x=json.loads(checks.read_text(encoding="utf-8")); failed=[k for k in GATES if not isinstance(x,dict) or x.get(k)!="PASS"]
    if failed:raise PGF1Error("PRODUCTION_NOT_ELIGIBLE: "+", ".join(failed))
    r={"adapter_conformance":{k:"PASS" for k in GATES},"production_eligibility":"PASS"};json_dump(out,r);return r
def main(argv:list[str]|None=None)->int:
    p=argparse.ArgumentParser();p.add_argument("--checks",required=True,type=Path);p.add_argument("--output",required=True,type=Path);a=p.parse_args(argv)
    try:validate(a.checks,a.output)
    except (OSError,json.JSONDecodeError,PGF1Error) as e:print(f"validate_adapter_conformance: {e}",file=sys.stderr);return 2
    return 0
if __name__=="__main__":raise SystemExit(main())
