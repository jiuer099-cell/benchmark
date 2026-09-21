"""Canonical-unit PG-F1 scoring primitives; no caller-specific code."""
from __future__ import annotations
import csv, gzip, hashlib, io, json, re
from collections import Counter
from collections.abc import Iterable, Mapping
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator, TextIO
EVALUATORS=("truvari","aardvark","vcfdist"); RELATIONS=frozenset({"NONE","LM","AM"})
FAILURE_STATUSES=frozenset({"UNSUPPORTED_FOR_SELECTED_TRACK","CONFIG_REJECTED","ENVIRONMENT_FAILURE","TOOL_RUNTIME_FAILURE","INVALID_NATIVE_OUTPUT","CANONICAL_MAPPING_FAILURE","AMBIGUOUS_COMPLEX_MAPPING","EVALUATOR_TRACE_FAILURE","CONSENSUS_EVIDENCE_INCOMPLETE","NORMALIZER_RECONCILIATION_FAILURE","UNIT_SET_MISMATCH","PGF1_CONTRACT_VIOLATION","PROVENANCE_FAILURE","LEADERBOARD_ADMISSION_FAILED","PRODUCTION_NOT_ELIGIBLE"})
class PGF1Error(ValueError): pass
@contextmanager
def open_text(path:Path,mode:str)->Iterator[TextIO]:
    if path.suffix!=".gz":
        with path.open(mode,encoding="utf-8",newline="") as h: yield h
    elif "w" in mode:
        with path.open("wb") as raw:
            with gzip.GzipFile(filename="",fileobj=raw,mode="wb",mtime=0) as z:
                with io.TextIOWrapper(z,encoding="utf-8",newline="") as h: yield h
    else:
        with gzip.open(path,mode,encoding="utf-8",newline="") as h: yield h
def sha256_file(path:Path)->str:
    d=hashlib.sha256()
    with path.open("rb") as h:
        for b in iter(lambda:h.read(1048576),b""): d.update(b)
    return d.hexdigest()
def canonical_alt_dosage(gt:str|None,alt_index:int=1)->int|None:
    if gt in {None,"",".","./.",".|."}: return None
    p=gt.replace("|","/").split("/")
    return sum(int(x)==alt_index for x in p) if len(p)==2 and all(x.isdigit() for x in p) else None
def is_positive(x:int|None)->bool: return x is not None and x>0
def exact_unphased_gt(a:int|None,b:int|None)->bool: return a is not None and b is not None and a==b
def vote_consensus(v:Mapping[str,int|None])->bool:
    missing=[e for e in EVALUATORS if v.get(e) not in {0,1}]
    if missing: raise PGF1Error("CONSENSUS_EVIDENCE_INCOMPLETE: "+", ".join(missing))
    return sum(int(v[e]) for e in EVALUATORS)>=2
def accounting(*,truth:int|None,query:int|None,consensus:bool)->tuple[int,int,int]:
    if is_positive(truth):
        if not is_positive(query): return 0,0,1
        return (1,0,0) if consensus and exact_unphased_gt(truth,query) else (0,1,1)
    return (0,1,0) if is_positive(query) else (0,0,0)
def f1_counts(tp:int,fp:int,fn:int)->dict[str,float|int]:
    p=tp/(tp+fp) if tp+fp else 0.; r=tp/(tp+fn) if tp+fn else 0.; f=2*p*r/(p+r) if p+r else 0.
    return {"tp":tp,"fp":fp,"fn":fn,"precision":p,"recall":r,"pg_f1":f*100}
UNIT_FIELDS=("unit_id","candidate_id","allele_id","chrom","pos","end","svtype","canonical_ref","canonical_alt","truth_gt","scorable_status")
EVIDENCE_FIELDS=("release_id","run_id","sample_id","tool_id","unit_id","evaluator","truth_record_id","query_record_id","relation","binary_vote","mapping_cardinality","mapping_status","match_group_id","superlocus_id","native_decision","native_reason","raw_evidence_sha256","normalizer_version","normalizer_sha256")
PGF1_LEDGER_FIELDS=("unit_id","truth_gt_raw","query_gt_raw","truth_gt_canonical","query_gt_canonical","truth_positive","query_positive","query_no_call","truvari_vote","aardvark_vote","vcfdist_vote","vote_sum","allele_consensus","exact_unphased_gt","count_tp","count_fp","count_fn","mapping_status")
def read_tsv(path:Path,required:Iterable[str])->list[dict[str,str]]:
    with open_text(path,"rt") as h:
        r=csv.DictReader(h,delimiter="\t")
        if r.fieldnames is None or not set(required).issubset(r.fieldnames): raise PGF1Error(f"invalid TSV header: {path}")
        return [dict(x) for x in r]
def write_tsv(path:Path,fields:Iterable[str],rows:Iterable[Mapping[str,Any]])->None:
    path.parent.mkdir(parents=True,exist_ok=True)
    with open_text(path,"wt") as h:
        w=csv.DictWriter(h,fieldnames=list(fields),delimiter="\t",lineterminator="\n"); w.writeheader()
        for x in rows:w.writerow({k:x.get(k,"") for k in w.fieldnames})
def load_units(path:Path)->dict[str,dict[str,str]]:
    out={}
    for x in read_tsv(path,UNIT_FIELDS):
        if not x["unit_id"] or x["unit_id"] in out or x["scorable_status"] not in {"SCORABLE","UNSCORABLE"}:raise PGF1Error("UNIT_SET_MISMATCH")
        out[x["unit_id"]]=x
    if not out:raise PGF1Error("UNIT_SET_MISMATCH: canonical unit set is empty")
    return out
def load_evidence(path:Path,evaluator:str,expected:set[str])->dict[str,dict[str,str]]:
    out={}
    for x in read_tsv(path,EVIDENCE_FIELDS):
        u=x["unit_id"]
        if x["evaluator"]!=evaluator or u in out or x["relation"] not in RELATIONS or x["binary_vote"] not in {"0","1"}:raise PGF1Error("NORMALIZER_RECONCILIATION_FAILURE")
        if int(x["binary_vote"]) != (x["relation"]=="AM") or x["mapping_cardinality"] not in {"one_to_one","not_applicable"}:raise PGF1Error("AMBIGUOUS_COMPLEX_MAPPING")
        if len(x["raw_evidence_sha256"])!=64:raise PGF1Error("EVALUATOR_TRACE_FAILURE")
        out[u]=x
    if set(out)!=expected:raise PGF1Error("CONSENSUS_EVIDENCE_INCOMPLETE")
    return out
def build_ledger(*,units:Mapping[str,Mapping[str,str]],query_gt:Mapping[str,str],evidence:Mapping[str,Mapping[str,Mapping[str,str]]])->tuple[list[dict[str,Any]],dict[str,Any]]:
    rows=[]; total=Counter()
    for u,x in sorted(units.items()):
        if x["scorable_status"]!="SCORABLE":continue
        match=re.search(r":A([1-9][0-9]*)$",u)
        if not match: raise PGF1Error("UNIT_SET_MISMATCH: invalid ALT-specific unit ID")
        alt_index=int(match.group(1)); t=canonical_alt_dosage(x["truth_gt"],alt_index); q=canonical_alt_dosage(query_gt.get(u,"./."),alt_index); votes={e:int(evidence[e][u]["binary_vote"]) for e in EVALUATORS}; c=vote_consensus(votes);tp,fp,fn=accounting(truth=t,query=q,consensus=c);total.update(tp=tp,fp=fp,fn=fn)
        rows.append({"unit_id":u,"truth_gt_raw":x["truth_gt"],"query_gt_raw":query_gt.get(u,"./."),"truth_gt_canonical":"NA" if t is None else t,"query_gt_canonical":"NA" if q is None else q,"truth_positive":str(is_positive(t)).lower(),"query_positive":str(is_positive(q)).lower(),"query_no_call":str(q is None).lower(),**{f"{e}_vote":votes[e] for e in EVALUATORS},"vote_sum":sum(votes.values()),"allele_consensus":str(c).lower(),"exact_unphased_gt":str(exact_unphased_gt(t,q)).lower(),"count_tp":tp,"count_fp":fp,"count_fn":fn,"mapping_status":"resolved"})
    score=f1_counts(total["tp"],total["fp"],total["fn"]);score.update({"scored_units":len(rows),"score_status":"valid"});return rows,score
def json_dump(path:Path,payload:Mapping[str,Any])->None:
    path.parent.mkdir(parents=True,exist_ok=True); tmp=path.with_name("."+path.name+".tmp");tmp.write_text(json.dumps(payload,indent=2,sort_keys=True,allow_nan=False)+"\n",encoding="utf-8");tmp.replace(path)
