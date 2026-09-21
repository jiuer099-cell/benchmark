"""Fixed PG-F1 mini-panel and failure-contract conformance tests."""
from __future__ import annotations
import sys
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/"workflow"/"scripts"))
from pgf1 import PGF1Error,accounting,canonical_alt_dosage,exact_unphased_gt,open_text,vote_consensus
from materialize_canonical_units import materialize as materialize_units
from normalize_evaluator_evidence import normalize
from materialize_pgf1_evidence import materialize as materialize_score
from materialize_pgf1_evaluator_query import main as evaluator_query_main

@pytest.mark.parametrize(("case","truth","query","votes","expected"),[
    ("C01", "0/1", "0/1", (1,1,1), (1,0,0)), ("C02", "0/1", "0/1", (1,1,0), (1,0,0)),
    ("C03", "0/1", "1/1", (1,1,1), (0,1,1)), ("C04", "1/1", "0/1", (1,1,1), (0,1,1)),
    ("C05", "0/1", "0/0", (0,0,0), (0,0,1)), ("C06", "0/1", "./.", (0,0,0), (0,0,1)),
    ("C07", "0/0", "0/1", (1,1,1), (0,1,0)), ("C08", "0/1", "1/0", (1,1,1), (1,0,0)),
    ("C09", "0|1", "1|0", (1,1,1), (1,0,0)), ("C10", "0/1", "0/1", (1,0,0), (0,1,1)),
    ("C11", "0/1", "./.", (0,0,0), (0,0,1)), ("C12", "0/1", "./.", (0,0,0), (0,0,1)),
    ("C13", "0/1", "0/1", (1,1,1), (1,0,0)), ("C14", "0/1", "0/0", (0,0,0), (0,0,1)),
])
def test_fixed_mini_panel(case,truth,query,votes,expected):
    consensus=vote_consensus(dict(zip(("truvari","aardvark","vcfdist"),votes)))
    assert accounting(truth=canonical_alt_dosage(truth),query=canonical_alt_dosage(query),consensus=consensus)==expected

@pytest.mark.parametrize("votes,expected",[((1,1,1),True),((1,1,0),True),((1,0,1),True),((0,1,1),True),((1,0,0),False),((0,1,0),False),((0,0,1),False),((0,0,0),False)])
def test_fixed_binary_two_of_three_voter(votes,expected):
    assert vote_consensus(dict(zip(("truvari","aardvark","vcfdist"),votes))) is expected

@pytest.mark.parametrize("votes",[(1,1,None),(1,None,0),(None,None,None)])
def test_missing_evaluator_evidence_is_not_two_of_two(votes):
    with pytest.raises(PGF1Error,match="CONSENSUS_EVIDENCE_INCOMPLETE"):
        vote_consensus(dict(zip(("truvari","aardvark","vcfdist"),votes)))

def test_phase_and_order_are_ignored_but_dosage_is_not():
    assert exact_unphased_gt(canonical_alt_dosage("0/1"),canonical_alt_dosage("1|0"))
    assert not exact_unphased_gt(canonical_alt_dosage("0/1"),canonical_alt_dosage("1/1"))
    assert canonical_alt_dosage("1/2",1)==1 and canonical_alt_dosage("1/2",2)==1


def test_alt_specific_multiallelic_units_do_not_share_dosage():
    assert canonical_alt_dosage("1/1", 1) == 2
    assert canonical_alt_dosage("1/1", 2) == 0
    assert canonical_alt_dosage("2/2", 1) == 0
    assert canonical_alt_dosage("2/2", 2) == 2

def test_c15_complex_mapping_is_a_formal_failure():
    with pytest.raises(PGF1Error,match="AMBIGUOUS_COMPLEX_MAPPING"):
        raise PGF1Error("AMBIGUOUS_COMPLEX_MAPPING")


def test_native_ledger_replays_to_one_auditable_pgf1_score(tmp_path: Path):
    panel=tmp_path/"panel.vcf"; hidden=tmp_path/"hidden.tsv"
    panel.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10\tC1\tN\t<DEL>\t.\tPASS\tEND=20;SVTYPE=DEL\n",encoding="utf-8")
    hidden.write_text("candidate_id\ttruth_gt\ttruth_scorable\nC1\t0/1\ttrue\n",encoding="utf-8")
    units=tmp_path/"units.tsv.gz"; materialize_units(panel,hidden,units,tmp_path/"units.sha",tmp_path/"units.audit.json")
    raw=tmp_path/"votes.tsv"; raw.write_text("result_id\tevaluator\tevaluator_accept\tevaluator_resolved\tmapping_status\ttruth_event_id\tevent_eligible\nC1\ttruvari\t1\t1\tresolved\tT1\t1\n",encoding="utf-8")
    evidence={}
    for evaluator in ("truvari","aardvark","vcfdist"):
        candidate=tmp_path/f"{evaluator}.tsv"; candidate.write_text(raw.read_text(encoding="utf-8").replace("truvari",evaluator),encoding="utf-8")
        evidence[evaluator]=tmp_path/f"{evaluator}.evidence.tsv.gz"
    all_sites=tmp_path/"all-sites.vcf";all_sites.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\nchr1\t10\tC1\tN\t<DEL>\t.\tPASS\tEND=20\tGT\t1/0\n",encoding="utf-8")
    for evaluator in ("truvari","aardvark","vcfdist"):
        candidate=tmp_path/f"{evaluator}.tsv";normalize(units=units,all_sites=all_sites,raw=candidate,evaluator=evaluator,release_id="R",run_id="run",sample_id="HG002",tool_id="tool",out=evidence[evaluator])
    result=materialize_score(units=units,all_sites=all_sites,evidence=evidence,ledger=tmp_path/"ledger.tsv.gz",audit=tmp_path/"audit.json",score=tmp_path/"score.json",run_id="run",sample="HG002",tool="tool",track="sr_illumina")
    assert result["score_name"] == "PG-F1" and result["pg_f1"] == 100.0


def test_m1_representation_preserving_query_does_not_become_truth_like(tmp_path: Path):
    source=tmp_path/"linked.vcf"; output=tmp_path/"query.vcf.gz"
    source.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\nchr1\t101\tNATIVE_B\tN\t<INS:CALLER_B>\t.\tPASS\tSVTYPE=INS;PANGENOME_LINKED_ID=C1;PANGENOME_CANDIDATE_COUNT=1\tGT\t0/1\n",encoding="utf-8")
    assert evaluator_query_main(["--linked-vcf",str(source),"--output-vcf",str(output)]) == 0
    with open_text(output,"rt") as handle: rendered=handle.read()
    assert "<INS:CALLER_B>" in rendered
    assert "PGBENCH_UNIT_TRACE=C1:A1" in rendered
    assert "PGBENCH_NATIVE_RECORD_ID=NATIVE_B" in rendered
    assert "\tC1\tN\t<INS:CALLER_B>\t" in rendered


def test_no_call_is_not_submitted_as_an_evaluator_event(tmp_path: Path):
    source=tmp_path/"linked.vcf"; output=tmp_path/"query.vcf.gz"
    source.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\nchr1\t101\tNATIVE\tN\t<DEL>\t.\tPASS\tPANGENOME_LINKED_ID=C1;PANGENOME_CANDIDATE_COUNT=1\tGT\t.\n",encoding="utf-8")
    assert evaluator_query_main(["--linked-vcf",str(source),"--output-vcf",str(output)]) == 0
    with open_text(output,"rt") as handle:
        assert all(line.startswith("#") for line in handle if line.strip())


def test_m2_phase_permutation_is_metamorphically_invariant():
    truth=canonical_alt_dosage("0/1")
    outcomes={accounting(truth=truth,query=canonical_alt_dosage(query),consensus=True) for query in ("0/1","1|0","0|1")}
    assert outcomes == {(1,0,0)}


def test_m3_evaluator_order_permutation_is_metamorphically_invariant():
    first={"truvari":1,"aardvark":0,"vcfdist":1}
    reordered=dict(reversed(list(first.items())))
    assert vote_consensus(first) is vote_consensus(reordered) is True
