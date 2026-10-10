import importlib.util
import json
from pathlib import Path
import pytest

path = Path(__file__).resolve().parents[2] / "plugins/vg/projection.py"
spec = importlib.util.spec_from_file_location("projection", path)
p = importlib.util.module_from_spec(spec)
spec.loader.exec_module(p)


def fixture(tmp_path, records):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nAAAAAAAAAA\n")
    Path(str(ref) + ".fai").write_text("chr1\t10\t6\t10\t11\n")
    panel = tmp_path / "panel.vcf"
    panel.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t3\tC1\tA\tAA\t60\tPASS\tAT=>1>2,>1>3>2;PANGENOME_ALLELE_ID=C1\n")
    native = tmp_path / "native.vcf"
    native.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n" + records)
    return ref, panel, native


def test_shifted_equivalent_path_preserves_native_quality_and_alleles(tmp_path):
    ref, panel, native = fixture(tmp_path, "chr1\t5\t>1>2\tA\tAA\t8.61\tlowad\tAT=>1>2,>1>3>2\tGT:GQ:DP:AD\t0/1:1:8:7,1\n")
    output = tmp_path / "calls.vcf"
    assert p.project(native, panel, output, "HG002", reference=ref)[:2] == (1, 0)
    f = next(l for l in output.read_text().splitlines() if not l.startswith("#")).split("\t")
    assert f[1:7] == ["5", "C1", "A", "AA", "8.61", "lowad"]
    assert f[8:] == ["GT:GQ:DP", "0/1:1:8"]
    audit = json.loads(output.with_suffix(".projection.json").read_text())
    assert audit["truth_used"] is False and audit["caller_alleles_preserved"]


def test_matching_graph_paths_do_not_override_sequence_mismatch(tmp_path):
    ref, panel, native = fixture(tmp_path, "chr1\t5\t>1>2\tA\tAC\t50\tPASS\tAT=>1>2,>1>3>2\tGT\t0/1\n")
    out = tmp_path / "calls.vcf"
    assert p.project(native, panel, out, "HG002", reference=ref)[:2] == (0, 1)
    assert "PGBENCH_NATIVE_STATUS=not_observed" in out.read_text()


def test_multiallelic_dosage_and_trace(tmp_path):
    ref, panel, native = fixture(tmp_path, "chr1\t3\t>1>2\tA\tAC,AA\t50\tPASS\tAT=>1>2,>1>4>2,>1>3>2\tGT\t1/2\n")
    out = tmp_path / "calls.vcf"
    assert p.project(native, panel, out, "HG002", reference=ref)[0] == 1
    assert "\tA\tAA\t50\tPASS\t" in out.read_text()
    assert out.read_text().rstrip().endswith("\tGT\t0/1")


def test_conflicting_native_records_fail_closed(tmp_path):
    ref, panel, native = fixture(tmp_path, "chr1\t3\tfirst\tA\tAA\t50\tPASS\t.\tGT\t0/1\nchr1\t3\tsecond\tA\tAA\t50\tPASS\t.\tGT\t1/1\n")
    out = tmp_path / "calls.vcf"
    assert p.project(native, panel, out, "HG002", reference=ref)[0] == 0
    assert "PGBENCH_NATIVE_STATUS=ambiguous_mapping" in out.read_text()


def test_wrong_sample_is_rejected(tmp_path):
    ref, panel, native = fixture(tmp_path, "")
    with pytest.raises(ValueError, match="single sample"):
        p.project(native, panel, tmp_path / "calls.vcf", "OTHER", reference=ref)


@pytest.mark.parametrize("gt,index,count,expected", [("1/2", 2, 2, "0/1"), ("2|2", 2, 2, "1/1"), ("./1", 1, 1, "./."), ("2/2", 1, 1, "./.")])
def test_dosage(gt, index, count, expected):
    assert p.dosage_gt(gt, index, count) == expected
