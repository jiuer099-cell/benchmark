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


@pytest.mark.parametrize("gt,expected", [("2/2", "./."), ("1/2", "./."), ("0/2", "./."), ("0/1", "0/1"), ("1/1", "1/1"), ("0/0", "0/0")])
def test_called_placeholder_is_no_call_but_uncalled_placeholder_does_not_mask(tmp_path, gt, expected):
    ref, panel, native = fixture(tmp_path, f"chr1\t3\tsite\tA\tAA,*\t50\tPASS\t.\tGT:GQ:DP\t{gt}:20:10\n")
    out = tmp_path / "calls.vcf"
    assert p.project(native, panel, out, "HG002", reference=ref)[0] == 1
    row = next(l for l in out.read_text().splitlines() if not l.startswith("#")).split("\t")
    assert row[9] == expected + ":20:10"
    import gzip
    with gzip.open(out.with_suffix(".native-trace.jsonl.gz"), "rt") as f:
        trace = json.loads(next(f))
    assert trace["native_record"][9] == gt + ":20:10"
    assert trace["projected_gt"] == expected


@pytest.mark.parametrize("content", ["", "##fileformat=VCFv4.2\n", "chr1\t3\tsite\tA\tAA\t.\t.\t.\tGT\t0/1\n"])
def test_empty_or_headerless_native_fails_without_publishing_output(tmp_path, content):
    ref, panel, native = fixture(tmp_path, "")
    native.write_text(content)
    out = tmp_path / "calls.vcf"
    with pytest.raises(ValueError, match="#CHROM"):
        p.project(native, panel, out, "HG002", reference=ref)
    assert not out.exists()
    assert not out.with_suffix(".projection.json").exists()


def test_valid_header_with_no_events_is_explicit_no_call(tmp_path):
    ref, panel, native = fixture(tmp_path, "")
    out = tmp_path / "calls.vcf"
    assert p.project(native, panel, out, "HG002", reference=ref)[:2] == (0, 1)
    assert out.read_text().rstrip().endswith("\tGT\t./.")


@pytest.mark.parametrize("sample_header", ["OTHER\tHG002", "HG002\tOTHER"])
def test_multi_sample_native_is_rejected(tmp_path, sample_header):
    ref, panel, native = fixture(tmp_path, "")
    native.write_text(native.read_text().replace("FORMAT\tHG002", "FORMAT\t" + sample_header))
    with pytest.raises(ValueError, match="single sample"):
        p.project(native, panel, tmp_path / "calls.vcf", "HG002", reference=ref)
