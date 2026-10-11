from __future__ import annotations

import importlib.util
import gzip
import json
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
TEMPLATE_PATH = ROOT / "plugins" / "adapter_template" / "run.py"
SPEC = importlib.util.spec_from_file_location("adapter_template_run", TEMPLATE_PATH)
assert SPEC is not None and SPEC.loader is not None
adapter_template = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = adapter_template
SPEC.loader.exec_module(adapter_template)


def _vcf(path: Path, records: list[str]) -> Path:
    path.write_text(
        "##fileformat=VCFv4.2\n"
        '##INFO=<ID=PANGENOME_ALLELE_ID,Number=1,Type=String,Description="x">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        + "".join(records),
        encoding="utf-8",
    )
    return path


def test_template_projection_writes_auditable_native_trace(tmp_path: Path) -> None:
    candidates = _vcf(
        tmp_path / "candidates.vcf",
        ["chr1\t10\tPGSV_1\tA\tT\t.\t.\tPANGENOME_ALLELE_ID=a1\tGT\t./.\n"],
    )
    native = _vcf(
        tmp_path / "native.vcf",
        ["chr1\t10\t.\tA\tT\t.\t.\tPANGENOME_ALLELE_ID=a1\tGT\t0/1\n"],
    )
    output = tmp_path / "raw" / "calls.vcf"

    assert adapter_template.project_all_sites(native, candidates, output, "HG002") == (1, 0)
    assert output.read_text(encoding="utf-8").endswith("GT\t0/1\n")
    trace = output.with_suffix(".native-trace.jsonl.gz")
    with gzip.open(trace, "rt") as stream:
        row = json.loads(stream.readline())
    assert row["candidate_id"] == "PGSV_1"
    assert row["match_strategy"] == "verified_identity"
    assert row["native_gt"] == "0/1"


def test_template_projection_traces_native_record_outside_scoring_universe(tmp_path: Path) -> None:
    candidates = _vcf(
        tmp_path / "candidates.vcf",
        ["chr1\t10\tPGSV_1\tA\tT\t.\t.\t.\tGT\t./.\n"],
    )
    native = _vcf(
        tmp_path / "native.vcf",
        ["chr1\t11\t.\tA\tC\t.\t.\t.\tGT\t0/1\n"],
    )
    output = tmp_path / "raw" / "calls.vcf"
    assert adapter_template.project_all_sites(native, candidates, output, "HG002") == (0, 1)
    trace = output.with_suffix(".native-trace.jsonl.gz")
    with gzip.open(trace, "rt") as stream:
        row = json.loads(stream.readline())
    assert row["candidate_id"] is None
    assert row["match_strategy"] == "outside_canonical_universe"


def test_template_rejects_multiple_primary_read_modes(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("PGBENCH_INPUT_FASTQ_R1", "/reads/r1.fq.gz")
    monkeypatch.setenv("PGBENCH_INPUT_LONG_READS_FASTQ", "/reads/lr.fq.gz")
    with pytest.raises(RuntimeError, match="more than one primary read evidence"):
        adapter_template.detect_input_mode()


def test_template_exposes_only_resolved_declared_inputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    panel = tmp_path / "panel.vcf"
    panel.write_text("x", encoding="utf-8")
    asset = tmp_path / "index"
    asset.mkdir()
    resolved = tmp_path / "resolved-inputs.json"
    resolved.write_text(
        json.dumps(
            {
                "inputs": [
                    {
                        "contract_name": "candidate_panel",
                        "environment_variable": "PGBENCH_CANDIDATE_VCF",
                    },
                    {
                        "contract_name": "adapter_asset.my_index",
                        "environment_variable": "PGBENCH_ADAPTER_ASSET_MY_INDEX",
                    },
                ]
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setenv("PGBENCH_RESOLVED_INPUTS", str(resolved))
    monkeypatch.setenv("PGBENCH_CANDIDATE_VCF", str(panel))
    monkeypatch.setenv("PGBENCH_ADAPTER_ASSET_MY_INDEX", str(asset))

    assert adapter_template.injected_inputs() == {
        "candidate_panel": panel,
        "adapter_asset.my_index": asset,
    }


def _project(tmp_path: Path, native_record: str, *, native_header: str | None = None,
             candidate_record: str = "chr1\t10\tPGSV_1\tA\tT\t.\tPASS\t.\tGT\t./.\n",
             reference: bool = False):
    candidates = _vcf(tmp_path / "candidates.vcf", [candidate_record])
    native = _vcf(tmp_path / "native.vcf", [native_record])
    if native_header is not None:
        native.write_text("##fileformat=VCFv4.2\n" + native_header + native_record)
    output = tmp_path / "raw" / "calls.vcf"
    kwargs = {}
    if reference:
        ref = tmp_path / "ref.fa"
        ref.write_text(">chr1\n" + "A" * 100 + "\n")
        fai = tmp_path / "ref.fa.fai"
        fai.write_text("chr1\t100\t6\t100\t101\n")
        kwargs = {"reference": ref, "reference_index": fai}
    count = adapter_template.project_all_sites(native, candidates, output, "HG002", **kwargs)
    records = [line.split("\t") for line in output.read_text().splitlines() if not line.startswith("#")]
    return count, records, output


def test_template_selects_target_sample_by_name_and_preserves_quality(tmp_path):
    _, records, output = _project(tmp_path,
        "chr1\t10\tPGSV_1\tA\tT\t12\tLowQual\tUK=0\tGT:GQ:DP:AD:PL\t1/1:2:8:0,8:90,30,0\t0/0:99:30:30,0:0,90,200\n",
        native_header="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG001\tHG002\n")
    assert records[0][5:7] == ["12", "LowQual"]
    assert "UK=0" in records[0][7]
    assert records[0][8:] == ["GT:GQ:DP:AD:PL", "0/0:99:30:30,0:0,90,200"]
    assert "\tFORMAT\tHG002\n" in output.read_text()


@pytest.mark.parametrize("header", [
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG001\n",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\tHG002\n",
])
def test_template_rejects_missing_or_duplicate_target_sample(tmp_path, header):
    with pytest.raises(ValueError, match="requested sample"):
        _project(tmp_path, "chr1\t10\tPGSV_1\tA\tT\t.\tPASS\t.\tGT\t0/1\n", native_header=header)


def test_template_rejects_claimed_id_with_wrong_allele(tmp_path):
    with pytest.raises(ValueError, match="conflicts"):
        _project(tmp_path, "chr2\t999\tPGSV_1\tG\tC\t.\tPASS\t.\tGT\t0/1\n")


@pytest.mark.parametrize("alt,gt,error", [
    ("T,C", "1/2", "multiallelic"), ("T", "1/2", "single ALT"),
    ("T", "0/1/1", "diploid"), ("<INS>", "0/1", "symbolic"),
    ("T", "garbage", "diploid"),
])
def test_template_unsupported_representation_fails_instead_of_no_call(tmp_path, alt, gt, error):
    with pytest.raises(ValueError, match=error):
        _project(tmp_path, f"chr1\t10\tPGSV_1\tA\t{alt}\t.\tPASS\t.\tGT\t{gt}\n")


def test_template_preserves_equivalent_native_indel_position(tmp_path):
    _, records, output = _project(tmp_path,
        "chr1\t11\tPGSV_1\tA\tAA\t23\tLowQual\t.\tGT:GQ\t0|1:8\n",
        candidate_record="chr1\t10\tPGSV_1\tA\tAA\t.\tPASS\t.\tGT\t./.\n", reference=True)
    assert records[0][:7] == ["chr1", "11", "PGSV_1", "A", "AA", "23", "LowQual"]
    audit = json.loads(output.with_suffix(".projection.json").read_text())
    assert audit["reference_equivalence_checked"] and audit["caller_alleles_preserved"]


@pytest.mark.parametrize("ref,alt", [("AA", "A"), ("A", "AA")])
def test_reference_mapping_finds_shifted_indels_without_candidate_ids(tmp_path, ref, alt):
    count, records, output = _project(tmp_path,
        f"chr1\t11\tnative_id\t{ref}\t{alt}\t23\tLowQual\tUK=0\tGT:GQ\t0|1:8\n",
        candidate_record=f"chr1\t10\tPGSV_1\t{ref}\t{alt}\t.\tPASS\t.\tGT\t./.\n", reference=True)
    assert count == (1, 0)
    assert records[0][:7] == ["chr1", "11", "PGSV_1", ref, alt, "23", "LowQual"]
    assert records[0][8:] == ["GT:GQ", "0|1:8"]
    with gzip.open(output.with_suffix(".native-trace.jsonl.gz"), "rt") as stream:
        trace = json.loads(stream.readline())
    assert trace["match_strategy"] == "reference_normalized_sequence"
    assert trace["native_record"][2] == "native_id"


def test_reference_mapping_does_not_match_non_equivalent_indels(tmp_path):
    count, records, _ = _project(tmp_path,
        "chr1\t11\tnative\tA\tAT\t23\tPASS\t.\tGT\t0/1\n",
        candidate_record="chr1\t10\tPGSV_1\tA\tAA\t.\tPASS\t.\tGT\t./.\n", reference=True)
    assert count == (0, 1) and records[0][9] == "./."


def test_reference_mapping_rejects_multiple_equivalent_candidates_even_with_claim(tmp_path):
    with pytest.raises(ValueError, match="ambiguous native-to-candidate"):
        _project(tmp_path, "chr1\t11\tPGSV_1\tAA\tA\t23\tPASS\t.\tGT\t0/1\n",
            candidate_record="chr1\t10\tPGSV_1\tAA\tA\t.\tPASS\t.\tGT\t./.\n"
                             "chr1\t12\tPGSV_2\tAA\tA\t.\tPASS\t.\tGT\t./.\n", reference=True)


def test_reference_normalization_rotates_tandem_repeat_and_handles_contig_start(tmp_path):
    from adapter_projection import Reference
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1\nATATATC\n")
    fai = tmp_path / "reference.fa.fai"
    fai.write_text("chr1\t7\t6\t7\t8\n")
    reference = Reference(fasta, fai)
    try:
        # Different anchors and payload rotations describe the same haplotype.
        deletion1 = ["chr1", "1", ".", "ATA", "A"]
        deletion2 = ["chr1", "2", ".", "TAT", "T"]
        insertion1 = ["chr1", "1", ".", "A", "ATA"]
        insertion2 = ["chr1", "2", ".", "T", "TAT"]
        for first, second in ((deletion1, deletion2), (insertion1, insertion2)):
            assert reference.equivalent(first, second)
            assert reference.normalized_edit(first) == reference.normalized_edit(second)
            assert reference.normalized_edit(first)[1] == 0
        assert reference.normalized_edit(["chr1", "7", ".", "C", "CT"])[1] == 7
    finally:
        reference.stream.close()


def test_reference_normalization_uses_bounded_blocks_for_long_repeats(tmp_path, monkeypatch):
    from adapter_projection import Reference
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1\n" + "A" * 20000 + "\n")
    fai = tmp_path / "reference.fa.fai"
    fai.write_text("chr1\t20000\t6\t20000\t20001\n")
    reference = Reference(fasta, fai)
    fetch = reference.fetch
    intervals = []
    def counted(chrom, start, end):
        intervals.append(end - start)
        return fetch(chrom, start, end)
    monkeypatch.setattr(reference, "fetch", counted)
    try:
        assert reference.normalized_edit(["chr1", "19999", ".", "AA", "A"]) == ("chr1", 0, "A", "")
        assert len(intervals) <= 5 and max(intervals) <= 4096
    finally:
        reference.stream.close()


@pytest.mark.parametrize("sequence", ["ATATATCAAA", "ACACAC", "AAAAAAAAAAAA"])
def test_reference_edit_keys_agree_with_exact_small_haplotype_equivalence(tmp_path, sequence):
    from adapter_projection import Reference
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1\n" + sequence + "\n")
    fai = tmp_path / "reference.fa.fai"
    fai.write_text(f"chr1\t{len(sequence)}\t6\t{len(sequence)}\t{len(sequence)+1}\n")
    reference = Reference(fasta, fai)
    by_haplotype, by_key = {}, {}
    try:
        for start in range(len(sequence)):
            for end in range(start, min(len(sequence), start + 3) + 1):
                for inserted in ("", "A", "T", "C", "AA", "AT", "TA", "ACA"):
                    haplotype = sequence[:start] + inserted + sequence[end:]
                    if haplotype == sequence:
                        continue
                    if start:
                        fields = ["chr1", str(start), ".", sequence[start-1:end], sequence[start-1] + inserted]
                    elif end < len(sequence):
                        fields = ["chr1", "1", ".", sequence[:end+1], inserted + sequence[end]]
                    else:
                        continue
                    normalized = reference.normalized_edit(fields)
                    assert by_haplotype.setdefault(haplotype, normalized) == normalized
                    assert by_key.setdefault(normalized, haplotype) == haplotype
    finally:
        reference.stream.close()


def test_template_rejects_ref_disagreeing_with_frozen_fasta(tmp_path):
    with pytest.raises(ValueError, match="REF disagrees"):
        _project(tmp_path, "chr1\t10\tPGSV_1\tG\tGT\t.\tPASS\t.\tGT\t0/1\n",
            candidate_record="chr1\t10\tPGSV_1\tA\tAT\t.\tPASS\t.\tGT\t./.\n", reference=True)


def test_template_counts_observed_no_calls_and_clears_placeholder_quality(tmp_path):
    count, _, _ = _project(tmp_path, "chr1\t10\tPGSV_1\tA\tT\t.\tPASS\t.\tGT\t./.\n")
    assert count == (1, 1)
    other = tmp_path / "other"
    other.mkdir()
    count, records, _ = _project(other, "chr1\t11\tnative\tA\tC\t.\tPASS\t.\tGT\t0/1\n")
    assert count == (0, 1)
    assert records[0][5:7] == [".", "."]
    assert records[0][9] == "./."


def test_template_rejects_duplicate_native_candidate_mapping(tmp_path):
    line = "chr1\t10\tPGSV_1\tA\tT\t.\tPASS\t.\tGT\t0/1\n"
    with pytest.raises(ValueError, match="duplicates candidate"):
        _project(tmp_path, line + line)


def test_template_trace_streams_large_outside_universe_output(tmp_path):
    import tracemalloc
    panel = _vcf(tmp_path / "panel.vcf", ["chr1\t10\tPGSV_1\tA\tT\t.\tPASS\t.\tGT\t./.\n"])
    native = _vcf(tmp_path / "native.vcf", [])
    with native.open("a") as stream:
        for i in range(15000):
            stream.write(f"chr1\t{i+100}\tN{i}\tA\tC\t.\tPASS\t.\tGT\t0/1\n")
    output = tmp_path / "raw" / "calls.vcf"
    tracemalloc.start()
    adapter_template.project_all_sites(native, panel, output, "HG002")
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    assert peak < 4 * 1024 * 1024
    with gzip.open(output.with_suffix(".native-trace.jsonl.gz"), "rt") as stream:
        assert sum(1 for _ in stream) == 15000


def test_template_rejects_plain_gzip_output_shortcut(tmp_path):
    panel = _vcf(tmp_path / "panel.vcf", ["chr1\t10\tPGSV_1\tA\tT\t.\tPASS\t.\tGT\t./.\n"])
    with pytest.raises(ValueError, match="uncompressed"):
        adapter_template.project_all_sites(panel, panel, tmp_path / "calls.vcf.gz", "HG002")


def test_native_subprocess_failure_is_fatal_not_a_no_call(tmp_path):
    import subprocess
    context = adapter_template.AdapterContext("short_fastq", "HG002", 1,
        tmp_path / "ref", tmp_path / "panel", tmp_path, tmp_path / "calls.vcf", tmp_path)
    with pytest.raises(subprocess.CalledProcessError):
        adapter_template.run_native(context, [sys.executable, "-c", "raise SystemExit(7)"])
    rows = [json.loads(line) for line in (tmp_path / "native-commands.jsonl").read_text().splitlines()]
    assert rows[-1]["returncode"] == 7
    assert not context.output_vcf.exists()
