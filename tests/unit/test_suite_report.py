from __future__ import annotations

import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from render_suite_report import SuiteReportError, render_suite  # noqa: E402


def _score(tool: str, mode: str, value: float) -> dict:
    return {
        "tuple_key": {
            "run_id": f"run-{tool}",
            "sample": "HG002",
            "tool": tool,
            "official_score_mode": mode,
            "primary_truth_profile": "giab_hg002_grch38_v5_0q",
        },
        "score_profile_sha256": "a" * 64,
        "comparable_score": value,
        "consensus_score": value + 1,
        "truth_eligible_count": 100,
        "total_evaluated": 90,
        "score_status": "valid",
    }


def _manifest(tool: str, paradigm: str, technology: str) -> dict:
    return {
        "id": tool,
        "paradigm": paradigm,
        "capabilities": {"technology": [technology]},
    }


def _asset_hashes() -> dict:
    return {
        "reference": "4" * 64,
        "truth_vcf": "5" * 64,
        "benchmark_bed": "6" * 64,
        "pangenome_manifest": "7" * 64,
        "challenge_hidden_ledger": "8" * 64,
        "graph_asset_lock": None,
    }


def _with_evaluator_profile(
    score: dict,
    profile_hash: str = "b" * 64,
    *,
    transitional_field_names: bool = False,
) -> dict:
    analysis_field = "analysis" if transitional_field_names else "formal_analysis"
    hash_field = (
        "evaluator_profile_hash"
        if transitional_field_names
        else "evaluator_profile_sha256"
    )
    score[analysis_field] = {
        hash_field: profile_hash,
        "evaluator_versions": {
            "truvari": {"sha256": "1" * 64},
            "aardvark": {"sha256": "2" * 64},
            "vcfdist": {"sha256": "3" * 64},
        },
        "asset_hashes": _asset_hashes(),
    }
    return score


def _with_evidence(
    score: dict,
    *,
    technology: str = "illumina_short_read",
    source: str = "HG002_ILLUMINA_2X151",
    digest: str = "c" * 64,
) -> dict:
    if "formal_analysis" not in score:
        _with_evaluator_profile(score)
    analysis = score["formal_analysis"]
    analysis["evidence_profile"] = {
        "contract": "pgbench_evidence_v1",
        "actual_technology": technology,
        "library_id": "HG002_LIB",
        "source_evidence_id": source,
        "evidence_kind": "paired_fastq",
        "official_score_mode": score["tuple_key"]["official_score_mode"],
        "input_assets": {
            "short_fastq_r1": {
                "sha256": digest,
                "size_bytes": 100,
                "path_type": "file",
            },
            "short_fastq_r2": {
                "sha256": "d" * 64,
                "size_bytes": 100,
                "path_type": "file",
            },
        },
        "resolved_inputs_sha256": "e" * 64,
        "coverage_x": 30.0,
        "read_count": 2,
        "read_bases": 302,
        "downsampling_seed": 0,
        "quantitative_metadata_complete": True,
    }
    return score


def test_unified_suite_accepts_same_sample_truth_and_profile() -> None:
    html = render_suite(
        [
            (
                _score("pangenie", "end_to_end_from_reads", 75.0),
                _manifest(
                    "pangenie",
                    "short_read_pangenome_genotyping",
                    "illumina_short_read",
                ),
            ),
            (
                _score("kanpig", "caller_only_shared_alignment", 80.0),
                _manifest("kanpig", "genotyping_only", "pacbio_clr"),
            ),
        ]
    )
    assert "PanGenie".casefold() in html.casefold()
    assert "kanpig" in html
    assert "ComparableScore" in html
    assert "illumina_short_read" in html
    assert "pacbio_clr" in html


def test_unified_suite_rejects_different_truth_profile() -> None:
    first = _score("a", "end_to_end_from_reads", 70.0)
    second = _score("b", "caller_only_shared_alignment", 80.0)
    second["tuple_key"]["primary_truth_profile"] = "different"
    with pytest.raises(SuiteReportError, match="share sample, truth profile"):
        render_suite(
            [
                (first, _manifest("a", "mapping_based", "illumina_short_read")),
                (second, _manifest("b", "mapping_based", "pacbio_clr")),
            ]
        )


def test_suite_accepts_legacy_scores_when_all_lack_evaluator_hash() -> None:
    html = render_suite(
        [
            (
                _score("legacy-a", "end_to_end_from_reads", 70.0),
                _manifest(
                    "legacy-a",
                    "short_read_pangenome_genotyping",
                    "illumina_short_read",
                ),
            ),
            (
                _score("legacy-b", "caller_only_shared_alignment", 80.0),
                _manifest("legacy-b", "genotyping_only", "pacbio_clr"),
            ),
        ]
    )
    assert "legacy-a" in html
    assert "legacy-b" in html
    assert "未声明正式评测器 profile" in html


def test_suite_report_contains_readable_utf8_chinese() -> None:
    html = render_suite(
        [
            (
                _score("vg", "end_to_end_from_reads", 75.0),
                _manifest("vg", "graph_pangenome", "pacbio_clr"),
            )
        ]
    )
    assert "泛基因组结构变异综合比较" in html
    assert "测序技术" in html
    assert "�" not in html


def test_suite_accepts_matching_formal_evaluator_hashes() -> None:
    html = render_suite(
        [
            (
                _with_evaluator_profile(
                    _score("pangenie", "end_to_end_from_reads", 75.0)
                ),
                _manifest(
                    "pangenie",
                    "short_read_pangenome_genotyping",
                    "illumina_short_read",
                ),
            ),
            (
                _with_evaluator_profile(
                    _score("kanpig", "caller_only_shared_alignment", 80.0),
                    transitional_field_names=True,
                ),
                _manifest("kanpig", "genotyping_only", "pacbio_clr"),
            ),
        ]
    )
    assert "ComparableScore" in html
    assert "b" * 64 in html


def test_suite_rejects_mixed_presence_of_evaluator_hash() -> None:
    with pytest.raises(SuiteReportError, match="cannot mix formal scores"):
        render_suite(
            [
                (
                    _with_evaluator_profile(
                        _score("formal", "end_to_end_from_reads", 75.0)
                    ),
                    _manifest(
                        "formal",
                        "short_read_pangenome_genotyping",
                        "illumina_short_read",
                    ),
                ),
                (
                    _score("legacy", "caller_only_shared_alignment", 80.0),
                    _manifest("legacy", "genotyping_only", "pacbio_clr"),
                ),
            ]
        )


def test_suite_rejects_different_evaluator_hashes() -> None:
    with pytest.raises(SuiteReportError, match="same evaluator profile hash"):
        render_suite(
            [
                (
                    _with_evaluator_profile(
                        _score("a", "end_to_end_from_reads", 75.0),
                        "a" * 64,
                    ),
                    _manifest(
                        "a",
                        "short_read_pangenome_genotyping",
                        "illumina_short_read",
                    ),
                ),
                (
                    _with_evaluator_profile(
                        _score("b", "caller_only_shared_alignment", 80.0),
                        "b" * 64,
                    ),
                    _manifest("b", "genotyping_only", "pacbio_clr"),
                ),
            ]
        )


def test_suite_rejects_different_evaluator_versions() -> None:
    first = _with_evaluator_profile(
        _score("a", "end_to_end_from_reads", 75.0)
    )
    second = _with_evaluator_profile(
        _score("b", "caller_only_shared_alignment", 80.0)
    )
    second["formal_analysis"]["evaluator_versions"]["truvari"]["sha256"] = (
        "f" * 64
    )
    with pytest.raises(SuiteReportError, match="version fingerprints"):
        render_suite(
            [
                (
                    first,
                    _manifest("a", "mapping_based", "illumina_short_read"),
                ),
                (
                    second,
                    _manifest("b", "mapping_based", "pacbio_clr"),
                ),
            ]
        )


def test_suite_rejects_different_benchmark_asset_hashes() -> None:
    first = _with_evaluator_profile(
        _score("a", "end_to_end_from_reads", 75.0)
    )
    second = _with_evaluator_profile(
        _score("b", "caller_only_shared_alignment", 80.0)
    )
    second["formal_analysis"]["asset_hashes"]["reference"] = "9" * 64
    with pytest.raises(
        SuiteReportError,
        match="same reference, truth, BED, pangenome manifest",
    ):
        render_suite(
            [
                (
                    first,
                    _manifest("a", "mapping_based", "illumina_short_read"),
                ),
                (
                    second,
                    _manifest("b", "mapping_based", "pacbio_clr"),
                ),
            ]
        )


@pytest.mark.parametrize(
    "field",
    ["pangenome_manifest", "challenge_hidden_ledger"],
)
def test_suite_rejects_different_shared_pangenome_assets(field: str) -> None:
    first = _with_evaluator_profile(
        _score("a", "end_to_end_from_reads", 75.0)
    )
    second = _with_evaluator_profile(
        _score("b", "caller_only_shared_alignment", 80.0)
    )
    second["formal_analysis"]["asset_hashes"][field] = "9" * 64
    with pytest.raises(SuiteReportError, match="pangenome manifest"):
        render_suite(
            [
                (
                    first,
                    _manifest("a", "mapping_based", "illumina_short_read"),
                ),
                (
                    second,
                    _manifest("b", "mapping_based", "pacbio_clr"),
                ),
            ]
        )


def test_suite_rejects_formal_analysis_without_complete_asset_hashes() -> None:
    score = _with_evaluator_profile(
        _score("a", "end_to_end_from_reads", 75.0)
    )
    score["formal_analysis"].pop("asset_hashes")
    with pytest.raises(SuiteReportError, match="asset_hashes is required"):
        render_suite(
            [
                (
                    score,
                    _manifest("a", "mapping_based", "illumina_short_read"),
                )
            ]
        )

    score = _with_evaluator_profile(
        _score("a", "end_to_end_from_reads", 75.0)
    )
    score["formal_analysis"]["asset_hashes"]["reference"] = None
    with pytest.raises(SuiteReportError, match="reference must be a non-empty"):
        render_suite(
            [
                (
                    score,
                    _manifest("a", "mapping_based", "illumina_short_read"),
                )
            ]
        )


def test_suite_allows_distinct_graph_locks_for_distinct_tools() -> None:
    first = _with_evaluator_profile(
        _score("vg", "end_to_end_from_reads", 75.0)
    )
    second = _with_evaluator_profile(
        _score("pangenie", "end_to_end_from_reads", 80.0)
    )
    first["formal_analysis"]["asset_hashes"]["graph_asset_lock"] = "a" * 64
    second["formal_analysis"]["asset_hashes"]["graph_asset_lock"] = None
    html = render_suite(
        [
            (first, _manifest("vg", "graph_based", "illumina_short_read")),
            (
                second,
                _manifest(
                    "pangenie",
                    "short_read_pangenome_genotyping",
                    "illumina_short_read",
                ),
            ),
        ]
    )
    assert "vg" in html
    assert "pangenie" in html


def test_suite_rejects_graph_lock_change_for_same_tool_and_mode() -> None:
    first = _with_evaluator_profile(
        _score("vg", "end_to_end_from_reads", 75.0)
    )
    second = _with_evaluator_profile(
        _score("vg", "end_to_end_from_reads", 80.0)
    )
    first["formal_analysis"]["asset_hashes"]["graph_asset_lock"] = "a" * 64
    second["formal_analysis"]["asset_hashes"]["graph_asset_lock"] = "b" * 64
    with pytest.raises(SuiteReportError, match="same graph asset lock hash"):
        render_suite(
            [
                (first, _manifest("vg", "graph_based", "illumina_short_read")),
                (second, _manifest("vg", "graph_based", "illumina_short_read")),
            ]
        )


def test_suite_uses_actual_frozen_evidence_instead_of_capability() -> None:
    score = _with_evidence(
        _score("tool", "end_to_end_from_reads", 75.0)
    )
    html = render_suite(
        [(score, _manifest("tool", "graph_based", "pacbio_clr"))]
    )
    assert "illumina_short_read" in html
    assert "HG002_ILLUMINA_2X151" in html
    assert "冻结输入证据" in html


def test_suite_rejects_different_evidence_for_same_technology_and_mode() -> None:
    first = _with_evidence(
        _score("a", "end_to_end_from_reads", 75.0), digest="a" * 64
    )
    second = _with_evidence(
        _score("b", "end_to_end_from_reads", 80.0), digest="b" * 64
    )
    with pytest.raises(SuiteReportError, match="identical frozen"):
        render_suite(
            [
                (first, _manifest("a", "graph_based", "illumina_short_read")),
                (second, _manifest("b", "graph_based", "illumina_short_read")),
            ]
        )


def test_suite_allows_different_evidence_across_actual_technologies() -> None:
    short = _with_evidence(
        _score("pangenie", "end_to_end_from_reads", 75.0),
        technology="illumina_short_read",
        source="ILLUMINA",
    )
    long = _with_evidence(
        _score("vg", "end_to_end_from_reads", 80.0),
        technology="pacbio_clr",
        source="PACBIO",
        digest="f" * 64,
    )
    html = render_suite(
        [
            (short, _manifest("pangenie", "panel", "illumina_short_read")),
            (long, _manifest("vg", "graph", "pacbio_clr")),
        ]
    )
    assert "ILLUMINA" in html
    assert "PACBIO" in html
