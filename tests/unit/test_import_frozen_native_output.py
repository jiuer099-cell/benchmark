from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from import_frozen_native_output import (  # noqa: E402
    FrozenNativeImportError,
    _assert_source_output,
    import_frozen_native_output,
)
from pgbench_provenance import sha256_file  # noqa: E402


def _source_manifest(raw: Path, resolved: Path, attempt: Path) -> dict[str, object]:
    return {
        "status": "success",
        "rule_name": "tool__paragraph__execute",
        "run_id": "paragraph-source-v3",
        "manifest_id": "source-manifest-id",
        "git_head": "a" * 40,
        "output_sha256": {
            "results/source/HG002/end_to_end_from_reads/paragraph/raw/calls.vcf": sha256_file(raw),
            "results/source/HG002/end_to_end_from_reads/paragraph/meta/resolved_inputs.json": sha256_file(resolved),
            "results/source/HG002/end_to_end_from_reads/paragraph/meta/attempt.json": sha256_file(attempt),
        },
    }


def test_imports_only_manifest_authenticated_native_artifacts(tmp_path: Path) -> None:
    source = tmp_path / "source"
    source.mkdir()
    raw = source / "calls.vcf"
    resolved = source / "resolved_inputs.json"
    attempt = source / "attempt.json"
    raw.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
    resolved.write_text(json.dumps({"inputs": []}), encoding="utf-8")
    attempt.write_text(json.dumps({"attempt": "source"}), encoding="utf-8")
    manifest = source / "tool-manifest.json"
    manifest.write_text(json.dumps(_source_manifest(raw, resolved, attempt)), encoding="utf-8")
    destination = tmp_path / "replay"
    args = type("Args", (), {
        "tool_id": "paragraph",
        "source_vcf": str(raw),
        "source_resolved_inputs": str(resolved),
        "source_attempt_record": str(attempt),
        "source_tool_manifest": str(manifest),
        "output_vcf": str(destination / "raw" / "calls.vcf"),
        "output_resolved_inputs": str(destination / "meta" / "resolved_inputs.json"),
        "output_attempt_record": str(destination / "meta" / "attempt.json"),
        "import_audit": str(destination / "meta" / "frozen-native-import.json"),
    })()

    import_frozen_native_output(args)

    assert sha256_file(destination / "raw" / "calls.vcf") == sha256_file(raw)
    audit = json.loads(
        (destination / "meta" / "frozen-native-import.json").read_text(encoding="utf-8")
    )
    assert audit["status"] == "valid"
    assert audit["source_run_id"] == "paragraph-source-v3"


def test_rejects_source_artifact_not_declared_by_manifest(tmp_path: Path) -> None:
    source = tmp_path / "source.vcf"
    source.write_text("changed", encoding="utf-8")
    with pytest.raises(FrozenNativeImportError, match="SHA-256"):
        _assert_source_output(
            {"output_sha256": {"results/raw/calls.vcf": "0" * 64}},
            source,
            "/raw/calls.vcf",
            "VCF",
        )
