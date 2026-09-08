from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import jsonschema
import yaml

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from build_coverage_fastqs import build  # noqa: E402
from freeze_information_contract import freeze  # noqa: E402


def _yaml(path: Path) -> dict:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert isinstance(value, dict)
    return value


def _fastq(path: Path, mate: int, count: int = 40) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for index in range(count):
            handle.write(f"@read{index}/{mate}\nACGT\n+\nIIII\n")


def test_every_adapter_uses_one_short_read_all_sites_contract() -> None:
    schema = _yaml(ROOT / "workflow" / "schemas" / "tool.schema.yaml")
    manifests = sorted((ROOT / "plugins").glob("*/tool.yaml"))
    assert {path.parent.name for path in manifests} == {
        "bayestyper",
        "example_genotyper",
        "graphtyper2",
        "pangenie",
        "paragraph",
        "varigraph",
        "vg",
    }
    for path in manifests:
        manifest = _yaml(path)
        jsonschema.Draft202012Validator(schema).validate(manifest)
        mode = manifest["supported_modes"]["end_to_end_from_reads"]
        assert {
            "short_fastq_r1",
            "short_fastq_r2",
            "pangenome_panel",
            "candidate_panel",
        }.issubset(
            mode["required_inputs"]
        )
        assert manifest["outputs"] == {
            "vcf": manifest["outputs"]["vcf"],
            "candidate_output_contract": "all_sites",
            "absence_semantics": "no_call",
        }
        information = manifest["information_contract"]
        assert information["target_truth_used_for_calling"] is False
        assert information["target_assembly_used"] is False
        assert information["target_family_genotypes_used"] is False


def test_bundled_and_community_sets_are_explicit_and_unified() -> None:
    manifests = {
        path.parent.name: _yaml(path)
        for path in sorted((ROOT / "plugins").glob("*/tool.yaml"))
        if path.parent.name != "example_genotyper"
    }
    assert {name for name, manifest in manifests.items() if manifest["source"] == "builtin"} == {
        "pangenie",
        "vg",
        "paragraph",
    }
    assert {name for name, manifest in manifests.items() if manifest["source"] == "external"} == {
        "graphtyper2",
        "varigraph",
        "bayestyper",
    }

    config = _yaml(ROOT / "config" / "config.unified-tools.example.yaml")
    registrations = config["external_plugins"]
    assert [entry["id"] for entry in registrations] == [
        "pangenie",
        "vg",
        "paragraph",
        "graphtyper2",
        "varigraph",
        "bayestyper",
    ]
    assert config["pangenome"]["id"] == "HG002_LOO_HPRC_GRCh38_SV_v1"
    assert len({config["sample"]["fastq_r1"], config["sample"]["fastq_r2"]}) == 2


def test_active_tree_contains_no_removed_benchmark_tracks() -> None:
    forbidden = ("kan" + "pig", "sva" + "rp", "pac" + "bio_clr", "caller" + "_only", "variant" + "_sites", "novel" + "_truth", "canonical" + "_fastq")
    roots = [ROOT / "Snakefile", ROOT / "config", ROOT / "plugins", ROOT / "workflow"]
    hits: list[str] = []
    for root in roots:
        paths = [root] if root.is_file() else root.rglob("*")
        for path in paths:
            if not path.is_file() or path.suffix in {".pyc", ".pyo"}:
                continue
            text = path.read_text(encoding="utf-8", errors="ignore").casefold()
            for term in forbidden:
                if term in text:
                    hits.append(f"{path.relative_to(ROOT)}:{term}")
    assert hits == []


def test_coverage_subsets_are_paired_deterministic_and_shared(tmp_path: Path) -> None:
    r1, r2 = tmp_path / "input.R1.fastq", tmp_path / "input.R2.fastq"
    _fastq(r1, 1)
    _fastq(r2, 2)
    first = build(
        r1=r1,
        r2=r2,
        output_dir=tmp_path / "first",
        source_coverage=40,
        coverages=[10, 20, 30],
        seeds=[1701, 1702, 1703],
    )
    second = build(
        r1=r1,
        r2=r2,
        output_dir=tmp_path / "second",
        source_coverage=40,
        coverages=[10, 20, 30],
        seeds=[1701, 1702, 1703],
    )
    assert len(first["replicates"]) == 9
    assert [item["fastq_r1_sha256"] for item in first["replicates"]] == [
        item["fastq_r1_sha256"] for item in second["replicates"]
    ]
    for item in first["replicates"]:
        with gzip.open(item["fastq_r1"], "rt", encoding="utf-8") as left, gzip.open(
            item["fastq_r2"], "rt", encoding="utf-8"
        ) as right:
            left_names = [line.split("/")[0] for line in left if line.startswith("@")]
            right_names = [line.split("/")[0] for line in right if line.startswith("@")]
        assert left_names == right_names


def test_information_contract_is_hash_frozen(tmp_path: Path) -> None:
    resolved = tmp_path / "resolved.json"
    names = ["short_fastq_r1", "short_fastq_r2", "candidate_panel"]
    resolved.write_text(
        json.dumps(
            {
                "inputs": [
                    {"name": name, "sha256": str(index) * 64, "size_bytes": 1, "path_type": "file"}
                    for index, name in enumerate(names, start=1)
                ]
            }
        ),
        encoding="utf-8",
    )
    documents = freeze(ROOT / "plugins" / "pangenie" / "tool.yaml", resolved)
    assert documents["allowed_inputs"]["status"] == "valid"
    assert documents["training_or_tuning_status"]["target_truth_inspected"] is False
    assert len(documents["parameter_manifest"]["parameter_sha256"]) == 64


def test_main_configuration_freezes_unified_panel_and_me_f1() -> None:
    config = _yaml(ROOT / "config" / "config.example.yaml")
    contract = config["benchmark_contract"]
    assert contract["chromosomes"] == [f"chr{i}" for i in range(1, 23)]
    assert contract["svtypes"] == ["DEL", "INS"]
    assert contract["minimum_sv_size"] == 50
    assert contract["maximum_sv_size"] == 10000
    assert contract["coverage_levels"] == [10, 20, 30, "full"]
    assert len(contract["downsampling_seeds"]) >= 3
    assert config["score"]["profile"] == "pgbench_me_f1_v1"
    assert config["catalogs"]["score_weights"] == "config/me_f1_scoring.yaml"
