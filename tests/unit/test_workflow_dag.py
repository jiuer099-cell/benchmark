from __future__ import annotations

import importlib.util
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]


def test_rules_resolve_python_from_the_activated_environment() -> None:
    for path in (
        ROOT / "Snakefile",
        ROOT / "workflow" / "modules" / "generic_external" / "Snakefile",
    ):
        source = path.read_text(encoding="utf-8")
        assert 'PYTHON_EXECUTABLE = "python"' in source
        assert "PYTHON_EXECUTABLE = sys.executable" not in source


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="Snakemake is required for workflow DAG integration tests",
)
def test_external_plugin_declared_vcf_and_helper_are_dag_inputs(
    tmp_path: Path,
) -> None:
    plugin_root = tmp_path / "example_genotyper"
    shutil.copytree(ROOT / "plugins" / "example_genotyper", plugin_root)
    helper = plugin_root / "helper.py"
    helper.write_text("VALUE = 1\n", encoding="utf-8")

    manifest_path = plugin_root / "tool.yaml"
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    manifest["outputs"]["vcf"] = "custom/nested/calls.vcf.gz"
    manifest_path.write_text(
        yaml.safe_dump(manifest, sort_keys=False),
        encoding="utf-8",
    )

    config = yaml.safe_load(
        (ROOT / "tests" / "fixtures" / "synthetic" / "config.yaml").read_text(
            encoding="utf-8"
        )
    )
    config["external_plugins"][0]["manifest"] = str(manifest_path)
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        yaml.safe_dump(config, sort_keys=False),
        encoding="utf-8",
    )
    runtime_cache = tmp_path / "snakemake-runtime-cache"
    runtime_cache.mkdir()

    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "snakemake",
            "--runtime-source-cache-path",
            str(runtime_cache),
            "--snakefile",
            "Snakefile",
            "--configfile",
            str(config_path),
            "--cores",
            "1",
            "--dry-run",
            "--printshellcmds",
        ],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 0, completed.stderr
    output = completed.stdout + completed.stderr
    assert (
        "results/synthetic_smoke/HG002/end_to_end_from_reads/example_genotyper/"
        "custom/nested/calls.vcf.gz"
    ) in output
    assert str(helper) in output


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="Snakemake is required for workflow DAG integration tests",
)
def test_graph_directory_and_lock_manifest_reach_external_plugin(
    tmp_path: Path,
) -> None:
    plugin_root = tmp_path / "graph_genotyper"
    shutil.copytree(ROOT / "plugins" / "example_genotyper", plugin_root)
    manifest_path = plugin_root / "tool.yaml"
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    manifest["supported_modes"]["end_to_end_from_reads"]["optional_inputs"] = [
        "graph_assets"
    ]
    manifest_path.write_text(
        yaml.safe_dump(manifest, sort_keys=False),
        encoding="utf-8",
    )

    graph_root = tmp_path / "graph"
    graph_root.mkdir()
    graph_manifest = graph_root / "graph-assets.lock.yaml"
    graph_manifest.write_text("producer: synthetic-vg\n", encoding="utf-8")
    for filename in ("graph.gbz", "graph.xg", "graph.min", "graph.dist"):
        (graph_root / filename).write_bytes(filename.encode())
    (graph_root / "samples.txt").write_text("HG001\nHG003\n", encoding="utf-8")

    config = yaml.safe_load(
        (ROOT / "tests" / "fixtures" / "synthetic" / "config.yaml").read_text(
            encoding="utf-8"
        )
    )
    config["external_plugins"][0]["manifest"] = str(manifest_path)
    config["pangenome"]["build_graph_assets"] = True
    config["pangenome"]["graph_assets"] = {
        "manifest": str(graph_manifest),
        "gbz": str(graph_root / "graph.gbz"),
        "xg": str(graph_root / "graph.xg"),
        "min": str(graph_root / "graph.min"),
        "dist": str(graph_root / "graph.dist"),
        "sample_list": str(graph_root / "samples.txt"),
        "reference_path": "GRCh38",
    }
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        yaml.safe_dump(config, sort_keys=False),
        encoding="utf-8",
    )
    runtime_cache = tmp_path / "snakemake-runtime-cache"
    runtime_cache.mkdir()

    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "snakemake",
            "--runtime-source-cache-path",
            str(runtime_cache),
            "--snakefile",
            "Snakefile",
            "--configfile",
            str(config_path),
            "--cores",
            "1",
            "--dry-run",
            "--printshellcmds",
        ],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 0, completed.stderr
    output = completed.stdout + completed.stderr
    assert f"graph_assets={graph_root}" in output
    assert (
        "results/synthetic_smoke/pangenome/synthetic_grch38_pg_v1/"
        "graph-assets.lock.yaml"
    ) in output
    assert str(graph_manifest) in output


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="Snakemake is required for workflow DAG integration tests",
)
def test_formal_dag_reaches_all_three_evaluators(tmp_path: Path) -> None:
    resources = tmp_path / "resources"
    resources.mkdir()

    def resource(name: str, content: bytes = b"fixture\n") -> str:
        path = resources / name
        path.write_bytes(content)
        return str(path)

    reference = resource("reference.fa", b">chr1\nA\n")
    reference_fai = resource("reference.fa.fai")
    reference_dict = resource("reference.dict")
    bam = resource("HG002.bam")
    resource("HG002.bam.bai")
    population = resource("population.vcf.gz")
    truth = resource("truth.vcf.gz")
    resource("truth.vcf.gz.tbi")
    regions = resource("truth.bed", b"chr1\t0\t1\n")
    graph_manifest = resource("graph-assets.lock.yaml", b"producer: fixture\n")
    graph_gbz = resource("graph.gbz")
    graph_xg = resource("graph.xg")
    graph_min = resource("graph.min")
    graph_dist = resource("graph.dist")
    graph_samples = resource("samples.txt", b"HG001\nHG003\n")

    config = yaml.safe_load(
        (ROOT / "config" / "config.example.yaml").read_text(encoding="utf-8")
    )
    config["sample"]["bam"] = bam
    config["reference"].update(
        fasta=reference,
        fai=reference_fai,
        dict=reference_dict,
    )
    config["caller_only"]["shared_bam"] = bam
    config["pangenome"]["population_vcf"] = population
    config["pangenome"]["graph_assets"].update(
        manifest=graph_manifest,
        gbz=graph_gbz,
        xg=graph_xg,
        min=graph_min,
        dist=graph_dist,
        sample_list=graph_samples,
    )
    config["evaluation"].update(
        truth_vcf=truth,
        benchmark_bed=regions,
    )
    config_path = tmp_path / "formal-config.yaml"
    config_path.write_text(
        yaml.safe_dump(config, sort_keys=False),
        encoding="utf-8",
    )
    runtime_cache = tmp_path / "snakemake-runtime-cache"
    runtime_cache.mkdir()

    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "snakemake",
            "--runtime-source-cache-path",
            str(runtime_cache),
            "--snakefile",
            "Snakefile",
            "--configfile",
            str(config_path),
            "--cores",
            "1",
            "--dry-run",
        ],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 0, completed.stdout + completed.stderr
    output = completed.stdout + completed.stderr
    for evaluator in ("truvari", "aardvark", "vcfdist"):
        assert (
            "results/hg002_clr_graph_panel_kanpig/evaluation/"
            f"kanpig/{evaluator}/votes.tsv"
        ) in output
    assert "materialize_formal_consensus_metrics.py" in output
