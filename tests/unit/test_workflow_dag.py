from __future__ import annotations

import importlib.util
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]


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

    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "snakemake",
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
        "results/HG002/end_to_end_from_reads/example_genotyper/"
        "custom/nested/calls.vcf.gz"
    ) in output
    assert str(helper) in output
