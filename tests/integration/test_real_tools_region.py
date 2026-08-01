"""Opt-in production-tool checks against a bounded real HG002 region.

The preflight gate only checks installed binaries, sandbox declarations, and
the fixture contract.  A separate, more explicit switch executes the formal
KanPIG adapter against the bounded regional inputs.
"""

from __future__ import annotations

import gzip
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import pytest
import yaml


ROOT = Path(__file__).resolve().parents[2]
PREFLIGHT_ENABLED = os.environ.get("PGBENCH_REAL_TOOL_TESTS") == "1"
EXECUTION_ENABLED = os.environ.get("PGBENCH_REAL_TOOL_EXECUTION") == "1"
PANGENIE_EXECUTION_ENABLED = (
    os.environ.get("PGBENCH_REAL_PANGENIE_EXECUTION") == "1"
)
VG_EXECUTION_ENABLED = os.environ.get("PGBENCH_REAL_VG_EXECUTION") == "1"
REGION_MANIFEST_NAME = "region.yaml"
MAX_REGION_SPAN_BP = 2_000_000
MAX_TOTAL_RESOURCE_BYTES = 2 * 1024**3
KANPIG_TIMEOUT_SECONDS = 10 * 60
PANGENIE_TIMEOUT_SECONDS = 20 * 60
VG_TIMEOUT_SECONDS = 20 * 60
VERSION_RE = re.compile(r"(?<!\d)(\d+(?:\.\d+){1,3})(?!\d)")
CLI_HELP_RE = re.compile(r"\b(usage|options?|version|commands?)\b", re.IGNORECASE)

pytestmark = [
    pytest.mark.integration,
    pytest.mark.real_tools,
]

PREFLIGHT_ONLY = pytest.mark.skipif(
    not PREFLIGHT_ENABLED,
    reason="set PGBENCH_REAL_TOOL_TESTS=1 to run real-tool preflight checks",
)
EXECUTION_ONLY = pytest.mark.skipif(
    not EXECUTION_ENABLED,
    reason=(
        "set PGBENCH_REAL_TOOL_EXECUTION=1 to run the formal KanPIG "
        "small-region execution test"
    ),
)
PANGENIE_EXECUTION_ONLY = pytest.mark.skipif(
    not PANGENIE_EXECUTION_ENABLED,
    reason=(
        "set PGBENCH_REAL_PANGENIE_EXECUTION=1 to run the formal PanGenie "
        "small-region execution test"
    ),
)
VG_EXECUTION_ONLY = pytest.mark.skipif(
    not VG_EXECUTION_ENABLED,
    reason=(
        "set PGBENCH_REAL_VG_EXECUTION=1 to run the formal vg small-region "
        "execution test"
    ),
)


@dataclass(frozen=True)
class ToolProbe:
    name: str
    plugin: str
    executable: str
    executable_env: str
    arguments: tuple[tuple[str, ...], ...]
    exact_version: str | None = None
    minimum_version: tuple[int, ...] | None = None


TOOL_PROBES = (
    ToolProbe(
        name="KanPIG",
        plugin="kanpig",
        executable="kanpig",
        executable_env="PGBENCH_KANPIG_BIN",
        arguments=(("--version",), ("version",), ("--help",), ("-h",)),
        exact_version="2.0.2",
    ),
    ToolProbe(
        name="PanGenie",
        plugin="pangenie",
        executable="PanGenie",
        executable_env="PGBENCH_PANGENIE_BIN",
        arguments=(("--version",), ("--help",), ("-h",)),
        exact_version="4.2.1",
    ),
    ToolProbe(
        name="PanGenie-index",
        plugin="pangenie",
        executable="PanGenie-index",
        executable_env="PGBENCH_PANGENIE_INDEX_BIN",
        arguments=(("--version",), ("--help",), ("-h",)),
        exact_version="4.2.1",
    ),
    ToolProbe(
        name="vg",
        plugin="vg",
        executable="vg",
        executable_env="PGBENCH_VG_BIN",
        arguments=(("version",), ("--version",)),
        minimum_version=(1, 63),
    ),
    ToolProbe(
        name="GraphAligner",
        plugin="vg",
        executable="GraphAligner",
        executable_env="PGBENCH_GRAPHALIGNER_BIN",
        arguments=(("--version",), ("-v",), ("--help",)),
    ),
)

PLUGIN_MANIFESTS = (
    ROOT / "plugins" / "kanpig" / "tool.yaml",
    ROOT / "plugins" / "pangenie" / "tool.yaml",
    ROOT / "plugins" / "vg" / "tool.yaml",
)

REQUIRED_RESOURCE_KEYS = {
    "reference_fasta",
    "reference_fai",
    "reference_dict",
    "truth_vcf",
    "truth_vcf_tbi",
    "benchmark_bed",
    "population_vcf",
    "population_vcf_tbi",
    "candidate_vcf",
    "pangenome_manifest",
    "shared_bam",
    "shared_bai",
    "long_fastq",
    "short_fastq_r1",
    "short_fastq_r2",
    "pangenie_panel_vcf",
    "pangenie_panel_vcf_tbi",
    "graph_manifest",
    "graph_gbz",
    "graph_xg",
    "graph_min",
    "graph_dist",
    "graph_samples",
}


def _load_yaml(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        value = yaml.safe_load(handle)
    assert isinstance(value, dict), f"{path} must contain a YAML mapping"
    return value


def _resolve_executable(environment_name: str, default: str) -> str:
    configured = os.environ.get(environment_name, default)
    candidate = Path(configured).expanduser()
    if candidate.is_absolute() or candidate.parent != Path("."):
        assert candidate.is_file(), (
            f"{environment_name} does not name a file: {candidate}"
        )
        assert os.access(candidate, os.X_OK), (
            f"{environment_name} is not executable: {candidate}"
        )
        return str(candidate)
    resolved = shutil.which(configured)
    assert resolved is not None, (
        f"cannot find {configured!r}; put it on PATH or set {environment_name}"
    )
    return resolved


def _probe_cli(executable: str, argument_sets: tuple[tuple[str, ...], ...]) -> str:
    attempts: list[str] = []
    for arguments in argument_sets:
        try:
            completed = subprocess.run(
                [executable, *arguments],
                check=False,
                capture_output=True,
                text=True,
                timeout=20,
            )
        except (OSError, subprocess.TimeoutExpired) as error:
            attempts.append(f"{arguments!r}: {error}")
            continue
        output = "\n".join(
            part.strip()
            for part in (completed.stdout, completed.stderr)
            if part.strip()
        )
        if output and (
            completed.returncode == 0
            or (
                completed.returncode in {1, 2}
                and CLI_HELP_RE.search(output) is not None
            )
        ):
            return output
        attempts.append(
            f"{arguments!r}: exit={completed.returncode}, output={output[:200]!r}"
        )
    pytest.fail(
        f"no usable version/help probe for {executable}; "
        + "; ".join(attempts)
    )


def _version_tuple(value: str) -> tuple[int, ...]:
    return tuple(int(part) for part in value.split("."))


def _assert_bwrap_isolates(executable: str) -> None:
    assert os.name == "posix", "the production bwrap smoke test requires Linux"
    completed = subprocess.run(
        [
            executable,
            "--die-with-parent",
            "--new-session",
            "--unshare-all",
            "--unshare-net",
            "--ro-bind",
            "/",
            "/",
            "--proc",
            "/proc",
            "--dev",
            "/dev",
            "--",
            "/bin/true",
        ],
        check=False,
        capture_output=True,
        text=True,
        timeout=20,
    )
    assert completed.returncode == 0, (
        "bwrap is installed but cannot create the network-isolated namespace "
        f"required by formal execution: {completed.stderr.strip()}"
    )


def _sample_paired_fastq_ids(
    path: Path,
    *,
    mate: int,
    limit: int = 128,
) -> list[str]:
    assert path.name.lower().endswith((".fastq.gz", ".fq.gz")), (
        f"short_fastq_r{mate} must be gzip-compressed FASTQ: {path}"
    )
    mate_token = re.compile(
        rf"(?:^|[._-])R?{mate}(?:[._-]|$)",
        re.IGNORECASE,
    )
    assert mate_token.search(path.name), (
        f"short_fastq_r{mate} filename does not identify mate {mate}: {path.name}"
    )

    identifiers: list[str] = []
    with gzip.open(path, mode="rt", encoding="utf-8") as handle:
        for record_number in range(1, limit + 1):
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline()
            separator = handle.readline()
            quality = handle.readline()
            assert sequence and separator and quality, (
                f"{path} ends inside FASTQ record {record_number}"
            )
            assert header.startswith("@"), (
                f"{path} record {record_number} has no @ header"
            )
            assert separator.startswith("+"), (
                f"{path} record {record_number} has no + separator"
            )
            assert len(sequence.rstrip("\r\n")) == len(quality.rstrip("\r\n")), (
                f"{path} record {record_number} sequence/quality lengths differ"
            )
            identifier = header[1:].split(maxsplit=1)[0]
            identifiers.append(re.sub(r"/[12]$", "", identifier))
    assert identifiers, f"{path} contains no FASTQ records"
    return identifiers


def _validated_region_resources() -> tuple[dict, dict[str, Path]]:
    configured_root = os.environ.get("PGBENCH_REAL_REGION_DIR")
    assert configured_root, (
        "PGBENCH_REAL_REGION_DIR must point to a prepared small-region directory"
    )
    region_root = Path(configured_root).expanduser().resolve(strict=True)
    assert region_root.is_dir(), f"not a directory: {region_root}"
    assert not list(region_root.rglob("*.aria2")), (
        "small-region resources contain an aria2 sidecar and are incomplete"
    )

    manifest_path = region_root / REGION_MANIFEST_NAME
    manifest = _load_yaml(manifest_path)
    assert manifest.get("schema_version") == 1

    region = manifest.get("region")
    assert isinstance(region, dict), "region.yaml requires a region mapping"
    assert isinstance(region.get("id"), str) and region["id"]
    assert region.get("sample_id") == "HG002"
    assert region.get("reference_id") == "grch38"
    assert isinstance(region.get("contig"), str) and region["contig"]
    start = region.get("start")
    end = region.get("end")
    assert isinstance(start, int) and start >= 1
    assert isinstance(end, int) and end >= start
    assert end - start + 1 <= MAX_REGION_SPAN_BP, (
        f"region spans {end - start + 1:,} bp; maximum is "
        f"{MAX_REGION_SPAN_BP:,} bp"
    )

    resources = manifest.get("resources")
    assert isinstance(resources, dict), "region.yaml requires a resources mapping"
    assert set(resources) == REQUIRED_RESOURCE_KEYS, (
        "resource keys differ; "
        f"missing={sorted(REQUIRED_RESOURCE_KEYS - set(resources))}, "
        f"extra={sorted(set(resources) - REQUIRED_RESOURCE_KEYS)}"
    )

    resolved_root = region_root.resolve()
    resolved: dict[str, Path] = {}
    for resource_id, raw_path in sorted(resources.items()):
        assert isinstance(raw_path, str) and raw_path, (
            f"{resource_id} must be a non-empty relative path"
        )
        relative_path = Path(raw_path)
        assert not relative_path.is_absolute(), (
            f"{resource_id} must be relative to PGBENCH_REAL_REGION_DIR"
        )
        path = (region_root / relative_path).resolve(strict=True)
        try:
            path.relative_to(resolved_root)
        except ValueError:
            pytest.fail(f"{resource_id} escapes the region directory: {path}")
        assert path.is_file(), f"{resource_id} is not a regular file: {path}"
        assert path.stat().st_size > 0, f"{resource_id} is empty: {path}"
        resolved[resource_id] = path

    unique_files = set(resolved.values())
    assert len(unique_files) == len(resources), (
        "each resource key must resolve to a distinct regional artifact"
    )
    total_bytes = sum(path.stat().st_size for path in unique_files)
    assert total_bytes <= MAX_TOTAL_RESOURCE_BYTES, (
        f"small-region resources total {total_bytes:,} bytes; maximum is "
        f"{MAX_TOTAL_RESOURCE_BYTES:,} bytes"
    )

    r1_ids = _sample_paired_fastq_ids(
        resolved["short_fastq_r1"],
        mate=1,
    )
    r2_ids = _sample_paired_fastq_ids(
        resolved["short_fastq_r2"],
        mate=2,
    )
    assert len(r1_ids) == len(r2_ids), (
        "sampled R1/R2 FASTQ record counts differ"
    )
    assert r1_ids == r2_ids, (
        "sampled R1/R2 FASTQ read identifiers are not paired in the same order"
    )
    return region, resolved


@PREFLIGHT_ONLY
@pytest.mark.parametrize("probe", TOOL_PROBES, ids=lambda probe: probe.name)
def test_real_tool_binary_and_version(probe: ToolProbe) -> None:
    """Every bundled production executable must be installed and identifiable."""

    executable = _resolve_executable(probe.executable_env, probe.executable)
    output = _probe_cli(executable, probe.arguments)
    versions = [match.group(1) for match in VERSION_RE.finditer(output)]
    assert versions, f"{probe.name} produced no parseable version: {output[:500]}"

    if probe.exact_version is not None:
        manifest = _load_yaml(ROOT / "plugins" / probe.plugin / "tool.yaml")
        assert str(manifest["version"]) == probe.exact_version
        assert probe.exact_version in versions, (
            f"{probe.name} version {versions} does not match the plugin contract "
            f"{probe.exact_version}"
        )
    if probe.minimum_version is not None:
        assert any(
            _version_tuple(version) >= probe.minimum_version
            for version in versions
        ), (
            f"{probe.name} version {versions} is below "
            f"{'.'.join(map(str, probe.minimum_version))}"
        )


@PREFLIGHT_ONLY
def test_bundled_plugins_have_usable_formal_sandboxes() -> None:
    """Formal plugins may not fall back to an unsandboxed execution backend."""

    backends: set[str] = set()
    for manifest_path in PLUGIN_MANIFESTS:
        manifest = _load_yaml(manifest_path)
        execution = manifest.get("execution")
        assert isinstance(execution, dict)
        backend = execution.get("sandbox_backend")
        assert backend in {"bwrap", "apptainer"}, (
            f"{manifest['id']} must use bwrap or apptainer for formal execution"
        )
        if backend == "apptainer":
            assert execution.get("container"), (
                f"{manifest['id']} declares apptainer without a container"
            )
        backends.add(backend)

    for backend in sorted(backends):
        environment_name = f"PGBENCH_{backend.upper()}_BIN"
        executable = _resolve_executable(environment_name, backend)
        output = _probe_cli(executable, (("--version",), ("version",)))
        assert VERSION_RE.search(output), (
            f"{backend} produced no parseable version: {output[:500]}"
        )
        if backend == "bwrap":
            _assert_bwrap_isolates(executable)


@PREFLIGHT_ONLY
def test_small_real_region_resource_contract() -> None:
    """The opt-in fixture must be complete, bounded, and safe to resolve."""

    _validated_region_resources()


@EXECUTION_ONLY
def test_kanpig_formal_small_region_execution(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Run the real KanPIG adapter once through the formal bwrap executor."""

    assert sys.platform.startswith("linux"), (
        "PGBENCH_REAL_TOOL_EXECUTION=1 requires a Linux server"
    )
    configured_binary = os.environ.get("PGBENCH_KANPIG_BIN")
    assert configured_binary, (
        "set PGBENCH_KANPIG_BIN to the KanPIG executable inside its Conda "
        "environment"
    )
    executable = Path(
        _resolve_executable("PGBENCH_KANPIG_BIN", "kanpig")
    ).resolve(strict=True)
    assert executable.parent.name == "bin", (
        "PGBENCH_KANPIG_BIN must resolve to <conda-prefix>/bin/kanpig so the "
        "complete runtime prefix can be mapped read-only into bwrap"
    )
    conda_prefix = executable.parent.parent.resolve(strict=True)
    assert conda_prefix.is_dir()

    bwrap = _resolve_executable("PGBENCH_BWRAP_BIN", "bwrap")
    _assert_bwrap_isolates(bwrap)
    region, resources = _validated_region_resources()

    monkeypatch.setenv("CONDA_PREFIX", str(conda_prefix))
    monkeypatch.setenv(
        "PATH",
        os.pathsep.join((str(executable.parent), os.environ.get("PATH", ""))),
    )

    scripts = ROOT / "workflow" / "scripts"
    monkeypatch.syspath_prepend(str(scripts))
    from pgbench_exec import execute_tool

    output_dir = tmp_path / "kanpig"
    output_dir.mkdir()
    result = execute_tool(
        tool_manifest_path=ROOT / "plugins" / "kanpig" / "tool.yaml",
        schema_path=ROOT / "workflow" / "schemas" / "tool.schema.yaml",
        mode="caller_only_shared_alignment",
        run_id=f"real-region-{region['id']}",
        sample_id="HG002",
        supplied_inputs={
            "shared_alignment": resources["shared_bam"],
            "shared_alignment_index": resources["shared_bai"],
            "reference": resources["reference_fasta"],
            "reference_index": resources["reference_fai"],
            "pangenome_manifest": resources["pangenome_manifest"],
            "pangenome_panel": resources["population_vcf"],
            "candidate_panel": resources["candidate_vcf"],
        },
        output_dir=output_dir,
        resolved_inputs_path=output_dir / "meta" / "resolved-inputs.json",
        attempt_record_path=output_dir / "meta" / "attempt.json",
        log_path=output_dir / "logs" / "kanpig.log",
        threads=2,
        memory_mb=4096,
        cache_policy="isolated_empty_tool_cache",
        alignment_kind="bam",
        execution_purpose="formal",
        timeout_seconds=KANPIG_TIMEOUT_SECONDS,
        random_seed=0,
    )

    attempt = _load_yaml(Path(result.attempt_record_path))
    assert attempt["attempt_id"] == result.attempt_id
    assert attempt["status"] == "success"
    assert attempt["mode"] == "caller_only_shared_alignment"
    assert attempt["execution_purpose"] == "formal"
    assert attempt["sandbox_backend"] == "bwrap"
    assert attempt["isolation_status"] == "sandboxed_bwrap"
    assert attempt["formal_score_eligible"] is True
    assert attempt["timed_out"] is False
    assert attempt["timeout_seconds"] == KANPIG_TIMEOUT_SECONDS
    assert attempt["output"] == result.output.to_dict()

    output_vcf = Path(result.output.path)
    assert output_vcf.is_file() and output_vcf.stat().st_size > 0
    assert result.output.record_count > 0
    assert Path(result.resolved_inputs_path).is_file()
    assert Path(result.attempt_archive_path).is_file()
    assert Path(result.log_path).is_file()
    assert Path(result.log_archive_path).is_file()

    command = list(result.command)
    assert "--unshare-net" in command
    assert "--ro-bind" in command
    assert str(conda_prefix) in command
    assert str(resources["shared_bai"]) in command
    assert str(resources["reference_fai"]) in command


def _activate_single_conda_runtime(
    monkeypatch: pytest.MonkeyPatch,
    executables: list[tuple[str, str]],
) -> Path:
    """Expose co-installed executables from one immutable runtime prefix."""

    resolved = [
        Path(_resolve_executable(environment, default)).resolve(strict=True)
        for environment, default in executables
    ]
    assert all(path.parent.name == "bin" for path in resolved), (
        "formal real-tool executables must live under <runtime-prefix>/bin"
    )
    prefixes = {path.parent.parent for path in resolved}
    assert len(prefixes) == 1, (
        "cooperating executables must be installed in the same runtime prefix "
        "for one read-only bwrap mapping"
    )
    prefix = prefixes.pop().resolve(strict=True)
    monkeypatch.setenv("CONDA_PREFIX", str(prefix))
    monkeypatch.setenv(
        "PATH",
        os.pathsep.join((str(resolved[0].parent), os.environ.get("PATH", ""))),
    )
    return prefix


def _assert_formal_execution_result(result, *, mode: str) -> None:
    attempt = _load_yaml(Path(result.attempt_record_path))
    assert attempt["status"] == "success"
    assert attempt["mode"] == mode
    assert attempt["execution_purpose"] == "formal"
    assert attempt["sandbox_backend"] == "bwrap"
    assert attempt["isolation_status"] == "sandboxed_bwrap"
    assert attempt["formal_score_eligible"] is True
    assert attempt["timed_out"] is False
    output_vcf = Path(result.output.path)
    assert output_vcf.is_file() and output_vcf.stat().st_size > 0
    assert result.output.record_count > 0
    assert Path(result.resolved_inputs_path).is_file()
    assert Path(result.attempt_archive_path).is_file()
    assert Path(result.log_path).is_file()
    assert "--unshare-net" in result.command


@PANGENIE_EXECUTION_ONLY
def test_pangenie_formal_small_region_execution(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Run the real PanGenie adapter once on paired regional Illumina reads."""

    assert sys.platform.startswith("linux")
    for variable in ("PGBENCH_PANGENIE_BIN", "PGBENCH_PANGENIE_INDEX_BIN"):
        assert os.environ.get(variable), f"set {variable} to an absolute executable"
    _activate_single_conda_runtime(
        monkeypatch,
        [
            ("PGBENCH_PANGENIE_BIN", "PanGenie"),
            ("PGBENCH_PANGENIE_INDEX_BIN", "PanGenie-index"),
        ],
    )
    _assert_bwrap_isolates(_resolve_executable("PGBENCH_BWRAP_BIN", "bwrap"))
    region, resources = _validated_region_resources()
    monkeypatch.syspath_prepend(str(ROOT / "workflow" / "scripts"))
    from pgbench_exec import execute_tool

    output_dir = tmp_path / "pangenie"
    output_dir.mkdir()
    result = execute_tool(
        tool_manifest_path=ROOT / "plugins" / "pangenie" / "tool.yaml",
        schema_path=ROOT / "workflow" / "schemas" / "tool.schema.yaml",
        mode="end_to_end_from_reads",
        run_id=f"real-region-{region['id']}",
        sample_id="HG002",
        supplied_inputs={
            "short_fastq_r1": resources["short_fastq_r1"],
            "short_fastq_r2": resources["short_fastq_r2"],
            "reference": resources["reference_fasta"],
            "pangenome_manifest": resources["pangenome_manifest"],
            "pangenome_panel": resources["pangenie_panel_vcf"],
            "candidate_panel": resources["candidate_vcf"],
        },
        output_dir=output_dir,
        resolved_inputs_path=output_dir / "meta" / "resolved-inputs.json",
        attempt_record_path=output_dir / "meta" / "attempt.json",
        log_path=output_dir / "logs" / "pangenie.log",
        threads=2,
        memory_mb=8192,
        execution_purpose="formal",
        timeout_seconds=PANGENIE_TIMEOUT_SECONDS,
        random_seed=0,
    )
    _assert_formal_execution_result(result, mode="end_to_end_from_reads")


@VG_EXECUTION_ONLY
def test_vg_formal_small_region_execution(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Run the real GraphAligner plus vg adapter once on regional long reads."""

    assert sys.platform.startswith("linux")
    for variable in ("PGBENCH_VG_BIN", "PGBENCH_GRAPHALIGNER_BIN"):
        assert os.environ.get(variable), f"set {variable} to an absolute executable"
    _activate_single_conda_runtime(
        monkeypatch,
        [
            ("PGBENCH_VG_BIN", "vg"),
            ("PGBENCH_GRAPHALIGNER_BIN", "GraphAligner"),
        ],
    )
    _assert_bwrap_isolates(_resolve_executable("PGBENCH_BWRAP_BIN", "bwrap"))
    region, resources = _validated_region_resources()
    monkeypatch.syspath_prepend(str(ROOT / "workflow" / "scripts"))
    from pgbench_exec import execute_tool

    output_dir = tmp_path / "vg"
    output_dir.mkdir()
    result = execute_tool(
        tool_manifest_path=ROOT / "plugins" / "vg" / "tool.yaml",
        schema_path=ROOT / "workflow" / "schemas" / "tool.schema.yaml",
        mode="end_to_end_from_reads",
        run_id=f"real-region-{region['id']}",
        sample_id="HG002",
        supplied_inputs={
            "canonical_fastq": resources["long_fastq"],
            "reference": resources["reference_fasta"],
            "pangenome_manifest": resources["pangenome_manifest"],
            "pangenome_panel": resources["population_vcf"],
            "candidate_panel": resources["candidate_vcf"],
            "graph_assets": resources["graph_manifest"].parent,
        },
        output_dir=output_dir,
        resolved_inputs_path=output_dir / "meta" / "resolved-inputs.json",
        attempt_record_path=output_dir / "meta" / "attempt.json",
        log_path=output_dir / "logs" / "vg.log",
        threads=2,
        memory_mb=8192,
        execution_purpose="formal",
        timeout_seconds=VG_TIMEOUT_SECONDS,
        random_seed=0,
    )
    _assert_formal_execution_result(result, mode="end_to_end_from_reads")
