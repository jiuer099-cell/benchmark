#!/usr/bin/env python3
"""BayesTyper 1.5 fixed-panel adapter; declared inputs and native defaults."""
from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from projection import digest, project
from candidates import prepare_candidates
from stability import analyze_file, RHAT_FLAG


def required(name: str) -> str:
    value = os.environ.get(name)
    if not value:
        raise RuntimeError(f"required environment variable is absent: {name}")
    return value


def run(command: list[str], *, cwd: Path | None = None) -> None:
    print("+ " + " ".join(command), flush=True)
    subprocess.run(command, check=True, cwd=cwd)


def kmc_command(binary: Path, reads: Path, prefix: Path, tmp: Path,
                threads: int, memory_mb: int) -> list[str]:
    if threads < 1 or memory_mb < 3072:
        raise ValueError("BayesTyper requires positive threads and at least 3072 MB")
    # BayesTyper 1.5 stores sample counts as unsigned bytes (maximum 255).
    # Saturate high counts; filtering them out with -cx would change the k-mer set.
    return [str(binary), "-k55", "-ci1", "-cs255", "-fq", "-sm",
            f"-t{threads}", f"-m{(memory_mb - 1024) // 1024}",
            f"@{reads}", str(prefix), str(tmp)]


def validate_sample(sample: str) -> None:
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", sample):
        raise ValueError("sample ID must be a safe single filename/TSV token")


def validate_reference_index(reference: Path, index: Path) -> None:
    if index.resolve() != Path(str(reference) + ".fai").resolve() or not index.is_file():
        raise ValueError("reference_index must be the declared .fai adjacent to the reference FASTA")


def declared_inputs() -> dict[str, Path]:
    payload = json.loads(Path(required("PGBENCH_RESOLVED_INPUTS")).read_text())
    result = {}
    for entry in payload["inputs"]:
        name = entry["contract_name"]
        if name in result:
            raise ValueError(f"duplicate resolved input: {name}")
        path = Path(required(entry["environment_variable"])).resolve()
        if not path.exists():
            raise ValueError(f"missing Core-resolved input: {name}")
        result[name] = path
    return result


def genotype_summary(path: Path) -> dict[str, int]:
    counts = {"total": 0, "called": 0, "no_call": 0}
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("#") or not line.strip():
                continue
            row = line.rstrip("\r\n").split("\t")
            if len(row) != 10:
                raise ValueError("Malformed projected single-sample VCF")
            values = dict(zip(row[8].split(":"), row[9].split(":")))
            gt = values.get("GT")
            if gt not in {"0/0", "0/1", "1/1", "./."}:
                raise ValueError("Unexpected projected diploid dosage")
            counts["total"] += 1
            counts["no_call" if gt == "./." else "called"] += 1
    return counts


def genotype_command(binary, unit, work, common, prefix, *, noise_genotyping=False):
    command = [str(binary), "genotype", "-v", str(unit), "-c",
               str(work / "bayestyper_cluster_data"), *common, "-z", "-o", str(prefix)]
    if noise_genotyping:
        command.append("--noise-genotyping")
    return command


def main(*, noise_genotyping: bool = False) -> int:
    inputs = declared_inputs()
    r1, r2 = inputs["short_fastq_r1"], inputs["short_fastq_r2"]
    panel, reference = inputs["candidate_panel"], inputs["reference"]
    validate_reference_index(reference, inputs["reference_index"])
    index, binaries = inputs["tool_index"], inputs["adapter_asset.binaries"]
    sample = required("PGBENCH_SAMPLE_ID")
    validate_sample(sample)
    threads, memory = int(required("PGBENCH_THREADS")), int(required("PGBENCH_MEMORY_MB"))
    seed = str(int(required("PGBENCH_RANDOM_SEED")))
    output_dir = Path(required("PGBENCH_OUTPUT_DIR")).resolve()
    output = Path(required("PGBENCH_OUTPUT_VCF")).resolve()
    output.relative_to(output_dir)
    if output.exists():
        raise RuntimeError("refusing to overwrite an existing adapter output")
    canon, decoy = index / "canon.fa", index / "decoy.fa"
    tools = {name: binaries / name for name in ("bayesTyper", "bayesTyperTools", "kmc")}
    for path in [canon, decoy, *tools.values()]:
        if not path.is_file() or path.is_symlink():
            raise RuntimeError(f"input asset must be a self-contained regular file: {path}")
    work = output_dir / "bayestyper-work"
    work.mkdir(parents=True, exist_ok=False)
    candidates = work / "candidates.vcf.gz"
    preparation = prepare_candidates(panel, candidates)
    kmers, temporary = work / "kmers", work / "kmc-tmp"
    kmers.mkdir(); temporary.mkdir()
    read_list = work / "reads.list"
    read_list.write_text(f"{r1}\n{r2}\n")
    prefix = kmers / sample
    run(kmc_command(tools["kmc"], read_list, prefix, temporary, threads, memory))
    run([str(tools["bayesTyperTools"]), "makeBloom", "-k", str(prefix), "-p", str(threads)])
    sample_file = work / "samples.tsv"
    # Current benchmark is autosomal/diploid: M/F produce identical ploidy.
    # The upstream parser requires M/Male or F/Female, not lower-case male.
    sample_file.write_text(f"{sample}\tM\t{prefix}\n")
    common = ["-s", str(sample_file), "-g", str(canon), "-d", str(decoy),
              "-p", str(threads), "-r", seed]
    run([str(tools["bayesTyper"]), "cluster", "-v", str(candidates), *common], cwd=work)
    units = sorted(work.glob("bayestyper_unit_*/variant_clusters.bin"),
                   key=lambda p: int(p.parent.name.rsplit("_", 1)[1]))
    if not units:
        raise RuntimeError("BayesTyper cluster produced no variant units")
    vcfs = []
    noise_artifacts = {}
    noise_diagnostics = []
    policy = ("official_v1.5_documented_joint_noise_fixed_panel" if noise_genotyping
              else "official_v1.5_defaults_seeded_fixed_panel")
    for unit in units:
        prefix = unit.parent / "calls"
        run(genotype_command(tools["bayesTyper"], unit, work, common, prefix,
                             noise_genotyping=noise_genotyping), cwd=work)
        vcfs.append(Path(f"{prefix}.vcf.gz"))
        noise = unit.parent / "calls_noise_parameters.txt"
        diagnostic = analyze_file(noise, noise_genotyping=noise_genotyping)
        if set(diagnostic.get("samples", {})) != {sample}:
            raise RuntimeError("Incomplete noise output or unexpected noise sample")
        noise_diagnostics.append(diagnostic)
        noise_artifacts["noise_parameters_" + unit.parent.name] = noise
    native = work / "native.vcf.gz"
    run(["bcftools", "concat", "-O", "z", "-o", str(native), *map(str, vcfs)])
    counts = project(native, panel, output, sample, reference=reference,
                     native_policy=policy)
    gt_summary = genotype_summary(output)
    noise_status = "REVIEW_REQUIRED" if any(d["status"] == "REVIEW_REQUIRED" for d in noise_diagnostics) else "NO_SCREEN_FLAG"
    noise_report = output.with_suffix(".noise-stability.json")
    noise_report.write_text(json.dumps({
        "schema": "pgbench-bayestyper-native-noise-screen-v1", "status": noise_status,
        "noise_genotyping": noise_genotyping, "unit_diagnostics": noise_diagnostics,
        "diagnostic_program_sha256": digest(Path(__file__).with_name("stability.py")),
        "local_review_threshold": RHAT_FLAG, "threshold_is_upstream_requirement": False,
        "limitations": "Classical split R-hat only; no ESS, full convergence or genotype accuracy proof.",
    }, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    artifacts = {"native_vcf": native, "samples": sample_file, "native_candidates": candidates,
                 "noise_stability_diagnostics": noise_report, **noise_artifacts}
    output.with_suffix(".native-audit.json").write_text(json.dumps({
        "contract": "pgbench-native-output-audit-v1", "tool_id": "bayestyper",
        "truth_used": False, "execution_policy": policy,
        "output_sha256": digest(output),
        "artifacts": {k: {"path": str(p.resolve()), "sha256": digest(p)} for k,p in artifacts.items()},
        "parameter_policy": {"k": 55, "ci": 1, "cs": 255, "seed": int(seed),
                             "native_genotype_filters": "unchanged", "population_prior_added": False,
                             "noise_genotyping": noise_genotyping},
        "candidate_preparation": preparation,
        "genotype_summary": gt_summary, "noise_stability": noise_status,
        "parameter_estimation_caveat": "v1.5 accepts few or no SNVs and uses non-nested variant types for sequential noise estimation; fewer than 100000 eligible variants triggers a warning, not rejection. Joint mode also uses nested variants. No external variation prior is added.",
    }, indent=2) + "\n")
    print(f"BayesTyper projection: matched={counts[0]} unmapped_candidates={counts[1]} "
          f"outside={counts[2]} gt_called={gt_summary['called']} gt_no_call={gt_summary['no_call']} "
          f"noise_stability={noise_status}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
