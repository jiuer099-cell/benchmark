"""Exercise isolated VCFdist scheduling without running native VCFdist."""

from __future__ import annotations

import hashlib
import shutil
import sys
import tempfile
import threading
import time
import unittest
from pathlib import Path
from unittest.mock import patch

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

import run_formal_evaluator as evaluator  # noqa: E402


class IsolatedVcfdistParallelRecoveryTest(unittest.TestCase):
    def test_cache_reuse_parallelism_and_stable_order(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            query = root / "source.vcf"
            ids = ["event-a", "event-b", "event-c"]
            query.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                + "".join(
                    f"chr1\t{index + 1}\t{event_id}\tA\tT\t.\tPASS\t.\n"
                    for index, event_id in enumerate(ids)
                ),
                encoding="utf-8",
            )
            tokens = {
                event_id: hashlib.sha256(event_id.encode()).hexdigest()[:20]
                for event_id in ids
            }
            by_token = {token: event_id for event_id, token in tokens.items()}
            output_dir = root / "evaluation" / "vcfdist"
            cache = (
                output_dir.parent
                / ".vcfdist.native-cache"
                / "isolated-unresolved-v1"
                / tokens["event-a"]
            )
            (cache / "artifacts").mkdir(parents=True)
            (cache / "native-run.json").write_text("{}", encoding="utf-8")
            (cache / "evaluator.log").write_text("cached", encoding="utf-8")
            (cache / "artifacts" / "summary.vcf").write_text("cached", encoding="utf-8")

            active = 0
            peak = 0
            lock = threading.Lock()
            native_calls: list[str] = []

            def write_variant_query(source, destination, wanted, **kwargs):
                header = [line for line in query.read_text().splitlines(True) if line.startswith("#")]
                record = [
                    line for line in query.read_text().splitlines(True)
                    if not line.startswith("#") and line.split("\t")[2] in wanted
                ]
                destination.write_text("".join(header + record), encoding="utf-8")

            def run(command, **kwargs):
                nonlocal active, peak
                if command[0] == "bcftools":
                    shutil.copyfile(command[-1], command[command.index("-o") + 1])
                elif command[0] == "vcfdist":
                    self.assertEqual(command[command.index("--max-threads") + 1], "4")
                    destination = Path(command[command.index("-p") + 1])
                    destination.mkdir(parents=True)
                    with lock:
                        active += 1
                        peak = max(peak, active)
                        native_calls.append(destination.name)
                    time.sleep(0.1)
                    (destination / "summary.vcf").write_text("new", encoding="utf-8")
                    with lock:
                        active -= 1
                return None

            def cache_key(**kwargs):
                token = kwargs["prepared_query"].name.removesuffix(".vcf.gz")
                return token, {}

            def parse(prepared, artifacts, **kwargs):
                event_id = by_token[prepared.name.removesuffix(".vcf.gz")]
                votes = {} if event_id == "event-c" else {event_id: event_id == "event-a"}
                return [event_id], votes

            with (
                patch.object(evaluator, "write_variant_query", side_effect=write_variant_query),
                patch.object(evaluator.subprocess, "run", side_effect=run),
                patch.object(evaluator, "native_evaluator_cache_key", side_effect=cache_key),
                patch.object(evaluator, "parse_vcfdist_isolated_event", side_effect=parse),
            ):
                recovered, details = evaluator.recover_vcfdist_unresolved_events(
                    source_query=query,
                    unresolved_ids=set(ids),
                    work=root / "work",
                    artifacts=root / "artifacts",
                    output_dir=output_dir,
                    bcftools_prefix=["bcftools"],
                    tabix_prefix=["tabix"],
                    vcfdist_prefix=["vcfdist"],
                    vcfdist_extra_args=[],
                    truth=root / "truth.vcf",
                    reference=root / "reference.fa",
                    regions=root / "regions.bed",
                    threads=4,
                    version_sha256="version",
                    credit_threshold=0.5,
                    isolated_workers=2,
                )

            self.assertEqual(recovered, {"event-a": True, "event-b": False})
            self.assertEqual([item["result_id"] for item in details], ids)
            self.assertEqual([item["native_cache_hit"] for item in details], [True, False, False])
            self.assertEqual(set(native_calls), {tokens["event-b"], tokens["event-c"]})
            self.assertEqual(peak, 2)


if __name__ == "__main__":
    unittest.main()
