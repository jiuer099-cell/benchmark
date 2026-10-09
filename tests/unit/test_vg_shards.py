import gzip
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

runner_path = Path(__file__).with_name("run.py")
if not runner_path.is_file():
    runner_path = Path(__file__).resolve().parents[2] / "plugins" / "vg" / "run.py"
spec = importlib.util.spec_from_file_location("vg_runner", runner_path)
assert spec and spec.loader
run = importlib.util.module_from_spec(spec)
spec.loader.exec_module(run)


class ShardCheckpointTests(unittest.TestCase):
    def test_split_keeps_pairs_and_rejects_changed_shard(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            inputs = [root / "reads.R1.fastq.gz", root / "reads.R2.fastq.gz"]
            for mate, path in enumerate(inputs, 1):
                with gzip.open(path, "wb") as handle:
                    for index in range(24):
                        handle.write(
                            f"@read{index} {mate}:N:0:TEST\nACGT\n+\n!!!!\n".encode()
                        )
            split_root = root / "shards"
            parts = run.split_pairs(inputs[0], inputs[1], split_root, "inputs-v1")
            self.assertEqual(len(parts), run.SHARD_COUNT)
            self.assertEqual(
                [item["pairs"] for item in json.loads((split_root / "manifest.json").read_text())["parts"]],
                [3] * run.SHARD_COUNT,
            )
            self.assertIsNotNone(run.valid_shards(split_root, "inputs-v1"))
            self.assertIsNone(run.valid_shards(split_root, "different-inputs"))
            with parts[0][0].open("ab") as handle:
                handle.write(b"damage")
            self.assertIsNone(run.valid_shards(split_root, "inputs-v1"))

    def test_gam_checkpoint_requires_matching_identity_and_digest(self):
        with tempfile.TemporaryDirectory() as directory:
            work = Path(directory) / "work"
            checkpoint = work / "checkpoints"
            checkpoint.mkdir(parents=True)
            gam = checkpoint / "part-000.gam"
            gam.write_bytes(b"valid gam")
            run.write_json_atomic(checkpoint / "part-000.json", {
                "mapping_key": "mapping-v1",
                "index": 0,
                "size": gam.stat().st_size,
                "sha256": run.sha256_file(gam),
            })
            self.assertEqual(run.cached_gam([work], 0, "mapping-v1"), gam)
            self.assertIsNone(run.cached_gam([work], 0, "mapping-v2"))
            gam.write_bytes(b"wrong gam")
            self.assertIsNone(run.cached_gam([work], 0, "mapping-v1"))

    def test_final_pack_checkpoint_requires_matching_digest(self):
        with tempfile.TemporaryDirectory() as directory:
            work = Path(directory) / "work"
            checkpoint = work / "final-pack"
            checkpoint.mkdir(parents=True)
            pack = checkpoint / "HG002.pack"
            pack.write_bytes(b"valid pack")
            run.write_json_atomic(checkpoint / "HG002.pack.json", {
                "pack_key": "pack-v1",
                "size": pack.stat().st_size,
                "sha256": run.sha256_file(pack),
            })
            self.assertEqual(run.cached_final_pack([work], "pack-v1"), pack)
            pack.write_bytes(b"wrong pack")
            self.assertIsNone(run.cached_final_pack([work], "pack-v1"))


if __name__ == "__main__":
    unittest.main()
