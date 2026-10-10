#!/usr/bin/env python3
"""Copy verified VCFdist v1 caches to stable query-content v2 keys.

The original key included a bcftools command/date header in the compressed
query.  This one-time migration keeps native artifacts and their logs intact.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import tempfile
from pathlib import Path

from run_formal_evaluator import sha256_file, vcf_records, write_variant_query


def canonical_hash(payload: dict) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def load_legacy(cache: Path) -> dict:
    manifest_path = cache / "native-run.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != 1 or canonical_hash(manifest) != cache.name:
        raise ValueError(f"legacy manifest/key mismatch: {cache}")
    for path in (
        cache / "evaluator.log",
        cache / "artifacts" / "query.tsv",
        cache / "artifacts" / "summary.vcf",
    ):
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"incomplete native cache: {path}")
    return manifest


def fixed_fields(manifest: dict) -> dict:
    return {
        key: value
        for key, value in manifest.items()
        if key not in {"schema_version", "prepared_query_sha256"}
    }


def mapped_inputs(root: Path) -> dict[str, Path]:
    matches: dict[str, Path] = {}
    for prepared in root.glob(".vcfdist.*/input/**/*.vcf.gz"):
        matches[sha256_file(prepared)] = prepared
    for prepared in root.glob(".vcfdist.*/input/*.vcf.gz"):
        matches[sha256_file(prepared)] = prepared
    return matches


def verified_log_token(cache: Path, root: Path) -> str:
    first = (cache / "evaluator.log").open(encoding="utf-8", errors="replace").readline()
    pattern = re.escape(str(root)) + r"/\.vcfdist\.[^/]+/input/isolated-unresolved-v1/([0-9a-f]{20})\.vcf\.gz"
    match = re.search(pattern, first)
    if not match:
        raise ValueError(f"isolated cache log does not bind this result root: {cache}")
    return match.group(1)


def migrate_one(
    cache: Path,
    manifest: dict,
    plain: Path,
    *,
    verification: str,
    apply: bool,
) -> dict:
    stable = fixed_fields(manifest)
    stable["schema_version"] = 2
    stable["query_content_sha256"] = sha256_file(plain)
    key = canonical_hash(stable)
    destination = cache.parent / key
    record = {
        "source_cache": str(cache),
        "destination_cache": str(destination),
        "query_content_sha256": stable["query_content_sha256"],
        "verification": verification,
        "existing": destination.is_dir(),
    }
    if not apply or destination.is_dir():
        return record
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{key}.migrating-", dir=destination.parent))
    shutil.copytree(cache / "artifacts", temporary / "artifacts")
    shutil.copy2(cache / "evaluator.log", temporary / "evaluator.log")
    migrated_manifest = {
        **stable,
        "prepared_query_sha256": manifest["prepared_query_sha256"],
        "migrated_from_cache_key": cache.name,
        "migration_verification": verification,
    }
    (temporary / "native-run.json").write_text(
        json.dumps(migrated_manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(destination)
    record["copied"] = True
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--source-query", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--bulk-plain", required=True, type=Path)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()

    # Results may be symlinked into acceptance storage; native logs retain the
    # path spelling used by the evaluator command.
    root = args.output_dir.parent.absolute()
    cache_root = root / ".vcfdist.native-cache"
    inputs = mapped_inputs(root)
    # The native summary VCF embeds command paths, but query.tsv is the vote
    # evidence.  Repeated v1 caches for one event must agree on that evidence.
    query_hashes: dict[str, str] = {}
    for cache in (cache_root / "isolated-unresolved-v1").iterdir():
        if not cache.is_dir() or not (cache / "native-run.json").is_file():
            continue
        token = verified_log_token(cache, root)
        digest = sha256_file(cache / "artifacts" / "query.tsv")
        previous = query_hashes.setdefault(token, digest)
        if previous != digest:
            raise ValueError(f"conflicting native query evidence for {token}")
    bulk_sha = sha256_file(args.bulk_plain)
    candidates = []
    for cache in cache_root.iterdir():
        if not cache.is_dir() or cache.name == "isolated-unresolved-v1":
            continue
        manifest = load_legacy(cache)
        prepared = inputs.get(manifest["prepared_query_sha256"])
        if prepared is not None and sha256_file(prepared.with_suffix("")) == bulk_sha:
            candidates.append((cache, manifest, prepared))
    if not candidates:
        raise ValueError("no v1 bulk cache has a verified matching plaintext query")
    if len({sha256_file(cache / "artifacts" / "query.tsv") for cache, _, _ in candidates}) != 1:
        raise ValueError("conflicting bulk native query evidence")
    candidates.sort(key=lambda item: (item[0] / "evaluator.log").stat().st_mtime)
    bulk_cache, bulk_manifest, bulk_prepared = candidates[-1]
    records = [
        migrate_one(
            bulk_cache,
            bulk_manifest,
            args.bulk_plain,
            verification=f"legacy_gzip_sha256_and_plaintext_sha256:{bulk_prepared}",
            apply=args.apply,
        )
    ]

    source_ids = {fields[2] for fields in vcf_records(args.source_query)}
    token_to_id = {
        hashlib.sha256(result_id.encode()).hexdigest()[:20]: result_id
        for result_id in source_ids
    }
    if len(token_to_id) != len(source_ids):
        raise ValueError("source query has a colliding recovery token")
    with tempfile.TemporaryDirectory(prefix="vcfdist-cache-migration-") as scratch:
        for cache in sorted((cache_root / "isolated-unresolved-v1").iterdir()):
            if not cache.is_dir() or not (cache / "native-run.json").is_file():
                continue
            manifest = load_legacy(cache)
            if fixed_fields(manifest) != fixed_fields(bulk_manifest):
                raise ValueError(f"isolated cache parameters differ from bulk: {cache}")
            token = verified_log_token(cache, root)
            result_id = token_to_id.get(token)
            if result_id is None:
                raise ValueError(f"isolated cache token is absent from query: {token}")
            plain = Path(scratch) / f"{token}.vcf"
            write_variant_query(
                args.source_query,
                plain,
                {result_id},
                reference=args.reference,
                phase_unphased_genotypes=True,
            )
            if len(vcf_records(plain)) != 1:
                raise ValueError(f"single-event input has unexpected count: {result_id}")
            prepared = inputs.get(manifest["prepared_query_sha256"])
            verification = "frozen_source_query_and_log_event_token"
            if prepared is not None:
                if sha256_file(prepared.with_suffix("")) != sha256_file(plain):
                    raise ValueError(f"legacy isolated plaintext differs: {result_id}")
                verification = f"legacy_gzip_sha256_and_plaintext_sha256:{prepared}"
            record = migrate_one(
                cache,
                manifest,
                plain,
                verification=verification,
                apply=args.apply,
            )
            record["result_id"] = result_id
            records.append(record)

    report = {
        "schema_version": 1,
        "output_dir": str(args.output_dir),
        "source_query_sha256": sha256_file(args.source_query),
        "bulk_plain_sha256": bulk_sha,
        "applied": args.apply,
        "records": records,
    }
    if args.apply and args.report is not None:
        if args.report.exists():
            raise FileExistsError(args.report)
        args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"applied": args.apply, "bulk": 1, "isolated": len(records) - 1, "records": records}, indent=2))


if __name__ == "__main__":
    main()
