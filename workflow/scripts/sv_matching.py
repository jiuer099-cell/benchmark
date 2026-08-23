#!/usr/bin/env python3
"""Frozen, deterministic structural-variant matching and genotype semantics."""

from __future__ import annotations

import bisect
import gzip
import hashlib
import heapq
import json
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path
from typing import Any, TextIO

import yaml  # type: ignore[import-untyped]

try:
    import edlib  # type: ignore[import-not-found]
except ImportError:  # Windows development hosts may lack a compiled wheel.
    edlib = None  # type: ignore[assignment]

from build_pangenome_manifest import default_variant_end, infer_svtype, parse_info


class SvMatchError(ValueError):
    """Raised when VCF records or the frozen matching profile are invalid."""


@dataclass(frozen=True)
class SvRecord:
    record_id: str
    chrom: str
    pos: int
    end: int
    svtype: str
    svlen: int
    ref: str
    alt: str
    gt: str
    filter_status: str = "PASS"

    @property
    def stable_key(self) -> tuple[str, int, int, str, int, str, str]:
        return (
            self.chrom,
            self.pos,
            self.end,
            self.svtype,
            self.svlen,
            self.ref,
            self.alt,
        )


@dataclass(frozen=True)
class SvMatch:
    query_index: int
    truth_index: int
    score: float
    start_distance: int
    end_distance: int
    size_similarity: float
    sequence_similarity: float | None


@dataclass
class _ResidualEdge:
    """Mutable unit-capacity edge used by the internal min-cost-flow solver."""

    to: int
    reverse_index: int
    capacity: int
    cost: Decimal
    match: SvMatch | None = None


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_evaluator_profile(path: Path) -> dict[str, Any]:
    try:
        loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise SvMatchError(f"cannot load evaluator profile {path}: {exc}") from exc
    if not isinstance(loaded, dict) or loaded.get("schema_version") != 1:
        raise SvMatchError("evaluator profile schema_version must be 1")
    for section in (
        "profile",
        "universe",
        "sv_match",
        "semantics",
        "evaluators",
        "bootstrap",
    ):
        if not isinstance(loaded.get(section), dict):
            raise SvMatchError(f"evaluator profile is missing {section}")
    universe = loaded["universe"]
    allowed_svtypes = universe.get("allowed_svtypes")
    if (
        not isinstance(allowed_svtypes, list)
        or not allowed_svtypes
        or any(not isinstance(value, str) or not value for value in allowed_svtypes)
    ):
        raise SvMatchError("universe.allowed_svtypes must be a non-empty string list")
    if set(allowed_svtypes) != {"DEL", "INS"}:
        raise SvMatchError("formal universe currently supports exactly DEL and INS")
    minimum_sv_size = universe.get("minimum_sv_size")
    if not isinstance(minimum_sv_size, int) or minimum_sv_size < 1:
        raise SvMatchError("universe.minimum_sv_size must be a positive integer")
    maximum_sv_size = universe.get("maximum_sv_size")
    if (
        not isinstance(maximum_sv_size, int)
        or maximum_sv_size < minimum_sv_size
    ):
        raise SvMatchError(
            "universe.maximum_sv_size must be an integer greater than or "
            "equal to minimum_sv_size"
        )
    if universe.get("truth_filter_policy") != "pass_only":
        raise SvMatchError("formal truth_filter_policy must be pass_only")
    if universe.get("query_filter_policy") != "pass_or_unfiltered":
        raise SvMatchError(
            "formal query_filter_policy must be pass_or_unfiltered"
        )
    if universe.get("biallelic_only") is not True:
        raise SvMatchError("formal universe requires biallelic_only: true")
    if universe.get("multiallelic_policy") != "exclude":
        raise SvMatchError("formal multiallelic_policy must be exclude")
    if universe.get("region_policy") != "fully_contained":
        raise SvMatchError("formal region_policy must be fully_contained")
    match = loaded["sv_match"]
    if match.get("one_to_one") is not True:
        raise SvMatchError("formal matching requires one_to_one: true")
    if int(match.get("max_breakpoint_distance", 0)) < 0:
        raise SvMatchError("max_breakpoint_distance must be non-negative")
    for key in ("min_size_similarity", "min_sequence_similarity"):
        value = float(match.get(key, -1))
        if not 0.0 <= value <= 1.0:
            raise SvMatchError(f"{key} must be in [0,1]")
    if match.get("sequence_similarity_algorithm") != (
        "truvari_v5.4.0_best_seqsim_edlib_1.3.9.post1"
    ):
        raise SvMatchError("sequence_similarity_algorithm is not frozen")
    if match.get("sequence_roll_policy") != "enabled":
        raise SvMatchError("sequence_roll_policy must be enabled")
    for evaluator, option in (
        ("truvari", "--sizemax"),
        ("vcfdist", "--largest-variant"),
    ):
        evaluator_settings = loaded["evaluators"].get(evaluator)
        if not isinstance(evaluator_settings, dict):
            raise SvMatchError(f"evaluator profile is missing {evaluator}")
        extra_args = evaluator_settings.get("extra_args")
        if not isinstance(extra_args, list) or not all(
            isinstance(item, str) for item in extra_args
        ):
            raise SvMatchError(f"{evaluator}.extra_args must be a string list")
        positions = [
            index for index, value in enumerate(extra_args) if value == option
        ]
        if (
            len(positions) != 1
            or positions[0] + 1 >= len(extra_args)
            or extra_args[positions[0] + 1] != str(maximum_sv_size)
        ):
            raise SvMatchError(
                f"{evaluator}.extra_args must freeze {option} "
                f"{maximum_sv_size}"
            )
    expected_evaluator_fingerprints = {
        "truvari": (
            "6587f29d09e453ae50a398c79804f127"
            "a615e76563f4c179ff7ae45fe750fad4"
        ),
        "aardvark": (
            "fb6fda42c6d7b8f3baabdf027e9bde9"
            "5a0114798673f84e676b1eff9b4c5c2ac"
        ),
        "vcfdist": (
            "714eb9f5def97655ca3cd6704c2f5dc1"
            "46604b2066f872e8f0a81d1239978c53"
        ),
    }
    for evaluator, expected_sha256 in expected_evaluator_fingerprints.items():
        observed = loaded["evaluators"].get(evaluator)
        if not isinstance(observed, dict):
            raise SvMatchError(f"evaluator profile is missing {evaluator}")
        if observed.get("expected_version_sha256") != expected_sha256:
            raise SvMatchError(
                f"{evaluator}.expected_version_sha256 is not frozen"
            )
    vcfdist_args = loaded["evaluators"]["vcfdist"]["extra_args"]
    for option, expected in (
        ("--cluster", "size"),
        ("--max-supercluster-size", "20000"),
    ):
        if option not in vcfdist_args:
            raise SvMatchError(f"vcfdist.extra_args must freeze {option}")
        position = vcfdist_args.index(option)
        if position + 1 >= len(vcfdist_args) or vcfdist_args[position + 1] != expected:
            raise SvMatchError(
                f"vcfdist.extra_args must freeze {option} {expected}"
            )
    bootstrap = loaded["bootstrap"]
    if int(bootstrap.get("replicates", 0)) < 100:
        raise SvMatchError("bootstrap.replicates must be at least 100")
    if float(bootstrap.get("confidence_level", 0)) != 0.95:
        raise SvMatchError("bootstrap confidence_level must be 0.95")
    if bootstrap.get("unit") != "fixed_genomic_block":
        raise SvMatchError("bootstrap.unit must be fixed_genomic_block")
    if bootstrap.get("method") != "paired_genomic_block_bootstrap":
        raise SvMatchError(
            "bootstrap.method must be paired_genomic_block_bootstrap"
        )
    if int(bootstrap.get("block_size_bp", 0)) < 1_000_000:
        raise SvMatchError("bootstrap.block_size_bp must be at least 1 Mb")
    loaded["_sha256"] = sha256_file(path)
    return loaded


def parse_gt(fields: list[str]) -> str:
    if len(fields) < 10:
        return "./."
    keys = fields[8].split(":")
    values = fields[9].split(":")
    if "GT" not in keys:
        return "./."
    index = keys.index("GT")
    return values[index] if index < len(values) else "./."


def gt_state(gt: str) -> str:
    if not gt or gt in {".", "./.", ".|."}:
        return "no_call"
    alleles = gt.replace("|", "/").split("/")
    if not alleles or any(allele in {"", "."} for allele in alleles):
        return "no_call"
    try:
        numeric = [int(allele) for allele in alleles]
    except ValueError:
        return "no_call"
    return "hom_ref" if all(allele == 0 for allele in numeric) else "variant"


def genotypes_equal(query_gt: str, truth_gt: str, *, require_phase: bool) -> bool:
    if gt_state(query_gt) == "no_call" or gt_state(truth_gt) == "no_call":
        return False
    if require_phase:
        return query_gt == truth_gt
    query = sorted(query_gt.replace("|", "/").split("/"))
    truth = sorted(truth_gt.replace("|", "/").split("/"))
    return query == truth


def _integer_info(info: dict[str, str | bool], key: str, default: int) -> int:
    value = info.get(key)
    if value is None or value is True:
        return default
    try:
        return int(str(value).split(",", 1)[0])
    except ValueError as exc:
        raise SvMatchError(f"{key} must be integer, found {value!r}") from exc


def record_from_fields(fields: list[str], *, prefix: str, index: int) -> SvRecord:
    if len(fields) < 8:
        raise SvMatchError("VCF record must contain at least eight columns")
    try:
        pos = int(fields[1])
    except ValueError as exc:
        raise SvMatchError(f"invalid VCF POS {fields[1]!r}") from exc
    info = parse_info(fields[7])
    svtype = infer_svtype(fields[4], info, ref=fields[3])
    end = _integer_info(
        info,
        "END",
        default_variant_end(
            pos=pos,
            ref=fields[3],
            alt=fields[4],
            svtype=svtype,
        ),
    )
    resolved_length = len(fields[4]) - len(fields[3])
    if fields[4].startswith("<") and fields[4].endswith(">"):
        span = max(0, end - pos)
        if svtype == "DEL":
            resolved_length = -span
        elif svtype in {"DUP", "INV"}:
            resolved_length = span
        else:
            # A symbolic insertion without SVLEN has no auditable event size.
            resolved_length = 0
    svlen = _integer_info(info, "SVLEN", resolved_length)
    record_id = fields[2]
    if not record_id or record_id == ".":
        payload = json.dumps(
            [fields[0], pos, end, svtype, svlen, fields[3], fields[4], index],
            separators=(",", ":"),
        )
        record_id = f"{prefix}_{hashlib.sha256(payload.encode()).hexdigest()[:20]}"
    return SvRecord(
        record_id=record_id,
        chrom=fields[0],
        pos=pos,
        end=end,
        svtype=svtype,
        svlen=svlen,
        ref=fields[3],
        alt=fields[4],
        gt=parse_gt(fields),
        filter_status=fields[6],
    )


def record_shape_in_universe(record: SvRecord, profile: dict[str, Any]) -> bool:
    """Return whether an SV has the frozen type, size, and allele shape."""

    universe = profile["universe"]
    return (
        record.svtype in set(universe["allowed_svtypes"])
        and abs(record.svlen) >= int(universe["minimum_sv_size"])
        and abs(record.svlen) <= int(universe["maximum_sv_size"])
        and (not universe["biallelic_only"] or "," not in record.alt)
    )


def query_record_in_universe(record: SvRecord, profile: dict[str, Any]) -> bool:
    """Apply the frozen submitted-query policy without inspecting genotype."""

    if not record_shape_in_universe(record, profile):
        return False
    policy = profile["universe"]["query_filter_policy"]
    if policy == "pass_or_unfiltered":
        return record.filter_status in {"PASS", "."}
    raise SvMatchError(f"unsupported query filter policy: {policy!r}")


def truth_record_in_universe(record: SvRecord, profile: dict[str, Any]) -> bool:
    """Apply the frozen truth policy to a materialized truth record."""

    if not record_shape_in_universe(record, profile):
        return False
    policy = profile["universe"]["truth_filter_policy"]
    if policy == "pass_only":
        return record.filter_status == "PASS"
    raise SvMatchError(f"unsupported truth filter policy: {policy!r}")


def load_vcf(
    path: Path,
    *,
    prefix: str,
    allow_empty: bool = False,
) -> list[SvRecord]:
    records: list[SvRecord] = []
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            records.append(
                record_from_fields(
                    line.rstrip("\n").split("\t"),
                    prefix=prefix,
                    index=len(records),
                )
            )
    if not records and not allow_empty:
        raise SvMatchError(f"VCF contains no records: {path}")
    return records


def _resolved_sequence(record: SvRecord) -> str | None:
    alt = record.alt
    if not alt or alt == "." or alt.startswith("<") or "[" in alt or "]" in alt:
        return None
    if record.svtype == "INS":
        if len(alt) <= len(record.ref):
            return None
        return alt
    if record.svtype == "DEL":
        if len(record.ref) <= len(alt):
            return None
        return record.ref
    return alt


def _edlib_sequence_similarity(left: str, right: str) -> float:
    """Match Truvari 5.4.0's public seqsim implementation exactly."""

    if edlib is None:
        raise SvMatchError(
            "edlib 1.3.9.post1 is required for formal sequence matching"
        )
    left = left.upper()
    right = right.upper()
    total_length = len(left) + len(right)
    if total_length == 0:
        return 1.0
    result = edlib.align(left, right, mode="NW", task="distance")
    edit_distance = result.get("editDistance")
    if not isinstance(edit_distance, int) or edit_distance < 0:
        raise SvMatchError("edlib could not calculate global edit distance")
    return (total_length - edit_distance) / total_length


def _smallest_rotation(sequence: str) -> str:
    doubled = sequence + sequence
    return min(
        doubled[index : index + len(sequence)]
        for index in range(len(sequence))
    )


def _unrolled_sequence_similarity(left: str, right: str, distance: int) -> float:
    offset = distance % len(right)
    unrolled = right[-offset:] + right[:-offset] if offset else right
    return _edlib_sequence_similarity(left, unrolled)


def _best_sequence_similarity(left: str, right: str, distance: int) -> float:
    """Match Truvari 5.4.0 best_seqsim, including default sequence rolling."""

    rolled = 0.0
    if len(left) < 500 and len(right) < 500:
        rolled = _edlib_sequence_similarity(
            _smallest_rotation(left),
            _smallest_rotation(right),
        )
    return max(
        rolled,
        _unrolled_sequence_similarity(left, right, distance),
        _unrolled_sequence_similarity(left, right, -distance),
        _unrolled_sequence_similarity(right, left, distance),
        _unrolled_sequence_similarity(right, left, -distance),
        _edlib_sequence_similarity(left, right),
    )


def _sequence_similarity(query: SvRecord, truth: SvRecord) -> float | None:
    if query.ref == truth.ref and query.alt == truth.alt:
        return 1.0
    query_sequence = _resolved_sequence(query)
    truth_sequence = _resolved_sequence(truth)
    if query_sequence is None or truth_sequence is None:
        return None
    start_distance = query.pos - truth.pos
    end_distance = query.end - truth.end
    if start_distance == 0 or end_distance == 0:
        return _edlib_sequence_similarity(query_sequence, truth_sequence)
    return _best_sequence_similarity(
        query_sequence,
        truth_sequence,
        start_distance,
    )


def compatibility(
    query: SvRecord,
    truth: SvRecord,
    match_profile: dict[str, Any],
) -> SvMatch | None:
    settings = match_profile["sv_match"]
    if query.chrom != truth.chrom:
        return None
    if settings["require_same_svtype"] and query.svtype != truth.svtype:
        return None
    distance = int(settings["max_breakpoint_distance"])
    start_distance = abs(query.pos - truth.pos)
    end_distance = abs(query.end - truth.end)
    if start_distance > distance or end_distance > distance:
        return None
    query_length = max(1, abs(query.svlen))
    truth_length = max(1, abs(truth.svlen))
    size_similarity = min(query_length, truth_length) / max(
        query_length, truth_length
    )
    if size_similarity < float(settings["min_size_similarity"]):
        return None
    sequence_similarity = _sequence_similarity(query, truth)
    if (
        sequence_similarity is not None
        and sequence_similarity < float(settings["min_sequence_similarity"])
    ):
        return None
    if (
        settings.get("sequence_required_for_resolved_insertions", False)
        and query.svtype == "INS"
        and (_resolved_sequence(query) is None) != (_resolved_sequence(truth) is None)
    ):
        return None
    distance_scale = max(1, distance)
    start_score = 1.0 - start_distance / distance_scale
    end_score = 1.0 - end_distance / distance_scale
    sequence_score = sequence_similarity if sequence_similarity is not None else 0.5
    score = (
        0.30 * start_score
        + 0.25 * end_score
        + 0.30 * size_similarity
        + 0.15 * sequence_score
    )
    return SvMatch(
        query_index=-1,
        truth_index=-1,
        score=score,
        start_distance=start_distance,
        end_distance=end_distance,
        size_similarity=size_similarity,
        sequence_similarity=sequence_similarity,
    )


def _record_order_key(
    record: SvRecord,
    original_index: int,
) -> tuple[Any, ...]:
    """Return an input-order-independent key, with index only for exact duplicates."""

    return (
        *record.stable_key,
        record.record_id,
        record.gt,
        original_index,
    )


def _solve_component(
    query_indices: set[int],
    truth_indices: set[int],
    edges: list[SvMatch],
    queries: list[SvRecord],
    truths: list[SvRecord],
) -> dict[int, SvMatch]:
    """Solve one sparse connected component by exact min-cost maximum flow."""

    query_order = sorted(
        query_indices,
        key=lambda index: _record_order_key(queries[index], index),
    )
    truth_order = sorted(
        truth_indices,
        key=lambda index: _record_order_key(truths[index], index),
    )
    query_rank = {original: rank for rank, original in enumerate(query_order)}
    truth_rank = {original: rank for rank, original in enumerate(truth_order)}

    source = 0
    query_offset = 1
    truth_offset = query_offset + len(query_order)
    sink = truth_offset + len(truth_order)
    graph: list[list[_ResidualEdge]] = [[] for _ in range(sink + 1)]

    def add_edge(
        start: int,
        end: int,
        cost: Decimal,
        match: SvMatch | None = None,
    ) -> _ResidualEdge:
        forward = _ResidualEdge(
            to=end,
            reverse_index=len(graph[end]),
            capacity=1,
            cost=cost,
            match=match,
        )
        reverse = _ResidualEdge(
            to=start,
            reverse_index=len(graph[start]),
            capacity=0,
            cost=-cost,
        )
        graph[start].append(forward)
        graph[end].append(reverse)
        return forward

    for rank in range(len(query_order)):
        add_edge(source, query_offset + rank, Decimal(0))
    for rank in range(len(truth_order)):
        add_edge(truth_offset + rank, sink, Decimal(0))

    maximum_score = max(Decimal.from_float(edge.score) for edge in edges)
    tracked_edges: list[tuple[SvMatch, _ResidualEdge]] = []
    for match in sorted(
        edges,
        key=lambda item: (
            query_rank[item.query_index],
            truth_rank[item.truth_index],
        ),
    ):
        # For a fixed cardinality, minimizing max_score - score exactly
        # maximizes the sum of the original binary floating-point scores.
        cost = maximum_score - Decimal.from_float(match.score)
        residual = add_edge(
            query_offset + query_rank[match.query_index],
            truth_offset + truth_rank[match.truth_index],
            cost,
            match,
        )
        tracked_edges.append((match, residual))

    node_count = len(graph)
    potentials = [Decimal(0)] * node_count
    while True:
        distances: list[Decimal | None] = [None] * node_count
        previous: list[tuple[int, int] | None] = [None] * node_count
        distances[source] = Decimal(0)
        queue: list[tuple[Decimal, int]] = [(Decimal(0), source)]

        while queue:
            distance_so_far, node = heapq.heappop(queue)
            if distances[node] != distance_so_far:
                continue
            for edge_index, edge in enumerate(graph[node]):
                if edge.capacity == 0:
                    continue
                candidate = (
                    distance_so_far
                    + edge.cost
                    + potentials[node]
                    - potentials[edge.to]
                )
                known = distances[edge.to]
                if known is None or candidate < known:
                    distances[edge.to] = candidate
                    previous[edge.to] = (node, edge_index)
                    heapq.heappush(queue, (candidate, edge.to))

        if distances[sink] is None:
            break
        for node, distance_value in enumerate(distances):
            if distance_value is not None:
                potentials[node] += distance_value

        node = sink
        while node != source:
            predecessor = previous[node]
            if predecessor is None:  # pragma: no cover - protects solver invariant
                raise SvMatchError("internal matching solver produced a broken path")
            previous_node, edge_index = predecessor
            edge = graph[previous_node][edge_index]
            edge.capacity = 0
            graph[node][edge.reverse_index].capacity = 1
            node = previous_node

    return {
        match.query_index: match
        for match, residual in tracked_edges
        if residual.capacity == 0
    }


def one_to_one_match(
    queries: list[SvRecord],
    truths: list[SvRecord],
    profile: dict[str, Any],
) -> dict[int, SvMatch]:
    """Return an exact, deterministic maximum-cardinality weighted matching.

    Objectives are applied lexicographically: maximize the number of matched
    events, then maximize the total compatibility score. Stable record keys
    and stable residual-graph traversal resolve any remaining ties.
    """

    distance = int(profile["sv_match"]["max_breakpoint_distance"])
    require_same_svtype = bool(profile["sv_match"]["require_same_svtype"])
    truth_index: dict[tuple[str, str], tuple[list[int], list[int]]] = {}
    grouped: dict[tuple[str, str], list[tuple[int, int]]] = {}
    for index, truth in enumerate(truths):
        group_type = truth.svtype if require_same_svtype else ""
        grouped.setdefault((truth.chrom, group_type), []).append((truth.pos, index))
    for key, rows in grouped.items():
        rows.sort(
            key=lambda row: (
                row[0],
                _record_order_key(truths[row[1]], row[1]),
            )
        )
        truth_index[key] = ([row[0] for row in rows], [row[1] for row in rows])

    edges: list[SvMatch] = []
    for query_index, query in enumerate(queries):
        group_type = query.svtype if require_same_svtype else ""
        group = truth_index.get((query.chrom, group_type))
        if group is None:
            continue
        positions, indices = group
        lower = bisect.bisect_left(positions, query.pos - distance)
        upper = bisect.bisect_right(positions, query.pos + distance)
        for truth_index_value in indices[lower:upper]:
            raw = compatibility(query, truths[truth_index_value], profile)
            if raw is not None:
                edges.append(
                    SvMatch(
                        query_index=query_index,
                        truth_index=truth_index_value,
                        score=raw.score,
                        start_distance=raw.start_distance,
                        end_distance=raw.end_distance,
                        size_similarity=raw.size_similarity,
                        sequence_similarity=raw.sequence_similarity,
                    )
                )
    if not edges:
        return {}

    query_edges: dict[int, list[SvMatch]] = {}
    truth_edges: dict[int, list[SvMatch]] = {}
    for match in edges:
        query_edges.setdefault(match.query_index, []).append(match)
        truth_edges.setdefault(match.truth_index, []).append(match)

    visited_queries: set[int] = set()
    visited_truths: set[int] = set()
    assignments: dict[int, SvMatch] = {}
    for seed in sorted(
        query_edges,
        key=lambda index: _record_order_key(queries[index], index),
    ):
        if seed in visited_queries:
            continue
        component_queries: set[int] = set()
        component_truths: set[int] = set()
        stack: list[tuple[bool, int]] = [(True, seed)]
        while stack:
            is_query, index = stack.pop()
            if is_query:
                if index in visited_queries:
                    continue
                visited_queries.add(index)
                component_queries.add(index)
                stack.extend(
                    (False, match.truth_index)
                    for match in query_edges.get(index, [])
                )
            else:
                if index in visited_truths:
                    continue
                visited_truths.add(index)
                component_truths.add(index)
                stack.extend(
                    (True, match.query_index)
                    for match in truth_edges.get(index, [])
                )
        component_edges = [
            match
            for query_index in component_queries
            for match in query_edges[query_index]
            if match.truth_index in component_truths
        ]
        assignments.update(
            _solve_component(
                component_queries,
                component_truths,
                component_edges,
                queries,
                truths,
            )
        )
    return assignments
