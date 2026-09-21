# Paragraph adapter

This external adapter consumes only the benchmark-supplied, frozen RG-v2
Illumina BAM/BAI plus the reference and canonical candidate panel. It does not
receive FASTQ files, does not map reads, and does not receive pangenome assets.
The target truth, target assembly, relatives' genotypes, and undeclared
callsets are not exposed to the process.

The BAM/BAI pair is a Core-managed immutable shared-alignment input. Its lock
records the source FASTQ and reference identities, while Paragraph's own
resolved-input manifest freezes the BAM/BAI content hashes. This is a generic
short-read adapter evidence option; it is not a Paragraph-specific Core path.

The adapter preserves every canonical candidate. Any site missing from the
native output is subsequently materialized as `./.` with an explicit output
status; it is never inferred to be `0/0`.

## Chunked genotyping (memory safety)

`multigrmpy.py` builds one pangenome graph per invocation and its peak RSS
scales with the number of graph sites it receives. Genotyping the whole
18,164-candidate panel in one invocation peaked at ~389 GB and was killed by
the kernel (2026-09-21 incident). The adapter therefore:

1. splits the candidate panel into disjoint chunks (one contig per chunk,
   further split by a candidate-count budget), and
2. runs the chunks with bounded concurrency.

Chunking never changes results: chunks are disjoint by contig and the final
projection deterministically covers every canonical candidate. Knobs:

- `PGBENCH_PARAGRAPH_MAX_CANDIDATES_PER_CHUNK` (default 1000) — budget per
  chunk; lower it if a chr1-scale chunk still exhausts memory.
- `PGBENCH_PARAGRAPH_CHUNK_PARALLELISM` (default 2, capped at
  `PGBENCH_THREADS`) — concurrent multigrmpy invocations; per-chunk threads
  are `PGBENCH_THREADS / parallelism`.
- `PGBENCH_PARAGRAPH_CONTIGS` (e.g. `chr22`) — pilot mode: genotype only the
  listed contigs. The run writes `paragraph-pilot.json`, stamps a
  `##PGBENCH_Paragraph_PilotContigs=` header, and warns that the output must
  not be used for a formal score. Use a pilot to measure real peak RSS before
  choosing the production chunk budget.
