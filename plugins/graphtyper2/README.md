# GraphTyper2 external adapter

This external plugin evaluates GraphTyper2 2.7.7 from the Core-managed frozen
RG-v2 shared BAM/BAI. It does not map FASTQs itself and receives only that
alignment pair, the reference/index, and the benchmark-owned blinded SV
candidate panel. It projects GraphTyper2 calls back onto immutable candidate
IDs; missing or conflicting candidates are emitted as explicit no-calls.

The reference FASTA and its adjacent `.fai` are both declared as immutable
sandbox inputs because GraphTyper2 opens the index by appending `.fai` to the
FASTA path.

The plugin runs in the formal network-disabled `bwrap` sandbox and never receives
the hidden HG002 truth ledger. The shared alignment remains Core-owned provenance,
not a GraphTyper2 mapping stage.

For a single input BAM, GraphTyper2 cannot usefully apply all requested threads
to SAM reading. The adapter therefore distributes the frozen, non-overlapping
10 Mb region list round-robin across at most `PGBENCH_THREADS` independent
single-thread GraphTyper2 processes. Each process writes to its own shard
directory; the adapter then performs one deterministic candidate-ID projection.
The aggregate worker count remains bounded by the benchmark resource contract.

The narrowly scoped immutable-output recovery utility is adapter-owned as
`recover_projection.py`; it is not part of Core orchestration or scoring.
