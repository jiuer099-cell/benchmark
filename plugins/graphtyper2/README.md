# GraphTyper2 external adapter

This external plugin evaluates GraphTyper2 2.7.7 as a complete paired-short-read
pipeline. It maps the registered R1/R2 reads with BWA-MEM2, creates a sorted BAM,
genotypes the benchmark-owned blinded SV candidate panel, and projects GraphTyper2
calls back onto the immutable candidate IDs. Missing or conflicting candidates are
emitted as explicit no-calls.

The reference FASTA and its adjacent `.fai` are both declared as immutable
sandbox inputs because GraphTyper2 opens the index by appending `.fai` to the
FASTA path.

The plugin runs in the formal network-disabled `bwrap` sandbox and never receives
the hidden HG002 truth ledger. Mapping is part of the measured tool-owned middle.

For a single input BAM, GraphTyper2 cannot usefully apply all requested threads
to SAM reading. The adapter therefore distributes the frozen, non-overlapping
10 Mb region list round-robin across at most `PGBENCH_THREADS` independent
single-thread GraphTyper2 processes. Each process writes to its own shard
directory; the adapter then performs one deterministic candidate-ID projection.
The aggregate worker count remains bounded by the benchmark resource contract.
