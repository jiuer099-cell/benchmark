# vg Giraffe graph-calling adapter

This end-to-end adapter maps the registered paired HG002 Illumina FASTQs to
the exact frozen GBZ graph with `vg giraffe`, accumulates support with
`vg pack`, computes snarls, and emits graph-supported variant sites with
`vg call`.

The adapter divides paired FASTQs into eight deterministic, pair-preserving
round-robin shards. Each shard creates a complete GAM checkpoint, so an
interrupted attempt loses at most the unfinished shard instead of the entire
sample's mapping. Completed shard FASTQs, GAMs, and the final pack are reused
from the adapter's read-only attempt history only after their input identity,
sizes, and SHA-256 checks pass. Concatenated GAM streams are fed into one
`vg pack`, followed by one whole-sample `vg call`; this avoids loading the
whole-genome graph separately for every shard's pack. The paired-fragment
parameters are frozen to the full-input Giraffe estimate (596.974 +/- 153.29)
so every shard uses the same mapping distribution. With 32 assigned threads,
Giraffe gets 32 and the final pack gets 4. The adapter waits for 240 GiB
available memory before mapping or packing and abandons only the current stage
if available memory falls below 96 GiB.

The `vg_giraffe_shortread` profile freezes `graph.gbz`, `graph.dist`,
`graph.shortread.withzip.min`, and `graph.shortread.zipcodes` together with
the graph manifest and HG002-free sample list. Explicit paths prevent Giraffe
from silently building or selecting mutable indexes at execution time. GAMs,
the final pack, shard FASTQs, and snarls files are attempt-local intermediates.

Formal attempts run without network access in the configured `bwrap` sandbox.
The runner receives no truth asset and emits a fresh VCF for the common formal
evaluators.
