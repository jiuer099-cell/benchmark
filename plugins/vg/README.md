# vg Giraffe graph-calling adapter

This end-to-end adapter maps the registered paired HG002 Illumina FASTQs to
the exact frozen GBZ graph with `vg giraffe`, accumulates support with
`vg pack`, computes snarls, and emits graph-supported variant sites with
`vg call`.

The `vg_giraffe_shortread` profile freezes `graph.gbz`, `graph.dist`,
`graph.shortread.withzip.min`, and `graph.shortread.zipcodes` together with
the graph manifest and HG002-free sample list. Explicit paths prevent Giraffe
from silently building or selecting mutable indexes at execution time. GAM,
pack, and snarls files are run-local intermediates.

Formal attempts run without network access in the configured `bwrap` sandbox.
The runner receives no truth asset and emits a fresh VCF for the common formal
evaluators.
