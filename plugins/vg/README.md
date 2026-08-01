# GraphAligner + vg graph-calling adapter

This end-to-end adapter exports the exact frozen GBZ graph to GFA, maps the
canonical HG002 PacBio CLR FASTQ with `GraphAligner`, creates a support pack,
computes snarls, and emits graph-supported variant sites with `vg call`.

The modern `vg_gbz_min_dist` profile contains validated `graph.gbz`,
`graph.min`, and `graph.dist` assets plus its manifest/sample list. `graph.xg`
is required only by the explicit legacy `vg_legacy_xg` profile. The exported
GFA, GAM, pack, and snarls are run-local intermediates, so node identity
remains tied to the locked GBZ.

Formal attempts run without network access in the configured `bwrap` sandbox.
