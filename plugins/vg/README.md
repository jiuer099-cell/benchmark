# GraphAligner + vg graph-calling adapter

This end-to-end adapter exports the exact frozen GBZ graph to GFA, maps the
canonical HG002 PacBio CLR FASTQ with `GraphAligner`, creates a support pack,
computes snarls, and emits graph-supported variant sites with `vg call`.

The graph directory must contain the validated `graph.gbz`, `graph.xg`,
`graph.min`, and `graph.dist` bundle. The exported GFA, GAM, pack, and snarls
are run-local intermediates, so node identity remains tied to the locked GBZ.
