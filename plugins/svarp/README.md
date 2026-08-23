# SVarp external adapter

This adapter evaluates SVarp 1.2.0 as an end-to-end PacBio CLR pangenome SV
discovery tool.  It converts only the registered input FASTQ to bgzip FASTA,
maps it to a frozen, leave-HG002-out HPRC minigraph rGFA, and runs SVarp on
CLR input. SVarp v1.2.0 has no read-type option, so the adapter does not invent
one. Its `wtdbg2`/`wtpoa-cns` local-assembly dependencies are frozen
in the plugin environment. No truth data, truth-derived panel, or linear-reference BAM is
exposed to the tool sandbox.

SVarp natively emits sequence-resolved local assemblies (`svtigs`), not a VCF.
For an auditable benchmark interface, the adapter aligns those svtigs to the
registered GRCh38 reference with frozen minimap2 2.31 and calls sequence-
resolved insertions/deletions from reference-oriented `cs` strings. The
converter requires one primary linear projection per svtig and verifies every
deleted sequence against the registered reference.
The resulting VCF uses the `variant_sites` contract: record presence is a
detection claim and the genotype is deliberately `./.`.  Consequently the
formal report separates SVarp's detection score from its non-scorable GT/no-
call diagnostics rather than inventing genotypes that SVarp did not produce.

The graph profile is `svarp_minigraph_longread`.  Its source manifest must
explicitly declare the HG002/NA24385 leave-out condition, and the rGFA is
content-locked before execution.  This profile must not be replaced by a GBZ
or a graph from a different HPRC release: the GAF must be generated from the
same rGFA that SVarp receives.

The v1.2.0 command contract explicitly passes support 5, breakpoint distance
100, wtdbg2, map ratio 0.90, precise clipping 0.97, alignment score 5000, and
the configured thread count. Phasing, skip-untagged, no-remap, assembly mode,
and debug remain explicitly disabled in `tool.yaml`. Unsealed cross-attempt
work checkpoints are never reused.
