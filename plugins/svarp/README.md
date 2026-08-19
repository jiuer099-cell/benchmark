# SVarp external adapter

This adapter evaluates SVarp 1.2.0 as an end-to-end PacBio CLR pangenome SV
discovery tool.  It converts only the registered input FASTQ to bgzip FASTA,
maps it to a frozen, leave-HG002-out HPRC minigraph rGFA, and runs SVarp in
`clr` mode.  Its `wtdbg2`/`wtpoa-cns` local-assembly dependencies are frozen
in the plugin environment. No truth data, truth-derived panel, or linear-reference BAM is
exposed to the tool sandbox.

SVarp natively emits sequence-resolved local assemblies (`svtigs`), not a VCF.
For an auditable benchmark interface, the adapter aligns those svtigs to the
registered GRCh38 reference with frozen minimap2 2.31 and calls sequence-
resolved insertions/deletions with the bundled `paftools.js call` procedure.
The resulting VCF uses the `variant_sites` contract: record presence is a
detection claim and the genotype is deliberately `./.`.  Consequently the
formal report separates SVarp's detection score from its non-scorable GT/no-
call diagnostics rather than inventing genotypes that SVarp did not produce.

The graph profile is `svarp_minigraph_longread`.  Its source manifest must
explicitly declare the HG002/NA24385 leave-out condition, and the rGFA is
content-locked before execution.  This profile must not be replaced by a GBZ
or a graph from a different HPRC release: the GAF must be generated from the
same rGFA that SVarp receives.
