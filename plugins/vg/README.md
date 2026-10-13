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

Native projection v2 preserves caller coordinates, REF/ALT, QUAL and FILTER.
Candidate linkage uses normalized sequence or exact graph allele paths, then
independent reference-haplotype equivalence checks. Truth and evaluator results
never participate. Multiallelic GT becomes ALT-specific unphased dosage; full
original FORMAT/sample data remains in the hashed native-record trace, and
GT/GQ/DP remain in the biallelic VCF. Conflicts fail closed. Unmatched sites are
`not_observed`, not proven addressable or reference.

Fresh calling uses pack `-Q 5` and call `-a -t N`, with memory monitoring for
calling too. Pack identities bind the recipe and source GAM digests. Verified
legacy GAMs can be reused only for the identical mapping recipe; old Q0 packs
cannot satisfy Q5 identities. Existing native-VCF replay is labelled legacy-Q0,
not a Q5 rerun. Default nested calling is retained pending separate validation.

VG sealing requires native VCF, full native trace and projection audit hashes.
A seal validates lineage and scoring consistency, not biological performance.
Experimental/single-sample qualifications remain separate from score validity.
# Native VCF validation

Native projection requires a valid `#CHROM` header with the requested single
sample and GT-bearing records. Empty/headerless VCFs and records preceding the
header fail before a projected VCF or successful audit is published. A valid
header-only VCF retains the explicit no-call semantics for unobserved sites.

Called non-sequence ALT alleles (including `*` and symbolic alleles) produce
`./.` rather than numeric dosage. Uncalled placeholders do not mask ordinary
sequence calls. The original native GT, quality and filters remain in the
trace. The pack/call parameters, mapping caches and existing run outputs are
not changed by this validation repair.
