# PanGenie adapter

This adapter implements the short-read, k-mer pangenome-genotyping family.
It indexes the benchmark-owned phased, sequence-resolved pangenome panel VCF
with `PanGenie-index`, genotypes one diploid sample with `PanGenie`, and returns
an all-sites VCF to the common normalization and evaluation stages.

PanGenie does **not** accept PacBio CLR reads. The formal adapter requires the
HG002 Illumina R1 and R2 FASTQ files as an explicit, complete pair. It validates
the pair at the plugin boundary, then combines and decompresses them only inside
the sandboxed attempt directory because PanGenie 4.2.1 consumes one
uncompressed k-mer stream.

Formal attempts run without network access in the configured `bwrap` sandbox.

The population VCF must be PanGenie-ready: fully phased from chromosome start
to end, multi-sample, non-overlapping/multi-allelic where necessary, and
sequence-resolved. HG002/NA24385 must not be a panel sample.
