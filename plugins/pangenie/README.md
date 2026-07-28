# PanGenie adapter

This adapter implements the short-read, k-mer pangenome-genotyping family.
It indexes the benchmark-owned phased, sequence-resolved pangenome panel VCF
with `PanGenie-index`, genotypes one diploid sample with `PanGenie`, and returns
an all-sites VCF to the common normalization and evaluation stages.

PanGenie does **not** accept PacBio CLR reads. Use an HG002 Illumina paired-end
dataset combined into one FASTQ stream. Gzip inputs are decompressed into the
attempt-local temporary directory because PanGenie 4.2.1 requires uncompressed
FASTA/FASTQ and VCF inputs.

The population VCF must be PanGenie-ready: fully phased from chromosome start
to end, multi-sample, non-overlapping/multi-allelic where necessary, and
sequence-resolved. HG002/NA24385 must not be a panel sample.
