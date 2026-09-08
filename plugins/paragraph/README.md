# Paragraph adapter

This bundled adapter begins from the benchmark-supplied paired Illumina FASTQ,
maps those reads inside the billed tool run, and genotypes the frozen canonical
candidate panel with Paragraph. The target truth, target assembly, relatives'
genotypes, and undeclared callsets are not exposed to the process.

The adapter preserves every canonical candidate. Any site missing from the
native output is subsequently materialized as `./.` with an explicit output
status; it is never inferred to be `0/0`.
