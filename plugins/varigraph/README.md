# Varigraph external adapter

This plugin evaluates Varigraph 1.0.8 as a complete paired-short-read
pangenome SV genotyper. It constructs its private graph index from the frozen
leave-HG002-out population panel, counts k-mers from both registered Illumina
mates, and projects only the blinded benchmark candidates into the scored VCF.

The adapter never receives truth data or the HG002 linear-reference BAM. Panel
records outside the blinded universe are discarded at the plugin boundary;
candidate records without a Varigraph genotype are emitted as explicit
`./.` no-calls. Execution uses the same network-disabled `bwrap` sandbox and
provenance contract as every other formal plugin.

The Varigraph conda package is distributed through its authors' `duzezhen`
channel. The benchmark freezes version 1.0.8 in `envs/environment.yaml`.
