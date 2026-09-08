# BayesTyper community adapter

This adapter runs BayesTyper 1.5 from the benchmark-supplied paired Illumina
FASTQs and the frozen canonical candidate panel. The plugin registration must
set `tool_index` to a frozen directory containing `canon.fa` and `decoy.fa`
from the GRCh38 BayesTyper reference bundle. Those files are hashed into the
resolved-input and allowed-information manifests.

The adapter uses the upstream project's documented k=55 KMC, `makeBloom`,
`cluster`, and `genotype` workflow. It then projects native calls onto every
canonical candidate; candidates without a resolved native call are emitted as
`./.` instead of being dropped.
