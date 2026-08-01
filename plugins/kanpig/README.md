# KanPIG adapter

This adapter runs KanPIG 2.x in the shared-alignment execution mode. It
genotypes the benchmark-owned blinded population SV panel from the common
HG002 GRCh38 BAM. KanPIG constructs local variant graphs internally; it does
not consume the whole-genome GBZ graph. Formal attempts run without network
access in the configured `bwrap` sandbox.
