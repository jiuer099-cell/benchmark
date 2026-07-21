# HG002 GRCh38 pangenome SV benchmark

The benchmark uses GRCh38 for the reference, truth set, stratifications,
alignments, and pangenome backbone. HG002 truth alleles must be excluded from
population/pangenome inputs. Production assets are content-addressed and the
workflow fails closed until their SHA-256 values are frozen.

Truvari, Aardvark, and vcfdist independently classify every normalized query
result in one shared universe. Results are counted by how many evaluators accept
them: three, two, one, or zero. Each evaluator has one equal vote. The score is
`(3*n3 + 2*n2 + n1) / (3*N) * 100`; no weighted layers exist. Provenance only
controls validity status and does not add points. Tool rankings are forbidden.
