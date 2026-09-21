# Frozen PG-F1 fixed mini-panel v1

This fixture names the release-blocking cases C01–C15 from the PG-F1 contract.
`tests/unit/test_pgf1_conformance.py` freezes their expected accounting,
2-of-3 behavior, phase-insensitive genotype behavior, missing-evidence failure,
and the fail-closed complex-mapping path. It is a conformance suite, not an
accuracy dataset.
