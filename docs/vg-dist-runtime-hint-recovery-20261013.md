# Recovery of VG distance-index runtime metadata — 2026-10-13

The initial Q5 run `vg_sr_pgf1_q5_v2_20261013` failed in `lock_graph_assets`
before native Q5 calling. The shared distance index had SHA-256
`995b2884fd67e50ac5faa3d5b781819ecddd70b1dd86b1ac05ed3052399f53ab`,
whereas the existing source lock required
`3f414317e6f22137719bf50a52e0b704056574740245188df3ecc18ad746259a`.

## Exact byte difference and its meaning

An HTTP Range request to the URL already recorded in `SOURCE.tsv` returned
the official first 4,096 bytes. Only offset 140 differs: the shared index has
`1`, the official index has `0`. Hashing the entire local index with only this
byte restored gives the exact frozen source hash. This proves that no other
byte differs from the previously hash-frozen artifact; the source lock is
not edited to accept an unknown file.

VG 1.73.0 uses libbdsg commit
`e74fb663a5f85bc1f76d159b2b3a3691ed85862f`. In that source, `yomo::Pointer`
stores a relative pointer offset followed by a mutable atomic `local` boolean.
`Pointer::get()` can cache `local=true` when the pointer resolves within the
same writable memory segment. The manager computes the same target when
the hint is false and suppresses this write for read-only segments.

The file's first allocated `MappedIntVector` object is at offset 92; its data
pointer is at offset 132 and this pointer's local hint is at offset 140.
The changed byte is the pointer-resolution shortcut, not the packed distance
data. File modification time alone does not identify which historical process
performed this update.

Primary sources:

- [VG Giraffe distance-index recommendations](https://github.com/vgteam/vg/wiki/Giraffe-best-practices)
- [Pinned libbdsg pointer layout and get implementation](https://github.com/vgteam/libbdsg/blob/e74fb663a5f85bc1f76d159b2b3a3691ed85862f/bdsg/include/bdsg/internal/mapped_structs.hpp)
- [Pinned mapped-memory manager](https://github.com/vgteam/libbdsg/blob/e74fb663a5f85bc1f76d159b2b3a3691ed85862f/bdsg/src/mapped_structs.cpp)

## Restoration and mapping verification

A separate bundle at
`/home/luzhiting/pgbench-workspace/shared/pangenome/hprc-v2.1-mc-grch38-eval/derived/vg-giraffe-d46-canonical-20261013`
contains the restored, read-only distance index. Other graph components are
unchanged hardlinks, and original source metadata is copied unchanged. The
shared historical bundle is not overwritten. Evidence is saved under
`/home/luzhiting/benchmark/recovery-assets/vg-dist-recovery-20261013`.

The targeted verification runs real VG 1.73.0 Giraffe on 200 paired reads with
each index variant, both mounted read-only, with the same 16 threads and frozen
fragment parameters. It compares all 400 alignment records after sorting,
including read identity/sequence/quality, paths, alignment scores, mapping
quality and pairing fields. Runtime annotations and GAM serialization order
are excluded. The restored index's full hash is checked again afterward.
The preparation script refuses to register cache transfer unless this test
records a successful comparison.

## Old GAM cache identity

The old graph directory's recorded Merkle hash
`01d7677b880db23b4588e97e23a2f5af4852f1d97fec037b2d6f12b89bea0325`
includes the runtime-hint byte. The canonical directory is
`9f37ffcc8113bb622d79584f5622463b5c26574c2d8a5e32ee492ea624bbbb1e`.
The two hashes must not be declared byte-identical.

After the exact-byte/source-code and real mapping checks, new checkpoint
metadata can explicitly certify this known hint equivalence. It retains the
original mapping key, source checkpoint manifest/hash and GAM digest, adds
the proof hashes, and derives a new mapping key for the canonical bundle and
unchanged mapper recipe. Original checkpoint files/metadata are untouched.
The existing adapter still validates the entire GAM SHA-256 and size before
using a transferred checkpoint; if the current input tree or recipe differs,
the key will not match. No Q0 pack, native VCF or previous score is reused.

## New recovery run

`vg_sr_pgf1_q5_v3_20261013` is preselected before calling or observing its score.
It uses the canonical read-only bundle, repaired adapter, Q5 pack and `call -a`,
and freshly produces all evaluator inputs and score/seal. The failed v2 run
is retained. Its input configuration/context is not rewritten into a successful
run. No genotype, scorer rule, source manifest checksum or Core input hash
check is weakened by this recovery.

This recovery does not invalidate the old Q0 score merely because a pointer
hint changed. Its historical version and input hashes remain distinct and
documented. Mapping time from reused GAMs must not be reported as a fresh
end-to-end mapping measurement.
