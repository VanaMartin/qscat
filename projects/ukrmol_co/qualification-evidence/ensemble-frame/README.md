# Physical CO ensemble-frame comparison

[Result](../../ENSEMBLE_FRAME_COMPARISON.md): all eight numerical comparisons
complete, but both A1 retained-five subspaces fail the frozen0.5 dominance screen.
Original command/container/owner **1/1/1** are the scientific rejection, not a
solver failure. New target/property qualifications remain valid; the original
seed rejection stays separate. Cost63.503s /.49293GiB, four reserved cores.

```sh
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/ensemble-frame
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/target-state-overlap-roots12
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/fifty-import
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/fifty-frame
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/fifty-properties
```

All **11 new payloads** fresh-verify; **18 prerequisite payloads** reauthenticate.
`results.tar.gz` preserves **85 files /108,765,937 uncompressed bytes**: eight new
CI arrays/logs/root records, ten AO/MO/self metrics,24 overlap arrays across three
routes, numerical/diagnostic records, progress/final reports, owner/source freeze,
commands/exits and resources. Every generated/owner file has a byte count/hash in
`collection-sha256.json`. Collection SHA256:
`278a93c6cdcbb03d10d4f0f616e5a92dd4b744b5222ab76aac5ba6fa14e2f5f8`.

The four prerequisite artifact pointers preserve exact prior publication bytes.
`prerequisite-input-map.json` maps **all32 original input paths** to the appropriate
prerequisite archive and exact member, SHA256 and byte count. In particular, the
eight old CI arrays come from the existing twelve-root publication, not a new
generation or duplicate publication. Extract each mapped member into a fresh
`prerequisites/` tree at its map key and verify all32 digests before replay.
Never extract over original producer runs.

Extract the new collection into a fresh directory, its `source.tar.gz` into
`source/` and the separate `verification-support.tar.gz` onto that source tree.
Authenticate all12 producer source hashes and the one reader-support hash against
their manifests. The reader overlay, frozen at **`e066879`**, adds a module without
changing the producer freeze **`869a588`**. Then, from the extraction root:

```sh
PYTHONPATH=source python -m projects.ukrmol_co.ensemble_frame_reading --root prerequisites --comparison outputs/comparison --contract source/projects/ukrmol_co/ensemble-frame-overlap-contract-20261021.json --output fresh-reading.json
```

This replay uses full occupied Slater minors from authenticated old/new CI arrays
and saved physical metrics, not the producer Schur contraction. Hermitian Gram
eigenvalues independently reconstruct the squared-SVD screen. All eight overlap
blocks reconstruct within2.55351e-14, and the original **false** decision remains.
The reader exits successfully when reconstruction succeeds even though the
scientific dominance verdict is false; do not reinterpret that as a scientific
pass. CI residual actions/AO integral generation are observed producer evidence,
not recomputed in this saved-only reading. No CI/native/scattering solve is run.

[Verification record](../../verification-ensemble-frame-public-20261021.json)
also authenticates every new/prerequisite payload,32 input pins, the frozen
sources and original exits. It reconstructs old-cache and new-root metadata gates
and preserves2 weak /0 ambiguous /3 retained-to-extra diagnostic flags across80
directed assignments. The extreme triplet-A1 minimum is numerically near zero;
its exact tiny decimal is not a resolved physical quantity.

No automatic root/channel ladder, retained-set change, new property export,
Hamiltonian or production sweep follows this rejection.
