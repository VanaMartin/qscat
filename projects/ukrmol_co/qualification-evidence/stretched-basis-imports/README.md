# Stretched TZ/aug-TZ imports and the seed-orbital discrepancy

At **R=2.5 bohr**, both CAS(10,11) native target calculations pass the
forty required roots and current-QC dipole checks. The aug-TZ import also
passes all **64 roots** against its independently covered seed. The original
TZ seed-spectrum verifier retains its **exit one**: triplet A1 root eight
differs by **1.563961e-7 Hartree**, exceeding the unchanged **1e-7** gate.

The discrepancy is diagnosed as the small orbital reoptimization during
import. Independent fixed-orbital solves on each import's **actual checkpoint**
agree with every native root within **4.984e-11 / 5.009e-11 Hartree** for
TZ/aug-TZ. The eighth triplet-A1 energy moves between the seed and imported
checkpoints; it is not a missing root or native eigensolver error. The original
strict seed-gate rejection remains part of the evidence.

| Basis | Native-import wall | Kernel peak | Required 40-root error against current QC | All-64-root error against seed |
|---|---:|---:|---:|---:|
| cc-pVTZ | 2986.70 s / 49.78 min | 5.265 GiB | 4.868e-10 Hartree | 1.564e-7 Hartree — reject |
| aug-cc-pVTZ | 3121.02 s / 52.02 min | 5.261 GiB | 5.421e-10 Hartree | 4.755e-8 Hartree — pass |

Both seed checkpoints pass **48-probe** coverage scans at five/eight roots,
spaces 40/80/160 and 600 cycles. Both imported checkpoints subsequently pass
**16-probe** checks at eight roots and spaces 80/160. Together these are
**128 raw probes / 880 eigenpair evaluations**. Physical/spin-penalized
residual maxima are **7.216e-10 / 9.990e-10 Hartree**, below the unchanged
1e-9 gates. The seed scan takes **524.56 s / 0.377 GiB**; the imported-orbital
diagnostic takes **307.42 s / 0.448 GiB**. Hamiltonian-action diagnostics
retain the independently tested Hubbard-dimer and perturbed-vector controls.

All native sectors use the unique newly solved CIDATA set. The saved
`%16.9f` table passes the exact-decimal half-last-unit bound **5e-10 Hartree**,
separately from physical root checks. Native/current-QC dipole errors are
below 1e-10 a.u.; same-basis core/active subspace overlaps agree to roundoff.
These are solver/import diagnostics. Basis, competing-start, state-identity
and geometry-continuity qualification remain open.

## Public evidence and reconstruction

The **1296 payloads** retain both native batches/runs, both seed coverage
scans, both imported-orbital diagnostics, exact checkpoints and sources,
the original failure log and completed controller exits **[0,0,1,0,0]**.
The verifier reconstructs every spectrum/spin flag, native/seed/checkpoint
comparison and **34607 resource samples**, and verifies **478 retained source
hashes**. Original checkout hash lists also identify incidental copies of
previously published archives and ignore rules; their digests are retained
without rebundling those parent archives. All payload digests, byte-identical
repackaging and byte-for-byte public fetch are checked.

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/stretched-basis-imports
EVIDENCE_DIR=$(mktemp -d)
tar -xzf projects/ukrmol_co/qualification-evidence/stretched-basis-imports/stretched-basis-imports.tar.gz \
  -C "$EVIDENCE_DIR"
PYTHONPATH=. uv run --with pyscf==2.11.0 --with h5py==3.15.1 \
  python "$EVIDENCE_DIR/state-averaged-evidence/verify-stretched-basis-imports.py" \
  "$EVIDENCE_DIR/state-averaged-evidence"
```

The [basis-QC parent](../stretched-basis-qc/README.md),
[native import recipe](../../calibration-sa11-stretched-basis-imports.json)
and [live handoff](../../CONTINUATION.md) specify the original seed lineage
and remaining gates. Publication and executable reconstruction persist on
Sadaharu under `prepared/publication-stretched-basis-imports-20261007/`.
