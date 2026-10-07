# Equilibrium aug-TZ → TZ competing-start qualification

All four new-name downward starts pass tight CAS(10,11) QC at **R=2.1323 bohr**.
Each projects one of the two reproducible aug-cc-pVTZ branches onto cc-pVTZ,
using numerical CI spaces 80/160 with the same forty-component orbital ensemble.
The repaired projection supplies thirteen core+active columns from the larger
92-column MO matrix; PySCF completes destination virtuals to a 60-column matrix.

| Source aug-TZ branch | CI space | Wall | Kernel memory peak |
|---|---:|---:|---:|
| Repaired fresh RHF lineage | 80 | 2085.26 s / 34.75 min | 0.352 GiB |
| Repaired fresh RHF lineage | 160 | 2049.96 s / 34.17 min | 0.409 GiB |
| Upward fresh-TZ lineage | 80 | 1386.34 s / 23.11 min | 0.353 GiB |
| Upward fresh-TZ lineage | 160 | 1326.20 s / 22.10 min | 0.409 GiB |

Original engine/batch/controller exits are zero. The sequential batch takes
**6849.90 seconds / 1.903 hours**. All orbital/fresh-CI, spin, MO/Pi and
ensemble checks pass. Two within-lineage space comparisons and four comparisons
against the qualified fresh-TZ target pass the unchanged same-basis restart gates:

- Maximum root difference: **4.156e-8 Hartree** (gate 1e-7).
- Maximum dipole difference: **6.916e-8 a.u.** (gate 1e-5).
- Minimum active-subspace singular value: **0.9999999999995466**.
- Maximum ensemble-objective difference: **8.527e-14 Hartree**.

Both aug-TZ lineages thus recover the existing TZ solution on downward
projection. This does not make their distinct aug-TZ solutions identical,
establish a unique orbital optimum, or qualify the full geometry curve.
The new aug-TZ branch retains its separate independent-import gate.

## Preserved failures and reproducibility

The [original competing-start companion](../competing-starts/README.md)
preserves both rejected 92→60-column interface attempts and the two passing
upward starts. This new companion retains four fresh attempts, six comparison
records, all source/seed hashes and the waiting-owner pause/resume records.
The exact projection repair was separately checked with four genuine PySCF
regressions and the CO test suite.

All **296 payloads** verify by size/digest with byte-identical repackaging and
a byte-for-byte public fetch.
The executable reconstruction checks **160 raw final QC states**, all source
and checkpoint identities, six spectrum/objective/dipole comparisons, six
checkpoint-based subspace comparisons and **34132 resource samples**. CSV
elapsed-time rounding is checked at its half-millisecond serialization bound.

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/downward-projections
EVIDENCE_DIR=$(mktemp -d)
tar -xzf projects/ukrmol_co/qualification-evidence/downward-projections/downward-projections.tar.gz \
  -C "$EVIDENCE_DIR"
PYTHONPATH=. uv run --with pyscf==2.11.0 --with h5py==3.15.1 \
  python "$EVIDENCE_DIR/state-averaged-evidence/verify-downward-projections.py" \
  "$EVIDENCE_DIR/state-averaged-evidence"
```

The [input recipe](../../calibration-sa11-tz-projection-coreactive-retries.json)
and [continuation](../../CONTINUATION.md) record the finite worker and scheduling
handoff. Its completion resumes the staged owner on CPUs 4–7; the covered
competing aug-TZ import follows the existing owner-completion contract.
