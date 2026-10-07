# Repaired fresh CAS(10,11) aug-TZ all-root import

The repaired fresh aug-cc-pVTZ target passes its eight-root-per-sector SLEPc
UKRmol import. All **64 native roots** agree with independent eight-root,
space-160, 600-cycle controls on the space-80 seed within **5.969e-11 Hartree**.
The separate space-160 repaired seed agrees with the same native roots within
**5.764e-11 Hartree**. Both pass the unchanged **1e-7-Hartree** physical gate.

There are forty equal-weight orbital-ensemble components, 64 computed roots
and forty configured retained target states. This is target-only execution;
it produces no scattering observable.

| Check | Measured difference |
|---|---:|
| All 64 native roots vs covered space-80 seed | ≤5.969e-11 Hartree |
| All 64 native roots vs covered space-160 seed | ≤5.764e-11 Hartree |
| All forty ensemble roots vs import-run QC | ≤4.672e-10 Hartree |
| DENPROP ground z dipole vs import-run QC | −5.169e-12 a.u. |
| Import-run QC roots vs original repaired seed | ≤3.621e-11 Hartree |
| DENPROP ground z dipole vs original repaired seed | +7.474e-11 a.u. |

The seed is slightly reoptimized in the import run; the first two comparisons
include that reoptimization. Minimum seed/import core and active-subspace
singular values are **0.9999999999999984 / 0.9999999999999990**. Spins, Pi
degeneracy and all QC flags pass. Same-basis CI-space repair agreement and
independent import are established for this fresh lineage; competing starts
and electronic-model convergence remain open.

Original engine/batch, supervisor and final verifier exits are all **zero**.
Wall is **3360.53 seconds / 56.01 minutes**, kernel peak **5656363008 bytes /
5.268 GiB**. Stage totals are:

| Stage | Seconds |
|---|---:|
| QC adapter | 643.11 |
| Integral preparation | 1.25 |
| Eight CONGEN sectors | 2.81 |
| Eight SCATCI sectors | 394.14 |
| DENPROP | 2319.03 |

The verifier uses the final CIDATA-set parser and exact decimal-token contract
validated by the [fresh TZ import](../fresh-tz-import/README.md). Native spectra
retain ten decimals, while the saved table uses `%16.9f`; format consistency
passes its half-last-unit **5e-10-Hartree** bound. This is distinct from physical
root/dipole gates of **1e-7 Hartree / 1e-5 a.u.**. Analytic controls reject corrupt
tokens, inherited-only spectra and missing final sets. The source, original
verification and acknowledged CPU-slot lease/resumption are preserved.

All **181 payloads** verify by size/digest, with byte-identical repackaging and
a byte-for-byte public fetch. The executable reconstruction reruns the ordinary
target/import analysis, reconstructs all 64 final-set native roots and checks
all **96 raw reference spectra / 624 eigenpair evaluations** across both repaired
seeds. Orbital overlaps are retained measurements with the executable PySCF
verifier. Parent evidence is the
[repair companion](../fresh-augtz-repairs/README.md) and
[coverage companion](../fresh-augtz-coverage/README.md); the original rejected
fresh checkpoint and nonzero exit remain exact in those earlier records.

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/fresh-augtz-import
EVIDENCE_DIR=$(mktemp -d)
tar -xzf projects/ukrmol_co/qualification-evidence/fresh-augtz-import/fresh-augtz-import.tar.gz \
  -C "$EVIDENCE_DIR"
PYTHONPATH=. uv run python "$EVIDENCE_DIR/state-averaged-evidence/verify-fresh-augtz-import.py" \
  "$EVIDENCE_DIR/state-averaged-evidence"
```

The [input recipe](../../calibration-sa11-fresh-augtz-import.json) preserves
the import contract. The [continuation](../../CONTINUATION.md) records four
same-geometry TZ↔aug-TZ projections at spaces 80/160, checking the
[competing starts](../../calibration-sa11-fresh-basis-competing-starts.json)
against each other and the covered/imported fresh solutions.

Those [competing-start results](../competing-starts/README.md) find a distinct,
numerically reproducible aug-TZ orbital solution with objective **0.538388 eV
lower** and minimum active-subspace overlap **0.141085**. This import remains
valid for its recorded fixed orbitals; it does not establish a unique optimum.
The new branch's coverage/import gates are independent.
