# Stretched CAS(10,11) DZ all-root import

At **R=2.5 bohr**, the covered cc-pVDZ target passes its eight-root-per-sector
SLEPc import. All **64 native roots** agree with independent space-160,
600-cycle CASCI controls on the space-80 and space-160 seeds within
**3.598e-10 / 3.836e-10 Hartree**. Ground-dipole and core/active-subspace gates
also pass against both seeds.

The calculation keeps forty equal-weight orbital-ensemble components,
64 computed target roots and forty configured retained target states.
This target-only run produces no scattering observable.

| Check | Measured difference |
|---|---:|
| All 64 native roots vs covered space-80 seed | ≤3.598e-10 Hartree |
| All 64 native roots vs covered space-160 seed | ≤3.836e-10 Hartree |
| Forty native ensemble roots vs import-run QC | ≤4.939e-11 Hartree |
| Saved ensemble roots vs import-run QC | ≤5.459e-10 Hartree |
| DENPROP ground z dipole vs import-run QC | +1.691e-11 a.u. |
| DENPROP ground z dipole vs space-80/160 seeds | +8.856e-11 / +1.367e-10 a.u. |

Minimum seed/import active-subspace overlaps are
**0.9999999999999990 / 0.9999999999999987**; core overlaps agree to roundoff.
The numerical import and refined fixed-orbital coverage pass. Basis refinement,
competing orbital starts and continuity across geometries remain open.

## Original failure and retained-data recheck

The original engine and batch exit **zero**, but supervisor **933694** and
its post-import verifier exit **one**. The verifier compares native ten-decimal
energies with the nine-decimal saved table at `1e-10 Hartree`; singlet A1 alone
contains valid rounding differences up to `5e-10 Hartree`.

The new recheck uses the previously validated unique-final-CIDATA-set parser
and exact decimal-token comparison. The `%16.9f` writer permits a half-last-unit
format error of **5e-10 Hartree**, independently of the unchanged physical
root/dipole gates of **1e-7 Hartree / 1e-5 a.u.** All 64 saved/native token pairs
pass. Analytic controls retain both half-unit ties and reject corrupt tokens,
inherited-only spectra and missing final sets. The retained calculation data
are reused for the recheck, which exits **zero**.

Original execution records, source hashes and failure logs remain preserved.
The parent [iteration-refinement evidence](../stretched-iteration/README.md)
retains the original 200-cycle coverage failures and the separately passing
600-cycle scans. This companion includes that parent's diagnostic bytes,
checking matching copied files against its immutable archive.

## Measured cost and reconstruction

Wall is **3254.85 seconds / 54.25 minutes**, kernel memory peak
**5657108480 bytes / 5.269 GiB**. Stage totals are:

| Stage | Seconds |
|---|---:|
| QC adapter | 371.96 |
| Integral preparation | 0.10 |
| Eight CONGEN sectors | 2.90 |
| Eight SCATCI sectors | 360.39 |
| DENPROP | 2519.30 |

All **488 payloads** verify by size/digest, with byte-identical repackaging and
a byte-for-byte public fetch.
The executable reconstruction checks ordinary target/import analysis,
**40 raw current QC states**, all 64 native roots against both seeds,
**96 raw reference spectra / 624 eigenpair evaluations**, both checkpoint-based
subspace comparisons and **16229 resource samples**. It also reproduces the
original saved-table rejection and checks the passing exact-decimal recheck.

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/stretched-import
EVIDENCE_DIR=$(mktemp -d)
tar -xzf projects/ukrmol_co/qualification-evidence/stretched-import/stretched-import.tar.gz \
  -C "$EVIDENCE_DIR"
PYTHONPATH=. uv run --with pyscf==2.11.0 --with h5py==3.15.1 \
  python "$EVIDENCE_DIR/state-averaged-evidence/verify-stretched-import.py" \
  "$EVIDENCE_DIR/state-averaged-evidence"
```

The [input recipe](../../calibration-sa11-stretched-import.json) records the
import settings. The passing recheck releases the existing same-geometry
[TZ/aug-TZ QC queue](../../calibration-sa11-stretched-basis-qc.json) on Sadaharu.
Its finite successor and independent follow-on gates are recorded in the
[continuation](../../CONTINUATION.md).
