# Fresh CAS(10,11) TZ all-root import

The covered fresh cc-pVTZ target completes its eight-root-per-sector SLEPc
UKRmol import. All **64 roots** agree with the independent eight-root,
space-160, 600-cycle controls on the seed checkpoint within
**3.196e-8 Hartree**, passing the existing **1e-7-Hartree** physical gate.
There are forty orbital-ensemble components and 64 computed target roots;
target-only execution has no scattering observable.

| Check | Measured difference |
|---|---:|
| All 64 native roots vs independent seed-orbital controls | ≤3.196e-8 Hartree |
| All forty ensemble roots vs import-run QC | ≤5.253e-10 Hartree |
| DENPROP ground z dipole vs import-run QC | −3.366e-11 a.u. |
| Import-run QC roots vs original fresh seed | ≤3.109e-8 Hartree |
| DENPROP ground z dipole vs original fresh seed | −5.119e-8 a.u. |

The import runner reoptimizes from the seed; the first comparison therefore
includes that small orbital change, rather than measuring only eigensolver
error. Minimum seed/import core and active-subspace singular values are
**0.9999999999999983 / 0.9999999999997468**. Spins, Pi degeneracy and QC
convergence pass the usual gates. Competing orbital starts and electronic-model
convergence remain open.

Original engine/batch exit is **zero**, wall **3262.90 seconds / 54.38 minutes**,
kernel peak **5658492928 bytes / 5.270 GiB**. Stage totals are:

| Stage | Seconds |
|---|---:|
| QC adapter | 537.91 |
| Integral preparation | 0.34 |
| Eight CONGEN sectors | 2.86 |
| Eight SCATCI sectors | 371.32 |
| DENPROP | 2350.20 |

The final supervisor verifier originally exits **one** after the successful
engine run. It compares ten-decimal native energies to the nine-decimal saved
table with an incorrectly tight 1e-10-Hartree format bound. The retained Perl
writer, `lib/ukrmollib.pm::save_target_energies`, uses `%16.9f`; half the last
printed unit is **5e-10 Hartree**. Exact decimal-token checks now verify that
bound. This format-consistency contract is separate from the unchanged
physical root/dipole gates.

A first read-only recheck also exits **one** because later SCATCI outputs echo
earlier CIDATA spectra. The final passing recheck selects the unique spectrum
following `CI data will be stored as set number`, verifies the expected sector
number, and ignores those earlier sets. Analytic controls accept exact/half-unit
rounding and reject corrupt tokens, inherited-only spectra and missing final
sets. Original source, logs and nonzero exits remain preserved for both failed
verifier attempts. The passing recheck takes **0.654 seconds** in a small
single-core analysis container; its brief overlap with augmented-TZ coverage is
explicitly recorded.

All **181 payloads** verify by size/digest, with byte-identical repackaging and
a byte-for-byte public fetch. The archive includes the successful import,
frozen batch/controller sources, raw sectors/dipoles, seed/import checkpoints,
reference QC and coverage inputs, both verifier failures, the successful
recheck and executable reconstruction. The ordinary import analyzer reruns,
all 64 final-set spectra reconstruct, and 48 raw reference spectra/spins are
checked. The orbital overlaps are retained measurements with their executable
PySCF verifier. Parent evidence is the
[fresh QC companion](../fresh-basis/README.md) and
[fixed-orbital coverage companion](../fresh-coverage/README.md).

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/fresh-tz-import
EVIDENCE_DIR=$(mktemp -d)
tar -xzf projects/ukrmol_co/qualification-evidence/fresh-tz-import/fresh-tz-import.tar.gz \
  -C "$EVIDENCE_DIR"
PYTHONPATH=. uv run python "$EVIDENCE_DIR/state-averaged-evidence/verify-fresh-tz-import.py" \
  "$EVIDENCE_DIR/state-averaged-evidence"
```

The [continuation](../../CONTINUATION.md) records the independent augmented-TZ
coverage worker and the coverage-gated
[augmented-TZ import recipe](../../calibration-sa11-fresh-augtz-import.json).
