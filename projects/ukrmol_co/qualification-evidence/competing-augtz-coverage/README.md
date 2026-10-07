# Competing aug-TZ orbital branch: fixed-orbital coverage and residuals

Both TZ-projected CAS(10,11) aug-TZ checkpoints pass independent fixed-orbital
coverage: **96 probes / 624 eigenpair evaluations**, including all forty ensemble
roots and eight roots in each singlet/triplet C2v sector. Original supervisor
and profiled calculation exits are **zero**. Wall is **610.35 seconds /
10.17 minutes**, kernel peak **407740416 bytes / 0.380 GiB**.

The [competing-start companion](../competing-starts/README.md) records their
space-80/160 numerical agreement and their failed comparison with the prior
fresh aug-TZ branch. The averaged objective is 0.538388 eV lower; coverage does
not erase that orbital disagreement or establish a unique global optimum.

Each fixed checkpoint has 5/8 roots × spaces 40/80/160 × eight sectors. The
600-cycle refinement retains energy/residual tolerances **1e-12 / 1e-9 Hartree**,
CI linear-dependence 1e-22, spin-penalty shift 1.0 and unchanged orbitals.
The ordinary target builder and `ci_probe.py` still default to **200 cycles**.

| Maximum across both checkpoints | Result |
|---|---:|
| Physical Hamiltonian-action residual | **6.630e-10 Hartree** |
| Spin-penalized residual | **9.977e-10 Hartree** |
| Difference from ensemble energies | **2.274e-13 Hartree** |
| Same-root-count trial-space difference | **1.706e-13 Hartree** |
| Common first-five root-count difference | **1.564e-13 Hartree** |
| CI-vector orthogonality error | **4.771e-14** |

Every probe passes solver convergence flags, finite energies, spin (1e-6),
ensemble consistency (1e-7 Hartree), CI orthogonality (1e-9), and physical/
penalized residuals (1e-9 Hartree). An analytic Hubbard-dimer benchmark verifies
energy offsets and exact-vector residuals; a perturbed eigenvector is rejected.

Input checkpoint SHA256s remain:

- Space 80: `2f0cc1e4cbd31ede3012cf19a9584ddefa88b9240b990ec908d5f2c100bc78fb`
- Space 160: `95ba854810510cb75580bc1cdc7eb058a3bd06ff983fdbb91bd99a77d24b9d7e`

All **335 payloads** verify by size/digest, with byte-identical repackaging and
public fetch. The executable verifier reconstructs every raw spectrum, spin
and convergence block, recomputes ensemble/space/root-count comparisons and
checks the recorded Hamiltonian-action residuals. Raw spectra are in aggregate
`run.log`; individual `casci.log` files are empty because PySCF holds stdout.
Sources, analytic controls, profiles, seed records, parent pair and acknowledged
CPU-slot lease/resumption are retained.

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/competing-augtz-coverage
EVIDENCE_DIR=$(mktemp -d)
tar -xzf projects/ukrmol_co/qualification-evidence/competing-augtz-coverage/competing-augtz-coverage.tar.gz \
  -C "$EVIDENCE_DIR"
PYTHONPATH=. uv run python "$EVIDENCE_DIR/state-averaged-evidence/verify-competing-augtz-coverage.py" \
  "$EVIDENCE_DIR/state-averaged-evidence"
```

The [downward-projection retries](../../calibration-sa11-tz-projection-coreactive-retries.json)
now run on the released CPU slot. A separately gated
[64-root SLEPc import](../../calibration-sa11-competing-augtz-import.json) waits
for that worker's completion, then checks native roots/dipole against the covered
seed and its independent spectra. Its physical gates do not depend on the
downward worker's restart-comparison verdict. The
[continuation](../../CONTINUATION.md) records the exact finite owner transitions.
