# Repaired fresh aug-TZ coverage and residuals

Both [fresh-lineage aug-TZ repairs](../fresh-augtz-repairs/README.md) pass their
independent final-orbital coverage/residual scans. Each checkpoint has 48
probes: both spins, all four C2v sectors, five/eight roots and CI trial spaces
40/80/160. The **600-cycle** contract keeps the tight tolerances, spin penalty,
orbitals and forty-component ensemble fixed. Ordinary QC's 200-cycle default
is unchanged.

| Repaired seed CI space | Eigenpair evaluations | Maximum physical residual (Hartree) | Maximum penalized residual (Hartree) | Maximum trial-space difference (Hartree) |
|---|---:|---:|---:|---:|
| 80 | 312 | 6.51464e-10 | 9.93623e-10 | 1.27898e-13 |
| 160 | 312 | 6.62042e-10 | 9.94099e-10 | 1.98952e-13 |

All convergence flags, spin, physical/spin-penalized residual and CI-vector
orthogonality gates pass across **96 probes / 624 eigenpair evaluations**.
Maximum ensemble-energy, common-first-five root-count and CI-vector
orthogonality differences are **1.85e-13 Hartree / 1.99e-13 Hartree / 1.23e-13**.
The analytic Hubbard-dimer residual and perturbed-vector rejection controls
pass before the scans. Residual tolerances remain **1e-9 Hartree**.

Original exit is **zero**, wall **600.90 seconds / 10.02 minutes**, kernel
peak **419938304 bytes / 0.391 GiB**, aggregate CPU time **1914.58 seconds**.
The frozen controller waits for the fresh TZ import to finish, independently
checks the repaired pair and exact seed hashes, safely leases the still-waiting
CPU-4–7 owner and resumes it afterward. The coverage-gated augmented-TZ import
then acquires that slot.

The archive preserves both unchanged input checkpoints and QC records, parent
repair/pair evidence, all probe records, raw aggregate CI logs, analytic
controls, resources and frozen source/lease records. The executable verifier
reconstructs all 96 raw spectra/convergence blocks and 624 energies/spins,
recomputes the numerical gates and checks the recorded Hamiltonian-action
residuals. Individual `casci.log` files are empty because PySCF retains its
own stdout; `run.log` contains the raw CI evidence.

All **253 payloads** verify by size/digest, with byte-identical repackaging and
a byte-for-byte public fetch. This establishes lowest-root coverage on these
fixed repaired orbitals. Competing starts, independent UKRmol import and
electronic-model convergence remain gates. The rejected fresh parent is
preserved in the repair companion with its original nonzero exit.

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/fresh-augtz-coverage
EVIDENCE_DIR=$(mktemp -d)
tar -xzf projects/ukrmol_co/qualification-evidence/fresh-augtz-coverage/fresh-augtz-coverage.tar.gz \
  -C "$EVIDENCE_DIR"
PYTHONPATH=. uv run python "$EVIDENCE_DIR/state-averaged-evidence/verify-fresh-augtz-coverage.py" \
  "$EVIDENCE_DIR/state-averaged-evidence"
```

The [continuation](../../CONTINUATION.md) records the now-active
[64-root augmented-TZ import](../../calibration-sa11-fresh-augtz-import.json).
