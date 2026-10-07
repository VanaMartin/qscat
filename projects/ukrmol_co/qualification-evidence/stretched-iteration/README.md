# Stretched-DZ CI iteration refinement

This follow-on repeats all nine failed settings from the
[original coverage companion](../stretched-coverage/README.md), increasing the
CI iteration limit from **200 to 600**. Orbitals, requested roots, trial spaces,
spin penalty and energy/residual tolerances remain the same. All nine failures
reproduce at 200 cycles and pass at 600. Complete **48-probe scans on both
checkpoints** then pass every requested root, spin, convergence flag, energy and
residual gate. The original 200-cycle failures retain their immutable archive.

The diagnostic explicitly computes
`||(H - E) c|| / ||c||` in Hartree for the spin-penalized Hamiltonian and
the physical Hamiltonian, using the latter's Rayleigh energy. Across the full
96-probe scans, maxima are **9.94e-10 / 7.07e-10 Hartree**, respectively; both
meet the existing 1e-9 threshold. CI-vector orthogonality also passes at 1e-9,
and ensemble energies agree with the recorded targets within **1.28e-13
Hartree**. The scans check 624 eigenpairs. Nine paired retries add eighteen
probe attempts, for **114** attempts in this companion.

The residual calculation is gated by an analytic two-site/two-electron Hubbard
dimer. For hopping `t=1` and on-site repulsion `U=4`, its singlet ground energy is
`(U - sqrt(U**2 + 16*t**2))/2`; a constant energy offset of 3 also checks the
inactive/core-energy subtraction. The exact vector's physical residual is
1.74e-15, while a deliberately perturbed vector gives 0.0504. Raw CO spectra,
spins and whole-probe convergence messages are reconstructed from `run.log`;
the source and report retain the actual Hamiltonian-action residuals.

The worker completes in **789.89 seconds / 0.285 GiB**, with original exit zero
and recorded resumption of its waiting CPU-slot owner. Its 281 payloads verify
by size and digest, and repackaging is byte-identical. The public client fetched
and verified the exact archive. The two initial checkpoints and target inputs
travel with it; the frozen numerical controller can repeat the residual solves
in the pinned CO image.

This supports fixed-orbital coverage under the **refined 600-cycle contract**.
Independent UKRmol import, basis, competing-start and electronic-model gates
remain required. The full target import and same-geometry TZ/aug-TZ follow-ons
are queued on Sadaharu. The main 46-attempt aggregate and original companion
retain their recorded outcomes.

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/stretched-iteration
EVIDENCE_DIR=$(mktemp -d)
tar -xzf projects/ukrmol_co/qualification-evidence/stretched-iteration/stretched-ci-iteration.tar.gz \
  -C "$EVIDENCE_DIR"
PYTHONPATH=. uv run python "$EVIDENCE_DIR/state-averaged-evidence/verify-stretched-iteration.py" \
  "$EVIDENCE_DIR/state-averaged-evidence"
```

See [the live continuation](../../CONTINUATION.md) for the dependency chain,
CPU groups, source snapshots and resource gates.
