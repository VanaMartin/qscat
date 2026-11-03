# Oxygen affinity: testing the open-shell triples orbital convention

## Question and scope

The frozen aug-5Z calculation gives EA(O)=1.426606459 eV, following
1.403072709 eV at aug-QZ. Correct charge, high-spin 2s2/2p4 and 2s2/2p5
occupations, one frozen 1s pair, convergence and energy subtraction reconstruct
from original logs. No R-matrix boundary is used in either atom.

A remaining implementation question is whether the perturbative triples energy
was evaluated in an appropriate orbital representation. PySCF 2.11's
ROHF-to-UHF conversion preserves the determinant; it explicitly does not produce
canonical UHF orbitals. Its UCCSD integral initializer recomputes the spin Fock
matrices but takes their diagonals as denominators without diagonalizing their
occupied/virtual blocks. Nonzero off-diagonal blocks therefore warrant an
explicit semicanonical control. This is a hypothesis, not a diagnosed error yet.

## Mathematics and invariants

Keep the ROHF determinant, frozen core, basis and Hamiltonian fixed. For each
spin diagonalize the active occupied and virtual Fock blocks separately, within
each spatial irrep. Leave the frozen core unchanged. With new coefficients
C'=CU, transform the real amplitudes covariantly:

$$
t'_{IA}=\sum_{ia}U_{iI}t_{ia}U_{aA},\qquad
t'_{IJAB}=\sum_{ijab}U_{iI}U_{jJ}t_{ijab}U_{aA}U_{bB}.
$$

The spin density and reference energy must be unchanged. CCSD correlation energy
must also be invariant when integrals and amplitudes are transformed consistently.
The diagonal-denominator triples approximation need not be invariant under an
arbitrary noncanonical occupied/virtual rotation; its change is the observable.
This does not reoptimize UHF or rerun SCF/CCSD iterations.

`oxygen_affinity_orbitals.py` exposes `block_rotation(fock, nocc, labels)` and
`rotate_amplitudes(t1, t2, rotations, nocc)`, pure NumPy real-orbital helpers.
Tests cover block/irrep/core separation, orthogonality, density invariance and
one-/two-body tensor-contraction and inverse-transformation invariance.
All energies are Hartree internally; only differences are converted through
`qscat.units`.

## Frozen diagnostic gates

Reuse only the independently reconstructed aug-5Z O/O- ROHF and CCSD checkpoints.
Authenticate all six original input files before and after execution. Rebuild
integrals and replay original CCSD energy/triples before the changed convention.
Do not modify the original capsules or endpoint totals.

- Original reference energy: 1e-11 Ha; original CCSD correlation/total replay:
  1e-8 Ha; original triples replay: 1e-9 Ha, all absolute/zero relative tolerance.
- Rotations/metric/spin density: 1e-10 absolute; active occupied/virtual Fock
  off-diagonal residual after rotation: 1e-9 Ha.
- Rotated CCSD correlation energy invariance: 1e-8 Ha. A failure stops before
  interpreting the new triples energy; no iterative CCSD rescue is released.
- Report each neutral/anion triples shift and their EA difference separately.
  An effect exceeding 0.1 meV is numerically resolved, not automatically an
  explanation of all physical error or an accuracy pass against experiment.

Four physical Sadaharu cores0–3, unchanged 16-GiB container/12-GB PySCF allowance,
24-GiB available-RAM and 10-GiB disk floors. Two sequential saved-state workers;
900 s each, 1800 s science, 1950 s container, 2100 s owner. Conditional 1–5 min
from the prior 5Z timings; no speedup claim. At most one hour implementation.
Durable ownership, 600-s active status observations plus immediate completion;
preserve original failures, resources and partial evidence. No fresh SCF/CCSD
iterations, 6Z, molecule, scattering, fit, holdout use or automatic successor.

## Decision interpretation

A significant semicanonical effect identifies an orbital-convention contribution,
not necessarily a complete remedy. A negligible effect rules it out at this
resolution. Remaining basis radial/diffuse/angular incompleteness, correlation
truncation, core correlation, relativistic/fine-structure conventions and external
reference uncertainty stay separate. Do not tune any of them to an observed EA.

The user requests diagnosing and correcting the affinity rather than another
unexplained basis ladder. This is one bounded reference-correctness stage,
separate from the PR under review and the neutral-dissociation calculation.
