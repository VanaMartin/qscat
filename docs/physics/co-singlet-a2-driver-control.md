# Bounded singlet-A2 driver control

## Decision served

The compressed CAS12/TZ QC attempt fails singlet-A2 CI despite a converged
orbital gradient; its final nominal singlets include S²≈2. The
[failed record](../../projects/ukrmol_co/qualification-evidence/compressed-cas12-tz-qc/README.md)
remains rejected. Neither a frozen-orbital diagnostic nor a new initial guess
can retroactively qualify that orbital objective. No physical basis verdict exists.

User-directed continuation selects **one separately frozen driver pilot**, not
another QC attempt. Review of the actual pinned PySCF2.11.0 source finds
`direct_spin0_symm.get_init_guess` builds alpha/beta-exchange-symmetric vectors
for a supplied point-group irrep, not only A1. Its action symmetrizes `C` with
`C.T`. The current CO selector intentionally enables it only for singlet A1;
singlet A2 still uses `direct_spin1_symm`. This supports testing a narrow A2
candidate, not changing the production selector or declaring the driver exact.

The source-review container only reads installed code: zero CI/QC/native/maps.
The thirty-minute design allocation is distinct from the new prospective pilot.
No package patch, dependency install, optimizer change, larger Davidson space,
iteration rescue, target import or additional electronic model is selected.

## Mathematics and exact control

All energies are Hartree. Use a real, spin-independent two-electron Hamiltonian

\[
H=\sum_{pq\sigma}h_{pq}a^\dagger_{p\sigma}a_{q\sigma}
+\tfrac12\sum_{pqrs\sigma\tau}(pq|rs)
a^\dagger_{p\sigma}a^\dagger_{r\tau}a_{s\tau}a_{q\sigma}.
\]

An independent Python oracle constructs its M_s=0 determinant matrix using
fermionic creation/annihilation signs, not a PySCF Hamiltonian contraction.
For one alpha/one beta electron, S² acts on coefficient matrices as `C−C.T`.
Project the point-group block before separating symmetric singlets and
antisymmetric triplets. These are ordinary Hermitian inner products, not ECS
c-products. Verify Hermiticity, eightfold ERI symmetry, `[H,S²]=0`, point-group
selection and normalized orthogonal eigenvectors.

The synthetic four-orbital irreps are A1,B1,B2,A2 (XOR labels0,2,3,1).
One-body energies are−1.2,−.6,.1,.7. Set J03=.8,K03=.12,J12=.65,K12=.08;
couplings `(01|32)=.09`, `(02|31)=.04`, with eightfold ERI completion.
On-site integrals2.7,2.2,1.9,1.4 and additional positive pair diagonals make
the pair-integral metric positive semidefinite. This is a declared synthetic
Hamiltonian, not a fitted CO/Coulomb-basis model.

In the A2 singlet pair basis `(03+30)/sqrt(2),(12+21)/sqrt(2)`,
the exact block is `[[.42,.13],[.13,.23]]`; the triplet block is
`[[.18,.05],[.05,.07]]`. Each has eigenvalues
`(a+d)/2 ± sqrt(((a−d)/2)^2+b^2)`. Compare both PySCF drivers with these analytic
singlets, retaining the original spin penalty and full-Hamiltonian/spin guards.
Include a triplet negative control that satisfies its own physical Hamiltonian
but must fail a singlet guard, plus a perturbed-vector residual control.
Exact energies/action tolerance1e-11; no CO input is touched by this gate.

Also exercise a small **orbital-moving mixed-spin A2 CASSCF**: four H atoms
at±1bohr on x and±1.6bohr on y, STO-3G/C2v, one A1 core, B1/B2 active pair,
two active electrons, equal singlet/triplet-A2 weights. Compare the optimized
two-state energies/objective with the independent two-electron oracle on its
final integral tensors. This is a shared-integral differential check, not an
independent molecular-basis test. Require nontrivial orbital movement, optimizer
and CI convergence, correct spins and full action. It tests a path the previous
A1-only Hubbard startup control did not cover.

## Conditional molecular diagnostic

Only after all exact/mixed-optimizer controls pass, read the unchanged rejected
CAS12/TZ checkpoint as a **diagnostic Hamiltonian**, not a qualified target.
Run direct_spin0_symm for singlet A2 with five and eight roots, space160,
ordinary200 CI cycles and unchanged1e-12 energy/1e-9 residual/1e-22 linear-
dependence settings. Do not rerun the original optimizer or the other seven
sectors; their saved flags/energies remain evidence, not newly qualified roots.

For every returned vector require convergence, finite values, S² error≤1e-6,
CI orthogonality≤1e-9, physical/penalized/full-Hamiltonian residual and spin-
penalty/full-action differences≤1e-9Ha. The independent action uses unrestricted
`direct_spin1.contract_2e`, without point-group compression or transpose
symmetrization. Require alpha/beta exchange error≤1e-9 and first-five energy
agreement between requests≤1e-7Ha. Persist vectors, integrals and actions; do not
infer complete lowest-root coverage or electronic accuracy from two requests.
No larger-space/cycle/seed retry follows a failure.

The API is project-local: `exact_controls()`, `mixed_optimizer_control()` and
`fixed_a2(checkpoint, roots, output)` return saved diagnostics; `run` owns their
order and failure gates. The existing `target.py` selection/defaults are unchanged.
No Rust optimization or qscat promotion is needed: this is a bounded validation
pilot around an existing external solver, not a new production CI implementation.

## Allocation and release boundary

Freeze/push source/tests/exact owner, package/failed-checkpoint/prerequisite pins
before CO use. Four physical cores0–3,8GiB, one BLAS environment thread/four
PySCF threads, read-only old inputs, no network. At most two new **CO** CI
calls, plus the declared synthetic/H4 validation controls; zero CO QC/native/
outer/pole calls. Conditional2–15minutes,1800s scientific /
1920s owner caps,16GiB available RAM/2GiB disk floors. At most90minutes narrow
implementation/freeze effort after this design decision; no general solver API.

Any failed exact control prevents molecular use. Any molecular flag/spin/action/
count/cap/integrity failure stops inconclusive without repair, full QC or model
successor. A verified pass would justify **considering** an A2-only experimental
selector plus mixed-optimizer regression tests and a separately frozen fresh
forty-state optimization from the qualified DZ seed. It does not release that
optimization, validate the failed TZ checkpoint or permit its properties/import.

All original failures, B1 model-position rejection,114 unanalyzed holdout S
matrices, equilibrium exclusion, nine working reference,2048 baseline and parked
new50 remain. B2, coverage/import/scattering, production curves and fits stay gated.
