# State-averaged CO targets for fixed-nuclei scattering

The experimental target builder in `projects/ukrmol_co/target.py` minimizes

\[
 E_{\rm SA}(\mathbf C)=\sum_{S,\Gamma,i}w_{S\Gamma i}
 E_{S\Gamma i}(\mathbf C),\qquad \sum w_{S\Gamma i}=1.
\]

All states share one orthonormal spatial-orbital set. The default ensemble
contains five roots in each of A1, B1, B2 and A2 for each of singlet and triplet
spin, with equal weight per C2v component (40 components, including both Pi
partners). Geometry is in bohr; energies and dipoles are in atomic units.
Two A1 core orbitals remain doubly occupied in every CI configuration, with
ten active electrons. The inactive orbitals can relax during CASSCF; freezing
their occupations is distinct from freezing their coefficients.
The molecule lies along z about the upstream center of mass (C/O masses
12.0110/15.9994), matching the origin of the R-matrix sphere.

PySCF's `state_average_mix_` combines eight symmetry-resolved CI solvers.
Each solver uses an appropriate spin projection and a spin penalty to exclude
higher-spin solutions. Every final root must have the requested S(S+1) within
1e-6, a ten-electron active density, and converged CI coefficients. CASSCF
requires an energy change below 1e-9 Hartree and an orbital-gradient tolerance
of 1e-5 by default; these controls are exposed for tightening checks.
The CI energy tolerance is 1e-10 Hartree. The augmented-Hessian controls are
also exposed: `--target-ah-lindep` (default 1e-14) and
`--target-ah-start-tolerance` (2.5). A small orbital trial vector can be dropped
by the optimizer's metric cutoff even when the requested gradient tolerance
has not been reached. Tight retries use 1e-20/1e-8 for these controls, retaining
the original convergence gates. `target-diagnostics.json` preserves optimizer
and CI convergence flags, final state energies and macroiteration history,
including failed target stages. The energy, gradient and maximum-rotation
sequence distinguishes convergence from stagnation.
`--target-ci-residual-tolerance` optionally sets the CI residual norm separately
from its energy tolerance. `--target-optimizer newton` selects PySCF's coupled
second-order orbital/CI optimizer before constructing the same state-average
ensemble. Its diagnostic gradient includes both orbital and CI variables;
the one-step gradient is orbital-only. Both routes retain the independent
spin, Pi, import and dipole gates. Compare a ground-only control and the full
ensemble when changing optimizer, rather than accepting an optimizer flag as
evidence of target quality.
The CI trial-space cutoff is separately exposed as `--target-ci-lindep`
(default 1e-14). PySCF's Davidson solver drops residual trial vectors when
their squared norm is at or below this cutoff; a stricter residual tolerance
therefore requires a compatible cutoff. The compressed-geometry retry uses
1e-9 residual tolerance and 1e-22 CI cutoff, with 1e-24 for the orbital/CI
augmented-Hessian cutoff. A unit-weight single-component ensemble uses the
ordinary single-state solver, preserving its objective while giving the
Newton implementation the expected CI array rather than a mixed-solver list.
The compatible-cutoff one-step R=1.9 target passes the strict gradient, spin,
Pi and UKRmol import checks. Recorded coupled-Newton retries still fail their
gradient checks, including a CI-refined single-component control. `newton`
remains a diagnostic option; it is not the qualified optimizer for these cases.
`--target-ah-tolerance` controls the augmented-Hessian eigensolver accuracy
(default 1e-12), independently of its metric cutoff and the outer orbital
gradient tolerance. Its residual stopping scale includes the square root of
this tolerance. `--target-verbosity 6` preserves the corresponding trial-space
and step diagnostics in the QC log. The stretched-geometry accuracy ladder
uses 1e-16 and 1e-20 while retaining the tight CI and orbital gates. Both
target-only pipelines pass; their maximum root-energy difference is 1.17e-11
Hartree and dipole difference is 4.35e-12 a.u. This two-point stability check
resolves the recorded tight-CI optimizer stall without qualifying the
electronic model.

The default active space for a state-averaged run should be supplied explicitly,
for example `[4,3,3,0]` for CAS(10,10). Symmetry labels, rather than numeric
irrep IDs, connect PySCF's ordering to UKRmol's A1/B1/B2/A2 ordering.
Core, active and external orbitals are exported in separate Molden blocks.
UKRmol selects them by Molden order, preserving the optimized active subspace
even when generalized-Fock energies overlap between blocks. The JSON inventory
records the actual orbital energies and the subspace of every exported MO.

An optional previous CASSCF checkpoint supplies projected starting orbitals.
The prior core/active sizes and active irrep counts must match; a new run copies
and hashes the checkpoint. The singular values of the initial/final active-space
overlap C_initial^T S C_final quantify subspace drift in the **new** AO metric.
They are orbital diagnostics, not an excited-state identity assignment.
Computed target roots can differ from orbital-ensemble roots. The uniform
`--target-roots` convenience option is overridden by four-component
`--target-singlet-roots` / `--target-triplet-roots` arrays. Increasing those
arrays while holding the `--sa-*` ensemble fixed isolates added channels;
increasing the ensemble itself changes the variational orbital objective.

## Verification contract

- RHF, every CI solver and the orbital optimization must converge.
- Exported orbitals must satisfy C^T S C = I within 1e-9.
- Every equal-weight B1/B2 target pair must agree within 1e-7 Hartree.
- UKRmol must independently reproduce **every averaged root's absolute
  energy** within 1e-7 Hartree. This checks AO normalization, symmetry and
  core/active selection as well as the target diagonalization. The imported
  Molden checksum must match the exported file.
- The ground-state dipole uses the ground-state CI density in the common
  orbitals, rather than the ensemble density. It includes the nuclear term.
  UKRmol's DENPROP must agree within 1e-5 atomic units. DENPROP records
  electronic position minus nuclear charge position for L=1,M=0, so its
  value is negated before comparison with the physical PySCF dipole.

`--target-only` runs both quantum chemistry and the UKRmol target pipeline,
before paying for scattering. The upstream slot and retained output filename
are called `psi4`; on `--orbitals state-averaged` the adapter actually invokes
PySCF and reads its explicit JSON MO inventory. Successful records identify
`target_backend: pyscf`; the log is not parsed as a Psi4 calculation.

These gates establish solver/import consistency, not electronic convergence.
The next checks vary basis, active space, ensemble and channel count, tighten
optimization tolerances, and inspect compressed/equilibrium/stretched target
continuity. Scattering additionally needs continuum and fit-window checks.
The state-averaged neutral energy is not the ground-state energy, and the
ground root in these common orbitals is not an independently correlated
neutral potential curve.
Different starting active spaces can yield different stationary solutions.
Compare the ensemble objective as well as each physical root: minimizing the
former can trade away ground-state accuracy, especially with diffuse orbitals.
Neither the lowest ensemble energy nor apparent agreement with one excitation
is a sufficient qualification criterion for the low-energy scattering target.

The published model motivating the 40-component ensemble is described in
[`dora-2016-epjd70-197.md`](../../reference/literature/dora-2016-epjd70-197.md),
p. 3 (orbital ensemble and Table 1 active spaces). This implementation uses
PySCF rather than the paper's quantum-chemistry backend; equivalence to its
optimized orbitals must be measured, not inferred from matching state counts.
