# CO electronic-model qualification

The experimental workflow qualifies the orbital model before extending its
fixed-nuclei scattering campaign. Geometry is in bohr, energies in Hartree and
dipoles in atomic units. The state-average objective and import contract are
defined in [co-state-averaged-target.md](co-state-averaged-target.md).

## Inexpensive orbital checks

`projects.ukrmol_co.run --qc-only` executes the existing PySCF target builder
without UKRmol integrals, configuration generation or diagonalization. It keeps
the same optimizer, ensemble, fresh-CI, spin, density, Pi and MO gates. The
analyzer checks the recorded ensemble composition, energies, optimizer flags,
reported gradient, fresh-CI spectrum/spins and exported Molden digest. These
records are `qc_validated`, distinct from independently import-validated runs.
`--target-only` continues to run both QC and UKRmol target checks.

The first CAS(10,11) refinement holds the 40-component ensemble and DZ basis
fixed. It requests energy/gradient tolerances 1e-11/1e-7 and CI energy/residual
tolerances 1e-12/1e-9, with the established compatible CI/Hessian cutoffs.
Starts are the successful fresh-CI checkpoint, the original rejected checkpoint
and canonical RHF orbitals. Compare the ensemble objective, individual roots,
ground dipole and active-subspace overlaps; do not select solely by one energy.
Final-orbital root-count/trial-space probes retain that orbital objective.
Passing finalists receive the full independent UKRmol import/dipole checks.

The matched CAS(10,11) ladder exposes CI convergence failures at stretched DZ
and equilibrium TZ after orbital convergence. `--target-ci-max-space` refines
the PySCF Davidson trial space without changing the five-root-per-sector
orbital ensemble, its weights or convergence tolerances. The default remains
`max(40, 8*roots)` in each sector; an override must exceed every ensemble root
count. Actual sector limits are retained in `target-diagnostics.json` and
checked by the QC-only analyzer when an override is requested. The completed
nine-job ladder has two QC passes, three CI failures and four projection setup
failures. PySCF 2.11.0's `project_init_guess` rejects simultaneous changes in
geometry and basis. Stage those starts through a passing same-geometry DZ
checkpoint, then change basis at fixed geometry; retain the rejected attempts.

Validate the override at 80/160 trial vectors against the passing CAS(10,8)
dense root/dipole control, then reoptimize each rejected final checkpoint in
fresh directories with those spaces. Compare root coverage, objective, dipole
and active subspaces; preserve the original failures. A fixed-orbital CI repair
alone does not establish minimization of the intended lowest-root ensemble.

## Independently correlated neutral pilot

The neutral pilot uses conventional spherical-basis RHF followed by
frozen-core restricted CCSD(T), with the two lowest occupied spatial orbitals
excluded from correlation. All 14 electrons enter RHF; ten enter correlation.
The retained energies are absolute; a neutral potential is formed later as
`E(R) - E(R_reference)` using one explicitly named reference geometry.

The initial basis comparison is aug-cc-pVTZ/aug-cc-pVQZ at R=1.9, 2.1323 and
2.5 bohr. Require converged RHF, CCSD and lambda equations, finite perturbative
triples and a stable RHF reference. Record internal/external RHF stability,
`||t1||_F / sqrt(10)` and the largest singular value of `t1`, without treating
a diagnostic threshold as a physical error bar. The dipole is from the CCSD
lambda density, explicitly not a CCSD(T) energy derivative. Its AO density
must contain 14 electrons and include the nuclear contribution.

A cc-pVDZ CO control compares independent PySCF and Psi4 conventional RHF and
frozen-core CCSD(T) energies within 1e-7 Hartree before the basis pilot is used.
This establishes implementation agreement in that control, not basis or
single-reference convergence. Compare relative neutral energies and dipoles
across the pilot bases; three geometries alone do not qualify a full curve.
Stretched/dissociation regions require their own reference-stability and
correlation checks before selecting a single- or multireference treatment.

The cost-constrained continuation adds aug-cc-pV5Z at the same three geometries
to measure a TZ/QZ/5Z sequence using the same frozen-core CCSD(T) method and
CCSD density. Reference each basis's relative curve to its own equilibrium
energy; report both successive shifts before attempting extrapolation. This
tests basis refinement, while stretched correlation-treatment checks remain
a separate requirement.

The completed 5Z sentinels reduce QZ-to-5Z relative-energy changes to
17.24/15.82 meV at R=1.9/2.5, versus 76.04/68.47 meV for TZ-to-QZ. This is
basis-refinement evidence, not a full-curve qualification. Further aug-TZ/QZ
checks at R=3.0 and 4.0 converge RHF but fail **external reference stability**
in all four cases. The conventional restricted CCSD(T) route stops at that
gate; those geometries need a separately designed and validated correlation
treatment before their neutral energies are used.

All calculations retain finite CPU/memory limits, persistent scratch, input
configs, raw logs, checkpoints, stage timings, source/image hashes and failures.
The methods remain experimental under `projects/ukrmol_co`; no library
promotion or potential-fitting qualification follows from these pilot checks.

## CAS(10,12) dense-diagonalization memory audit

The installed engine's `SCALAPACKDiagonalizer_module` queries both `PDSYEVD`
and `PDORMTR`, using the larger of the eigenvector workspace and the orthogonal
transformation workspace plus `2*N`. It allocates full distributed Hamiltonian
and eigenvector matrices even when only five output eigenpairs are requested.
`projects/ukrmol_co/scalapack_workspace.f90` performs the same workspace queries
without allocating either matrix. Compile with `mpifort -fdefault-integer-8`
against the pinned ScaLAPACK/OpenBLAS libraries in the build/toolchain image,
then execute in the gated runtime. Its inputs are dimension and process-grid
rows/columns; block size 64 matches the recorded engine run.

The queried array floor includes `A`, `Z`, real/integer workspaces and eigenvalues.
It excludes build caches, MPI/library overhead and other live engine arrays.
The four-rank singlet control reproduces lwork/liwork 935713305/302376 and a
55.62-GiB array floor, close to the measured 56.10-GiB peak.

For the triplet dimension 70674, the four-rank eigensolver query returns only
2615578 doubles on rank zero, below the necessary `2*35346*35346` eigenvector
workspace alone. This is consistent with overflow at the 32-bit integer limit;
its positive return value must not be accepted as a usable memory estimate.
Eight ranks (2x4 or 4x2) return valid queries with 148.92-GiB array floors;
sixteen ranks (4x4) require 148.99 GiB. Both exceed the current host's 123.5 GiB.
The recorded raw queries, per-rank sums, source/licence copy and commands are
under `diagnostics/cas12-scalapack-workspace-audit/` in the calculation evidence.

A dense retry needs a larger-RAM host, a process grid with valid workspace
queries, and headroom above this array floor. A selected-spectrum eigensolver
would be a separate engine configuration requiring its own build/reference and
all-root/dipole gates; increasing only the current four-rank memory limit does
not establish a viable calculation.

## Selected-root target experiment

The pinned UKRmol-in 3.3.0 engine already includes a disk-backed Davidson
diagonalizer. The experiment's `--target-diagonalizer davidson-serial` selects
`igh=0, forse=1` in the **target** `&cinorn` template. `forse=1` is essential:
without SLEPc, the distributed dispatcher falls back to dense ScaLAPACK even
when `igh=0`. MPI still constructs the Hamiltonian, while the master performs
the disk-backed iteration. The default `auto` preserves upstream dispatch.

This is an external-engine configuration, not a new QSCAT eigensolver. It
retains the requested spin-adapted target roots and the same Hamiltonian,
orbitals, density/import checks and scattering diagonalizer. The selected
eigenvectors require O(N*k) memory; matrix construction and disk storage may
still scale quadratically and must be measured. The runner passes `crite`
and `maxiter` explicitly (defaults 1e-12 and 500). The pinned Davidson wrapper
hardcodes vector/residual/orthogonalization controls to 1e-10/1e-8/1e-7;
changing those namelist entries would not tighten that wrapper's controls.

Before using the route for larger targets, compare all eight sectors against
the independently gated dense CAS(10,8) and tight CAS(10,11) controls. Require
every requested root to agree within 1e-7 Hartree and the ground dipole within
1e-5 a.u.; retain the fresh-QC spin, Pi and orthogonality gates. Verify that
all eight raw logs identify Davidson, that iteration succeeds without fallback,
and that requested eigenpairs reach DENPROP and the usual analyzer. Probe
five/eight target roots and tighter `crite` at fixed QC ensemble, recording
root coverage and common-root differences. Record wall/CPU time, container
memory and disk peaks, including failed controls. Only after those gates may
CAS(10,12) receive a full local target-import trial. Scattering retains its
all-spectrum dense path and needs a separate actual-dimension memory audit.

The initial Davidson control is rejected: five/eight/sixteen/thirty-two roots
miss required roots or return inconsistent energies despite IERR=0. The first
singlet-A1 stored Hamiltonian independently diagonalizes to the expected
spectrum. Larger Davidson jobs remain gated out.

The next experiment enables the engine's existing SLEPc/PETSc support with
Docker build argument `WITH_SLEPC=ON`. Those shared libraries and Fortran
modules come from the same digest-pinned toolchain; record their byte digests
in build provenance and repeat the upstream/reference gates. The runner's
`--target-diagonalizer slepc` selects `igh=-1, forse=0` and the upstream
distributed Krylov–Schur backend. In this pinned **target** path,
`CI_Hamiltonian_module` initializes `MAT_DENSE`, and `SLEPCMatrix_module`
calls `MatCreateDense`. It retains an O(N²) distributed real Hamiltonian,
while avoiding the all-spectrum ScaLAPACK workspace and full eigenvector
matrix. The generic engine's sparse-storage support does not make these target
jobs sparse. Reject a
silent dense fallback in raw-log analysis. The backend warns about degeneracy;
all-root import, Pi and root-count refinement checks remain mandatory.

The completed five-root CAS(10,11) comparison passes all 40 imported roots
within 5.34e-10 Hartree and the ground dipole within 8.04e-11 a.u. Its complete
wall/peak are 3152.81 seconds / 5.25 GiB, versus 5865.15 seconds / 18.42 GiB
for the dense control. Eight SCATCI stages total 351.70 versus 3248.81 seconds
(89.17% lower), but DENPROP takes 2365.46 seconds, limiting the observed full
wall reduction to 46.25%. Different concurrent workloads/CPU placement make
these observed comparisons, not controlled MPI scaling.

The eight-root-per-sector CAS(10,11) refinement at `crite=1e-13` also passes
all 40 ensemble-root imports within 5.33e-10 Hartree and the dipole within
7.90e-11 a.u. Complete wall/peak are 3086.61 seconds / 5.27 GiB; the 299.95-second
SCATCI total and 2409.43-second DENPROP total again have different CPU placement
and concurrent workloads from the earlier controls. Common excitations match
the dense and five-root SLEPc references at printed precision.
Sixteen independent fixed-orbital PySCF probes, eight roots at trial spaces
80/160 in all sectors, fully converge and check **every one of the 64 requested
UKRmol roots** within 5.34e-10 Hartree. Raw CASCI energies, root ordering and
spins are verified. The repaired scan takes 128.22 seconds at a 0.280-GiB peak;
its failed pre-probe process-identity guard and repaired launch remain evidence.
This completes the CAS(10,11) extra-root solver gate for a larger target trial;
active-space/basis/ensemble convergence and the scattering path remain separate.

For the largest CAS(10,12) target dimension 70860, one dense Hamiltonian alone
requires `8*N² = 40169116800 bytes` (37.41 GiB), or about 9.35 GiB per rank
on four ranks. The local trial therefore uses a 64-GiB container and 16-GiB
internal matrix budgets, conditional on 80 GiB available host RAM, 20 GiB
free scratch and the full extra-root/import gates. This is a matrix floor,
not a measured total-job peak. The all-spectrum scattering route still needs
its own actual CONGEN dimensions and workspace audit.
