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
