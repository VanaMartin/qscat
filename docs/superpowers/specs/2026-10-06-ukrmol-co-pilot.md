# UKRmol+ CO fixed-nuclei pilot

## Scope and physical model

Use the installed UKRmol+ suite as an external numerical engine. Generate a
neutral CO RHF target with Psi4 and calculate its electron-scattering doublet
B1/B2 symmetries in C2v, containing the two components of the Pi sector and
higher odd angular projections. SEP permits one
valence target electron to excite into the virtual space. Freeze the two core
orbitals and retain the other five occupied orbitals. Begin at R = 2.1323 bohr
with aug-cc-pVDZ and 20 symmetry-balanced virtual orbitals.

Geometry and propagation distances use bohr, and internal target/resonance
records use Hartree. Explicitly label the external electron-energy grid in eV
and cross sections in bohr². The initial 491-point energy grid covers 0.1–5 eV.

## Execution interface

`python3 -m projects.ukrmol_co.run --workdir PATH` creates a new immutable-input
run directory, obtains the checksum-pinned UKRmol-scripts release, prepares the
CO deck, runs it and records resources. Options vary bond length, SE/SEP model,
virtual counts by symmetry, continuum radius/angular cutoff, precision,
deletion threshold, MPI ranks and energy grid.

`python3 -m projects.ukrmol_co.analyze PATH` validates completed outputs and
writes a small JSON summary with target energy, native resonance fits and
resource measurements. No calculation outputs are promoted into QSCAT during
this pilot.

## Success criteria

- All requested target and scattering stages complete successfully.
- The full requested electron-energy grid is present for B1 and B2.
- Cross sections and eigenphase sums are finite, and cross sections nonnegative.
- Pi components agree modulo pi to native text precision, and their cross
  sections agree to an explicit relative tolerance.
- Psi4 and UKRmol target energies agree within the tolerance justified by the
  SCF/integral approximations used.
- Resource records preserve wall time, per-stage time, aggregate container
  memory and allocated persistent working-directory disk usage.

An automatic resonance fit is useful evidence only when its position and width
are positive, finite and stable under model/continuum refinement. Pipeline
success alone does not establish convergence or validate a reduced local
potential. The subsequent three-geometry pilot uses R = 1.9, 2.1323 and 2.5
bohr after the equilibrium deck passes these checks.

## Convergence extension

The equilibrium SEP virtual-space check failed a useful production tolerance.
The calibration therefore varies frozen occupied space (two versus four A1
orbitals), atomic basis, continuum angular cutoff, sphere radius, deletion
threshold and propagation radius. Independent Docker jobs are CPU-pinned
and profiled in separate cgroups.

An upstream contracted CAS-A close-coupling diagnostic uses two frozen core
orbitals, ten active electrons in [4,2,2,0] or [4,3,3,0] active orbitals,
and equal singlet/triplet root counts in all four C2v irreps. HF orbitals give
a small first diagnostic; ground-state CASSCF orbitals are an additional
choice supported by the upstream Psi4 template. These are distinct from
published multi-irrep state-averaged CASSCF models. Explicit excited channels
and a correlated CAS target permit tests of model balance.

Provisional testbed acceptance criteria, assessed at compressed, equilibrium
and stretched sentinel geometries, are changes below 0.05 eV in resonance
position, 5% in width (with a 0.001 eV absolute floor), and 0.05 radians in
the symmetry eigenphase sum modulo pi over the fit window. Fits must also
survive energy-grid and window changes. Electronic model uncertainty is
assessed separately by active-space, atomic-basis and channel-count changes;
agreement with a single published energy is insufficient. Near the crossing,
pole/threshold diagnostics must replace an unconstrained broad-resonance fit.
The neutral curve requires a separate correlated electronic calculation.

The finite batch runner snapshots project source and records image/source
hashes, disjoint CPU groups, aggregate cgroup CPU time, sampled anonymous/page
cache memory, kernel peak memory and allocated run disk. Existing K-matrices
can be replayed through RESON with varied background order/detection threshold.
Calibration snapshots include failed attempts and distinguish engine failures
from later analyzer fixes.

The Docker packaging adds a checksum-pinned source build of UKRmol-in 3.3.0
(including its bundled GBTOlib) and UKRmol-out 3.3.0.1. It reuses the
digest-pinned upstream CPU toolchain/Psi4 layer; source-built double-precision
executables and libraries replace the reference engine after serial/parallel
upstream checks. The `pilot` target retains the measured reference engine.

The durable execution and measured-result descriptions are in
`docker/ukrmol-plus/README.md` and `projects/ukrmol_co/README.md`. Generic engine
targets are separate from the `co-pilot`/`co-source` experiment layers; the
default build target is `co-source`.
