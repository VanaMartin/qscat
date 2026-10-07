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

Both small 80/160-vector controls now pass all forty imported roots within
4.98e-10 Hartree and dipoles within 4.50e-11/3.07e-11 a.u. Their printed
excitations agree, with a maximum change of 2.10e-8 Hartree from the earlier
dense control. Both stretched-DZ CAS(10,11) 80/160-vector reoptimizations now
pass QC in 389.55/363.26 seconds. Both projected equilibrium TZ retries remain
CI-unconverged in singlet B1/B2 despite orbital convergence (807.23/806.98
seconds). Larger trial space alone does not repair that projected solution.
Independent lowest-root coverage and paired subspace checks of the passing
stretched targets are queued; fresh/projected orbital comparisons remain open.

## Completed CAS(10,11) scattering and extraction checks

The tight DZ/40-component l=4 baseline completes at contracted B1/B2 dimensions
27546 (raw CONGEN counts 344124), in 4.423 hours at a 25.20-GiB kernel peak.
Its default native candidate is at 2.520805 eV with full width 1.151002 eV.
Relative to the earlier CAS(10,10) l=4 run, position and phase changes of
64.84 meV / 0.109809 rad exceed the chosen 0.05-eV / 0.05-rad gates; width
changes by 2.80%. QC/fresh-CI controls also differ, so this is a recorded-model
comparison rather than an isolated active-space perturbation.

The numerically matched tight CAS(10,10) control completes in 1282.76 seconds
(21.38 minutes) at a 2.864-GiB kernel peak. All forty independently imported
roots pass within 4.69e-10 Hartree and the dipole within 3.28e-11 a.u. Every
ensemble root passes all 48 final-orbital coverage probes within 9.95e-14
Hartree; the extra eighth singlet B1/B2 roots fail at space 40 and repair at
80/160. The candidate is **2.455964890/1.119642705 eV** in position/full width.
Relative to the loose CAS(10,10) control, tightening changes position/full width
by only −0.544/+0.395 micro-eV and phases by at most 1.30e-6 rad. Matching the
numerical controls therefore leaves the CAS(10,10)→CAS(10,11) position shift
**+64.8401 meV** and phase difference **0.1098099 rad** resolved; both exceed
the chosen gates. Full width changes by +31.3589 meV (about +2.80%). The active
counts and same-active-space seeds differ; all QC/continuum/channel controls
match. Newly serialized automatic-diagonalizer defaults are inactive, and all
eight generated target SCATCI inputs are byte-identical. This qualifies the
numerical comparison, while orbital continuity and electronic/continuum model
selection require their remaining checks.

Twelve background/detection replays span background terms 1–4 and detection
thresholds 0.7/1.0/1.3. Position changes remain below 6.65 meV, but width
changes by 12.65% from default, exceeding the 5% gate. Twenty-four additional
replays cover both Pi components, terms 1–4, detection 1.0 and automatic or
clipped 1.8–3.3 / 2.0–3.15-eV intervals. At fixed background, window clipping
changes position by at most 6.03 meV and width by 3.62%, within the chosen gates.
This smaller window sensitivity does not qualify the unresolved background
dependence or certify a pole assignment.

The pinned `source/libouter/reson.f`, `RESONC` lines 868–915, intersects the
automatically detected interval with `ABVTHR`/`BELTHR` relative to the adjacent
thresholds. With `GETETA=.false.`, it uses saved grid points without interpolation.
Raw fit grids contain 39/30/23 points, approximately 1.6–3.5 / 1.8–3.25 /
2.0–3.1 eV; the requested endpoints need not coincide with retained points.
The engine source/manual/licence, exact native grids and unit conversion are
preserved. A standalone B2 replay initially binds the wrong K-matrix unit
(921 instead of 922) and fails with EOF; its original exit/log and a fresh
template-unit-aware retry remain evidence. Both default fits reproduce the
original scattering outputs. Continuum controls and electronic-model selection
remain pending.

### Independent phase-fit diagnostic

`projects.ukrmol_co.phase_fit` checks the native one-candidate fit independently
on the retained energy grid. Its input energy/width unit is Hartree; phases are
radians. Adequate sampling is assumed: adjacent physical phase changes must
be smaller than pi/2 for modulo-pi unwrapping to recover the correct branch.
Below the first excited threshold, the diagnostic model is

\[
 \delta(E)=\operatorname{atan2}(\Gamma/2,E_r-E)
             +\sum_{j=0}^{n_{\rm back}-1}b_j x^j,\qquad
 x=(E-E_{\rm mid})/s,\quad s=(E_{\max}-E_{\min})/2.
\]

The resonant factor has positive full width \(\Gamma\). On the real energy axis
\(S_{\rm res}=(E_r-E+i\Gamma/2)/(E_r-E-i\Gamma/2)\) has unit modulus, and its
phase rises by pi. Fit phases are unwrapped modulo pi. A constant phase branch
is absorbed by the background constant. Centering the polynomial changes its
coefficients, not its span, and avoids powers of small Hartree energies.

For each trial \((E_r,\log\Gamma)\), linear least squares eliminates the
background coefficients. A bounded nonlinear least-squares solve then minimizes
the phase residual. The center must lie in the fitted interval, and the width
is bounded between 1e-6 and twenty times its half-range. Boundary solutions are
explicit diagnostics. The fit reports residual RMS/maximum/sum of squares and
the singular values of the full, scaled parameter Jacobian. These are measures
of fit quality and local identifiability, not statistical uncertainty estimates
for deterministic electronic data.

Validation uses an independently constructed analytic unitary S factor with a
known position/width, several polynomial backgrounds, three energy meshes and
phase wrapping/branch shifts. Recovery must agree within 1e-9 Hartree, and
predicted phases within 1e-9 rad. Differential checks against native RESON use
the native input grid reconstructed from EIGENP, matching printed fit points
within their Rydberg rounding error; printed residues are accurate to 5e-8 rad.
Interleaved held-point predictions distinguish background inadequacy from an
in-sample improvement. They do not establish a scattering pole, a global
background model or a physical error bar. Existing position/width/phase gates
and independent electronic/continuum checks continue to apply.

The pinned `source/libouter/rsolve.f` line 142 and `eigenp.f` line 77 use
`RYD=0.073500` to convert input eV to Rydberg and back. Their labeled eV energy
grid therefore corresponds to physical energies larger by a factor
**1.0000184445400588**, using the project's Hartree-to-eV constant. This is
18.44 ppm, or approximately 46.5 micro-eV near the candidate. The diagnostic
verifies both source literals and uses `E_input_eV * 0.0735 / 2` Hartree as its
actual fit abscissa. Native fitted Rydberg positions/widths are converted using
`qscat.units`, as before. The supplied window clips and printed grid labels
must be distinguished from exact physical eV. Source and native values remain
pinned; this audit does not change the historical engine's constants.

The analytic unitary-S/mesh/branch and held-point-background controls pass
(14 tests). Independent fits to all 24 native replays agree in position/full
width within **4.74e-7 / 2.50e-7 eV**. Three separated initial guesses per fit
agree within 1.55e-9 eV in position and 5.90e-10 in relative width. The native
printed residues and their sum-of-squares goodness are independently checked.
Five interleaved held-point fits per replay give the following automatic-window
B1 diagnostic; B2 reproduces the comparison:

| Background terms | In-sample RMS (rad) | Held-point RMS (rad) | Held-point maximum (rad) | Scaled Jacobian condition |
|---|---:|---:|---:|---:|
| 1 (constant) | 0.016595 | 0.016926 | 0.040083 | 5.90 |
| 2 (linear) | 0.001095 | 0.001147 | 0.003584 | 12.44 |
| 3 (quadratic) | 0.000469 | 0.000510 | 0.001396 | 24.06 |
| 4 (cubic) | 0.000102 | 0.000128 | 0.000523 | 57.75 |

The constant background predicts the retained points substantially worse;
additional terms improve held-point predictions as well as the fitted points.
Terms 2–4 change the automatic-window width by about 1.03% from default.
The shorter windows improve residuals but increase Jacobian conditioning, up
to 166 for the cubic 2.0–3.15-eV clip. These diagnostics explain the background
sensitivity without selecting a new acceptance threshold after seeing the data,
discarding the constant-background evidence or supplying physical error bars.

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
