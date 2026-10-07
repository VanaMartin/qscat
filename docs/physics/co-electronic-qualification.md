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
The 96 independent stretched-target probes now complete: 87 fully converge,
including all 64 at spaces 80/160. Three five-root space-40 controls fail the
fifth ensemble-root flag, rejecting both targets under the existing strict
all-probe gate. Their first-five spectra agree within 9.95e-14 Hartree, but
energy agreement does not remove failed convergence flags. Pair roots/objective/
dipole and core/active subspaces agree. Fresh/projected orbital comparisons,
independent imports and qualification of these repairs remain open.

The iteration-refinement follow-on repeats all nine failed settings at the
original 200 and refined 600 cycles. All nine original failures reproduce and
repair at 600 without changing orbitals, roots, trial spaces or tolerances.
Full 48-probe scans on both checkpoints then pass every requested root,
convergence flag and spin. Explicit physical/spin-penalized residual norms
`||(H-E)c||/||c||` stay below 7.07e-10/9.94e-10 Hartree, respectively, and
ensemble energies agree within 1.28e-13 Hartree. An analytic Hubbard-dimer
eigenstate and perturbed-vector control validate the residual calculation.
The worker takes 789.89 seconds / 0.285 GiB. This qualifies fixed-orbital
coverage under the refined 600-cycle contract; the original 200-cycle gate
rejections remain preserved. The validated refinement releases an independent
64-root stretched import and same-geometry TZ/aug-TZ QC, with further
electronic-model and competing-start checks still required.

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
numerical comparison. The matched CAS(10,10) initial/final core/active subspaces
also pass continuity, with minimum overlaps 0.9999999999996424 /
0.9999999998858343. Electronic/continuum model selection requires its remaining
checks.

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

### Outer-only energy-grid refinement

Retained CAS(10,11) channels and R-matrix amplitudes permit new outer-region
solves without rebuilding the inner Hamiltonian. The 0.1–5.0-eV requested
interval now has three step sizes: **0.05/0.025/0.0125 eV**, or 99/197/393 points.
Both Pi components have complete, finite grids, nonnegative cross sections and
consistent final-state sums. Their phases at the 99 common points match the
baseline exactly at printed precision. The four finer-grid outer pipelines
take **381.99 seconds / 0.719 GiB** together on four physical cores.

Independent linear-background fits hold both endpoints fixed in labeled input
eV, converting those points through the audited native energy convention:

| Fixed interval (input eV) | Fit points, coarse→medium→fine | Medium−coarse position (meV) | Fine−coarse position (meV) | Medium−coarse width (meV) | Fine−coarse width (meV) |
|---|---|---:|---:|---:|---:|
| 1.6–3.5 | 39→77→153 | +0.09591 | +0.14459 | +0.19120 | +0.29127 |
| 2.0–3.1 | 23→45→89 | +0.06606 | +0.10025 | +0.17973 | +0.27739 |

Both components pass the chosen position/width/phase gates. The largest width
change is 0.0254%, and successive medium-to-fine changes are smaller. This
qualifies these energy meshes for this fit representation; it does not resolve
the background dependence or electronic/continuum model selection.

Native RESON retains `MAXFIT=100` (`source/libouter/reson.f` lines 32, 899–915).
The fine automatic fits truncate the saved-point interval at that limit and
print a warning. Both truncated records remain evidence, with native
position/width acceptance left unset. Independent fits consume the full fixed
intervals. The archived driver failures include both MPI input-protocol EOFs,
a pre-container pause-acknowledgement race and a completion-message mismatch
after four successful native stages. Fresh retries preserve their sources,
original exits and CPU-slot resumption records.

### Outer-only propagation refinement

At fixed CAS(10,11) inner amplitudes and the original 99-point grid, refine
radial propagation from 8→16→32 subranges at a 100-bohr matching radius.
Then increase the Gailitis matching radius to 150/200 bohr with 26/36
subranges, keeping interval lengths at most 5.125 bohr and Legendre order ten.
All forty native stages for both components succeed in **1669.54 seconds /
1.868 GiB**. Grid, cross-section and Pi checks pass.

The full-precision binary K-matrix audit resolves the identical printed
refined results: maximum baseline phase change is **2.56845e-6 rad modulo pi**,
and successive refined differences stay below **2.28e-10 rad**. Native phases
reconstruct as `sum(arctan(eigvalsh(K)))` within 4.87e-8 rad, their rounding
precision. Linear-background fits of the binary phases over the fixed
1.6–3.5/2.0–3.1-input-eV intervals change position/full width by at most
**2.63e-8 / 4.35e-8 eV**. These pass the existing numerical gates.

Every requested sector is verified against all four rank logs. The pinned
source audit follows `BPROP` through `RPROP1_MPI`, `CURLYR`/`GAILIT` and
`RPROPX` to the K-matrix solve, confirming an active refinement family.
The [propagation companion](../../projects/ukrmol_co/qualification-evidence/outer-propagation/README.md)
retains raw binary K data, all rank logs, sources and the reconstruction.
The plateau is resolved well below the chosen phase tolerance; electronic-model,
continuum/channel and fit-background qualification remain open.

## Matched MPI scaling and concurrent throughput

The tight equilibrium CAS(10,10)/cc-pVDZ/40-channel benchmark completes all
seven replicas with original exits zero. On CPUs 12–15, one/two/four ranks
take **3105.91 / 1830.62 / 1239.61 seconds**, with aggregate CPU **3105.86 /
3435.09 / 4268.46 seconds** and kernel peaks **2.849 / 2.743 / 2.822 GiB**.
Full-job speedups are **1.697× / 2.506×**, whereas the scattering SCATCI
stage speeds up **1.980× / 3.764×**. This measures the complete pipeline's
nonparallel costs rather than assuming eigensolver scaling applies everywhere.

Four one-rank replicas on CPUs 12/13/14/15 finish in **3444.36 seconds**, giving
**1.440× throughput** relative to repeating the measured four-rank runner.
Aggregate CPU is **13675.56 seconds**, versus **17073.86 seconds** for four
repetitions, and summed per-job memory peaks bound the concurrent envelope
at **11.342 GiB**. These are single-shot measurements on one geometry with
recorded concurrent host work; they do not establish larger-model scaling.

All seven target/import, grid, cross-section-sum and Pi gates pass. Both
contracted dimensions remain **8350**. Different initial-checkpoint lineage is
checked explicitly: QC roots/dipoles match the baseline within **6.119e-9
Hartree / 2.986e-8 a.u.**, and core/active subspaces agree to roundoff.
All forty import roots/dipoles agree with each replica's QC within **5.044e-10
Hartree / 4.295e-11 a.u.**. Maximum phase and native position/full-width
changes are **2e-7 rad / 0.137 micro-eV / 0.041 micro-eV**, well inside the
unchanged numerical gates. The
[benchmark companion](../../projects/ukrmol_co/qualification-evidence/mpi-throughput/README.md)
reanalyzes all seven replicas and the exact published baseline and verifies
98997 raw profile samples. Execution-layout equivalence does not resolve
active-space/basis/background-width sensitivity or certify a resonance pole.

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

## Fresh larger-basis QC starts

The fresh-RHF-start equilibrium batch now completes: CAS(10,11)/cc-pVTZ and
CAS(10,12)/cc-pVTZ/aug-cc-pVTZ pass tight QC in 39.75/112.87/91.74 minutes.
CAS(10,11)/aug-cc-pVTZ is rejected after 22.80 minutes: orbital optimization
and its optimized CI converge, but fresh singlet-A1 CI does not. Matching
energies within 1.42e-13 Hartree cannot replace that convergence flag.
The [fresh-basis companion](../../projects/ukrmol_co/qualification-evidence/fresh-basis/README.md)
retains all four attempts and reconstructs 160 raw final-state energies/spins.

For fresh cc-pVTZ starts, CAS(10,11)→CAS(10,12) changes ground energy by
−0.391782 eV and z dipole by +0.152827 a.u. For fresh CAS(10,12), TZ→aug-TZ
changes ground energy by +0.479725 eV, ensemble objective by −1.795712 eV
and z dipole by −0.041475 a.u. The ensemble objective is the variational
quantity; its reduction does not enforce a lower individual ground root.
Minimum initial-to-final active overlaps are 0.171594/0.152714/0.036790.
These significant changes require competing-start and orbital/state continuity
checks before assigning a basis-convergence interpretation.

The independent 600-cycle coverage/residual scans now pass all **144 probes /
936 eigenpair evaluations**. All flags, spins, CI-vector orthogonality and
physical/spin-penalized residual gates pass; residual maxima are **6.96e-10 /
9.99e-10 Hartree**, respectively. Maximum trial-space energy difference is
3.13e-13 Hartree and common first-five root-count difference is 3.98e-13
Hartree. Wall is **32.36 minutes / 0.728 GiB**. The
[coverage companion](../../projects/ukrmol_co/qualification-evidence/fresh-coverage/README.md)
reconstructs every raw spectrum/spin/convergence block and retains the unchanged
parent checkpoints and frozen measurement source. This validates lowest-root
coverage on those fixed orbitals; competing starts and all-root/dipole imports
remain required. The covered fresh CAS(10,11) TZ target's 64-root UKRmol import
is queued with independent eight-root/space-160 controls and the existing gates.

Both space-80/160 same-basis restarts of the rejected fresh CAS(10,11) aug-TZ
checkpoint now pass QC in **641.47/638.23 seconds**, with kernel memory peaks
**0.419/0.477 GiB**. All forty roots agree within **1.667e-11 Hartree**, dipoles
within **3.230e-12 a.u.**, and core/active subspaces to roundoff. The
[repair companion](../../projects/ukrmol_co/qualification-evidence/fresh-augtz-repairs/README.md)
reconstructs all 80 raw final states and preserves the rejected parent exactly.
The default 200 CI cycles and physical/ensemble/tolerance settings are unchanged.
These are numerical repairs of one fresh lineage; competing starts and all-root
imports remain gates. A 96-probe, 600-cycle fixed-orbital coverage/residual scan
is queued for both checkpoints after the active fresh TZ import.

The fresh CAS(10,11) TZ import now completes in **54.38 minutes / 5.270 GiB**.
All **64 native roots** pass within **3.196e-8 Hartree** against the independent
seed-orbital controls. Import-run QC changes the original seed roots by at most
3.109e-8 Hartree; the all-64 comparison includes this small reoptimization.
All forty imported ensemble roots agree with import-run QC within 5.253e-10
Hartree, and the DENPROP dipole agrees within 3.366e-11 a.u. Core/active
subspaces agree to roundoff. DENPROP takes 2350.20 seconds of the total
3262.90-second run, a larger cost than the eight SCATCI sectors' 371.32 seconds.

The [import companion](../../projects/ukrmol_co/qualification-evidence/fresh-tz-import/README.md)
retains a successful engine run with original supervisor exit one: its saved-table
consistency check was tighter than the writer's nine-decimal precision. A first
recheck also rejects echoed earlier CIDATA spectra. The passing recheck selects
the newly solved sector and uses an exact-decimal half-last-unit format bound
of 5e-10 Hartree, with corrupt-token/missing-sector controls. These are parser
and format checks; physical root/dipole gates remain 1e-7 Hartree / 1e-5 a.u.
Both verifier failures retain their original source and exits. The augmented-TZ
coverage worker advances independently, followed by a coverage-gated all-root
import. Competing-start and electronic-model gates remain open.

The augmented-TZ coverage worker subsequently passes all **96 probes / 624
eigenpair evaluations** on both repaired checkpoints in **10.02 minutes /
0.391 GiB**. Physical/spin-penalized residual maxima are **6.62e-10 / 9.94e-10
Hartree**; trial-space and common-first-five root-count differences stay below
1.99e-13 Hartree. The
[coverage companion](../../projects/ukrmol_co/qualification-evidence/fresh-augtz-coverage/README.md)
reconstructs every raw spectrum/spin/convergence block. Its unchanged repaired
seed and pair gates release the now-active 64-root aug-TZ import. This closes
fixed-orbital coverage for those repairs, while independent import, competing
starts and electronic-model convergence remain open.

The repaired aug-TZ import subsequently passes all **64 native roots** against
the covered space-80 seed within **5.969e-11 Hartree**, and against the separate
space-160 seed within **5.764e-11 Hartree**. All forty ensemble import roots
agree with import-run QC within **4.672e-10 Hartree**, and the DENPROP dipole
within **5.169e-12 a.u.**. Core/active subspaces match to roundoff; all original
engine/supervisor/verifier exits are zero. Cost is **56.01 minutes / 5.268 GiB**,
including **2319.03 seconds** for DENPROP and **394.14 seconds** for eight SCATCI
sectors. The [import companion](../../projects/ukrmol_co/qualification-evidence/fresh-augtz-import/README.md)
reconstructs every native root and both seeds' 96 raw reference spectra.

### Competing-start failure and the new aug-TZ branch

Four fixed-geometry TZ↔aug-TZ projections retain the same core/active sizes,
forty-component ensemble, tight tolerances and ordinary 200-cycle CI limit.
Both upward projections pass QC in **1737.00 / 1705.00 seconds**, with a
space-80/160 pair agreement of **5.941e-12 Hartree** in roots and
**1.901e-11 a.u.** in dipole. Both fail comparison with the earlier fresh
aug-TZ solution: averaged objective **−0.538388 eV**, ground energy **−1.115756
eV**, maximum root shift **2.126321 eV**, z dipole **−0.06670992 a.u.** and
minimum active overlap **0.141085**. These are differences between separately
converged orbital solutions in the same basis/ensemble, not basis convergence.
The prior coverage/import checks remain valid for their fixed orbitals but do
not establish a unique optimum. The new solution's coverage/import remain gates.

Both downward projections fail before QC because PySCF 2.11.0 rejects the full
92-column source matrix in a 60-orbital destination. The adapter now supplies
only its thirteen core+active columns in this case and lets PySCF construct
the destination virtual complement. The exact failed checkpoint control retains
the requested active irreps and has **1.932e-14** orthogonality error. Four genuine
PySCF regressions exercise the repaired interface; simultaneous geometry/basis
changes remain rejected. Original batch/controller exits stay one. All
[277 competing-start payloads](../../projects/ukrmol_co/qualification-evidence/competing-starts/README.md)
verify, including both exceptions, 80 raw final states and independent
checkpoint-based reconstruction of all three subspace comparisons.

Both new aug-TZ checkpoints subsequently pass independent 600-cycle fixed-orbital
coverage: **96 probes / 624 eigenpair evaluations**, maximum physical/penalized
residuals **6.630e-10 / 9.977e-10 Hartree**, with original exit zero. Ensemble,
trial-space and common-first-five root-count differences are at most **2.274e-13 /
1.706e-13 / 1.564e-13 Hartree**. The
[335-payload coverage companion](../../projects/ukrmol_co/qualification-evidence/competing-augtz-coverage/README.md)
reconstructs every raw spectrum/spin/convergence block. Cost is **10.17 minutes /
0.380 GiB**. Four new-name downward retries from the fresh and competing
aug-TZ branches now run, with a separately gated 64-root SLEPc import queued
behind their completion. Orbital-branch disagreement remains a scientific result.
The qualified fixed-DZ l=3/l=5 continuum controls continue independently of
CAS(10,12) import, which subsequently completes. Electronic-model convergence
remains open.

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
free scratch and the full extra-root/import gates. This is a matrix floor;
the completed local job now measures its total peak below. The all-spectrum
scattering route still needs its own actual contracted dimensions and workspace audit.

### Completed local CAS(10,12) target import

The tight-restart DZ/forty-component target completes with original engine,
batch and controller exits zero in **22324.47 seconds / 6.201 hours**, aggregate
CPU **9.437 hours**, kernel peak **40.405 GiB**. Five roots per singlet/triplet
sector all pass: native-current-QC error **4.801e-11 Hartree**, saved-current-QC
error **5.465e-10 Hartree**, covered-seed first-five error **6.815e-9 Hartree**,
independent RHF-start error **3.485e-8 Hartree**. The DENPROP ground dipole agrees
with current QC within **1.279e-11 a.u.**, with the seed within **3.810e-9 a.u.**
and the independent RHF start within **5.881e-8 a.u.** Root/dipole physical gates
remain 1e-7 Hartree / 1e-5 a.u.; final-CIDATA-set parsing and nine-decimal saved
rounding consistency remain separate checks. Core/active subspaces match the
seed to roundoff and the separate RHF-start active space to **0.9999999999995508**.

All forty ensemble roots pass the original **48** fixed-orbital reference probes;
two extra eight-root/space-40 failures remain failures and their larger-space
controls pass. This import computes forty native roots; it does not assert the
extra-root import gate at CAS(10,12). The
[335-payload companion](../../projects/ukrmol_co/qualification-evidence/cas12-import/README.md)
reconstructs all native spectra/dipole, forty raw current QC states, 48 raw
reference spectra/spin/convergence blocks and 111339 resource samples, with
independent checkpoint-based seed/RHF-start subspace comparisons.

Eight SCATCI stages total **2808.79 seconds / 46.81 minutes**, QC **1472.31
seconds / 24.54 minutes**, and serial DENPROP **18030.97 seconds / 5.009 hours**,
or **80.768%** of profiled wall. The largest dense PETSc matrix remains dimension
70860 / 37.41 GiB before overhead. Numerical target import is feasible on
Sadaharu; the all-spectrum scattering path, its contracted dimensions and
electronic-model convergence retain separate gates. The stretched CAS(10,11)
import/basis queue advances on released CPUs 8–11.
