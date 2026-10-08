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

The subsequent R=2.5-bohr/DZ import now passes all **64 native roots** against
the independently covered space-80/160 seeds within **3.598e-10 / 3.836e-10
Hartree**. Forty native ensemble roots match current QC within **4.939e-11
Hartree**, and the DENPROP ground dipole within **1.691e-11 a.u.** Core/active
subspaces agree to roundoff against both seeds. Cost is **54.25 minutes /
5.269 GiB**, including **2519.30 seconds** for DENPROP.

The original engine/batch pass, while the supervisor and post-import verifier
retain exit one: a `1e-10-Hartree` saved-table comparison rejects valid rounding
from the `%16.9f` writer. The passing retained-data recheck uses the established
final-CIDATA-set parser and exact-decimal **5e-10-Hartree** half-unit format
bound; physical root/dipole gates remain **1e-7 Hartree / 1e-5 a.u.** The
[488-payload companion](../../projects/ukrmol_co/qualification-evidence/stretched-import/README.md)
reconstructs all roots, 40 raw QC states, both seeds' 96 reference spectra/624
eigenpairs, both orbital-subspace comparisons and 16229 resource samples.
Both subsequent same-geometry TZ/aug-TZ QC targets pass in **14.57 / 20.33
minutes**, at **0.345 / 0.434 GiB**, with all ordinary/fresh-CI and orbital
flags accepted. The [170-payload QC companion](../../projects/ukrmol_co/qualification-evidence/stretched-basis-qc/README.md)
reconstructs 80 raw states and 10430 resource samples. DZ→TZ and TZ→aug-TZ
ground-energy changes are **−0.825431 / −0.0248054 eV**; maximum ranked-excitation
changes **0.0883686 / 0.0755466 eV**, z-dipole changes **−0.00374061 /
−0.0110101 a.u.** Cross-basis minimum active overlaps are **0.985481 / 0.998189**.
Ranked spectra do not establish state identities. Both 48-probe basis
coverage scans now pass, and aug-TZ passes its 64-root native import against
the covered seed within **4.755e-8 Hartree**. TZ's original comparison rejects
the eighth triplet-A1 root at **1.564e-7 Hartree**, exceeding the unchanged
**1e-7** gate. A separate fixed-orbital diagnostic on both imports' actual
checkpoints passes every native root within **4.984e-11 / 5.009e-11 Hartree**
(TZ/aug-TZ), identifying seed-to-import orbital drift rather than a native
eigensolver error. The original seed-gate failure is retained. The
[1296-payload companion](../../projects/ukrmol_co/qualification-evidence/stretched-basis-imports/README.md)
reconstructs **128 probes / 880 eigenpair evaluations** and **34607 resource
samples**. Native imports take **49.78 / 52.02 minutes** at **5.265 / 5.261
GiB**. Basis, competing-start and across-geometry model qualification remain
open.

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

The completed finite near-equilibrium pilot adds **23 geometries per basis**, using
aug-QZ and aug-5Z over **1.9–2.5 bohr** at **0.025-bohr** spacing and reusing
the three original anchors. Each basis has 26 points including R=2.1323.
The [recipe](../../projects/ukrmol_co/calibration-neutral-near-equilibrium.json)
and [predeclared contract](../../projects/ukrmol_co/neutral-near-equilibrium-contract.json)
define independent midpoint checks against a 0.05-bohr knot grid plus the
reference geometry. Cubic not-a-knot and PCHIP interpolation use no
extrapolation; diagnostic budgets are **1 meV** in held-point energy,
**0.001 a.u.** in held-point dipole and **20 meV** in QZ→5Z relative energy.
These budgets test this labelled pilot. They do not certify the correlation
treatment or extend the stable-reference domain toward dissociation.

All **52 points** now pass the electronic/stability/density checks. Independent
reconstruction of the completed **46-job** queue gives maximum held-energy
errors **0.11386/0.11330 meV** for cubic QZ/5Z and **0.95428/0.91453 meV** for
PCHIP. Maximum held z-dipole errors are **1.84e-7 a.u.** for cubic and
**7.21e-6 a.u.** for PCHIP. The maximum QZ→5Z relative-energy difference is
**17.2433 meV**. All predeclared pilot budgets pass. New-job walls sum to
**3.210 hours**; maximum kernel container/anonymous-memory charges are
**24.000/17.138 GiB**. The
[public neutral companion](../../projects/ukrmol_co/qualification-evidence/neutral-near-equilibrium/README.md)
verifies 959 payloads and 63867 profile samples, preserves the six anchors and
four externally unstable stretched references, and reconstructs every midpoint
and basis comparison. Correlation-model/dissociation qualification remains open.

## Equilibrium comparison with published calculations

At **R=2.1323 bohr**, the strongest completed CAS(10,11)/cc-pVDZ/40-channel
baseline gives the native fit candidate **2.520805 / 1.151002 eV**
(position/full width relative to its neutral target). The matched CAS(10,10)
baseline gives **2.455965 / 1.119643 eV**. Both pass numerical target import
and pipeline checks; active-space/phase and background-width dependence
still prevent an electronic-model or pole qualification.

| Published electronic model | Position / full width (eV) | CAS(10,11) position difference | CAS(10,11) width difference |
|---|---:|---:|---:|
| Laporta et al. 2012, six-electron SEP | 1.67 / 0.82 | +0.8508 eV | +40.37% |
| Dora et al. 2016, cc-pVTZ/CAS(10,10), CC50 | 1.73 / 0.84 | +0.7908 eV | +37.02% |
| Dora and Tennyson 2020, cc-pV6Z/CAS(10,10), 41 C2v components | 1.8744 / 1.2916 | +0.6464 eV | −10.89% |
| Dora et al. 2016, cc-pVDZ/CAS(10,11), CC40 | 2.20 / 0.95 | +0.3208 eV | +21.16% |

Sources: [Laporta 2012](../../reference/literature/laporta-2012-psst21-045005.md),
preprint p. 4; [Dora 2016](../../reference/literature/dora-2016-epjd70-197.md),
p. 6, Table 4; [Dora and Tennyson 2020](../../reference/literature/dora-2020-jpb53-195202.md),
p. 4, Table 2. The last row is the closest nominal active-space/basis/channel
comparison; continuum and implementation choices still differ. These are
published model values with substantial inter-model spread, not a shared
experimental error bar. Laporta's tabulated width is unadjusted: the later
nuclear calculation applies a 10% width increase and 0.035-eV level shift
(preprint pp. 4, 6).

The earlier aug-TZ/four-frozen-orbital/39-virtual SEP check gives
**1.623439 / 0.768383 eV**, only **−46.56 meV / −6.29%** from Laporta's
unadjusted pair. This is closer agreement in the initial calibration but
does not qualify that model: increasing to 60 virtuals moves it to
**1.3292 / 0.5743 eV**. Selecting a virtual space by agreement would hide
the measured sensitivity. Raw values remain in
`projects/ukrmol_co/calibration-results.json`.

Target construction fares better in the direct compact-model plausibility
check. Our tight CAS(10,10)/cc-pVDZ ground energy **−112.8947375 Hartree**,
lowest triplet-Pi excitation **6.3952 eV** and dipole **about 0.234 D** are
close to Dora 2016's **−112.89473 Hartree / 6.40 eV / 0.234 D**
(p. 4, Table 2). Agreement of these properties supports the implementation;
the accompanying resonance difference remains unresolved. Independently,
the aug-5Z CCSD lambda-density neutral dipole is **about 0.124 D**, compared
with the **0.122 D** experimental magnitude quoted in that table. This is
a single-property comparison, not a CCSD(T) dipole derivative or full-curve
validation. Neutral basis refinement still shifts the R=1.9/2.5 relative
energies by approximately **−17.24 / +15.82 meV** from QZ to 5Z.

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
aug-TZ branches subsequently pass QC. Their two space-80/160 comparisons and
four comparisons against qualified fresh TZ pass all same-basis restart gates:
maximum root/dipole differences **4.156e-8 Hartree / 6.916e-8 a.u.**, minimum
active overlap **0.9999999999995466**, maximum objective difference
**8.527e-14 Hartree**. Both source lineages recover the existing TZ solution,
without identifying their distinct aug-TZ solutions. Individual QC walls span
**22.10–34.75 minutes**, peaks **0.352–0.409 GiB**, sequential batch **1.903 hours**.
The [296-payload companion](../../projects/ukrmol_co/qualification-evidence/downward-projections/README.md)
reconstructs 160 raw QC states, six complete checkpoint-based restart comparisons
and 34132 resource samples. All original batch/controller/comparison exits are
zero; the earlier interface failures remain preserved. The covered competing
aug-TZ 64-root import subsequently passes within **3.8811e-8 Hartree**, with
native/seed dipole error **1.8756e-7 a.u.**, in **3259.47 seconds / 5.260 GiB**.
Orbital-branch disagreement at aug-TZ remains a scientific result. The
[continuum/import companion](../../projects/ukrmol_co/qualification-evidence/continuum-competing/README.md)
reconstructs its exact checkpoint, spectra, table serialization and subspace checks.

### Completed equilibrium angular controls

At matched CAS(10,11)/DZ targets, l=3/4/5 common-grid comparisons pass the
declared numerical gates. Relative to l=4, l=3 changes position/full width by
**+3.9602 meV / +0.7434%**, with maximum phase difference **0.0166575 rad**;
l=5 changes them by **+0.07987 meV / +0.2655%**, with phase difference
**0.0162174 rad**. Independent complete-grid fits on common windows and
background orders 1–4 remain within **3.982 meV / 0.7214%** for l=3 and
**0.3930 meV / 0.2279%** for l=5. The above public companion reconstructs all
36 native background/detection replays and raw profiles. Background-width spans
are **12.60% / 12.65% / 11.98%** at l=3/4/5, still failing the 5% extraction
gate. Angular convergence therefore does not establish a qualified width.
The earlier fine-grid `MAXFIT=100` truncation verdicts remain unset.

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
scattering route has the separate completed dimension/workspace audit below.

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
Sadaharu; electronic-model convergence retains separate gates.

### Completed CAS(10,12) native scattering preflight

The [814-payload preflight companion](../../projects/ukrmol_co/qualification-evidence/scattering-preflight/README.md)
reuses the qualified targets and executes only fresh native integral generation
and B1/B2 CONGEN. The known CAS(10,11) completed native SCATCI dimension is
reproduced from retained-target contraction: **26136 L² + 1410 continuum =
27546** configurations, versus **344124 raw** CONGEN configurations.
CAS(10,12) measures **84942 L² + 1410 continuum = 86352** contracted
configurations and **990990 raw** configurations in each Pi sector.
Native preparation takes **172.99 seconds / 12.730 GiB**.

The installed-library all-spectrum workspace query at dimension 86352 needs
**222.391 / 222.555 GiB** aggregate arrays on 4×4/4×8 process grids, before
other live engine arrays or MPI/library overhead. One dense matrix alone
requires **55.556 GiB**. The 2×2 query overflows and retains exit five with
no usable estimate. Raw rank-local rows/columns, query work sizes, total bytes
and maximum-rank bytes reconstruct independently. The current dense full-spectrum
CAS(10,12) scattering path exceeds Sadaharu's RAM; actual full-job peak, runtime
and MPI scaling remain unmeasured. Increasing target selected-root efficiency
does not itself replace that scattering expansion. This audit does not rule out
a separately validated sparse/iterative scattering route.

Preparation retains the original scheduling rejection and four native restart
failures. The passing successor preserves exact target/checkpoint/CI hashes,
sets `LNDO=NBMX=100000000` separately from enlarged CONGEN dimensions, and
initializes retained sectors from unique newly solved native CIDATA set numbers.
The generated inventory requires five roots in each of eight sectors; DENPROP
execution remains skipped. Missing L² groups and changed inventories reject.
The companion verifies 245 source hashes and 1099 raw resource samples.

The finite CAS(10,11) three-anchor pilot continues with an independently covered
compressed 64-root import and separately gated compressed/stretched scattering.
It uses the equilibrium l=4/radius-18 model. The broad compressed and finer
near-threshold stretched grids diagnose geometry/extraction behavior; an empty
automatic candidate does not establish zero width or a bound state. State,
orbital and pole continuity remain required before a full resonant sweep.

At **R=1.9 bohr**, the [compressed companion](../../projects/ukrmol_co/qualification-evidence/compressed-anchor/README.md)
now reconstructs all **48** fixed-orbital residual probes, the **64-root** import
within **9.050e-9 Hartree**, and both complete 27546-dimensional dense scattering
sectors. Their native candidate is **3.826384 / 2.322242 eV** position/full width;
cost is **11439.93 seconds / 25.277 GiB**. Coverage/import/pipeline consistency
passes; active-space/basis/ensemble selection and geometry/extraction continuity
remain open.

The [original CAS(10,12) repair companion](../../projects/ukrmol_co/qualification-evidence/cas12-base-repairs/README.md)
preserves two compressed QC rejections, two passing stretched QC retries and
**96** original fixed-orbital **200-cycle** probes. Fifth singlet-A2/triplet-Pi
roots fail at the smallest trial space despite accurate energies. The separately
owned **600-cycle** residual successor reuses the previously validated protocol
with unchanged saved orbitals and physical/penalized-residual gates. Original
flags continue to reject their original scans; successors require their own
complete verification.

The [staged-target diagnostics](../../projects/ukrmol_co/qualification-evidence/staged-targets/README.md)
now preserve four original CAS12 QC rejections and both passing fifty-component
equilibrium pilots. On their common forty-root objective, ensemble changes are
**+0.170147 / +0.004204 eV** for CAS11/CAS12; reported ground-energy changes
are **+0.456512 / +0.007930 eV**. Independent fifty-component lowest-root coverage
and native import remain required before scattering interpretation.

The [CAS12 residual successor](../../projects/ukrmol_co/qualification-evidence/cas12-residual/README.md)
now independently verifies all **96** 600-cycle probes on the two unchanged
stretched DZ repair checkpoints. Maximum physical/penalized residuals are
**7.638e-10 / 9.997e-10 Hartree**, and ensemble/root-space differentials are
below **2.701e-13 Hartree**. Public reconstruction passes, preserving original
200-cycle failures. A finite independently gated TZ/aug-TZ QC pair and DZ
64-root native import now follow that qualified checkpoint; electronic-model
selection and scattering remain separate gates.

The [stretched CAS12 basis-coverage companion](../../projects/ukrmol_co/qualification-evidence/cas12-stretched-basis-coverage/README.md)
then checks **96** fixed-orbital TZ/aug-TZ probes at 5/8 roots, spaces40/80/160
and 600 cycles. TZ passes every probe. Aug-TZ rejects only the eighth triplet-A2
root at space40: CI flag false, physical/penalized residuals **9.41093e-7 /
1.38254e-6 Hartree**, despite accurate energies and passing ensemble roots.
The original verifier/controller exits one and remains rejected. Scan cost is
**1726.94 seconds / 0.7026 GiB**, with **354 payloads / 8603 profiles** publicly
reconstructed. A separately declared
[80/160/240 supported-space sequence](../../projects/ukrmol_co/cas12-stretched-augtz-supported-spaces-contract.json)
uses fresh probes at unchanged orbitals, cycles and tolerances, requiring its
own full reconstruction; it does not automatically replace the rejected import
proof. The independently qualified TZ branch retains its native-import eligibility.

The [fresh supported-space companion](../../projects/ukrmol_co/qualification-evidence/augtz-supported-spaces/README.md)
subsequently passes **48/48 aug-TZ probes** at spaces80/160/240: maximum physical/
penalized residuals **7.084e-10 / 9.986e-10 Hartree**, ensemble/space/root-count
errors **1.990e-13 / 2.274e-13 / 1.421e-13 Hartree**. Cost is **871.86 seconds /
0.8383 GiB**, and **269 payloads / 4346 profiles** publicly reconstruct. The
original space40 rejection and independent outer-wrapper failures remain public.
A [separate fresh import release](../../projects/ukrmol_co/cas12-augtz-supported-import-contract.json)
uses this supported production-space80 family after a successful original TZ
import/verifier, preserving the original rejected aug-TZ proof.

The [complete stretched outer-window companion](../../projects/ukrmol_co/qualification-evidence/stretched-window-controls/README.md)
replays both Pi sectors at R=2.5 bohr on 0.01–2.0-requested-eV grids with 0.005/
0.0025-eV spacing. All 20 native stages complete in **390.39 seconds / 0.6342
GiB**; **310 payloads / 1947 profiles** reconstruct. Common phases differ from
the prior shorter window by at most **2.0e-7 rad modulo π**. All 48 complete-data
fits converge with multistart/five-fold held-point checks, and grid-pair position/
width differences stay below **2.205e-6 / 4.293e-5 eV**. The all-background width
spread nevertheless reaches **6.72–6.74%**, rejecting the 5% gate. Terms 2–4 alone
span about 0.131%; improved held-point residuals motivate a separately justified
representation choice, not retroactive exclusion of the constant control. Native
fine-grid **MAXFIT** leaves its fit verdict unset. Threshold/geometry identity and
full extraction remain open.

#### Compressed extraction and three-anchor orbital screens

The [finite low-memory companion](../../projects/ukrmol_co/qualification-evidence/small-anchor-controls/README.md)
uses the completed CAS11 anchors at **R=1.9/2.1323/2.5 bohr**. At the compressed
anchor, three below-threshold windows (1.5–6.5, 2–6 and 2.5–5.5 requested eV),
background terms1–4, three starts and five held-point folds give **24 cases /
192 data fits / 120 held folds**. Synthetic exact recovery and phase-branch
invariance pass; every declared case is numerically valid and multistart-stable.
Nevertheless, the fitted position/full-width spans are **0.08993 / 0.74645 eV**,
with widths spanning **32.206%** of their **2.31773-eV** median. Both position
and width gates reject. Constant-background held errors reach **0.120 rad**,
also exceeding the 0.05-rad gate. Higher-background residual improvement does
not erase these controls. Both Pi sectors reproduce the same findings.

The saved-orbital screen rechecks native/QC thresholds, all forty ensemble
roots, Pi degeneracy and AO-metric orthogonality. Adjacent physical core/active
minimum overlap singular values are **0.9132/0.9693** and **0.8072/0.9319**;
the end-to-end pair is **0.6001/0.8249**. Physical core follow-up triggers are
preserved. A separately declared AO-following Lowdin frame uses
`X(R) = S_AO(R)^(1/2) C(R)` with identical ordered atom/basis/AO labels.
Singular values of `X(R_left).T @ X(R_right)` give adjacent core/active minima
**0.9986/0.9873** and **0.9983/0.9873**, end-to-end **0.9940/0.9500**.
This contrast diagnoses sensitivity to motion of the atom-centred basis;
the AO-indexed transport is an explicit coefficient-frame identification,
not a physical wavefunction overlap. Identity, reverse-pair reciprocity and
subspace sign/rotation checks pass. Neither screen establishes many-electron
state identity or assigns a resonance pole.

The three screens cost **0.605/0.402/0.402 seconds**, kernel peaks
**0.0874/0.0931/0.0683 GiB**. Public reconstruction verifies **94 payloads /
10 profiles** and byte-identical repackaging. A separately declared
[compressed outer-grid replay](../../projects/ukrmol_co/cas11-compressed-outer-grid-contract.json)
tests the same saved channels/amplitudes at 0.025/0.0125-requested-eV spacing;
energy-grid stability remains separate from the all-background rejection.

The [completed compressed outer-grid companion](../../projects/ukrmol_co/qualification-evidence/compressed-outer-grid/README.md)
then verifies all **20 native stages**, 317/633-point grids, **48 complete-data
fits / 240 held folds** and original retained binary inputs. Cost is **634.56
seconds / 0.7989 GiB**. Common original-grid phases agree within **3.2e-6 rad
modulo π**; matched finest-pair fit position/full-width differences are at most
**0.3651 / 2.0296 meV**, passing the grid gates. All-background position spans
remain **0.08886/0.08832 eV**, width spans **32.096/32.043%**, rejecting extraction.
Native MAXFIT appears on both grids, with empty finest-grid candidate lists;
native fit acceptance and pole identity stay unset. **257 payloads / 3170
profiles** publicly reconstruct, including the parent screens and full
0.05/0.025/0.0125-requested-eV fit refinement sequence.

### Compressed CAS12 singlet-A1 driver diagnosis

At R=1.9 bohr, both preserved CAS(10,12) DZ QC checkpoints fail their original
singlet-A1 CI convergence flag. A [complete fixed-orbital diagnosis](../../projects/ukrmol_co/qualification-evidence/compressed-cas12-diagnostic/README.md)
requests 5/8 roots at spaces 40/80/160 with a 600-cycle cap. All 96 probes
execute, but only the 84 probes outside singlet A1 pass. The 12 singlet-A1
probes retain false flags, physical residuals up to 0.19803 Hartree and
spin-square values near 2 for some requested singlets. Common energies differ
by 1.6561 Hartree between five/eight-root requests despite search-space
differences below 3.695e-13 Hartree. These failures explain why stable printed
energies and a small orbital gradient cannot qualify the original objective.

The [separate eight-probe driver pilot](../../projects/ukrmol_co/qualification-evidence/compressed-singlet-driver-pilot/README.md)
uses `direct_spin0_symm` only for singlet A1, preserving the same physical
Hamiltonian, orbitals, spin penalty and tolerances. This driver restricts the
CI matrix to alpha/beta-exchange symmetric vectors. That restriction admits
even-spin states as well as singlets, so total-spin certification still requires
the physical spin-square and spin-penalty-action gates. An independent
`direct_spin1.contract_2e` action is compared with the restricted-driver action
on every returned vector; its residual also uses an independently computed
physical Rayleigh energy. Both actions use the same effective one-/two-electron
Hamiltonian in Hartree, with the same frozen-core offset.

All eight fixed-orbital probes pass. Maximum physical/penalized residuals are
8.708e-10 / 9.988e-10 Hartree; independent full-action difference/residual are
6.752e-10 / 5.514e-10 Hartree. Common-root/space differences are below
2.843e-13 Hartree. Closed-form two-site Hubbard controls and perturbed-vector
residual detection pass. The two checkpoints nevertheless differ from the
rejected optimizer's saved singlet-A1 energies by up to 1.670142 Hartree, so
they must be **reoptimized under the repaired CI objective** before any import.

The experimental Python interface is `--target-singlet-a1-driver spin0` in the
state-averaged CAS runner; `build_target` receives the equivalent
`target_singlet_a1_driver` configuration field. The default remains `spin1`.
The selected driver is recorded per spin/irrep in QC diagnostics. All residual,
spin, CI-convergence and fresh-start gates remain mandatory, and all other
sectors keep the established driver. Repository validation additionally uses
a mixed-singlet/triplet CASSCF toy system with a nonzero optimized active-space
rotation, comparing both sector energies at the final orbitals with the full
CI oracle. These analytic and moving-orbital checks pass in the pinned image.
The [fresh QC contract](../../projects/ukrmol_co/cas12-compressed-singlet-driver-qc-contract.json)
declares two finite new optimizations at trial spaces 80/160, with 200 CI cycles
and 150 orbital cycles; separate independent coverage and native-import gates
are required afterward. This remains an experimental project-stage capability.

The [fresh optimizer companion](../../projects/ukrmol_co/qualification-evidence/compressed-singlet-qc/README.md)
then records a mixed result. Space80 rejects its singlet-A1 CI flag despite
orbital convergence. Space160 passes all optimizer/fresh-CI/spin/Pi/MO checks
in **1724.00 seconds / 0.9610 GiB**, with objective **−112.23910084461343 Hartree**,
gradient **5.932e-8** and fresh-root error **1.279e-13 Hartree**. Saved roots
agree across the two new checkpoints within **9.813e-11 Hartree**, but numerical
agreement cannot override the space80 rejection. The public companion reconstructs
**228 payloads / 16864 profiles** and preserves both exits and all sources.
The [independent coverage contract](../../projects/ukrmol_co/cas12-compressed-singlet-qc-coverage-contract.json)
declares 96 fixed-orbital probes at spaces80/160/240, with physical/spin-penalized
and independent singlet Hamiltonian-action gates. Newly passing QC remains
unqualified for import/scattering until that separate coverage and native
root/dipole/orbital-identity checks pass.

That [coverage subsequently passes all 96 probes](../../projects/ukrmol_co/qualification-evidence/compressed-singlet-coverage/README.md)
in **1420.60 seconds / 0.8013 GiB**. Physical/penalized residual maxima are
**8.889e-10 / 9.968e-10 Hartree**; independent singlet full-action difference/
residual are **7.081e-10 / 5.690e-10 Hartree**. Ensemble-root differences are
below **3.837e-13 Hartree**, space/root-count differences below **2.558e-13
Hartree**. Both new checkpoints agree in roots within 9.814e-11 Hartree and
core/active subspaces to approximately machine precision. Successful fixed-
orbital coverage does not repair the rejected space80 optimizer record.
Only space160 jointly passes original QC and independent coverage. The
[separate native-import contract](../../projects/ukrmol_co/cas12-compressed-singlet-qualified-import-contract.json)
therefore uses that checkpoint, with all 64 native-root/dipole/subspace gates and
existing queue/resource priority. Neither native import nor scattering is yet
qualified by this electronic coverage alone.

### Sparse/iterative scattering qualification contract

**Status: fixed CAS(10,10) 2048-root implementation qualifies; CAS(10,11) qualification in progress.**
Qualify this local route before allocating a large-memory host. The goal is the
same fixed-nuclei scattering observables and physical target model, with a measured
lower-memory implementation and an independently checked spectral approximation
if fewer inner-region eigenstates are retained.

Distinguish three independent properties: sparse Hamiltonian storage, iterative
eigensolution, and selected eigenpairs. The successful CAS(10,12) target calculation
uses SLEPc/Krylov–Schur with a dense PETSc Hamiltonian and five roots per sector;
its 40.405-GiB peak does not measure sparse scattering. The 222.4-GiB scattering
array floor measures the installed full-spectrum ScaLAPACK implementation.

Source inspection of the pinned engine identifies a concrete alternative:
`Contracted_Hamiltonian_module.build_contracted_hamiltonian` initializes the
matrix structure as `mat_sparse`; `SLEPCMatrix_module` has a PETSc
`MatCreateSBAIJ` branch as well as `MatCreateDense`. The source-level sparse
preallocation heuristic is not a measured nonzero fraction or memory budget.
`Options_module.compute_expansions` interprets `nstat=0` as the complete contracted
spectrum. Selecting an iterative backend without changing that requested spectrum
does not establish a low-memory selected-root calculation. The existing runner's
`--target-diagonalizer` option changes only the target template.

Before deriving a scattering recipe, resolve those symbols in the checksum-pinned
upstream source matching the built image. Audit dispatch, actual PETSc matrix type,
eigenpair selection/coverage, convergence and boundary-amplitude output, and the
SCATCI-to-outer-region interface. Record source and PETSc/SLEPc library digests
alongside generated namelists. Supported generic engine branches are leads, not
proof that every downstream stage accepts a truncated scattering spectrum.

The [first native sparse control](../../projects/ukrmol_co/qualification-evidence/sparse-native-control/README.md)
now verifies B1 CAS(10,10) roots 1–128 against its 8350-dimensional dense oracle:
energy error **5.97e-13 Hartree**, absolute residual **1.34e-10 Hartree**, and
continuum-coefficient error **7.88e-8**, with `mpisbaij` storage at a **1.228-GiB**
kernel container peak. Its legacy SWINTERF continuation fails on a reduced
partitioned CI file. The pinned source additionally reads the partitioned
diagonal into a local array that is never returned to its consumer; full
coefficient output alone therefore does not qualify that branch.

The native MPI boundary exporter instead emits an **uncorrected truncated pole
sum**, with no omitted-state correction. Its first B1 export matches dense
channels/thresholds/multipoles byte-exactly and boundary amplitudes within
**3.46e-8**, but its 128-pole phase difference reaches **1.56 rad**, rejecting
the observable gate. The [finite replay guide](../../projects/ukrmol_co/SPARSE_SCATTERING.md)
specifies 128/512/2048 iterative controls and independent dense-spectrum omission
diagnostics. Those checks separate eigensolver/export correctness from omitted
spectral-background error; neither result qualifies CAS(10,11)/(10,12) yet.

The [completed boundary companion](../../projects/ukrmol_co/qualification-evidence/sparse-boundary/README.md)
now qualifies the fixed CAS10 **2048-root** implementation in both sectors:
complete phases agree within **2.0e-7 rad**, and maximum identical fixed-window
position/full-width differences are **5.847e-8 / 2.381e-7 eV**. Energy/residual/
boundary errors are **9.664e-13 Hartree / 1.336e-10 Hartree / 1.003e-6**.
Dense-omission refinement through **4096/6144/8192/8350** confirms stability;
the complete endpoint reproduces the dense boundary bytes and phase grid.
Two-sector replay cost is **1340.74 seconds / 1.861 GiB**; prior electronic/
integral/CONGEN generation is excluded. Public fetch and independent reconstruction
of all **306 payloads / 25504 profiles** pass.

The [finite CAS11 contract](../../projects/ukrmol_co/sparse-scattering-cas11-contract.json)
then declares dense omissions **2048/4096/8192/16384/24576/27546** and native
counts **2048/4096/8192/16384** against the unchanged forty-state target. Native
controls run only if dense omission establishes an eligible count and stable
higher-count refinement. Four-rank containers are limited to **24 GiB**, with
**48-GiB available-RAM / 30-GiB disk** floors. Each control receives the same
phase/position/width gates and identical backgrounds. CAS12 still requires its
own qualified smaller-model prerequisites and measured local resource headroom.

The [CAS11 refinement companion](../../projects/ukrmol_co/qualification-evidence/cas11-sparse-refinement/README.md)
independently reconstructs **2048/4096/8192** native controls in both Pi sectors.
All eigenpair, static-interface/boundary, complete-grid phase and fixed-window
position/width gates pass. Maximum phase errors are **5.1e-6 / 1.0e-7 / 1.0e-7
rad modulo π**, with wall/peak charges **2.454 h / 7.693 GiB**, **3.335 h /
7.365 GiB**, and **7.060 h / 13.160 GiB**. The dense omitted-spectrum controls
remain stable through 27546 poles. At 8192 native roots, energy/residual maxima
are **1.478e-12 / 1.330e-10 Hartree**, and fixed-window position/full-width
differences are at most **1.242e-8 / 2.083e-7 eV**.

The original 16384-root B1 solve instead reaches its **24-GiB cgroup ceiling**
and exits **137**, with a recorded Docker OOM event. Its original failed controller
stops before B2 or fit continuation. The separately declared
[resource-only successor](../../projects/ukrmol_co/sparse-scattering-cas11-memory-repair-contract.json)
raises only the container cap to **40 GiB**, retains the same scientific contract,
and requires **80-GiB available-RAM / 40-GiB disk** headroom plus completion of
the active CAS12 DZ target import/verifier. It permits one fresh attempt and
independent reconstruction; any failure stops it. No CAS12 scattering is
automatically released. The public companion also preserves the original
interrupted local monitoring record and its QC-reporting exception.

That 40-GiB successor later hits its declared **21600-second B1 stage timeout**,
with total wall time **21601.48 seconds** and kernel memory peak **31.994 GiB**.
It has no completed eigenpair/export/observable verdict; a wall-time interruption
does not establish numerical nonconvergence. The [public DZ-import/failure companion](../../projects/ukrmol_co/qualification-evidence/dz-import-queue-failures/README.md)
independently reconstructs this original failure and the successfully completed
stretched CAS12 DZ64-root import. The latter passes all native root/dipole/
subspace gates in **21336.38 seconds / 40.371 GiB**, maximum root difference
**5.460e-8 Hartree** and dipole difference **8.919e-8 a.u.**

A separately declared [time-budget successor](../../projects/ukrmol_co/sparse-scattering-cas11-time-budget-successor-contract.json)
uses a **24-hour per-stage limit**, the same 40-GiB cap, immutable scientific
inputs and independent eigenpair/boundary/phase/fixed-window gates. The original
24-GiB OOM and six-hour timeout remain separate failures. The later TZ import
owner fails before any native engine launch because its package snapshot lacks
`projects/__init__.py`; both dependent imports stop at their guards. A fresh
[package-complete target-import queue](../../projects/ukrmol_co/cas12-native-import-packaging-successor-contract.json)
retains unchanged scientific arguments/proofs/checkpoints and TZ → supported-space
aug-TZ → compressed priority, with fresh run names. Each native verifier remains
mandatory; the original aug-TZ coverage rejection is still explicitly recorded.
These execution repairs do not release a CAS12 scattering pilot or model sweep.

The user-directed [eight/twelve-core policy](../../projects/ukrmol_co/sadaharu-cpu-allocation-policy.json)
subsequently redeclares only the zero-step pending target queue as a fresh
[twelve-rank successor](../../projects/ukrmol_co/cas12-native-import-mpi12-successor-contract.json).
Scientific arguments/checkpoints/proofs and all import gates are unchanged;
the original waiter/source remain exact administrative supersedure evidence.
Native imports use twelve physical cores/ranks, with four cores reserved for
small controls and the existing memory/disk floors. Larger-model rank scaling
is unmeasured until a matched comparison exists; more cores do not qualify
the model or relax a rejected electronic proof.

For a spectral R-matrix construction, each retained inner-region eigenstate supplies
a pole and boundary amplitudes; omitted states can alter the scattering background.
A few accurately converged eigenpairs alone therefore do not qualify scattering.
Any omitted-spectrum treatment must be explicit and validated, or a retained-state
refinement must demonstrate that omission errors are below the observable gates.

#### Differential and convergence gates

1. Start with a qualified CAS(10,10) dense control. Freeze target checkpoints,
   retained neutral states, CI/model inputs, continuum, radius, deletion threshold,
   energy grid and extraction settings. Verify native solver residuals and
   symmetry, compare boundary amplitudes in a sign/degeneracy-aware manner, and
   confirm complete downstream execution without silent dense fallback.
2. Use the completed CAS(10,11) equilibrium l=4 dense calculation as the larger
   reference, dimension 27546. After the source/interface audit, predeclare a finite
   eigenpair-count and spectral-coverage sequence, fixed eigensolver tolerances and
   resource limits. Increasing root count changes the scattering spectral expansion;
   it must not change the forty-state neutral target or its orbital ensemble.
3. Compare both B1/B2 phase grids at common native energies over the declared
   energy domain. Require maximum phase differences **≤0.05 rad modulo π**,
   position differences **≤0.05 eV**, and full-width differences within the existing
   **5% gate with 0.001-eV absolute floor**. Use identical complete-data fit windows
   and backgrounds for each pair. Require agreement with the dense reference and
   stability on further spectral refinement; do not rely on literature closeness
   or adjust the background to compensate for an omitted-spectrum error.
4. Preserve residuals, requested/converged eigenpair counts, spectral coverage,
   matrix dimensions/nonzeros, boundary data and raw phases for every refinement.
   Empty native candidates stay unset; `MAXFIT=100` truncation cannot supply a
   passing fit verdict. Existing background-width sensitivity remains a separate
   extraction gate even if the low-memory solver reproduces the dense result.
5. Measure full-pipeline wall/CPU time, aggregate memory, eigenvector/workspace
   allocation and disk/scratch peaks, with exact CPU ownership and RAM/disk floors.
   Only a passing numerical comparison and a resource envelope with headroom
   release a same-model CAS(10,12) pilot on Sadaharu. Recheck spectral convergence
   at CAS(10,12); a smaller-model root count does not establish its required count.

#### Evidence for reproducing the route

Publish checksum-pinned build/configuration recipes, dense-reference pointers,
immutable target lineage, generated solver inputs, actual backend/storage logs,
eigenpair/refinement tables, boundary data, raw phases, independent fixed-window
fits and complete resource profiles. Include an executable verifier and preserve
original failed attempts. State the validated model/energy domain, numerical gates,
resource limits and any unresolved omitted-spectrum correction. A demonstrated
solver memory reduction does not establish active-space/basis/ensemble convergence
or qualify a full geometry sweep. If the local route fails a gate, record that
specific limitation before reconsidering the deferred large-host fallback.
