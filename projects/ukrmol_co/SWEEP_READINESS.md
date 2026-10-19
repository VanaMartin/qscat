# CO: gates before the production sweep

The goal is reproducible fixed-nuclei scattering and an independently correlated
neutral curve, followed by a physically grounded QSCAT two-dimensional potential
with held-out scattering validation. **A production sweep is not yet released.**
The anchor campaign is R=1.9/2.1323/2.5 bohr. A wider geometry range needs its own
electronic,continuity and neutral-reference qualification.

This roadmap is decision-driven: use the smallest retained eigenpair count that
passes observable and omission checks; every refinement must address a measured
discrepancy. The full execution record is in [CONTINUATION.md](CONTINUATION.md).
Completed/failed owners,raw outputs and dated overview releases remain evidence.

The [quarter evidence](qualification-evidence/quarter-state-tracking/README.md)
now independently reconstructs all672 states and112 matrices on the host and
from public downloads,with byte-identical repackaging.
Compressed proof replay passes all438 payloads/96 probes,and the fresh native
import passes. The [compressed matched A2 control](A2_COMPRESSED_CONTROL.md)
now passes independent/public reconstruction and completes the three-anchor
channel comparison. A [supplementary same-engine active-space reading](qualification-evidence/same-engine-active-space/README.md)
preserves the CAS11→CAS12 rejection. The
[prospective extraction study](PHASE_EXTRACTION.md) now passes its finite
30-case calibration but rejects compressed/equilibrium saved nine-A2 phases;
stretched provisionally passes. A three-anchor width curve remains unqualified.
The subsequent [analytic identifiability audit](PHASE_IDENTIFIABILITY.md) finds
94 wrong-width counterexamples despite 45/45 true-pole smooth controls passing;
local phase fidelity/unitarity alone cannot justify a stronger width claim.
The [saved-matrix audit](SAVED_MATRICES.md) now reconstructs all six nine-A2
anchor/sector K-to-S cases, but the full twelve-case gate stops with six unset
five-A2 interfaces: target-table rows are dense while channel labels retain
original property IDs. Existing nine-A2 data support the next physical benchmark
without a new inner solve; independent propagation/poles remain unqualified.
The [first physical outgoing-pole benchmark](PHYSICAL_POLE_BENCHMARK.md) now
passes all 36 real-axis profiles and 81 solves against independent exact Jost/
Lambert-W oracles, with three-rung convergence and no reseeding. This qualifies
the scalar-shell boundary convention only; coupled polar tails and actual
native amplitude normalization/threshold sheets remain open.
The [separate matrix tail/threshold controls](COUPLED_TAIL_BENCHMARK.md) now pass
513 S samples, 27 pole solves and twenty frozen checks. Tail truncation remains
visible despite unitary matrices; wrong closed sheets/flux factors are detected.
The [combined exterior boundary-map study](COMBINED_EXTERIOR_BENCHMARK.md)
completes 436 paired-method cases but **rejects spectral refinement**:
7.90991e-7 exceeds 1e-7 at ε=30,L=80,θ=20° for all three dipoles. Finest
agreement 2.654e-9 does not override it; original exits 1/1 remain. No successor
grid is automatic. A [native source trace](../../docs/physics/co-native-boundary-conventions.md)
resolves the already absorbed sqrt(2), total-Hartree length-valued zero-Bloch
sum and already-Rydberg multipoles. Full native K normalization remains open.

The separately frozen [bounded resolution control](COMBINED_EXTERIOR_RESOLUTION.md)
passes all sixteen checks with degrees128/256/512 and unchanged1e-7 budgets.
The final spectral change is4.28129e-8; all436 cases complete, exits0/0.
The parent rejection remains. Increased condition estimate4.692e8 and finest
domain/method differences5.207e-8/4.350e-8 limit the interpretation. Next is a
prospectively specified independent saved-native open-open K reconstruction,
including radius/step/tail, closed-channel and Gailitis-normalization controls.
No missing export justifies a new inner solve; native K/poles remain unset.

The [native sign/flux source trace and independent Numerov method](qualification-evidence/native-exterior-method/README.md)
now pass17 focused tests and fresh public replay. A separately costed
[21-case equilibrium-B1 pilot](NATIVE_EXTERIOR_PILOT.md) completes at three
energies with independent step/endpoint/tail-order controls, all298 channels and
all2048 saved poles. **It rejects S/phase reproduction at5eV:**1.49310e-6 and
2.72174e-6rad moduloπ exceed1e-6 despite K,unitarity,reciprocity and all
refinement checks passing. Original exits1/1 remain. Measured cost735.161s /
0.08564GiB sampled,below the975.5s estimate/2400s cap. Independent public replay
reproduces all21 cases,all ten decisions,false and exit1. A separately justified
[native-subrange attribution design](../../docs/physics/co-native-exterior-attribution.md)
is frozen at **`666c056`**: native8→16→32 subranges at unchanged100bohr plus two
finer independent5-eV cases. Estimate1667s,cap3600s,4GiB; no new inner solve.
This control now **passes all ten checks**,original exits0/0,and accounts for
99.72% of the measured5-eV S mismatch. Refined native32-to-independent K/S/phase
differences are2.68388e-9/4.08632e-9/1.15175e-8rad; native-only16→32 changes
over all99 energies are≤5.588e-9/3.139e-10/3.477e-10rad. Cost1339.551s /
0.56024GiB sampled. Fresh public replay verifies ten payloads,297 native records,
both independent cases and all ten checks,exit0. Original pilot
false/exits1/1 remain; full594-matrix independent reconstruction,CO poles and
production still need a prospective all-anchor/sector/energy control design.

That [six-case reconstruction design](NATIVE_FULL_RECONSTRUCTION.md) is now
specified: native baseline/16/32 controls over all99 energies, then18 independent
step/radius/order probes, then gated paired-step completion of594 matrices.
The stretched5-eV point lies only0.02524 eV below the first excited threshold.
The frozen scope forecasts4.5 h on eight physical cores, with16 GiB and a10 h cap.
All K/S/phase gates remain separate at1e-6; false/unset stops new submissions.
The study, frozen at **`fd7c871`**, now stops with **false / exits 1/1**.
All six native controls and all 18 independent probes pass, but the remainder
rejects absolute K changes at stretched 1 eV in both sectors. Coverage stops at
132/594 energies and 354 independent solves (112 remainder passes / 2 rejects).
S/phase and invariants pass at the rejected points; they cannot waive the K veto.
The [saved-matrix diagnosis](NATIVE_FULL_RECONSTRUCTION.md#saved-matrix-diagnosis)
resolves large-tangent coordinate sensitivity through the exact Cayley identity,
without identifying the remaining propagation/reference error or releasing a
complex-pole interpretation. Original cost is 1.91033 h / 3.14685 GiB; no new
inner solve is justified. Further numerical work needs a separately frozen,
costed attribution experiment rather than automatic full-grid continuation.
Fresh public replay verifies all saved decisions and reproduces the four
rejected-point solves exactly. All 31 focused tests pass. The original false
verdict/exits 1/1 and the 462 uncomputed energy identities remain.

The subsequent [two-identity attribution](NATIVE_K_ATTRIBUTION.md) is initially
**unset**,native32 timeouts at180 s, native exits−9/−9, runner/container exits1/1.
No K file or independent solve was produced. Fresh public-byte replay preserves
the unset decision and passes12 payloads/19 source files/29 tests. The native
cost forecast was wrong: energy-independent setup diagonalizes a2980-dimensional
Hamiltonian in every subrange before the energy loop; saved native32 calls cost
866–920 s. A separately frozen [resource successor](../../docs/physics/co-native-k-attribution-resources.md),
commit **`4a558b3`**, completes the same finite scientific contract: expected3400 s,
upper5400 s,batch cap9000 s,native-stage cap3600 s,four physical cores,4 GiB.
Measured1669.647 s /1.973831 GiB; all four native calls and ten independent
invariant checks pass, but the study is **false / exits1/1**: both sectors reject
step, endpoint and finest/native K. The smaller step worsens the K discrepancy
to about7.256e-5, while tail order6→8 changes only about4.4e-11. Native32→64
changes pass but are not absolute error bounds. Fresh public saved-decision
replay preserves all thirty decisions and passes29 tests. All ten numerical
controls reexecute from public bytes with exactly reproduced K/stored S
components and unchanged false decisions;15 final payloads verify. No scientific
container now runs. The original timeout
and six-case rejection remain. Further work needs a fixed-grid precision/algebra
or independent-propagation design, not another automatic spacing rung.

The subsequent [fixed-grid increment control](NUMEROV_INCREMENT_CONTROL.md)
passes16 analytic checks and24/24 two-point CO checks, exits0/0. It retains
the exact-arithmetic problem while changing cancellation-sensitive algebra;
finest/native K discrepancies reduce to1.49e-7/1.09e-7. All49 focused tests
pass. Fresh public numerical replay reproduces all ten K/stored S outputs
exactly with identical24 passing checks;20 final payloads verify. No scientific
container runs. This removes the observed two-point vetoes under the declared control,
not the original reconstruction rejection. At step0.0125/endpoint800 the
native K diagnostic remains1.59e-6/1.46e-6, so not every row is qualified.
At that stage, cross-anchor/energy and endpoint/precision validation needed a separate contract
before any full-grid release;462 original energies remain uncomputed.

The subsequent [cross-anchor/all-finest-endpoint study](NUMEROV_CORNER_CONTROL.md)
passes both stages,14 identities,294/294 case checks and36/36 reference checks,
exits0/0. All finest-step endpoints200/400/800 pass the unchanged1e-6 native gates;
the previously problematic800-bohr rows give K2.15388e-7/5.59271e-8 at h=.00625.
Measured78 new solves/six reused rows,31.3469 min /0.406113 GiB. Independent
saved-decision replay and57 public-byte tests pass. The predeclared14-row
numerical subset reproduces every K/stored S component exactly, exits0/0, and
all15 final payloads verify. No scientific container runs. These identities were
already visited: original132/594 coverage and462 uncomputed identities remain; untested
interior/large-K energies need a separate finite contract before full-grid release.

The subsequent [interior stress selection and analytic prerequisite](NUMEROV_INTERIOR_CONTROL.md)
examines all594 native matrices and freezes the largest interior K plus fixed
midpoints. Its24-solve high-momentum/large-tangent analytic prerequisite **rejects
four of68 checks / original exits1/1**: scalar discrete-oracle drift1.21e-7>1e-7
and rotated open/closed continuous K1.11–1.13e-6>1e-6. S/phase/invariant and
high-precision-oracle checks pass but cannot waive those gates. Cost72.95 s /
0.04408 GiB; independent saved decisions and62 fresh public tests pass. The
conditional60-solve CO scope is withheld: zero CO/native/inner solves, candidate
qualification remains14 and historical coverage132/594/462 uncomputed remains.
Next needs a separately frozen fixed-grid precision/algebra or independent
propagation attribution in this analytic regime, not a spacing ladder or relaxed
tolerance. The earlier corner result remains qualified in its finite scope.
Full public numerical replay reproduces all24 stored K matrices exactly, all68
decisions and scientific false/exits1/1; verification exits0/0. All16 final public
payloads verify. Replay costs72.22 s /0.04161 GiB. No scientific container runs.

The subsequent [fixed-grid precision attribution](NUMEROV_PRECISION_ATTRIBUTION.md)
passes46/46 diagnostic checks/exits0/0 with ten analytic solves. A matrix discrete
oracle separates roundoff from truncation: drift dominates by57–120× at rejected
200/400-bohr points. Equivalent open-frame controls reduce drift, but miss the
predeclared all-endpoint tenfold criterion at800 bohr. This is an attribution,
not a repair: interior false/exits1/1 and coverage132/594/462 uncomputed remain.
Next is a separately costed arithmetic-precision control on those same grids.
Public replay reproduces all ten new K matrices exactly and all diagnostic/
hypothesis values;68 tests and19 final payloads verify. Verification-only
dependency/preflight failures remain archived. No scientific container runs.

## 1. Finish trustworthy target inputs

**Current:** CAS11 DZ imports/scattering and CAS12 DZ imports pass at three
anchors. Supported stretched aug-TZ passes all64
root/dipole/subspace/seed-preservation checks. The original stretched TZ import
failed seed preservation; its reoptimized orbital frame is a separately
qualified numerical target. Compressed CAS12 space160 ordinary QC and96-probe
coverage pass;its private-cache native import completes with original gates
preserved. Both stretched and equilibrium nine-A2 fixed-orbital48-state imports
pass independent full-vector/property reconstruction. The equilibrium gate
completes51.817minutes/5.425GiB;public reading verifies44raw files/48vectors/
930multipoles. Original forty-component orbital averaging remains unchanged.

Compressed48-state inputs also pass independent/public44-file reconstruction,
with48fullvectors/915multipoles at51.178minutes/5.477GiB.

**Next:** publish/fetch remaining native-import evidence and qualify any new
40→50 ensemble/basis target before matched scattering. CAS11 16K/TZ public reconstruction already
passes. Preserve every original failed exit and the original TZ seed rejection.

**Pass:** ordinary QC and fresh-CI flags; roots≤1e-7 Hartree; physical/penalized/
spin-penalty actions≤1e-9 Hartree; dipoles≤1e-5 a.u.; spin-square errors≤1e-6;
MO orthogonality≤1e-9; Pi degeneracy and predeclared seed/subspace checks.
Successful import establishes engine/orbital consistency,not model convergence.

## 2. Establish an affordable CAS12 scattering calculation

**Current:** CAS11 selected-root scattering passes both-sector independent
eigenpair,boundary and observable checks;2048 roots is the smallest qualified
count. A further16K confirmation cost19.65 hours without changing that decision.
CAS12 dense workspace exceeds222 GiB and does not fit Sadaharu.

**Completed:** equilibrium CAS12 audited2048/4096 pair stability passes in both
sectors,with maximum phase difference1.67e-5rad and passing exports/residuals/Pi/
fixed-window fits. Public raw reading now reconstructs both resolutions. The
user accepts these controls as sufficient for now;use economical2048 for
matched experiments. The administratively stopped8192owner and partial outputs
remain preserved. No further resolution follows automatically.

**Pass for current controls:** physically relevant phases and fitted position/
width meet the declared spectral-sensitivity budgets. Use **0.05 eV
position,5% full width with0.001 eV floor,and0.05 rad phase moduloπ**. Independent
eigenpair/boundary checks must also pass. Successive differences are not absolute
omitted-spectrum bounds;that uncertainty remains explicit. Record solve/export/
outer cost separately. Escalate only for a named decision-changing discrepancy.

## 3. Select a converged electronic and channel model at the anchors

**Current:** both CAS10→CAS11 and CAS11→CAS12 equilibrium changes exceed
position/phase budgets. The saved CAS11/2048→CAS12/2048 change is0.1212798rad
and63.652–73.447meV in fitted position;the rejection persists withCAS12/4096.
Stretched five→nine-A2 rejects0.05rad at0.0679062rad;its declared nine-A2
2048→4096 refinement is1.96e-5rad and the refined channel rejection persists.
Basis,active-space and full retained-channel convergence remain open.
Fifty-component equilibrium orbital coverage passes,but native fifty-component
imports and matched scattering have not been performed.

**Completed:** equilibrium five-versus-nine-A2at2048passes all budgets at
0.0163665rad,centers9.427–10.924meV lower,widths2.672–3.645meV smaller.
Public reconstruction passes73raw files/58reader sources. Five/nine costs
1.656/1.817h and6.408/6.896GiB;dense baseline phase5.1e-6rad passes.

**Completed:** compressed five-versus-nine-A2 at 2048 passes at **0.0115944 rad**,
centers 13.241–18.210 meV lower and widths 2.218–4.404 meV smaller. Independent/
public reconstruction verifies 73 raw files / 58 sources; five/nine costs are
1.651/1.827 h and 6.011/7.246 GiB. Original native/resource/owner exits are zero.
The fresh dense-baseline phase change is 2.2e-6 rad. The matched 2–5 requested-eV
window retains orders 1–4 and does not replace the earlier 24-case extraction
rejection. Five-A2 is not uniformly insensitive over the domain. Nine-A2 is
the larger tested channel representation for subsequent matched controls;
full channel convergence remains open.

**Same-engine supplementary control:** the fresh equilibrium CAS11 baseline
exactly reproduces the older selected phases and preserves the CAS11→CAS12
phase/position rejection on the same audited engine and eight ranks. The
strict byte-deck preflight fails on two spacing/comment lines; that original
failure and unset verdict remain preserved. A separately declared exact-line
audit confirms physical equality without changing any numerical budget. Public
replay reconstructs 197 raw files and six owner-provenance records. Target
diagonalizers/imports and optimized orbital frames still differ.

Use the measured active-
space/channel effects to select remaining matched three-anchor and
DZ→TZ→aug-TZ controls with all other settings fixed.
Independently test explicit channel count and forty→fifty orbital averaging:
these change different parts of the model and need separate comparisons.
Check all three anchors so agreement at equilibrium cannot hide a stretched
failure. Any needed larger orbital/channel model requires its own QC/import
proof before scattering. Diagnose distinct orbital optima using more than the
mean state-average objective.

**Pass:** position,width and phase shifts fall within the above budgets,with
qualified numerical/extraction errors smaller than the model comparison.
Document remaining basis/correlation uncertainty. Literature values provide
context; agreement with a paper cannot substitute for successive-model
convergence or justify choosing a convenient orbital solution.

## 4. Qualify continuous target thresholds and resonance extraction

**Current:** public twelve-root anchor and two-start midpoint reconstruction
passes. Quarter controls greatly strengthen the named singlet-A2 root5 route,
but singlet-A2 root4 and triplet-A2 root5 still change endpoint labels under
path refinement. Separate quarter wavefunction replay passes on the host.
Earlier native background/window changes move widths by roughly12% at
equilibrium,32% compressed and6.7% stretched. The new multi-window/three-start/
five-fold control rejects extraction with width spans17.235%/16.635% for
equilibrium CAS11/CAS12 and6.999%/6.991% for stretchedfive/nine A2. These are
different declared control sets;all exceed5%. Constant background dominates
the measured sensitivity,including failure of equilibrium held-phase prediction.

**17 October prospective estimator:** model selection and unused-data rules
were frozen before calibration/native use. All 30 analytic cases meet their
declared calibration requirements (17 accepted; all smooth-exponential cases
reject). Saved nine-A2 stretched passes in both sectors with width span
**0.380 meV** and untouched-audit maximum **0.000298345 rad**. Compressed and
equilibrium reject required-window prediction at the fixed **0.001-rad**
representation threshold; best missing-window maxima are 0.00225936/0.00130390
and 0.00108757 rad, respectively. Original native-reader exit **1** and every
earlier extraction/channel rejection remain preserved. The calibrated estimator
does not qualify a native pole or arbitrary smooth backgrounds. See
[design, results and stopping decision](PHASE_EXTRACTION.md).

**17 October endpoint-isolated identifiability audit:** the fixed 405-cell matrix finds **94**
local analytic/unitary representations with **±10% known width errors** that
pass every unused phase audit within 0.001 rad. Counts are compressed-like
27, equilibrium-like 13, narrow-control 54, resolved/stretched-like 0. All 45
mandatory true-pole smooth-nonpolynomial controls pass. These background families
are not qualified finite-range causal CO models; zero counts do not establish
uniqueness. The result rules out treating added background flexibility and
better local prediction alone as width qualification. No native data was read.
See [physical limits and complete evidence](PHASE_IDENTIFIABILITY.md).
The original producer's 98 matches and exit 0 remain preserved, alongside its
failed endpoint-unused declaration (14/45 overlapping grids). Corrected source
was frozen before rerunning; all 45 corrected grid sets are independently disjoint.

**Next:** audit the two
remaining retained components: inspect competing row scores,energies and
complete nearly degenerate subspace projectors. Determine whether a larger
retained channel manifold is required. Test only geometries/model controls
that address those specific failures. Track relevant resonance features across
neighboring geometries with phase/K-matrix evidence and consistent thresholds;
electronic-state assignment alone does not establish resonance-pole identity.
The completed prospective estimator stops at its three-anchor rejection.
The saved-data inventory now finds boundary/channel/multipole data and native
open-open K samples already present. Six nine-A2 cases reproduce native phases
to 4.94e-8 rad; six five-A2 strict interface readings stay unset on target IDs,
preserving aggregate false and reader exit 1. Trace the five-A2 ID consumers and
freeze any mapped decoder separately. The scalar physical-shell gate now passes:
maximum finest S/q errors are 2.584e-10/6.495e-12. Separate exact dipole-tail and
unequal-threshold shell benchmarks also pass (finest matrix/q errors ≤6.895e-11/
3.029e-12), including sheet and flux controls. Next combine noncommuting thresholds
and polar tails under independent real/complex checks: the first combined
boundary-map study now rejects the fixed spectral refinement budget despite
passing finest-method/domain controls. Freeze any successor resolution family
against that measured discrepancy before proceeding. The native zero-Bloch
surface/energy/CF source conventions are derived; full native K and Gailitis
normalization remain unqualified before an independent outer
propagator uses saved nine-A2 data. Require full K-matrix reproduction and
radius/step/tail controls before CO complex poles. No missing export currently
justifies a new inner solve. A new estimator needs independently
justified background/threshold/causal constraints or direct pole evidence,
mandatory nonpolynomial known-truth controls, and frozen selection/unused-data
rules. Preserve every earlier all-order/window rejection. Small held-phase
prediction error alone does not establish identifiable pole parameters.
Revisit energy resolution and threshold effects on the selected scattering model;
no relaxed budget or automatic fitting-order ladder follows this rejection.

**Pass:** relevant channel/subspace coverage and thresholds remain stable under
declared path controls; no unaccounted retained-to-omitted exchange affects the
observable. Position/width/phase meet the extraction budgets under independent
controls. Empty fits and `MAXFIT=100` leave verdicts unset. Physical wavefunction
overlap and AO-following coefficient transport stay separately labelled.
If Breit–Wigner fits cannot meet the budget,declare and validate an appropriate
alternative diagnostic before building a pole curve; weak fits do not define
zero widths or bound anions.

The current finer-route audit finds that both remaining retained components
have physical squared overlaps above0.5 and unambiguous local row margins,but
temporarily enter root6,which lies outside each five-root sector cutoff. Singlet
A2 root4 follows4→4→5→6→6; triplet-A2 root5 follows5→6→6→5→5. This makes retained
channel/subspace coverage a specific next check. Increasing the A1 roots in a
fifty-component orbital ensemble does not by itself enlarge A2 channel coverage.
Strong local links do not prove path convergence.

## 5. Qualify the neutral reference over the chosen range

**Current:** the QZ/5Z frozen-core RHF/CCSD(T) near-equilibrium dataset has52
points over1.9–2.5 bohr. QZ→5Z relative energies differ by at most17.24 meV,
below the20-meV pilot budget. Independent held-point interpolation is below
1 meV for both spline and PCHIP. Stretched R=3/4 bohr fails external RHF
stability before correlation; its replacement remains undesigned.

**Next:** declare the production geometry domain and uncertainty requirements.
For the existing range,assess correlation-model uncertainty beyond basis
convergence. For a wider range,design and independently validate a stable
correlated treatment and its connection to the near-equilibrium curve before
extending production. Keep the energy zero at each basis's own equilibrium;
dipoles are all-electron CCSD lambda-density values.

**Pass:** every used geometry has stable reference/converged correlation and
density checks,plus basis/correlation and held-out interpolation errors within
predeclared production budgets. Near-equilibrium pilot success does not release
a dissociation-range neutral potential.

## 6. Freeze the sweep and prove restartable pilot coverage

**Depends on:** a selected model passing steps1–5 over its declared domain.
Some anchor tests and neutral work can proceed alongside each other.

**Next:** freeze basis,active/core spaces,orbital ensemble,channel manifold,
root-selection/count,continuum/angular/radial/outer settings,energy zero and
extraction rule. Define geometry and energy grids that resolve the resonance
and threshold structure; reserve independent geometries and energies as
held-out data. Run a small interior/boundary pilot with checksums,phase/spin/
Pi gates,resource forecasts,finite stops and checkpoint/restart verification.

**Pass:** the pilot and independent controls meet all declared tolerances;
geometry interpolation and resource extrapolation are supported by measured
points. Freeze immutable sources/configs and an explicit sweep release record.
Use8 physical cores for large jobs,12 for very large jobs,and reserve0–3 for
small checks; bound concurrency by measured memory/disk. Paid provisioning
remains deferred. This record is the trigger to launch the production sweep.

## 7. Sweep,fit and validate the QSCAT potential

Run the frozen sweep incrementally,with per-geometry gates and additive
publication. Fit the neutral and electron–CO interaction to the qualified
training data with physical asymptotic/threshold constraints. Independently
validate QSCAT discretization and scattering at the reserved geometries and
energies against UKRmol+ data not used in fitting. A potential matching only
the training pole positions/widths is insufficient. Failed held-out comparisons
return to the measured model/extraction/fit discrepancy before promotion.

## Immediate queue

1. Design and independently validate an extraction rule addressing the measured
   background failure using known-pole phases and unused-data predictions;
   preserve every rejected all-order gate.
2. Use the completed three-anchor channel effects and same-engine active-space
   rejection to declare independent ensemble/basis controls, including target
   QC/import gates and explicit measured-cost forecasts before native launch.
3. Complete the remaining matched model/basis/ensemble controls,continuity and
   neutral correlation/domain qualification;then freeze and pass a restartable
   production pilot. Full-sweep readiness follows the gates,not a fixed date.
