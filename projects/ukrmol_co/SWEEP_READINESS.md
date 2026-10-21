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

The separately frozen [extended-arithmetic control](NUMEROV_ARITHMETIC_CONTROL.md)
then passes35/35 control and10/10 qualification checks, original exits0/0:
14 analytic solves in64.80 s. On unchanged grids,64-bit-significand longdouble
reduces drift453–2458× versus the same typed binary64 implementation. All six
extended drift errors pass1e-8 and finest continuous K/S/phase pass1e-6. Both old
SciPy baselines reproduce exactly. The small free-exterior solver is qualified
only for constant diagonal potentials, not the coupled CO exterior. Next design
and separately cost a dtype-preserving coupled inverse-power formulation with
independent controls, before revisiting the withheld ten-identity interior scope.
All earlier rejections and coverage132/594/462 uncomputed remain.
Public replay reproduces all14 campaign matrices exactly, including extended
decimal components, plus every decision/hypothesis. Six short-grid preflight
controls pass remotely in both dtypes;74 local tests pass with one wider-format
skip. All14 final public payloads verify, exits0/0; no scientific container runs.

The subsequent [coupled precision and dense-cost control](NUMEROV_COUPLED_PRECISION.md)
passes89/89 checks/exits0/0,32 typed propagations in45.78 s. Noncommuting
short-grid discrete arithmetic and continuous references, exact inverse-square
physical convergence, archived stress endpoints and synthetic298-channel solves
pass. Complete public replay reproduces every decimal K component and89 decisions;
80 tests pass with two host-dependent wider-format skips, exercised remotely.
Dense extended arithmetic costs79× binary64. Profiling the two existing dense
rows attributes61.92% to Taylor setup and26.11% to mixed solves. Next separately
measure setup-versus-interval scaling and target the proven hot path before
scheduling interior CO. The eight-interval average, including setup, is not a
CO throughput measurement. All14 public payloads verify. Original interior
false/exits1/1 and coverage132/594/462 uncomputed remain; no CO successor runs.

The subsequent [dense setup/interval benchmark](NUMEROV_DENSE_SCALING.md) passes
166/166 checks/exits0/0 in207.11 s:24 unchanged-kernel298-channel propagations,
two potentials, N=8/32/128 and two repetitions. Extended marginal rates are
.135263/.271376 s per interval for free/varying W. Recurring solves dominate
at128 intervals; all timing hypotheses pass. Synthetic60-solve estimates are
23–46 single-worker days, not a native CO forecast. Next isolate mixed-solve
residual products versus LU, then optimize/requalify the proven primitive under
a separate contract before an interior CO pilot. The existing Python numerical
reference and backward-error gates carry forward; all earlier rejections remain.
Full public replay reproduces all24 decimal/binary64 matrices,166 decisions and
six oracle pairs exactly. All14 final payloads verify;85 local tests pass with
two wider-format skips, exercised remotely. No scientific container runs at the
13:56UTC capture.

The subsequent [mixed-solve source-line profile](NUMEROV_SOLVE_PROFILE.md) passes
76/76 checks/exits0/0 in49.63 s. Eight existing-row reexecutions reproduce K and
diagnostics exactly. Extended residual evaluation accounts for94.6% of recurring
solve time,27.5× factor+initial+correction cost; timing hypotheses all pass.
Six actual A/B/X captures support the next precision-preserving residual kernel
experiment without another capture propagation. Public replay reproduces all
eight K/diagnostic records,76 decisions and six numeric operand triples exactly;
recomputed backward errors pass. All14 payloads verify,89 local tests pass/two
wider-format skips exercised remotely. Next resolve compiled extended arithmetic
and NumPy format handling, then benchmark/differentially qualify before integration.
Original interior rejection and all production gates remain; no scientific
container runs at the14:11UTC capture.

## 1. Finish trustworthy target inputs

The [native residual primitive](NUMEROV_NATIVE_RESIDUAL.md) now passes69/69 gates
with2.46–2.49× end-to-end speedup on the three extended captures. Full matrices
match NumPy; independent paired references/cancellation/layout controls pass.
Complete public replay reproduces all13 numerical rows/23 output hashes/69
decisions exactly; all17 public payloads verify. Original artifact-layout
preflight exit1/unset is preserved alongside its passing successor. Next freeze
integrated mixed-solve/coupled stress and interval-performance requalification
before drawing CO cost conclusions. No propagation or new CO identity is added.

That [native integration qualification](NUMEROV_NATIVE_INTEGRATION.md) now passes
**538/538 checks / exits 0/0**, 50 propagations / 1,737,200 mixed solves. All
1,735,263 audited native residual evaluations match NumPy; zero gate disagreements
and at most one correction. Four-channel embeddings force the mixed path for
short/physical and all three finest stress endpoints. Separate dense timing
rows yield **2.26–2.27×** marginal interval gain and **1.84–1.88×** full N128 gain.
Complete public replay reproduces all numerical rows/audits/decisions exactly;
all 14 public payloads verify. Producer/replay cost 415.73/416.17 s, about .138 GiB,
on one worker within reserved physical cores 0–3 and 4 GiB. Conditional synthetic
recurring-only estimates fall to 10.09–20.34 single-worker days, within that
synthetic scope only. The **next** experiment is a separately costed actual
saved-CO-matrix conditioning/throughput canary, before an interior pilot or any
full 60-solve proposal. Dense exterior matrices are small; memory rules excluding
dense inner-region Hamiltonians remain. Original interior false/exits1/1 and
historical coverage 132/594 remain unchanged.

The subsequent [actual-CO prefix canary](NUMEROV_CO_CANARY.md) passes **303/303
checks / exits0/0**,36 prefixes /3,996 mixed solves across six cases. All180
audited residuals match NumPy; zero gate disagreements, one correction maximum.
Early-region rates are .12392–.12499 s per native interval,2.23–2.24× reference
gain, only3.42–4.32% above the synthetic varying rate. Minimum binary64 rcond
.00426065 is a preconditioner observation. Runtime468.93s / .179GiB; all17 remote
tests pass. The original public strict-estimate replay fails exits1/1 on rcond
last bits while science fields/operands agree. A separate explicit supplement
keeps science records exact, checking rcond estimates within32 binary64 epsilons
relative; six fresh prefixes /90 solves and303 saved decisions pass exits0/0.
All19 final public payloads verify, including the failed replay and diagnostic.
These prefixes end before18.2bohr, without tail matching or K. They justify the
separately costed long-range radial pilot below.

The [bounded radial pilot](NUMEROV_CO_RADIAL_PILOT.md) now passes **131/131
checks /five full propagations /282,400 mixed solves**, original exits0/0/owner0,
using one already corner-qualified stretched B1 identity. Finest200/400-bohr K
discrepancies from native64 are1.56617e-7 /1.81307e-7, below separate1e-6 gates;
S/phase, final step,200→400 endpoint and invariant gates pass. Full200 native
and NumPy-residual K/diagnostics/sampled A/B/X hashes are exact. Producer cost
4.31638h wall /6.04903 sampled CPUh /.27837GiB. Maximum recorded accepted error
3.24329e-18 is below but close to32-extended-epsilon3.46945e-18; minimum sampled
binary64 rcond9.45956e-4 is not a global condition or forward-error bound.
Every solve keeps its backward gate; only48 native solves /96 residual
evaluations receive independent in-process audits. All12 captured full residual
matrices and36 paired80/110-digit entries are separately reconstructed.

[Independent public saved verification](qualification-evidence/numerov-co-radial-pilot/README.md)
reconstructs all131 decisions, exits0/0, with26 remote tests passing. At15:03UTC
the separately frozen **single full400-bohr numerical replay is running**, after
fresh public-download authentication. Expected7800s /upper10800s /cap14400s on
cores4–11 /one single-threaded worker /4GiB. Scientific records and full captured
array entries must remain exact; auxiliary rcond alone has the prior rounding
budget. **Next:** collect its original verdict/exits/resources and publish/fetch
the numerical replay before assessing a bounded interior CO comparison. One
identity through400bohr does not qualify other identities,800bohr or global
conditioning. Interior false/exits1/1,132/594 historical coverage and zero new
unique identities remain; no automatic60-solve release or new inner solve.

**20 October,17:04UTC closure:** the full400-bohr propagation completes and its
scientific record and three full captured operand sets reproduce exactly. The
original wrapper fails exits1/1/owner1 because nine tuple coordinates compare
unequal to JSON lists. A separately frozen saved-only diagnostic reproduces this
schema failure, reconstructs131 original and22 fresh decisions, and checks all
captured entries exactly; exits0/0/owner0,12.29s, zero new propagations. Auxiliary
rcond drift6.77226e-16 passes its unchanged rounding budget. The original failed
wrapper and all output bytes remain public; no full-row retry is needed to
establish the saved result. All26 payloads fresh-verify. This finite numerical
milestone is closed; **next cost electronic-model and resonance-width decisions**
using existing native results first. Additional numerical work must address a
specific error capable of changing the intended physical calculation.

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

**21 October,18:20UTC:** [signed target-CSF/address correspondence](NATIVE_CSF_MAP.md)
passes all8 sectors /158,994 CSFs /all74 native roots; only16 singlet-A1 signs
change. Independent canonical determinant decoding and fresh public replay pass.
Native source confirms phase-aware import and source centre-count handling:
retain original CIDATA unchanged, no QC/target/property rerun. Producer16.496s /
.60803GiB, independent38.100s. Next qualify live native read-only consumer,
integral target subblock and actual scattering phase/property/channel interfaces.
Both A1 physical subspaces still reject; no scattering solve or sweep runs.

**21 October,17:40UTC:** [full-continuum integrals and8 target CSFs](IMMUTABLE_CONTINUUM.md)
are generated at fixed new-frame orbitals,154 MOs /141 continuum, no QC/target
diagonalization/property/scattering solve. The original CONGEN trailer-check
failure/exits1/1/1 stays beside a seven-sector completion0/0/0 with no repeated
native work. Total15.049s /peak2.46218GiB. Full-layout CSF hashes match baseline,
but signed native-vector/address conversion remains unqualified; prove it or cost
a finite native reimport before compatible target/property/scattering use. Both
A1 physical retained subspaces still reject. No scientific job or sweep runs.

**21 October,16:49UTC:** the [physical ensemble-frame screen](ENSEMBLE_FRAME_COMPARISON.md)
**rejects both A1 retained subspaces**, minima0.06930/approximately zero versus0.5,
while all eight numerical cases pass. Original scientific exits1/1/1; cost63.503s /
.49293GiB. Old A1 root5 maps to new6; new triplet-A1 root3 is poorly captured even
in the old twelve-root buffer. The proposed clean same-manifold contrast is not
released. Keep the baseline; explicitly distinguish combined prescription
sensitivity from a separately justified matched-manifold selection. Design and
cost immutable-frame continuum/native preparation before engine work, without
another QC optimization. New CI waves are saved; no automatic root/channel ladder,
changed property set, Hamiltonian or production sweep. No scientific job runs.

**21 October,15:38UTC closure:** [explicit40-state properties](FIFTY_SELECTED_PROPERTIES.md)
pass native/reader gates,40 IDs /676 moments, native/container/owner0/0/0,
37.203min /.43713GiB. No scientific job runs. Qualified new-frame data now support
physical cross-frame channel/subspace correspondence and compatible sparse
preparation/cost design without repeating QC. Matching labels/counts and the620
exact same-frame common moments do not establish cross-frame identity. No new
inner solve, model-convergence claim or production sweep is released.

**21 October,14:54UTC advance:** the [completed new frame](FIFTY_FRAME_COVERAGE.md)
passes48 independent fixed-orbital probes and74 native binary root/vector/EOF
checks, exits0/0/0,5.144min /.27065GiB. Original seed rejection stays unchanged.
An [explicit40-state property export](FIFTY_SELECTED_PROPERTIES.md) is running
because global-low50 retained only3 singlet-A2 states. Forecast30min /cap60min;
no new QC or target diagonalization. Then physical channel correspondence and
compatible scattering-preparation/cost design; no inner solve is yet released.

**21 October,13:47UTC result:** [fifty-component native target/analyzer pass,
seed-frame preservation rejects](FIFTY_COMPONENT_IMPORT.md), original verifier
1/1/owner1. Current native/QC agreement5.29468e-10 Hartree is separate from
six seed-root shifts up to3.38450e-7, exceeding1e-7. Do not inherit old-frame
coverage, relax that gate or repeat QC automatically. Next cost/freeze independent
coverage for the completed new frame, then map/cost matched ensemble scattering.
Actual target cost2.899h /5.518GiB. No successor is running.

**21 October decision update:** follow the [costed decision queue](PRODUCTION_DECISIONS.md).
The equilibrium CAS11 fifty-component **native target import is running** from
its already covered checkpoint, forecast1h expected /3h hard cap on cores4–11 /
24GiB. This advances the independent orbital-ensemble axis after a+0.456512eV QC
ground shift. A pass releases mapping/costing for matched scattering, not a sweep.
Use matched real-axis phases for model contrasts while widths remain unqualified.

1. Design and independently validate an extraction rule addressing the measured
   background failure using known-pole phases and unused-data predictions;
   preserve every rejected all-order gate.
2. Use the completed three-anchor channel effects and same-engine active-space
   rejection to declare independent ensemble/basis controls, including target
   QC/import gates and explicit measured-cost forecasts before native launch.
3. Complete the remaining matched model/basis/ensemble controls,continuity and
   neutral correlation/domain qualification;then freeze and pass a restartable
   production pilot. Full-sweep readiness follows the gates,not a fixed date.
