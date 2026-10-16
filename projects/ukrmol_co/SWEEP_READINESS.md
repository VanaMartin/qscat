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
import passes. The current experiment is the
[compressed matched A2 scattering control](a2-compressed-scattering-contract.json),
following [saved model/extraction controls](SAVED_MODEL_CONTROLS.md).

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

**Next:** follow the active matched compressed scattering control,then
publish/fetch remaining native-import evidence. CAS11 16K/TZ public reconstruction already
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

**Active:** compressed five-versus-nine-A2at2048passes dense replay and the
fresh five-A2baseline. Independent37-file reading confirms2.2e-6rad phase
agreement with dense,cost1.651h/6.011GiB. The declared nine-A2model is running.
Shared48-state target representation/original orbitals and audited8-rank
runtime isolate retention. The compressed matched2–5requested-eVwindow retains
orders1–4and does not replace the earlier24-case extraction rejection.
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

**Next:** audit the two
remaining retained components: inspect competing row scores,energies and
complete nearly degenerate subspace projectors. Determine whether a larger
retained channel manifold is required. Test only geometries/model controls
that address those specific failures. Track relevant resonance features across
neighboring geometries with phase/K-matrix evidence and consistent thresholds;
electronic-state assignment alone does not establish resonance-pole identity.
Revisit fit windows/background orders,energy resolution and threshold effects
at all anchors on the selected scattering model.

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

1. Read and reconstruct the active compressed nine-A2/2048model and its matched
   channel comparison. Equilibrium passes and stretched rejects;use the missing
   anchor effect/cost to select further model/basis/ensemble controls.
2. Design and independently validate an extraction rule addressing the measured
   background failure;preserve every rejected all-order gate.
3. Complete the remaining matched model/basis/ensemble controls,continuity and
   neutral correlation/domain qualification;then freeze and pass a restartable
   production pilot. Full-sweep readiness follows the gates,not a fixed date.
