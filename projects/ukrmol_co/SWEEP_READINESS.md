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
Compressed proof replay passes all438 payloads/96 probes; its fresh private-cache
native import is active.

## 1. Finish trustworthy target inputs

**Current:** CAS11 DZ imports/scattering pass at three anchors. Equilibrium and
stretched-DZ CAS12 imports pass; supported stretched aug-TZ now passes all64
root/dipole/subspace/seed-preservation checks. The original stretched TZ import
failed seed preservation; its reoptimized orbital frame is a separately
qualified numerical target. Compressed CAS12 space160 ordinary QC and96-probe
coverage pass; compressed native import is being continued with complete public
proof and a private script cache.

**Next:** finish that compressed import and independently verify its64 roots,
spin0/direct-spin1 actions,dipole and orbital subspaces against the original
covered seed. Publish/fetch the remaining expensive evidence,including CAS11
16K and the TZ-frame diagnostic. Preserve every original failed exit.

**Pass:** ordinary QC and fresh-CI flags; roots≤1e-7 Hartree; physical/penalized/
spin-penalty actions≤1e-9 Hartree; dipoles≤1e-5 a.u.; spin-square errors≤1e-6;
MO orthogonality≤1e-9; Pi degeneracy and predeclared seed/subspace checks.
Successful import establishes engine/orbital consistency,not model convergence.

## 2. Establish an affordable CAS12 scattering calculation

**Current:** CAS11 selected-root scattering passes both-sector independent
eigenpair,boundary and observable checks;2048 roots is the smallest qualified
count. A further16K confirmation cost19.65 hours without changing that decision.
CAS12 dense workspace exceeds222 GiB and does not fit Sadaharu.

**Next:** after CAS11 public review,declare one modest equilibrium CAS12 sparse
pilot,with explicit initial root count,8 physical cores,MPI ranks,energy/fit
window,resource/time limits and diagnostic outputs. Confirm native root selection,
residuals and boundary exports against independent matrix actions; compare a
small higher-count control and stable omission estimates. A passing pilot then
extends to compressed and stretched anchors. CAS11's passing count is an initial
cost guide,not a CAS12 convergence result.

**Pass:** physically relevant phases and fitted position/width stabilize under
the declared spectral control. Use the existing observable budgets: **0.05 eV
position,5% full width with0.001 eV floor,and0.05 rad phase moduloπ**. Independent
eigenpair/boundary checks must also pass. If no tractable spectral comparison
can bound truncation,the method remains unqualified. Record solve/export/outer
cost separately; choose the smallest passing count. Escalate only for a named
observable/truncation discrepancy.

## 3. Select a converged electronic and channel model at the anchors

**Current:** the matched CAS10→CAS11 equilibrium change exceeds position/phase
budgets. Basis,active-space and retained-channel convergence remain open.
Fifty-component equilibrium orbital coverage passes,but native fifty-component
imports and matched scattering have not been performed.

**Next:** compare CAS11→CAS12 and DZ→TZ→aug-TZ with all other settings fixed.
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
Native background/window changes move widths by roughly12% at equilibrium,
32% compressed and6.7% stretched,exceeding the5% width budget.

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

1. Complete checksum-pinned compressed proof replay,then its fresh native import
   under [the private-cache successor](cas12-compressed-private-cache-successor-contract.json).
2. Run [independent quarter reconstruction](cas11-quarter-state-tracking-replay-contract.json)
   on reserved cores and publish/fetch/reconstruct passing evidence.
3. Complete CAS11/TZ public evidence and review,then declare a modest CAS12 pilot.
4. Use those outcomes to choose the next specific electronic/channel/extraction
   experiment. Full-sweep readiness is determined by the gates,not a fixed date.
