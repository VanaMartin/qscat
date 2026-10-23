# CO production decisions

The goal is a geometry sweep supplying a physically grounded two-dimensional
potential, followed by held-out scattering validation. Engine/interface progress
is supporting work, not evidence that the electronic model or potential is ready.

## Goal-convergence checkpoint — every one or two days of campaign work

Review at least every48 hours during active work, before opening a new experimental
branch, and after a decision-changing result. On resuming after a gap, review
before launching work. An authorized long job need not be interrupted just because
a review is due. This is a working/handoff rule, not an unattended scheduled task.

- [ ] **Goal:** name the sweep/potential deliverable this cycle serves, including
  geometry domain, observables and required accuracy. Do not silently substitute
  an exploratory, model-conditional potential for the validated deliverable.
- [ ] **Critical path:** identify the production blocker: electronic-model choice,
  position/width extraction or direct scattering fitting, continuity/neutral
  reference, or restartable geometry pilot. Label optional/supporting work as such.
- [ ] **Decision before experiment:** what will a pass, rejection or inconclusive
  result change? Name the choice, eliminated candidate, reuse or stop decision.
  If no outcome changes a production decision, do not launch the experiment.
- [ ] **Reuse/minimum scope:** use saved data and accepted numerical controls.
  Explain why this is the smallest decision-changing check; no repeated QC or
  properties, and no automatic root/channel/refinement ladder.
- [ ] **Bound the cycle:** state compute cost/caps, implementation/diagnosis effort
  allowed and a stopping rule. A seconds-long calculation can still consume days
  of supporting development. Continuing that branch into another cycle needs a
  renewed decision rationale and effort bound, not merely its next technical task.
- [ ] **Evidence delta:** compare cycle start/end. Which physical uncertainty
  shrank, candidate was ruled out, production choice became possible or pilot
  prerequisite was removed? Explicitly list unchanged blockers. Test/payload/
  coefficient counts are not readiness metrics; a useful rejection is progress.
- [ ] **Drift:** are we solving an interesting subsidiary problem instead of the
  named blocker? Is an optional comparison becoming a universal sweep launch gate?
- [ ] **Disposition:** record **continue / redirect / stop**, why, one next
  deliverable and the next review trigger. Without a production-relevant evidence
  delta, default to reassessment/redirection rather than another refinement.
  Preserve failed gates; do not lower tolerances to manufacture convergence.

Each dated checkpoint records UTC date, starting blocker, bounded question,
evidence delta, remaining uncertainty, disposition and next deliverable. Report
goal-level progress to the user first; link technical evidence separately.
The checkpoint itself does not authorize an expensive calculation.

### Checkpoint — 23 October 2026: redirect toward production decisions

**Starting position:** reproducible/affordable scattering is substantially in
place. Electronic-model selection and reliable widths remain open; no qualified
potential or production-style geometry pilot exists. Neutral data are reusable
over1.9–2.5bohr, not a qualified wider-domain reference.

**Evidence delta:** saved target reuse and scattering phase/property mapping now
avoid repeating QC, CONGEN and properties for one optional comparison. These are
preparation obstacles removed, not electronic accuracy or width identifiability
settled. Both A1 physical subspaces reject a clean same-manifold ensemble contrast.

**Disposition: redirect.** The fifty-component branch is not a prerequisite for
every possible sweep. Retain the original forty-component frame and accepted
sparse/root controls. [The route decision](POTENTIAL_ROUTE_DECISION.md) parks that
optional comparison, including remaining native consumer closure, for this cycle.
Reopen only if a specified production decision needs it, with an explicit effort
bound and finite contract. No repair campaign, rebuild, refinement ladder or
multi-hour inner solve follows automatically.

**Completed next-step decision:** [model/extraction route triage](POTENTIAL_ROUTE_DECISION.md)
keeps the existing model-specific references without prematurely freezing a
production model; unqualified widths cannot yet supply potential-factory targets.
Next specify one finite feasibility test choosing physically restricted extraction
versus direct scattering fitting, or name the single missing ingredient. One
working-cycle design/review cap; no new numerical job is authorized by this note.
Do not open another propagation/pole-method campaign without that decision link.

**Next review:** after that deliverable or any closure failure, and by the next
one-to-two-day working-cycle boundary (25 October if work remains active).
Next readiness milestone: choose an explicitly qualified model and extraction/
fitting route, then a small restartable geometry pilot—not more interface passes.
No sweep or dense Hamiltonian is released.

**Bounded next test selected,23 October:** [direct-S representability](../../docs/physics/co-direct-scattering-representability.md)
uses existing six-channel CAS11/nine-A2 matrices to test whether the bare
isotropic potential family is already incompatible with full physical-channel
scattering. Rejecting its larger diagonal-unitary relaxation avoids building
an unsuitable fitter. It neither rejects tail-dressed scalar reductions nor
qualifies widths. Freeze/tests precede the one saved-only60s expected /300s cap
application;480 development samples,114 prospective route holdouts left out.
No new native solve, propagation, QC, new50 work or inner Hamiltonian. Any outcome
returns to the goal decision, rather than automatically expanding the family.

### Goal checkpoint update —23 October,10:12UTC: bare direct-S shortcut stopped

[The bounded screen](DIRECT_SCATTERING_SCREEN.md) rejects even arbitrary diagonal
unitary phases at all three anchors/both sectors:392/480 development matrices
fail the mixing bound;114 route holdouts remain unanalyzed. Maximum necessary
relative-eigenphase distances0.492 /0.523 /0.439rad exceed the newly declared
0.05rad full-S metric. This is not a retroactive failure of eigenphase-sum gates
or a rejection of tail-dressed scalar models.

**Evidence delta:** an incompatible fitter specification is eliminated before
implementation, in0.872s container time. **Unchanged:** no electronic model
chosen, no qualified three-anchor width curve, no potential or production pilot.
**Disposition:** stop the bare full-S fitter; return to the physical target choice.
Specify one justified channel-to-effective-model map, retained tail/background
and validation observable within the next bounded working cycle, or report the
assumption needing a decision. No coupled fitter, general pole campaign, inner
solve or new50 restart follows automatically. Review on that deliverable or next
working-cycle boundary.

## Supporting branch: equilibrium orbital/channel-prescription comparison

**Parked by the23 October goal review.** The historical technical next steps
below are conditional, not the current campaign queue or a universal sweep gate.

### Current decision — 21 October,21:55 UTC

**Saved scattering phases/root/property mapping now passes:**
[317,988 phase-defining terms](SCATTERING_INPUT_MAP.md) reconstruct exactly across
B1/B2. Saved CONGEN algebraic definitions need not be regenerated. The new
forty-state property file needs **25 changed idtarg positions**; all selected
energies match original CIDATA within4.91639e-10Ha. Properties/CIDATA stay
unchanged;1,410 expanded coefficients per symmetry. Parser failure1/1/1 is
preserved beside parser-only completion0/0/0,6.833s /.703GiB process RSS.

The remaining interface qualification is **native CIRMAT with the actual B1/B2
scattering destination phases**, plus downstream property/channel semantics,
not another CONGEN, QC, target solve or property export. The saved-only result
does not qualify those native consumer calls or overturn physical A1 rejection.
Compact public replay verifies all first-term witnesses, not the full1.4GB
original N+1 CSF files. Keep lowest-five sensitivity labelled combined orbital/
channel-prescription sensitivity; only after this interface and scientific
contrast are fixed cost the sparse2048-root comparison.

**Native target reuse is now supported:** [the original native consumer and saved
target integral blocks](NATIVE_TARGET_INTERFACE.md) pass. All74 native vectors /
relative energies are byte-exact under the proved signs; h/overlaps exact,
max ERI delta4.44089e-15Hartree, conservative common-basis operator-difference
bound6.26708e-13Hartree. No target/QC/property repetition is needed. Three build/
ABI failures remain1/1/1; a matched-double-header integral-only completion passes
0/0/0 without repeating the successful target read. Four-job total3.761s /
peak.12254GiB; no native diagonalization or scattering. Native header mismatch
is repaired for this client and guarded; future Docker header copying is corrected
but no image rebuild is claimed.

Next qualify actual **scattering CONGEN target-channel phase/root/property maps**
and declare combined-prescription sensitivity versus separately justified matched-
manifold selection before costing a sparse scattering solve. Do not describe
both A1 physical rejections as resolved by this technical target reuse pass.

**Earlier signed-basis obstacle resolved:** [all eight target-CSF pairs](NATIVE_CSF_MAP.md)
map exactly,158,994 CSFs /1,670,879 terms per layout; only16 singlet-A1 signs
change. All74 native target vectors preserve root IDs/energies and round-trip
coefficient bytes exactly. Independent canonical expansions and fresh public
replays pass. Native `CIRMAT/READCIP/CPAKPZ` already reads source centre metadata
and corrects phases: **retain original CIDATA unchanged**, not a converted or
padded file. No QC/target/property repeat is justified for this layout change.
Producer16.496s /.60803GiB; independent38.100s, no native calls. Next bound the
native read-only consumer, target integral-block and actual scattering-channel
phase/property interfaces. The physical A1 rejection still requires an explicit
combined-prescription or matched-manifold scientific choice before scattering.

**Earlier17:40UTC technical advance:** [full-continuum integrals and eight target-CSF files](IMMUTABLE_CONTINUUM.md)
are now generated on the unchanged new frame:154 MOs /141 continuum, exact
baseline continuum cards. One integral plus eight CONGEN calls total; no QC,
target eigensolver, properties or scattering Hamiltonian. An original CONGEN
trailer-check wrapper failure/exits1/1/1 is retained; a separately frozen seven-
sector completion reuses its integral/first-CSF outputs and passes0/0/0.
Total15.049s /peak2.46218GiB. All eight CSFs match recorded baseline full-layout
hashes. At that milestone, signed coefficient/address reuse was still unset;
the18:20UTC result above now settles the finite signed CSF correspondence.
The physical A1 rejection and the need to declare the scientific contrast remain.

The [physical ensemble-frame comparison](ENSEMBLE_FRAME_COMPARISON.md) now
**rejects retained-subspace dominance in both A1 sectors**, minima0.0692964 and
approximately zero versus0.5; six sectors pass. All eight numerical comparisons
complete, original scientific exits1/1/1, **63.503s /.49293GiB**. Old retained
singlet/triplet A1 root5 maps to new root6 at squared overlaps0.9075/0.9838.
New triplet-A1 root3 has only0.3886 total capture in the old twelve-root buffer.
This is not fixed by relabeling root5 alone. No scientific job remains.

The new orbital frame and explicit40-state properties remain qualified, but the
proposed **clean same-manifold orbital-average contrast is not released**.
Keep the old frame as the current baseline. Choose explicitly between a combined
lowest-five orbital/channel-prescription sensitivity test and a separately
designed matched-manifold selection/projection. Do not describe the former as an
isolated orbital-average effect, or silently tune the retained set into a pass.

For either route, use the completed **immutable-frame continuum/native preparation**
instead of `--target-initial-checkpoint`, which reoptimizes. Original phase-aware
target vectors are reusable under the proved representation map; actual native
consumers and the scattering target-channel map remain to qualify. Cost each
finite interface stage before submitting it, not a whole-job ladder. The
new CI arrays are saved; no QC, root/channel ladder, new properties or dense
Hamiltonian follows this rejection. [Public numerical replay](verification-ensemble-frame-public-20261021.json)
reproduces all eight overlap blocks within2.55351e-14, preserving false/exits1/1/1.

### Original recipe and subsequent milestones

The existing CAS11 fifty-component QC frame changes its ground energy by
**+0.456512 eV** and the common-forty objective by **+0.170147 eV**, relative to
the forty-component frame. These are diagnostics of different optimized frames,
not estimates of a resonance error. The independent fixed-orbital coverage already
passes48 probes on each of the CAS11/CAS12 frames. The missing step is native
import followed by matched scattering.

The [one-job recipe](calibration-sa11-fifty-import-20261021.json) and
[finite contract](fifty-import-contract-20261021.json) advance the equilibrium
CAS11 frame first:

- Same CAS(10,11), two frozen A1 spatial orbitals, cc-pVDZ and R=2.1323 bohr.
- Fifty equal-weight orbital-average components: `[10,5,5,5]` per spin.
- Native roots `[13,8,8,8]` per spin, **74 total**. All have independent
  fixed-orbital coverage. The DENPROP import setting uses50 states.
- Start from the existing converged fifty-component checkpoint. Preserve the
  ordinary200-cycle CI cap and original QC tolerances. Reoptimization must
  preserve every ensemble root within1e-7 Hartree, dipole within1e-5 a.u. and
  core/active singular values above1−1e-8.
- Verify all74 binary native eigenvalues/vectors and complete EOF, the native
  solver's requested/converged counts, all50 averaged-root imports and ground
  dipole. This qualifies the input frame, not ensemble convergence.

**Cost:** earlier same-basis64-root import3086.6s /5.266GiB on four ranks.
New forecast3600s expected /7200s upper /10800s hard cap; eight MPI ranks on
physical cores4–11, single-threaded native libraries and24GiB memory. The
existing PySCF target builder uses eight threads during QC. Host admission
requires32GiB available RAM /30GiB disk. Extra ranks and nonuniform74 roots have
not been benchmarked. All output, original exits and source/input hashes persist
under `/home/kooza/ukrmol/co-sa-20261006` on Sadaharu.

The subsequent scattering comparison must retain the **same physical channel
manifold** on both frames. Using50 DENPROP states here does not authorize a
50-versus40/48 channel comparison. Before a new inner solve: reconstruct explicit
property/dense-row maps, inspect physical subspace correspondence and retained/
omitted exchanges, and cost the resulting sparse2048-root solve. A2 root9 is not
covered by this74-root import. If a48-state nine-A2 comparison is selected, its
extra A2 input qualification needs a separately declared finite scope.

For orientation, the old equilibrium CAS11 full dense pipeline cost4.423h /
25.201GiB; this is **not** a forecast for a newly mapped sparse calculation.
The accepted CAS12/2048 baseline and2048→4096 control remain reusable.
No dense scattering Hamiltonian is authorized. No automatic channel/root ladder
follows a rejection.

**Completed21 October:** target and analyzer pass, but seed-frame preservation
rejects6/50 averaged roots, maximum3.38450e-7 Hartree against1e-7. Native energies
agree with the newly optimized frame to5.29468e-10 Hartree. Original verifier
exits1/1/owner1; overall import remains unset. Actual cost **2.899h /5.518GiB**,
dominated by135min QC and35.5min DENPROP, not the3.33min eigensolvers.
The [result/evidence](FIFTY_COMPONENT_IMPORT.md) preserves all153 original files.
The next candidate is independent coverage on the completed **new** fixed-orbital
frame, not another reoptimization or a weakened seed gate. Earlier96 probes on
two frames cost1280.58s; a one-frame forecast of10–20min requires its own finite
contract before submission. No successor is yet running.

**21 October,14:54UTC advance:** the [new fixed frame](FIFTY_FRAME_COVERAGE.md)
now passes all48 independent CI probes and all74 native binary root/vector/EOF
checks, exits0/0/0, in **5.144min /.27065GiB**. The old seed gate stays rejected.
The global-low50 property export contains only3 singlet-A2 states. An explicitly
selected40-state [property-only calculation](FIFTY_SELECTED_PROPERTIES.md) now
runs, forecast30min /upper45min /cap60min, using the qualified vectors and no new
QC/target diagonalization. Then establish physical cross-frame channel identity
and compatible scattering preparation/cost, rather than comparing unequal maps.

**21 October,15:38UTC closure:** the selected40-state property export and
independent reader pass **40 IDs /676 moments**, native/container/owner0/0/0,
37.203min /.43713GiB. All620 shared moment values across38 same-frame shared
states are exact. All8 public payloads and9 prerequisites authenticate. No job
remains. The next decision is now **physical cross-frame retained-channel/
subspace correspondence and compatible scattering preparation/cost**, using the
qualified frame without another QC optimization. No new inner solve follows
automatically; matching labels/counts are not state identity.

Eight candidate cached equilibrium12-root CI arrays and the validated physical
`ci_overlap.py` are available from the earlier state-overlap campaign. Authenticate
them against the actual scattering-frame checkpoint before reuse. Its full
three-anchor24-sector /48-matrix run cost4322.37s /1.177GiB. The new frame's
coverage saved energies/residuals but no CI arrays; cost only the missing new-frame
wavefunctions and physical comparisons, not another three-anchor campaign or QC
optimization. No physical-overlap successor is yet submitted.

## Decision queue and release conditions

Ordered by the23 October goal review. Rows are scientific decisions, not an
automatic sequence of calculations; parked work does not block every route.

| Decision | Existing evidence | Next deliverable / release condition | Cost basis |
|---|---|---|---|
| Fitting-information route | Width extraction remains unqualified; bare scalar full-S shortcut rejects392/480 development matrices | Define one channel-to-effective-model target map, retained tail/background and validation observable; no automatic coupled fitter | Completed screen0.872s; next target-definition cycle only, not a new numerical campaign |
| Active space / channels | CAS11→CAS120.1213rad and64–73meV candidate-position shifts; stretched five→nine A2 changes0.0679rad | Select model with explicit decision-changing checks and remaining uncertainty budget; no accuracy claim from a larger model alone | Reuse existing2048 results; no automatic active-space/channel ladder |
| Basis | Qualified CAS11 larger-basis target frames exist; ground-state accuracy alone is insufficient | One matched basis contrast on the selected ensemble/channel frame, phase first | Reuse imported targets; new inner solve requires its own measured forecast |
| Continuity | Retained/omitted exchanges and ambiguous coarse links remain | Neighboring-geometry physical subspace/resonance correspondence for the selected model | Bound the geometry pilot after the model is selected |
| Neutral reference |52 QZ/5Z points pass basis/interpolation over1.9–2.5bohr | Correlation-treatment uncertainty over this domain; widening to3–4bohr needs a stable reference | Reuse52 points; previous46 new points cost3.21h wall sum |
| Production launch | Model/extraction choices open | Frozen restartable small pilot, held-out geometries/energies, measured time/disk forecast | Cost from the selected model; no defensible total forecast yet |
| Optional orbital ensemble — parked | Both A1 subspaces reject; saved scattering phase/property map and target reuse pass | Reopen only if combined-prescription sensitivity could change a named model choice; native consumer work is not complete or a universal sweep prerequisite | Preserve all inputs; no QC, CONGEN, property or native-interface repeat now |

The position/full-width/phase scientific budgets remain **0.05eV /5% with
0.001eV floor /0.05rad**. A failed extraction cannot supply a certified width
comparison; phases and target properties can still reject model insensitivity.
A passing pair is finite sensitivity evidence, not an absolute electronic error.

## Width work that can proceed from existing data

The discrepancy is named: compressed best inner-held phase errors0.00225936 /
0.00130390rad, equilibrium broad-window0.00108757rad, against0.001rad. Adding
background flexibility can reduce those errors while destroying identifiability.
The94 counterexamples establish that a local analytic/unitary phase factor does
not qualify as a causal finite-range CO Hamiltonian. No replacement width is
published merely because a more flexible function predicts phases.

The next extraction design must specify the physical background class and
threshold dependence, the meaning of its pole, its known-truth generators and
mandatory negative controls **before** a new CO application. Reuse the untouched
audit points, original rejected windows and counterexamples. A failure to find
another width is not proof of uniqueness. In the meantime, the model decision
uses matched real-axis phases without waiting for a new numerical-pole campaign.

## Numerical scope supporting these decisions

The radial full400 scientific record and three full captured operand sets are
exact. The original wrapper schema failure/exits1/1/1 remain preserved beside
the passing saved-only diagnosis. One radial identity does not qualify every
energy,800bohr or complex propagation. Original interior false/exits1/1 and
historical132/594 coverage remain. A further numerical pilot must name the
specific accuracy question affecting the proposed physical calculation and its
finite cost; filling every matrix or extending every radius is not a launch gate
by itself.
