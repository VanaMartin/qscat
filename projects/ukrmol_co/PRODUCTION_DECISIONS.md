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

### Goal checkpoint —23 October,10:56UTC: target map defined, no numerical launch

**Starting blocker/question:** what part of coupled CO scattering can the scalar
potential represent without repeating the rejected bare full-S fit?
**Decision:** use a continuously assigned physical outgoing pole, with full width
Gamma = -2 Im z, as the prospective factory target. The
[specification](../../docs/physics/co-effective-resonance-target.md) retains
native polar/closed-channel shifts in that pole and explicitly omits residues,
channel branching and background. Do not add the polar shift twice.

**Evidence delta:** the target, energy/sheet conventions and claim boundary are
now explicit. Validation must predict an excluded geometry's native pole; a
resonant determinant factor reconstructed at a training geometry is not a test.
A reduced pass is not full-scattering qualification. This is a design decision,
not new physical evidence, an uncertainty reduction or a qualified CO width.

**Reuse/bounds:** current source and saved reference metadata only; zero new
scientific calculations, zero holdout S matrices read. The bounded design cycle
is complete. Factory T1 requires at least four actual tracked training geometries
plus a separate validation geometry; the current three anchors cannot supply
that test. No asymptote, width floor or extra geometry is invented to run it.

**Unchanged:** electronic-model choice, qualified pole/width curve, neutral
correlation/wider-domain reference, continuity, full-scattering embedding and
restartable pilot. **Disposition: stop before implementation.** Next decision:
authorize a one-anchor pole-reference feasibility contract with explicit effort/
compute caps and qualified prerequisites, or redirect to a coupled effective
model. No solver, coupled fitter, QC/inner solve or new50 work is launched.
Review before either branch is opened, or at the next active working-cycle boundary.

### Continuation decision —23 October: one-anchor outgoing-map feasibility

**Production link:** the pole-target map has a concrete missing input: a qualified
complex outgoing reference. Before any CO root search, test whether one small
general-channel outgoing-map implementation is numerically viable on saved data.
The [canary design](../../docs/physics/co-outgoing-map-canary.md) and
[contract](outgoing-map-canary-contract-20261023.json) select compressed/nine/B1,
indices 60/70/80 only, then one fixed non-pole complex point with finite contour
variations. All real full-S gates must pass before complex probes. Five analytic/
independent BVP, boundary-Numerov, flux and sheet/cap controls pass locally.

**Bound/reuse:** no new native calls or inner/QC/property solves. At most eleven
propagations, 90s each, 300s runtime allocation /900s upper /1800s hard cap, four
physical cores reserved, one BLAS thread, 4GiB. These are finite allocations,
not measured CO forecasts. One working-cycle implementation/diagnosis cap; no
retry, precision escalation or follow-on solver campaign. Preserve failures.
The equilibrium pole exclusion and all 114 route holdout S values stay protected.

**Disposition: continue only with this canary.** Pass supplies finite map
feasibility, not a qualified pole/width or full scattering model. Accuracy
failure rejects this implementation/scope; timeout or numerical failure is
inconclusive and returns to route selection. Neither releases extra geometries,
an automatic pole search, model choice, fitter or sweep. Review immediately on
the result or failure, before any next numerical branch.

### Canary goal review —23 October,13:14UTC: packaging stop, one narrow successor

Original `c00c09a` application stops before any CO propagation: the independent
test grid helper imports absent `mpmath`. Original command/container/owner exits
**2/2/2**, scientific verdict unset; all123 inputs and source hashes verify.
This is dependency packaging, not a physical or numerical rejection. Wall .972s;
no holdout, pole search or new native call. Preserve the original directory.

**Disposition: continue only through one separately frozen dependency-only
successor**, [contract](outgoing-map-dependency-successor-contract-20261023.json).
Replace the incidental test helper import by polynomial-coefficient
differentiation at identical nodes, retaining the same BVP and all tolerances.
The propagator, runner and scientific contract remain byte-identical and pinned.
No dependency install, image rebuild, numerical retuning or scope expansion.
This small change is the entire packaging-diagnosis allowance; another failure
ends the cycle with no further successor. Runtime remains90s per propagation /
1800s total, four reserved physical cores /one thread /4GiB. This is a
goal-reviewed exception, not relabelling the original stop or authorizing a
general repair campaign. Review immediately on the successor result.

### Goal checkpoint —23 October,13:21UTC: finite outgoing-map feasibility passes

**Question completed:** can the missing complex outgoing map be tested on one
saved anchor without a general solver campaign? [Yes, in a finite scope](OUTGOING_MAP_CANARY.md):
three compressed/nine/B1 development energies reproduce full S at two endpoints,
worst entry error1.51042e-7 versus1e-6. Five settings at one complex non-pole point
give largest radius-scaled map change2.06671e-11 versus1e-7. All40 row gates pass;
all123 inputs verify. Fresh public saved reconstruction independently rebuilds
every gate, preserving original dependency-stop2/2/2 and successor0/0/0.

**Evidence delta:** the route's298-channel outgoing-map feasibility and recurring
cost are now measured, not inferred from a two-channel toy. Container55.850916s,
3.53–8.27s/map, cgroup peak.32147GiB; original packaging stop.972009s retained.
No numerical source/tolerance change, engine call, inner solve, pole search or
holdout S analysis. This reduces feasibility uncertainty, not CO width/model
uncertainty. Complex invariance is local, not independent CO pole qualification.

**Disposition: continue to one bounded pole-reference proposal.** Define a finite
compressed-anchor search region, sheet, root-evaluation cap, independent
construction/residue witness and error gate before execution. If that cannot be
bounded, report the blocker; no automatic propagation/precision campaign.
Electronic model, neutral correlation/wider domain, continuity, extra training
nodes and full-scattering embedding remain open. All114 route holdouts and the
equilibrium pole exclusion stay protected; new50 stays parked. Review on the
proposal or next active working-cycle boundary, before any root search.

### Selected continuation —23 October: one finite-model pole and visible residue

**Production question:** can a numerically witnessed outgoing pole replace the
unidentifiable phase-fit width at the compressed development anchor? The
[design](../../docs/physics/co-compressed-pole-reference.md) and
[contract](compressed-pole-reference-contract-20261023.json) bound one box and
seed pair in incident Hartree, three prescribed map settings, a nested16/32-node
interior resolvent contour, and derivative/contour agreement of a nonzero6x6
open-channel residue. An invisible interior state or a matching-coordinate
artifact must not qualify. Eleven analytic/independent controls pass locally.

**Reuse/cost bound:** unchanged native2048-root/298-channel compressed/nine/B1
inputs and completed map-canary evidence; no engine/inner/QC/property solve,
geometry, new model, potential fit or holdout S analysis. At most77 maps,8 physical
cores /8 contour workers /one BLAS thread each /8GiB. Root stages sequential.
Measured3.53–8.27s/map informs600s target /1200s upper /1800s hard cap,90s per
map, without a guarantee of conditioning near a pole. Implementation/diagnosis
is limited to this active working cycle. Freeze/push precedes CO root use.

**Disposition: continue only with this experiment.** Pass establishes a pole
of the finite reference with local numerical/residue witnesses; the two pole
constructions share the exterior ODE, so no independent exterior, omitted-spectrum
bound or absolute electronic accuracy follows. A rejected witness/cap/box exit/
failure stops without a repair, extra seed, contour refinement or successor.
Review immediately on the outcome. Model, continuity, neutral, scalar/full-
scattering and sweep gates remain; all114 route holdouts and equilibrium pole
exclusion stay protected, new50 parked.

### Goal checkpoint —23 October,14:55UTC: one visible finite-reference pole

**Question completed:** can the compressed development reference supply an
outgoing pole rather than unidentifiable local phase-fit parameters? The
[finite experiment](COMPRESSED_POLE_REFERENCE.md) finds **3.7604eV position /
2.0846eV full width**, compressed/nine-A2/B1 CAS11,2048 roots/all298 channels.
The boundary zero is recovered through interior contour moments and has a
nonzero rank-one open-channel residue. Original exits0/0/0;47/77 allocated maps,
180.990698s wall /1.94329GiB cgroup peak, eight physical cores. Independent public
saved reconstruction verifies every witness with no new propagation. All123
original inputs and25 source entries verify; eleven controls pass locally/runtime.

**Evidence delta:** one pole/residue datum of the declared finite model now
replaces a missing local fitting-information ingredient. This resolves that
local pole-identification question, not electronic accuracy or a CO width curve.
Both stricter settings accept the same candidate without new iterates; raw
contour/algebra agreement is not an absolute uncertainty enclosure. Constructions
share the exterior ODE. Omitted spectrum, paired-sector assignment, model,
neutral, continuity, extra training nodes and full-scattering embedding stay open.

**Disposition: continue to one paired-sector/character consistency proposal at
the same compressed anchor**, separately bounded before any B2 work. B1/B2
degeneracy alone is not a Lambda assignment; require open-channel angular
character before a Pi target. Failure would block
the proposed single-Pi assignment; pass should return effort to a matched
electronic-model decision rather than a general numerical campaign. No next
sector/anchor, fit, geometry or root ladder is launched. All114 route holdouts
and equilibrium pole exclusion remain protected; new50 stays parked. Review on
that proposal or the next active working-cycle boundary.

### Goal checkpoint —23 October,17:21UTC: freeze character before its measurement

**Question:** does the compressed finite-reference B1 pole have the intended
open-Pi character, and is a B2 counterpart test worth running? The
[paired design](../../docs/physics/co-paired-pole-character.md) and
[prospective contract](paired-pole-character-contract-20261023.json) use verified
native `(channel,target,l,m)` labels, not channel order, degeneracy or a single l.
The decision thresholds and negative/invariance controls are fixed before any
actual CO residue-character score. Nine synthetic controls pass; those counts
are implementation evidence, not new sweep/model evidence.

**Minimum reuse and bounds:** first classify all four saved B1 residues, no new
map, native call, QC/property/inner solve or geometry. Four physical cores0–3,
one BLAS thread,1GiB,60s hard cap. A resolved rejection or inconclusive B1 screen
stops before B2. If it passes, the prospective B2 allocation is at most67 maps:
two real development-S checks, at most28 root evaluations and the original fixed
37 center/derivative/contour witnesses. Eight physical cores0–7 /8GiB /600s target,
1200s upper,1800s hard cap. The B2 adapter and exact owner must separately freeze
before propagation. No automatic retry or general solver development.

**Evidence delta so far:** label provenance/order is resolved and an economical
necessary-condition test is executable; actual B1/B2 character and pairing remain
unmeasured. **Disposition: continue only to the frozen saved B1 screen first.**
Implementation/diagnosis is bounded to this active cycle, at most two additional
hours after label reconnaissance; stop with the exact missing step if the B2
adapter cannot freeze within it. Review immediately after the B1 result. A full
pair pass would return priority to matched electronic-model pole evidence, not
another geometry or exterior campaign. Electronic accuracy, omitted spectrum,
continuity, neutral/factory/asymptote/full-scattering gates remain open. All114
route holdouts and equilibrium pole exclusion stay protected; new50 stays parked.

### Goal checkpoint —23 October,17:24UTC: saved B1 character passes; bounded B2 eligible

Freeze/push `cee9e3e` precedes all actual character scores. The saved-only B1
application passes at17:23:40UTC, original scientific/command/container/owner
0/0/0/0, container1.322034s, no propagation. Leakage outside the ground-target
`|m|=1` subspace is1.14749065e-9 against the prospectively frozen1e-4 pass gate.
All four residue estimates pass; maximum leakage change1.45e-19, profile TV
change4.10e-14 and relative norm change2.51e-7. The finite-model profile is mixed
in l (l=1 about.3271, l=2 about.6663); this is not a p-wave-only model or a
branching probability. All123 inputs /86 archived files /17 source entries
authenticate; nine controls pass in the pinned runtime. The paired verdict is unset.

**Evidence delta:** the B1 pole's proposed open-Pi support passes a necessary
character condition, without a new solve. B2 pairing remains unmeasured, as do
electronic accuracy, exact axial conservation, omitted spectrum and continuity.
**Disposition: continue to the single B2 counterpart already allocated by the
unchanged prospective policy**, after its adapter/source/owner freeze. The
[execution contract](paired-pole-execution-contract-20261023.json) pins the
original policy, base pole contract and saved B1 result/exits. Only the existing
pole runner's compressed-sector guard changes; its numerical routines/settings
do not. Two development-S checks precede any complex map; at most67 new maps,
eight physical cores /8GiB /1800s hard cap. Original witnesses and same-root
budgets remain, with no B1 repropagation, alternate seed/box, repair or successor.
Review on any B2 outcome; all114 holdouts and equilibrium pole exclusion remain
protected. This is still within the one-cycle/two-hour implementation effort cap.

### Goal checkpoint —23 October,17:36UTC: finite open-Pi pair; return to model evidence

**Question completed:** the [compressed finite-reference B1/B2 pair](PAIRED_POLE_CHARACTER.md)
passes its prospectively frozen angular-character/counterpart conditions, not
degeneracy alone. Leakage outside target1 `|m|=1` is1.14749e-9 /1.14744e-9,
profile TV2.581e-11, both visible open residues and interior contours pass. B2
accepts the B1 starting candidate once per setting; identical reported positions
are not independent root estimates or zero splitting uncertainty. Distinct
contour estimates agree within4.26e-11eV position /3.42e-11eV full width, not an
absolute CO error bound. Source/policy freeze precedes each stage as required.

**Reuse/cost:** saved B1 screen1.322s, zero maps; B2 counterpart42/67 new maps,
173.055525s container /1.83207GiB cgroup peak /eight physical cores0–7. Original
exits0/0/0/0, no timeout/OOM/repair, no new engine/QC/property/inner solve, geometry
or potential fit. Independent fresh public reconstruction verifies10 payloads,
all saved pole/real-S/character gates and123 original inputs without propagation.
Implementation remains inside this active cycle's two-hour effort cap.

**Evidence delta:** a local finite-model open-Pi pair is now supported; the
assignment blocker is removed for this datum. Electronic-model accuracy and
the known active-space/channel sensitivity are unchanged. This is not exact
Lambda conservation, a qualified CO curve, scalar fit or full-S validation.
**Disposition: redirect from supporting pole/character development to one
matched electronic-model pole-contrast proposal.** At a development geometry,
identify the smallest same-R/same-definition comparison from completed sparse
inputs. If no compatible saved contrast exists, name the missing artifact and
preparation cost; do not automatically open an interface repair or another solver
campaign. One design/review cycle, no numerical job yet. Do not mix old phase-fit
parameters with outgoing poles or consume the excluded equilibrium pole.

Keep2048 roots; no automatic root/channel ladder, extra sector/anchor, geometry,
fit or sweep. All114 route holdouts and equilibrium pole exclusion remain
protected, new50 parked. Omitted spectrum, independent exterior, continuity,
neutral correlation, training/validation nodes, asymptote and physical full-S
embedding stay open. Review on that proposal or next active cycle boundary.

### Goal checkpoint —25 October: resume on matched-model decision, not solver work

**Starting blocker:** the compressed finite-reference open-Pi pair is witnessed;
electronic-model accuracy and a qualified multi-geometry curve remain open.
**Question this cycle:** which compatible same-R/same-outgoing-pole contrast is
already supported by saved sparse inputs, and what exact prerequisite is missing?
Inventory only the completed model/control records and native interface metadata.
Keep the equilibrium pole reserved and all114 route holdout S values unanalyzed.

**Minimum scope/bound:** one design/evidence-review cycle, at most90 minutes of
additional work; metadata/source/hash inspection only, no S/pole evaluation,
propagation, QC/property/inner solve, new50 work or native executable. Reuse the
2048-root references. A viable contrast gets a finite staged proposal with the
missing reader/input gate and measured costs; an incompatible contrast is
declined, not silently repaired. No new numerical job is released by this review.

**Evidence delta at resume:** none yet. **Disposition: continue only this bounded
saved-input inventory and contrast proposal**, then review immediately. Passing
the local pair tests is not electronic convergence. No automatic root/channel/
basis ladder, new geometry, fit or sweep. If an input gap cannot be bounded in
this cycle, name it and stop rather than opening another interface campaign.

### Goal checkpoint —25 October,21:17UTC: saved channel contrast selected; reader gap bounded

**Decision completed:** the [compressed A2 outgoing-pole proposal](../../docs/physics/co-compressed-a2-pole-contrast.md)
selects CAS11 five→nine retained A2 states per spin at1.9bohr,2048 roots. This is
the cheapest compatible development-geometry control found: the original
40/250 and48/298 target/channel models already exist, target-table subsets match
exactly and recorded shared native inputs/outer decks agree. The completed
CAS11→CAS12 scattering pair is at reserved equilibrium; inspected development
CAS12 records are target imports, not a compatible completed scattering pair.
That scoped inventory is not a global absence proof.

**Production link:** distinguish whether the next expensive model allocation
should address retained-channel effects on the **pole target**, rather than
trust the compressed historical phase pass. A valid pole-sensitivity rejection
prioritizes channel-model evidence; a pass returns priority to a justified
basis/active-space contrast while retaining nine, not downgrading to five across
the domain. Neither proves nine-channel convergence or selects a production model.

**Exact gap/cost:** an explicit mapped250-channel compressed-five model view and
narrow existing-runner adapter remain unimplemented/unqualified.12 channels per
sector need property-ID→dense-row lookup; six IDs exceed40. The original strict
reader failure stays intact. The authenticated zero-Buttle/fixed-nuclei K source
path bypasses the problematic target-pointer dereferences and uses stored
thresholds; this does not qualify other native consumers. No engine patch or
new inner eigenpairs are indicated by this trace. Repeating five/nine would
duplicate their measured1.651263/1.827292h calculations, not fill the reader gap.

**Evidence delta:** a concrete saved contrast and exact prerequisite replace a
generic electronic-model proposal. Four scoped interfaces/26 input pins verify
against two fresh public payloads; target/threshold metadata reconstruct exactly.
This is availability/source evidence, **not new physical/model-insensitivity
evidence**. Zero propagation, S construction, pole evaluation or native call.
The preliminary metadata-only path typo exits1 and is preserved in the inventory.

**Disposition: continue only to the mapped-view/adapter implementation and
source/tests freeze**, one active cycle /at most two hours. No job is released
by the [prospective allocation](compressed-a2-pole-contrast-proposal-20261025.json).
Saved-view gate:4 cores /4GiB /60s; eventual conditional five/B1 then B2 at most
79+67=146 maps /8 physical cores /8GiB /3600s combined hard cap, reusing all nine
witnesses. Fresh resource/owner freeze is required. Stop if the view cannot be
qualified in that effort bound; no repair/rebuild/precision campaign. Review on
that freeze or first gate/result, and by27 October if work remains active.
All114 holdouts and equilibrium pole exclusion remain protected; new50 parked,
2048 baseline retained. Model accuracy, omitted spectrum, continuity, neutral/
factory/asymptote/full-S gates remain open. No fit, new geometry or sweep.

### Implementation checkpoint —25 October: explicit view and adapter; saved gate only

The bounded implementation now exposes original property IDs and dense target
rows separately, preserving strict default failure. A private compressed/five
view checks pilot/deck IDs and the48-state subset, then shares the native W/D/CF
assembly with the unchanged-default nine-model loader. Its own2048 eigenpairs
are retained, not sliced from nine. The narrow pole adapter translates the
existing allocation and requires a pinned, measured two-point real guard before
complex five-model use; ODE/root/contour/derivative routines are unchanged.
Synthetic mapping/CF/identity, scope, guard and full-width-band controls precede
source freeze; no actual five pole or S value has been evaluated.

**Disposition: continue only with the separately frozen saved-view gate**
([contract](compressed-five-view-gate-contract-20261025.json)): four physical
cores0–3 /one BLAS thread /4GiB /60s hard cap, zero maps. Authenticate all123
original inputs, pilot/deck/source pins, strict failure, own W/D bytes and
nine identity-view fields. A gate failure stops without a repair or retry;
a pass removes this reader prerequisite but does not release the numerical
contrast. Review immediately afterward and freeze the exact real-S/contrast
owner separately. This remains inside the one-cycle/two-hour implementation
bound; no new QC/properties/inner solve, geometry or fit, and no actual electronic-
model evidence yet.114 holdouts/equilibrium exclusion/new50 protections remain.

### Goal checkpoint —25 October,23:03UTC: mapped view qualified; end interface work

**Starting blocker/question:** can the existing compressed five-model eigenpairs
be consumed explicitly without rewriting native labels or replacing W/D with a
nine-model slice? This prerequisite serves the selected retained-A2 pole-target
contrast, not a production-model selection. Source/tests/contract froze/pushed
at `9deb5a1` before the single saved application.

**Completed evidence:** [both mapped views pass](progress-compressed-five-view-result-20261025.json),
all250 channels /2048 own roots per sector,12 remapped labels each, exact selected
target metadata and nine-model identity-view fields. Strict default failure and
the original false aggregate/six unset cases remain. Independent public-byte
reconstruction verifies six arrays per sector,53 files /39 sources /seven pins /
123 inputs;39 source hashes also match the frozen Git blobs.36 runtime controls
pass (39 focused local tests). Zero S construction, propagation, pole evaluation
or native call. No repeated QC/property/inner solve or nine repropagation.

**Cost/bound:** four physical cores0–3 /one thread /4GiB /60s hard cap;
container2.073691s, gate.147262s, cgroup peak.223804GiB /CPU1.486406s sampled at
final gate save including controls. Original exits0/0/0/0, no timeout/OOM or
scientific retry. A read-only public-index fetch403 before reconstruction is
recorded separately; first actual independent reconstruction exits0. This
deliverable completes within the active cycle's two-hour implementation bound.

**Evidence delta:** the explicit saved-reader prerequisite is removed. **No
physical uncertainty shrank:** five/nine pole sensitivity, electronic accuracy,
omitted spectrum, independent exterior, continuity, neutral correlation, factory
nodes/asymptote and physical full-S embedding remain open. No qualified curve,
production model or restartable geometry pilot exists.

**Disposition: redirect from completed interface work to the exact real-S/contrast
execution freeze.** Stop numerical work at this handoff: the saved gate is not
execution authorization. One next deliverable is a separately reviewed/pushed
owner pinning saved nine witnesses, proposal, inputs, two measured development-S
guards per evaluated sector and original pole/character gates. Do not repeat
this saved gate or open another general solver/interface campaign. The unchanged
prospective ceiling is79 B1 +67 conditional B2 maps /8 physical cores /8GiB /
3600s combined; resolved B1 sensitivity rejection stops before B2. Input, identity,
numerical or cap failures stop inconclusive without repair/retry/successor.
Owner work must fit the original active-cycle effort cap; otherwise stop with
the precise missing owner, not an extended development allocation. Review before
any numerical launch and on its first outcome, or next active-cycle boundary.
All114 holdouts, equilibrium exclusion, parked new50 and2048 baseline remain.

### Execution checkpoint —25 October,23:17UTC: continue to the decision-changing contrast

The user authorizes continued productive use of the available hours, not a stop
after every short supporting check. **Disposition: continue to the single
[frozen numerical contrast](compressed-a2-contrast-execution-20261025.json)**,
after source/tests/exact persistent owner freeze and push. The question remains
whether retained A2 channels materially change the coupled pole target at1.9bohr;
reuse both completed2048-root calculations and all nine witnesses. No actual
five-model S/pole score exists at this checkpoint. No scientific gate is relaxed.

The narrow owner calls the existing real-S and pole machinery: two real guards
persist and pin before complex use, B1 first, B2 only after a qualified B1 pass,
with all original character/identity/counterpart and boundary-band gates. Valid
B1 sensitivity rejection stops B2; numerical/input/identity failures stay
inconclusive, preserving nested original exits. Source W/D/CF, ODE, bordered
search, contour and derivative algorithms are unchanged. Routing/failed-guard/
after-integrity controls precede freeze. Stay within the original active-cycle
two-hour implementation cap (22:35:57UTC start), not a new supporting-work cycle.

Fresh Sadaharu resource inspection23:12UTC finds eight distinct physical cores
0–7 available to pin,114GiB available RAM /64GiB free persistent disk; unrelated
ComfyUI remains untouched. Allocation8 physical cores /one BLAS thread per worker /
8GiB, at most79 B1 +67 B2=146 maps,90s per map /1800s per sector /3600s combined.
No paid provisioning, image rebuild, QC/property/inner solve or root/channel ladder.

**Next review:** at the scientific outcome, independently reconstruct/publish it
and choose the next electronic-model allocation in the same continuation, rather
than ending solely because the compute took minutes. A rejection prioritizes
retained-channel evidence; a pass retains nine and prioritizes a justified basis/
active-space contrast. Any further expensive calculation needs that outcome-led
minimum-scope contract first; free capacity is not an accuracy argument. All114
holdouts, equilibrium exclusion, parked new50 and2048 baseline stay protected.
No production model/qualified curve/potential fit or sweep is released.

### Goal checkpoint —25 October: contrast passes; allocate missing matched CAS12 properties

**Physical question completed:** [compressed five/nine pole-target insensitivity
passes both sectors](qualification-evidence/compressed-a2-contrast/README.md).
Changing retained A2 targets lowers position10.876503meV /full width4.480322meV,
within50meV /104.230719meV budgets and frozen boundary bands. Own real-S,
root/contour/residue/Pi/identity and five-model counterpart gates pass. B1 root
search uses6/1/1 evaluations; B2 accepts the B1 candidate once per setting,
not independent iterates or zero splitting uncertainty. Independent public
reconstruction completes23:24:18UTC, authenticating168 files /42 sources /
24 prerequisites /123 inputs and all scientific gates without propagation.

**Cost/evidence delta:**89/146 maps,216.873008s /1.737797GiB, original exits0/0/0/0,
no timeout/OOM/retry. No QC/property/inner solve or nine repropagation. Retained
channels do not dominate this **compressed finite pole-target** contrast at the
chosen budgets. This is new model-sensitivity evidence, not absolute accuracy;
stretched rejection, full-domain channel convergence, omitted spectrum, independent
exterior, continuity, neutral/factory/asymptote/full-S gates remain open.
Keep nine as working reference, not a five-state production downgrade.

**Production decision: redirect to one matched compressed CAS11/CAS12 model
control**, not another channel/root/precision ladder. Both native target sets
already exist at1.9bohr/cc-pVDZ. The witnessed CAS11/five pair supplies only the
local firstfive-per-sector comparator; CAS12/nine and full-domain convergence
remain unqualified. This changes optimized orbital frames/imported target
representations too, so no pure active-orbital-count claim or equilibrium pole use.

**Exact prerequisite found:** original CAS12 DENPROP retained global lowest40,
not five per sector: A1 singlet/triplet6/4, B1=5/7, B2=5/7, A2=3/3. Missing A2
roots4/5 and tripletA1root5 already exist among the64 independently qualified
native vectors. Reuse all vectors, target integrals, CONGEN, checkpoint/Molden;
only a **different, matched property inventory** needs computing. No QC, target
diagonalization/reimport, new integrals or scattering solve is released now.

**Disposition: continue to the single multi-hour property export** after its
source/tests/exact owner freeze/push, under the
[property contract](compressed-cas12-matched-properties-contract-20261025.json).
Change only DENPROP NTGTL to5 in all eight sectors, preserving original set order.
30 saved-input hashes and the all64-root proof authenticate reuse. Original
DENPROP took4.43h; conditional forecast3–6h, native8h /owner9h hard caps,
64GiB /8 physical cores4–11 /one library thread, minimum80GiB available RAM /
10GiB disk. Whole old target-pipeline41.12GiB peak is a conservative anchor,
not isolated DENPROP cost. No paid provisioning or speculative speedup.
Implementation/diagnosis is bounded to60minutes after this model decision;
this is a new physical input allocation, not extending the completed view work.

**Next while it runs:** prepare the costed exact-input full-continuum integral/
CONGEN contract and one2048-root CAS12 B1 model contrast, reusing these properties.
Actual dimensions and input/target identity must freeze before an inner solve;
no dense build or automatic B2/root/basis ladder. Property failure stops without
repair/retry/successor; property success alone releases no Hamiltonian or pole job.
Review at that contract or property outcome, and next active-cycle boundary.
All114 holdout S matrices/equilibrium exclusion/new50/2048 protections remain.
No production model, qualified curve, potential fit or sweep is released.

The original pre-execution23:20 caption was prospective/approximate; its measured
Git freeze is23:17:52UTC before launch23:18:30UTC. Current caption is corrected;
the frozen capsule and original chronology remain preserved in the result record.

**Launch update23:40UTC:** `bfa3a88` freezes/pushes the exact export source/tests/
contract/owner before the durable service starts23:40:13UTC. One launch-time
observation finds native DENPROP active, three runtime controls passing,30 inputs /
23 sources verified;69 combined focused local tests preceded freeze. Property
verdict/terminal exits stay pending. Initial100.59% CPU /.133GiB is not a final
peak or claimed eight-core speedup. Original data/workloads remain untouched.
This uses the available hours for a necessary different model input, not a
repeat of successful QC, target eigenpairs or original property selection.

### Goal checkpoint —26 October,04:16UTC: matched properties verified; prepare the model input

**Decision served:** remove the missing matched-retention CAS12 property input
before the selected compressed CAS11/CAS12 model control. The
[single different DENPROP export passes](qualification-evidence/compressed-cas12-properties/README.md):
40 states/five per sector,725 unique multipoles, ground dipole agreement2.075e-9au;
native CIDATA energies agree within4.768e-10Ha. This qualifies input identity/
inventory, **not electronic convergence, a new pole or a production model**.
Nine remains the working reference; five is only a local matched comparator.

The source-frozen job finishes04:06UTC, original native/scientific/command/
container/owner/service exits0, no timeout/OOM/scientific retry. Native4.4296h,
container15948.199065s /15979.841532CPU seconds /1.903126GiB kernel/sample peak.
Eight physical cores were allocated but execution is roughly one-core; supervisor
23.9M/7.126CPU seconds are not scientific resources. Saved30 inputs,12 copies
and23 source hashes remain unchanged; no QC/target solve/reimport was repeated.
Public reconstruction reads native CIDATA and properties independently, nine
payloads /63 files. Original archive-reader exit1 and symlink-byte recollection
are preserved beside corrected reader0; no scientific calculation was rerun.

**Disposition: continue to the one bounded native preparation**, under the
[exact contract](compressed-cas12-preparation-contract-20261026.json), after
source/tests/owner freeze/push. One full-continuum integral call +two CONGEN
calls,4 physical cores0–3 /16GiB /one library thread;600s conditional target /
3600s combined owner cap.36 input pins include the completed property bytes/exits;
geometry-free CAS12 algebra is reused, but compressed Molden/checkpoint/own target
vectors/properties are the only physical inputs. All14 target/core orbitals
must survive; derive actual continua/L2/contracted counts, not force equilibrium
86352/990990 planning anchors. No target/QC/property repeat or inner solve.

**Implementation and stop bounds:** source/test/exact owner freeze≤60minutes
after04:06UTC property outcome. Reuse native stage/inventory/completion helpers,
with explicit CAS12 inventory rather than changing the default CAS11 gate.
Failed prerequisite/native/count/identity/cap stops inconclusive, preserving
exits/partial data, with no repair/retry/successor. A pass needs independent
saved-output authentication before the separate2048-root B1 allocation review.
That very large job should use12 physical cores; MPI-layout changes, actual
dimensions, RAM/disk and measured source costs must be explicit before release.
No B2 or root/channel/basis ladder is automatic.114 holdouts, equilibrium exclusion,
2048 baseline and parked new50 remain protected; no fit, curve or sweep released.

**Preparation launch:** `6bad622` freezes/pushes source/tests/exact owner before
04:21:57UTC; six runtime controls pass at the04:23 observation.76 focused local
tests passed before freeze. This stays within the60minute implementation cap.
The [B1 cost-only proposal](compressed-cas12-inner-cost-proposal-20261026.json)
records prior same-dimension equilibrium cost5.95h at8ranks and whole-run21.97GiB
peak; a12-rank forecast is conditional, not linear scaling evidence. It releases
no inner job before actual preparation dimensions and independent qualification.

### Goal checkpoint —26 October,04:31UTC: preparation qualifies; allocate one B1 model solve

[Independent preparation reconstruction passes](qualification-evidence/compressed-cas12-preparation/README.md):
all14 target/core MOs, continuum51/34/34/22,990990 raw CSFs and actual86352
contracted states in each sector. Three native exits0; original scientific/
command/container/owner/service exits0, no timeout/OOM/retry. Owner171.42s /
10.700GiB kernel peak /5.617GiB allocated disk.36 inputs/copies and34 sources
authenticate; public replay verifies six payloads/83 files without a native call.
This removes the exact same-frame input/dimension blocker, not CIRMAT phase,
scattering-interface, pole or electronic accuracy gates. No QC/target/property
calculation repeated. Nine remains working reference; five is a local comparator.

**Disposition: continue to one B1-only compressed CAS12/2048 inner/export job**,
after the [exact contract/source/tests/owner freeze](compressed-cas12-inner-contract-20261026.json)
and push. This is the smallest missing expensive calculation that supplies own
W/D/CF/K for the selected model contrast. Reuse all completed preparation and
target/property bytes; use the already qualified242-artifact audited engine,
not a rebuild. Preserve existing residual/CI-boundary/native metadata gates and
strict40-target/250-channel assembly. B1 alone leaves Pi pairing/model/pole
verdicts unset; a success does not launch B2 or a pole search automatically.

**Cost and limits:** same measured dimension as the saved equilibrium CAS12
baseline; B1 previously took21416.594209s (5.95h) at8ranks. Whole-run kernel
peak21.97GiB/disk5.10GiB are anchors, not new compressed measurements.
Conditional3–8h forecast, native12h /owner13h caps,64GiB /8GiB shared backing /
**12 physical cores0–11** /one library thread per rank.12-rank MPI layout differs
from the qualified8-rank run; no linear speedup or exact8-rank replay claim.
Recheck80GiB available RAM /16GiB disk and no other scientific workload before
launch. Preserve residual1e-7Ha and CI/boundary1e-10Ha gates unchanged.

**Implementation bound:** at most60minutes from04:31UTC review to exact narrow
adapter/owner freeze; existing default two-sector pilot stays unchanged. B1-only
is explicit, failed prerequisites block native work, and no B2-only route is
added. Freeze the actual solver deck and own property map before native use,
then require copied-input/deck integrity and strict own export. No general
mapping/solver framework or supporting precision cycle. Failure stops inconclusive
without repair/retry/successor; preserve every original nested native exit.

**Next at outcome:** independently qualify the own2048-root export, then freeze
the two real-S guards and outgoing-pole/character/model comparison separately.
No scalar/factory export, fit or sweep follows native success.114 holdouts,
equilibrium pole exclusion,2048 baseline and parked new50 remain protected.

**Launch04:42UTC:** source/tests/exact owner freeze/push `6ae5b79`04:41:34UTC
precedes service start04:42:17UTC. At04:43, six runtime controls/50 input pins/
44 source entries/242 installed artifacts and runtime numerical/telemetry gate
pass. The12-rank native matrix is86352 with7196 local rows; target reads are in
the compressed1.9bohr frame. Startup1199.39% CPU /13.39GiB is not a final peak or
scaling result.114 focused local tests preceded freeze. No final native/export/
Pi/model/pole verdict yet; no successor is released by this bounded launch.

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
| Fitting-information route | Compressed finite open-Pi pair passes; saved same-R CAS11 A2 five/nine contrast identified, no pole-insensitivity or electronic-accuracy result | Qualify/freeze one explicit250-channel mapped view and narrow adapter; no new inner solve; eventual pole contrast triages the next model allocation | Two-hour implementation cap; future at most146 maps /8 cores /8GiB, no job released;114 holdouts/equilibrium exclusion protected |
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
