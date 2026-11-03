# CO research plan: close physical references before more scattering

## Active direction — 3 November 2026

The user approved this reordering after reviewing the stalled campaign. This is
the assistant's restartable work queue. The current goal checkpoints and execution
contracts remain the authority for actual allocations; historical next steps do
not override the order below. Update the status table when a physical deliverable
lands, not for every passing test or packaging change.

**End product:** a neutral curve and lowest-doublet-Pi resonance/bound-anion
continuation, on a consistent atomic energy zero, constraining a QSCAT 2-D model
with excluded-geometry validation. A channel-specific vibrationally assisted DA
claim requires a separate nuclear-observable check. This does not promise total
thermal CO DA near10eV, where other resonance families contribute.

**Current verdict:** useful local data and calculated atomic endpoints, not a
dissociation-qualified potential.
The correlated neutral segment is1.9–2.5bohr only; molecular3/4bohr restricted
references failed stability before CCSD(T). The first four-atom aug-QZ pilot now
supplies independently reconstructed finite-method sums; absolute EA accuracy,
a consistent resonance curve,bound continuation and production model remain
unqualified. Qualified stretched DZ/TZ poles differ103.370meV/44.172meV in position/
full width. Aug-TZ preparation is preserved but its inner/pole work is paused.

## Rules preventing another open-ended supporting campaign

1. Every new job must answer a named physical decision below. Reuse accepted
   electronic data and2048-root controls; no automatic basis/root/channel ladder.
2. Keep numerical error distinct from model error. Preserve0.1meV position/full
   width and1e-6rad same-model numerical criteria, and existing50meV/5%-with-1meV-
   floor/0.05rad physical criteria in their specified metrics. Earlier rejections
   remain rejections. The discussed0.1eV/10–20% exploratory tier is **not adopted**
   by this plan; any lower-quality release needs an explicit new scope decision.
3. No model is selected because it has a larger basis or one convenient literature
   match. Published phase-fit widths are not exact outgoing-pole benchmarks.
4. Small jobs: four physical Sadaharu cores; larger/very large jobs: eight/twelve
   only when justified. No paid provisioning, image rebuild or unrelated workload
   changes by default. Preserve immutable inputs, exits, failures and measured cost.
5. Freeze and push exact inputs/source/tests/owner before a quantum job. A stage
   failure stops that release; do not turn it into an unbounded repair project.
6. Review at each deliverable and at most48h of active work. One working-cycle
   implementation cap per stage; exceeding its stated bound requires reassessment.
7. Keep114 S holdouts untouched, reserve the equilibrium outgoing pole for the
   specified validation, retain nine as working reference/new50 parked. No fit,
   sweep or automatic B2 rescue is released by this plan.
8. User requests automatic continuation and status checks every5–15minutes.
   During active runs,target600s checks plus immediate terminal observation;
   record actual gaps. Continue within an approved bounded stage and independently
   qualify its result. A new stage needs its goal-linked tested/frozen release,
   not another preference question. Stop on scientific/resource/scope blockers.
   This is active-run monitoring,not a perpetual scheduler while no job runs.

## Ordered work packages and decision gates

| Order | Deliverable | Completion / decision gate | Bound and stop |
|---|---|---|---|
| 0 | Reference/energy-zero specification and actual missing-data ledger | All three atomic limits and native-pole-to-neutral map explicit; unknown totals remain null, not literature placeholders | Saved source only; one planning session |
| 1 | Four-atom C(3P), O(3P), O-(2P), C-(4S) pilot | Correct charge/spin/occupation/core conventions; finite converged totals; internally consistent sums/EAs; external comparison labeled by convention | First one augmented basis matched to existing neutral data; four physical cores, no geometry/scattering job; exact costed contract before launch |
| 2 | One matched literature pass and existing-neutral shape assessment | Actual source data coverage/zeros/full-width convention authenticated; compare only overlapping R and estimator; use existing52 points for Re/curvature/low-v assessment with domain checks | At most one60minute acquisition/normalization pass; report inaccessible/missing data, do not replace it by high-precision figure digitization |
| 3 | Neutral dissociation/correlation decision | Select one justified stable spin-consistent treatment, or explicitly declared literature/hybrid reference; demonstrate overlap with accepted well data and a route to fragment consistency | One bounded method-selection cycle; a small diagnostic set, not a new dense curve; numerical release/caps specified separately |
| 4 | Single electronic-model decision | Use literature, atom/neutral evidence and existing sensitivity to choose a prescription with uncertainty, reject it, or freeze a limited model-conditional claim | Finish at most the already prepared aug-TZ contrast if it can change this choice; no mandatory augmented run or unlimited model ladder |
| 5 | Minimal compatible resonance/bound curve | At least four actual training geometries plus separate validation; assignment continuous; cover binding transition/escape region for DA; correct bound continuation and asymptotic residual | Choose nodes after prior gates; reuse compatible records; per-node costs and stopping conditions before allocation |
| 6 | Fit once, freeze, predict held-out observables | Excluded geometry passes the declared reduced target; independently specified channel/background treatment for scattering; nuclear observable before a DA claim | No fitting training poles back to themselves as validation; rejection/inconclusive if missing physics dominates |

Packages1 and2 are independent and may be interleaved without delaying the atomic
pilot behind a publisher-access problem. A second atomic basis is conditional on
the first actual result: estimate the endpoint/EA basis effect only if needed;
do not launch a QZ/5Z/CBS ladder in advance. The four-atom pilot alone does **not**
certify molecular dissociation or absolute EA accuracy.

## Energy and state contract

Use Hartree and bohr internally. Required free-fragment sums are

\[
A_0=E_C+E_O,\qquad A_O=E_C+E_{O^-},\qquad A_C=E_{C^-}+E_O.
\]

Define `V0(R)=E_CO(R)-A0`. On that zero, the anion limits are
`AO-A0=-EA_O` and `AC-A0=-EA_C`. The scattering pole is relative to its own
native neutral threshold: `z=E_res-i*Gamma/2`. The proposed effective map is
`Vminus=V0_ref+Re(z)`, `Gamma=-2*Im(z)`. Substituting an independently correlated
neutral is a declared hybrid reference construction, not a calculated correlated
absolute anion energy. Validate its joining/assignment assumptions.

Specify nonrelativistic spin-free term energies for calculated atoms, not
fine-structure ground-J values. External EA/D0 comparisons must label isotope,
term/ground-J convention and any unapplied corrections. Keep the electronic
minimum, fixedR=2.1323 reference, and initial rovibrational zero distinct;
`D0=De-ZPE`. Ground-fragment DA thresholds are `D0-EA`, not the low-energy
fixed-nuclei resonance position.

Match the existing valence correlation by freezing **one1s spatial orbital per
atom**, retaining its total core energy; molecular CO freezes two occupied
spatial orbitals. Use spherical functions and augmented bases for both charge
states. Open-shell atomic reference/correlation choices must be named explicitly;
do not label an automatically converted ROHF/UHF implementation as restricted
closed-shell CCSD(T). Multiplicity alone is not an orbital-term check.

Free-atom sums are necessary endpoints, not finite-R matching proof. Later check
spin-coupled neutral dissociation, charge localization, size-consistency/BSSE and
residual long-range interactions. For C+O-, the isotropic polarization term uses
neutral **carbon's** polarizability; open-shell anisotropic interactions may also
matter. Do not invent a CO matching radius or impose a spherical tail without
testing its applicability.

## Aug-TZ branch is parked, not rejected

Own properties/preparation qualify and must not be repeated. The4.94144eV
excited thresholds block the unchanged99-point/six-open-channel export view.
No new inner allocation or general variable-channel framework is on this queue.
If package4 genuinely needs this one finite comparison, first prescribe an
explicit subthreshold diagnostic reference, retain **all** closed-channel
couplings and qualify the modified domain. No silent grid truncation or new
claim about threshold-crossing/high-energy DA follows.

## Restart status

| Item | Status | Next action |
|---|---|---|
| User direction and work order | Approved; recorded here | Follow this queue instead of the old inner-export successor |
| Reference specification | Complete; [protocol](../../docs/physics/co-atomic-energy-references.md) and [unset ledger](atomic-asymptote-inputs-20261103.json) | Freeze source/tests/owner before quantum execution |
| Open-shell runtime | [Source-only inspection passes](atomic-runtime-preflight-20261103.json);PySCF2.11 conversion supported,external ROHF stability unavailable | Use supported internal checks;external status must remain unavailable,not passed |
| Actual atomic sums | [Four aug-QZ finite-method endpoints independently reconstructed](qualification-evidence/atomic-qz-endpoints/README.md),17.914s/1.943GiB | One bounded basis/reference-convention decision;no automatic second basis |
| Neutral QZ/5Z well data | Reuse; no dense rerun | Read existing records for shape/spectroscopic diagnostics |
| Dora2020 numerical supplement | One acquisition pass reached human-verification page; [receipt](literature-reference-acquisition-20261103.json) | No bypass/digitized precision benchmark; proceed with independent atomic work |
| Neutral stretched reference | Restricted approach stopped | Select spin-consistent correlation/dissociation strategy, not an SCF retry |
| Augmented scattering | Paused; no inner/pole result | Reopen only after the physical-model decision needs it |
| Potential/DA/production release | Not ready | Complete reference, curve and independent validation gates |

**Completed continuation:** four atoms execute under their separate frozen release;
all original exits0,public raw-log/source/resource reconstruction passes. Computed
EAs O1.403073/C1.235049eV have correct ordering but differ from rounded external
context by approximately58/27meV; reference-convention/correlation/basis accuracy
remains unresolved. Completion occurs before the first600s periodic check; the
immediate terminal observation is17.907s after start.

**Immediate next deliverable:** at most30minutes saved evidence/source review to
decide whether one matched5Z endpoint control is needed to isolate basis sensitivity
and connect the existing5Z molecular zero. Reconcile spin-free/ground-J conventions
before scoring precise EA errors. No second basis or molecular/scattering job
follows automatically; any selected quantum control needs a new exact freeze.
No active campaign job/monitor remains. Continue the approved research without
repeated preference questions,but do not waive scientific stops.

## Evidence and provenance

- [Neutral near-equilibrium evidence](qualification-evidence/neutral-near-equilibrium/README.md).
- [Finite pole sensitivity](qualification-evidence/stretched-pole-counterpart/README.md).
- [Own augmented preparation](qualification-evidence/stretched-augtz-preparation/README.md)
  and [threshold blocker](qualification-evidence/stretched-augtz-inner-preflight/README.md).
- [Effective target and held-out claim boundaries](../../docs/physics/co-effective-resonance-target.md).
- [Laporta2012 reference](../../reference/literature/laporta-2012-psst21-045005.md),
  preprint pp.4–6; [Dora2016](../../reference/literature/dora-2016-epjd70-197.md),
  pp.4–7; [Dora2020](../../reference/literature/dora-2020-jpb53-195202.md), pp.3–8.

Planning starts from current branch source at`74ae008`. Main retrieval still
errors on incompatible shared-service dependencies; last-known indexed commit
`423d5770d3de337571f51e62702ae919bad97e8c`, integrity/upstream currency unverified.
Current source and scoped local branch differences, not stale retrieval, govern
this queue. Prior research scratch notes are not a load-bearing dependency.
