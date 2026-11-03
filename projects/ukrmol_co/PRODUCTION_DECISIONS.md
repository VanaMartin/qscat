# CO production decisions

## Current user-selected affinity correctness checkpoint

The user asks to solve the atomic-affinity discrepancy and identify an incorrect
setting if present. The saved5Z orbital-convention control completed: EA changes
only +0.029095 meV,below its0.1meV criterion;density/core and CCSD energy invariant.
All original exits0,76.515963s owner,11.817177GiB kernel peak,zero SCF/CCSD
iterations. Thus this convention is not the main cause. Occupation/core/subtraction
audit finds no configuration bug. External atomic benchmarks support slow basis
convergence and missing post-CCSD(T) correlation,not an empirical correction.

**Stop the attempted diffuse-SCF diagnostic:** the fixed-QZ105-AO extension
passes integral-only parent-space controls;neutral O completes,but O− ROHF fails
100cycles before correlation. Worker exits0/1,outer exits2;no affinity or diffuse
shift. No gate relaxation/retry,timeout or OOM. Owner15.860873s,CPU44.236479s,
kernel peak5,896,859,648bytes. Original partial evidence and failure remain intact.
This does not establish negligible diffuseness or absence of a stationary ROHF
state;it does show that adding diffuse functions is not a reliable automatic fix.

**Redirect:** at most60minutes saved-source/method selection,zero quantum calls,
for one state/core-consistent correlation control with spin/residual gates and a
costed pilot. Small-basis FCI is a candidate,not an automatic successor or a5Z
correction. Absolute affinity remains unresolved;no6Z/CBS,empirical adjustment,
molecular/scattering or holdout work. See the
[completed diagnosis and limits](../../docs/physics/oxygen-affinity-diagnosis.md).

## Current disposition

**Reference-first; production model, potential fit and sweep remain stopped.**
The neutral well and finite atomic/pole calculations are reusable. Their
physical insensitivity and dissociation qualification do not follow from passing
execution or independent replay. The permanent
[results and rejected-routes note](../../docs/physics/co-results-and-rejected-routes.md)
is the compact scientific record; per-study evidence retains original failures.

The matched aug-QZ→5Z control rejects atomic basis insensitivity: O EA changes
23.533750 meV against 10 meV; matching-basis atom-to-fixed-molecular-reference
offsets change 52.916179/29.382428/46.965473 meV against 20 meV. These are not
De/D0 or initial-v attachment thresholds. Stop convergence/production inference;
do not launch 6Z or extrapolate automatically. Both finite datasets stay qualified.

## Decision queue

| Question | Settled status | What can change the decision |
|---|---|---|
| Atomic/neutral reference | Finite endpoints exist; basis sensitivity remains; restricted neutral continuation fails at 3/4 bohr | One stable spin-consistent correlation/reference treatment with well-overlap and fragment checks |
| Electronic model | Compressed CAS11/CAS12 and stretched DZ/TZ pole insensitivity reject | A justified model/reference choice with explicit residual uncertainty, not simply a larger basis |
| Channel/ensemble prescription | Five A2 insufficient at stretched anchor; nine working reference; new50 parked | Named observable-sensitive comparison, not channel/root escalation by default |
| Resonance extraction and identity | Local finite poles exist; phase-fit widths need physical constraints; no global continuous branch | Compatible model, thresholds, pole/residue assignment and bound continuation over a declared domain |
| Aug-TZ scattering | Properties/preparation qualify; unchanged export blocked at open-channel threshold crossing | Reopen only if needed for the physical-model choice, with explicit domain/export design |
| Production | No fit or sweep; 114 holdouts untouched | All [readiness gates](SWEEP_READINESS.md), then one frozen restartable pilot |

## Unchanged error budgets and protections

- Physical comparisons: **50 meV position**, **5% full width with a 1-meV floor**,
  **0.05 rad in the specified phase metric**. Failed extraction cannot certify
  a width comparison. A passing finite pair is not an absolute electronic error.
- Same-model numerical pole controls: **0.1 meV** position/full width and
  **1e-6 rad**. Other studies retain their distinct frozen K/S/flux/representation
  budgets; do not conflate them with physical tolerances.
- The discussed lower-accuracy exploratory tier was not adopted.
- Keep 114 S-matrix holdouts, the specified equilibrium outgoing-pole exclusion,
  nine working/new50 parked and the qualified 2048-root baseline.
- Preserve immutable input/source/image pins, original failed exits and actual
  resource measurements. Packaging failure, numerical rejection and physical
  model sensitivity are distinct outcomes.

## Goal-convergence checkpoint — every one or two days of campaign work

Review at most every 48 hours of active work, on resumption, before opening a
new experimental branch and after a decision-changing result. Record:

1. The physical/potential decision and observable/domain/accuracy it serves.
2. Reused data and the smallest check that could change that decision.
3. Compute/resource caps, implementation-effort cap and first-failure stop.
4. Actual evidence delta, unchanged uncertainty and measured cost.
5. **Continue / redirect / stop**, with one bounded next deliverable.

Passing tests, payload counts and interface fixes are supporting evidence, not
sweep progress. Without a decision-relevant delta, reassess rather than refine.
Historical launch entries and contract flags do not release new execution.

## Current checkpoint and next deliverable

The atomic-basis control reused QZ and all molecular data, calculating only
four 5Z atoms with identical calculation-source bytes and state/core/method
conventions. All original calculation exits zero; owner 114.673835 s,
CPU 302.461764 s, kernel peak 11.883816 GiB on four physical cores within 16 GiB.
The separate zero-quantum AST-control stop remains failed; exact-source-byte
equality changes only the control, not physics or tolerances. Independent public
reading reconstructs 48 raw files, 21 executed Git sources, 573 samples and QZ.
The result changes the physical decision: endpoint insensitivity cannot be used
to justify production alignment.

**Redirect:** at most 60 minutes saved evidence/source method selection, zero
quantum calls. Choose one stable neutral dissociation/correlation treatment or
declared hybrid/literature reference with overlap/fragment checks and a costed
small diagnostic set, or record a blocker. Review at that decision and before
any new quantum release. No active job/monitor or automatic basis/scattering
successor. The detailed remaining order is in [the research queue](RESEARCH_PLAN.md).

Earlier checkpoint chronology is optional audit material in the
[public history archive](qualification-evidence/campaign-history/README.md), not
a dependency for understanding these decisions.
