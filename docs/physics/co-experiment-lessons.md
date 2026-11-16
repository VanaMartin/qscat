# Lessons from the CO experiment campaign

These are reusable planning lessons, not a run ledger or a new qualification of
the calculations. Measurements, detailed review and original failures remain in
the [external scientific review](https://data.qscat.org/co-completion-review-20261103T1610/review.6b36f1e638e5.md)
and its immutable archive. Apply the [experiment-planning gate](../experiment-planning.md).

| Mistake or unsupported inference | Why it failed | Future rule |
|---|---|---|
| Treating successful replay as independent electronic validation | Reader and producer can share the same physical approximation or convention | Name the oracle's independence and shared errors before launch |
| Chasing an affinity setting before establishing expected method error | Finite-basis, valence CCSD(T) omits contributions present in a precision threshold | Build the component/convention error budget from primary references first |
| Taking a negligible orbital-convention effect as a complete correctness proof | It excludes one proposed mechanism, not every implementation or model error | State exactly which hypothesis a control rejects |
| Interpreting failed diffuse SCF as a diffuse-affinity result | The anion never reached a qualified correlated energy | Preserve the failure; keep the affinity difference unset |
| Extending a stable restricted well reference into open-shell dissociation | The reference and correlation balance change with physical regime | Test a spin- and fragment-consistent bridge at a sentinel before a dense curve |
| Equating inactive orbitals with frozen-core correlation | Those settings have implementation-specific meanings | Verify the actual runtime API and molecule/fragment core convention before release |
| Treating a spin label or inactive-orbital count as state identity | Near-degenerate iterations can mix spins or rotate nominal core orbitals into a different subspace | Measure final spin and orbital-space witnesses; do not accept a failed iterate because its labels look correct |
| Inferring that an engine is unavailable because another interpreter cannot import it | Packaged engines may have their own launcher and Python environment; existing controls may already use them | Inspect the supported entrypoint and reuse existing capability evidence before proposing a rebuild |
| Comparing different geometries, zeros or width estimators | Those changes can dominate the discrepancy being scored | Match definitions or explicitly quantify the conversion and uncertainty |
| Reporting tight numerical tolerances as molecular accuracy | Basis, active-space and ensemble sensitivity can be much larger | Carry separate numerical and physical error budgets into the observable |
| Selecting by a larger model or proximity to a familiar number | Error cancellation and state changes can mimic improvement | Use controlled contrasts, independent references and frozen predictions |
| Letting support work become the project | More controls and records need not close a physical decision | Give each support task an effort cap and a goal-linked exit |
| Requesting a residual smaller than the solver's correction cutoff permits | A cutoff on squared residuals can discard the next direction before the requested norm is reached | Inspect the stopping rules and cutoff units; test a near-eigenvector analytically without weakening acceptance gates |
| Treating failure of a restricted initializer as failure of a multireference method | An orbital initializer need not be the final stationary reference | When supported, project accepted core/active spaces into an orthonormal frame and qualify the final CASSCF state independently |
| Treating the end of one bounded attempt as exhaustion of the scientific task | A stopped release preserves its limits but does not rule out a diagnosed alternative | Retain the failure, reassess the goal and freeze a justified next attempt; stop the task only at resolution or a substantiated blocker |
| Treating a failed two-basis test as a permanent ban on extrapolation | Sensitivity rejection is not a statement that systematic convergence is impossible | Reopen only with a justified, separately frozen convergence design |
| Shrinking tracked experiment records instead of removing the data/code boundary violation | Small run ledgers and dated reports are still experiment data | Store immutable experiment objects externally; Git holds reusable code/guidance and minimal pointers |
| Making every later physical claim an entrance gate to a first local dataset | Global dissociation, exact affinities and full-S embedding answer different questions from finite E0/Er/Gamma rows | Put row acceptance, curve identity, fit quality and physical/global qualification in separate stages; retain failed physical diagnostics without discarding valid finite data |
| Waiting for an authoritative potential before fitting our own energies | The external curve was mistaken for the output we intend to compute | Use our accepted neutral energies as inputs; external curves are optional benchmarks |
| Letting a global fitter API dictate missing physics | Molecular-zero tables are not atom-zero EMO inputs; mandatory asymptotic arguments can tempt invented constants | Handle local offsets and deferred global constraints explicitly; store valid local data before adapting the fitter |

## What to retain

Preserving original failures, verifying provenance, testing invariants and
protecting held-outs were necessary. They were not substitutes for independent
physical validation. Continue those controls at the precision and effort needed
for the final scientific decision, rather than building a second campaign around
every supporting uncertainty.

## How this changes the route

Produce the two data outputs in the [CO goal table](../../projects/ukrmol_co/RESEARCH_PLAN.md)
using accepted records before demanding their global completion. An unavailable
matched implementation is a limitation, not a mandatory framework project or a
reason to relabel a different method as a matched oracle. A global extension
then needs a coherent molecular bridge: well overlap, stretched-state identity
and approach to correct fragments. Atomic refinement is useful only when its
uncertainty blocks that decision; an experimentally exact affinity is not the
end product.

Do not fit a global curve before branch identity, common zeros and bound joining
are demonstrated. A fitted pole does not supply channel/background embedding or
nuclear dynamics. Close the original physical goal only with its independently
validated observable; label a methods-only rescope or a blocked result honestly.

## Standing CO protections

Keep the 114 S-matrix holdouts and specified equilibrium pole exclusion untouched
during selection, the qualified 2048-root baseline, nine working/new50 parked,
and accepted finite calculations. Small Sadaharu pilots prefer four physical
cores with explicitly bounded RAM/time and 600-second active checks plus immediate
terminal observation. Do not provision paid resources or launch a sweep by
default. Historical execution flags are not renewed launch authority.
