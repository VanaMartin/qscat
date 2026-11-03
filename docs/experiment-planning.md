# Goal-directed experiment planning

Use this gate before a new calculation, a retry, a model extension or a change
of physical regime. The purpose is a scientific decision, not a larger dataset.
Reusable code and this guidance belong in Git. Filled plans and experiment
records belong in immutable external storage, including unsuccessful attempts.

## 1. Start from the end product

State the observable, physical domain and accuracy needed to finish the project.
Name the unresolved decision on the shortest dependency path to that product.
Explain why the proposed result can change that decision. If all plausible
outcomes merely say “calculate more,” redesign the experiment before allocating.

Distinguish three endings: a validated scoped prediction; a bounded scientific
rejection of the candidate; and BLOCKED/inconclusive because evidence or runtime
support is missing. A narrower methods study requires an explicit scope decision;
it is not completion of the original physical objective. Cleanup/merge is separate.

## 2. Use the literature before diagnosing a bug

Establish the expected error of the chosen approximation before interpreting a
discrepancy. Use primary sources with edition/page/table or equation locators.
Record inaccessible data as unavailable; do not substitute a plotted estimate
for a precision numerical benchmark or bypass an access challenge.

Match geometry, state, charge/spin, Hamiltonian, core treatment, basis, reference,
energy zero, units and estimator. Separate electronic/rovibrational thresholds,
spin-free/lowest-J energies, isotope/field conventions and experimental/theoretical
values. A literature correction is not a correction calculated for our state.
State whether the task tests implementation agreement or physical accuracy.

## 3. Declare what the evidence can establish

| Evidence | Establishes | Does not establish alone |
|---|---|---|
| Hashes, exits and saved-log replay | Execution provenance and arithmetic consistency | Correct electronic implementation or physical accuracy |
| Invariants and analytic controls | The tested identity/limiting case | Completeness of the physical model |
| Independent implementation with matched inputs | Agreement of implementations in that regime | Accuracy of a shared basis/approximation |
| Controlled refinement | Sensitivity/convergence of the specified family | A bound on omitted physics without further evidence |
| Independent physical benchmark | Accuracy for the matched observable/domain | Validity outside that domain |
| Withheld prediction | Generalization within the frozen validation design | Unmeasured observables or regimes |

Specify which errors are shared with the oracle. Partition numerical, model,
reference and joining/embedding uncertainties in units of the final observable;
do not add correlated uncertainties in quadrature without justification.
Choose numerical precision sufficient for the physical decision, not arbitrarily
tighter controls while the model error dominates.

## 4. Choose the smallest discriminating calculation

Reuse qualified inputs/results by digest. Start with a sentinel at the known
failure or regime boundary, then add only points needed to resolve its decision.
Changing basis, state average, root identity and core simultaneously does not
isolate a basis effect. Label unavoidable changes as a new protocol.

Check implementation support **before** release, including frozen-core semantics,
spin adaptation, state following and the actual solver/orbital convention. An
inactive CASSCF core does not automatically mean frozen-core dynamic correlation.
Plan for open-shell fragments and difficult SCF references; a convergence failure
does not prove that a physical state or a stationary solution is absent.

Systematic basis/CBS studies are allowed when they answer the named decision:
justify the refinement family, normally use at least three resolutions, separate
energy components, retain state identity and test the residual/extrapolation error.
A failed two-basis insensitivity test neither supplies a CBS bound nor permanently
prohibits a properly designed convergence study. Reopening a stopped route needs
a new frozen rationale; no automatic basis/root/channel/solver ladder.

## 5. Freeze one decision card before launch

Fill this template **outside the checkout**. A small repository catalogue may
reference its immutable URL and SHA256; do not commit the filled card.

```text
ID / parent evidence digests:
End product / observable / domain / accuracy / DONE condition:
Blocking decision and why this experiment is on its dependency path:
Competing explanations and expected result under each:
Decision on PASS / REJECT / INCONCLUSIVE (one next action each):
Reused data and smallest new calculation; what remains deliberately unmeasured:
Hamiltonian / state / geometry / core / basis / reference / zero / units / estimator:
Independent comparator, primary-source locator, shared errors and access limits:
Numerical / model / reference / joining uncertainty budgets and gate rationale:
State identity, convergence/residual and physical acceptance tests:
Runtime capability preflight and synthetic/analytic controls:
Implementation-effort cap / conditional runtime forecast / hard compute caps:
CPU topology / memory / disk / scratch / owner / observation and terminal checks:
First-failure stop; specifically authorized retry or extension, if any:
Protected validation data; parameters frozen before access:
Exact executed source / inputs / environment pins and external record destination:
```

The release is the frozen card plus tested source and a bounded owner—not this
template, an old contract flag or a historical “next” entry. Record forecasts as
conditional, not measured performance. Capability inspection is not permission
to run a scientific calculation.

## 6. Execute, assess and stop

Use persistent external scratch, a fresh attempt ID and capped resources. Preserve
original exits, failures, partial products and measured costs. Observe active
runs at the agreed cadence and immediately on completion; no perpetual idle
scheduler. Do not silently loosen gates, switch solvers or repeat successful peers
after failure. A changed attempt needs its own bounded release and a link to the
original failure. A numerical failure provides no missing physical observable.

Assess the frozen observable independently where possible. Record evidence delta,
remaining uncertainty and continue/redirect/stop. After one or two supporting
stages without a physical decision, reassess the route instead of adding controls.
Benchmark acquisition and implementation support each get an effort cap; an
unavailable oracle is not authorization for an infrastructure campaign.

Freeze training/validation separation before fitting. Do not inspect held-outs
for state/parameter/seed selection and later call them independent predictions.
Report a rejected prediction without tuning it away.

## Record boundary

Archive filled plans, input decks, source snapshots, runtime manifests, logs,
checkpoints, profiles, analyses, tables, figures and reports as write-once objects.
Corrections are new objects with explicit parent/supersedes links, not edits to
old observations. Authenticate a fresh download before citing an archive as
portable evidence. Keep code independent of the presence of data-heavy Git
history; archive exact executed bytes rather than requiring old data commits.
See [storage rules](artifacts.md) and [lessons from CO](physics/co-experiment-lessons.md).
