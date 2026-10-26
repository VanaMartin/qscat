# Next electronic-model control: compressed CAS12 DZ versus TZ

## Decision

**26October12:39UTC outcome:** the [QC prerequisite stops at singlet A2](../../projects/ukrmol_co/qualification-evidence/compressed-cas12-tz-qc/README.md).
Orbital gradient passes after25 macros, but the singlet-A2 CI flag fails; raw
printed nominal singlet roots include S²≈2. No successful target/Molden export
or fresh-CI audit exists. Original runner1, analyzer unrun, owner/scientific2;
no timeout/OOM,1.36h. The physical DZ/TZ comparison remains **inconclusive**.
This failed optimization is not a target input, even if a future fixed-orbital
CI calculation could recover valid singlets. Its objective was contaminated.

**Stop this route under the frozen contract.** The24-probe proposal below is
blocked, not released on the saved failed checkpoint. No retry, alternative
seed/larger space, native import, scattering, pole, B2 or other model successor.
Next is a separately reviewed singlet-A2 correctness/optimizer-objective decision
from saved evidence, not another full QC run or an automatic solver campaign.

The matched CAS11/CAS12 B1 outgoing-pole comparison rejects position
insensitivity:86.774662meV exceeds50meV, with passing numerical/character/local
identity gates. Its full-width difference16.066613meV passes. Neither model is
qualified; optimized orbital frames/imported representations also change.
[Completed evidence](../../projects/ukrmol_co/qualification-evidence/compressed-cas12-pole/README.md).

Select **CAS(10,12), cc-pVTZ, R=1.9bohr, the unchanged forty-component ensemble**
as the next physical model control. The intended eventual comparison is against
the saved CAS12/DZ B1 coupled pole, using five retained states in each sector,
own2048 roots and the same scattering/pole definition. It asks whether DZ
remains an adequate atomic-basis model at the tested active-space size.
Reoptimization is part of each basis model, not an isolated frozen-orbital AO
perturbation. No basis or active-space completeness claim follows from one pair.

The first and only released calculation is **one missing TZ QC target**, projected
from the already qualified compressed CAS12/DZ import checkpoint. No completed
compressed CAS12/TZ target is available in the inspected run inventory. Existing
compressed TZ/aug-TZ attempts are CAS11 and rejected; successful stretched or
equilibrium targets do not supply this geometry/model. Do not repeat them.

## Why this allocation rather than the alternatives

| Candidate | Decision value / reason |
|---|---|
| CAS12 B2 | Cannot reverse the resolved B1 model rejection; remains stopped. |
| More2048-root refinement or exterior controls | No numerical gate failure explains the observed shift; preserve existing economical baseline. |
| CAS13 by adding another A1 orbital | Extends only one symmetry sector and changes another optimized frame. No saved qualified target or measured scattering cost; not selected merely because13 is larger. |
| Common-frame CAS11/CAS12 attribution | Could split relaxation from correlation, but requires another deliberately unoptimized target/property/scattering model. It would not test the still-missing atomic-basis direction; defer rather than open a new interface/attribution campaign. |
| New50 orbital frame | Previously rejected physical A1 retained-subspace screens remain; parked. |
| CAS12 DZ→TZ, unchanged40-component average | Tests an independent required model direction without changing active counts/ensemble/geometry; reuses the qualified DZ seed and unchanged solver. Selected. |

The literature motivates checking basis/active/state-count dependence, not
choosing a target to reproduce published equilibrium numbers. Dora et al.2016,
printed p.6/Table4, reports substantial model dependence and non-monotonic active-
space improvement; p.3 describes differing DZ/TZ orbital ensembles. **We hold
our forty-component ensemble fixed**, rather than silently adopting its extra
A1 states for TZ. Dora and Tennyson2020, printed pp.3–4/Table2, does not make its
large basis an accuracy certificate. See the tracked
[2016](../../reference/literature/dora-2016-epjd70-197.md) and
[2020](../../reference/literature/dora-2020-jpb53-195202.md) reference notes.
Those are source-note claims, not newly verified full-text measurements here.

## Exact initial scope and stop rules

Use existing `run --qc-only` → `target.build_target` → `analyze.qc_properties`.
No numerical-source changes, engine build, package installation, new optimizer
or CI rescue. Active counts6/3/3/0, two frozen A1 cores, ten active electrons,
five roots in every spin/irrep with equal1/40 weights, singlet-A1 `spin0`, fresh
CI, trial space160,150 macroiterations/ordinary200-cycle CI cap. Retain existing
energy/gradient/CI/residual and MO/spin/Pi/fresh-root gates. Project the actual
qualified native DZ checkpoint at the same geometry; only the basis changes.

Four distinct physical cores0–3; one BLAS environment thread with the existing
PySCF four-thread setting,24GiB container /16000MB PySCF memory. Conditional
0.5–4h forecast: compressed DZ QC took1724s, while older projected equilibrium
CAS12/TZ failed after29193s. Neither is a controlled TZ timing prediction.
Eight-hour hard container cap, bounded owner teardown, minimum32GiB available
RAM /8GiB free disk. At most60minutes implementation/freeze effort for the exact
owner/recipe/positive-negative routing controls; no solver development.

Freeze and push source/tests/exact owner/seed hashes before launch. Keep inputs
read-only; copy the cached pinned scripts archive to fresh output cache rather
than modifying the old cache. Preserve source/archive/image/checkpoint hashes,
original runner/analyzer/container/owner/service exits and profiled resources.
Failure, timeout, wrong spin, optimizer/fresh-CI/gradient failure or identity
failure stops inconclusive **without repair/retry, competing-start queue,
automatic import or another model job**.

## Preserved prospective release policy — QC pass did not occur

**Launch observation26October11:19UTC:** the exact source/tests/owner freeze
`317ae15` precedes11:17:37UTC launch. All eight runtime controls, including the
analytic singlet/full-action control, pass without skips in PySCF2.11.0. Nine
root/seven source pins and23 Git-matched capsule sources authenticate. RHF
converges and the intended CAS(10,12)/40-state optimizer is active; no QC verdict
or terminal runner/analyzer exit exists at this observation.
[Launch record](../../projects/ukrmol_co/progress-compressed-cas12-tz-qc-launch-20261026.json).

QC pass supplies an orbital/self-consistency candidate only. Independently
authenticate its output and compare DZ/TZ ranked energies/excitations, dipole
and MO metadata as diagnostics—not identified-state continuity or a pole budget.
Then separately review the minimum fixed-orbital/full-action coverage and native
target import needed for this new TZ frame. Existing successful DZ inputs remain
unchanged. No integral, DENPROP, native target/scattering solve, pole map, B2 or
sweep is released by this QC contract.

The [conditional coverage proposal](../../projects/ukrmol_co/compressed-cas12-tz-coverage-proposal-20261026.json)
is **unreleased**: one new checkpoint,24 probes rather than repeating the old
96-probe/two-seed scan. Per spin/irrep use5roots/space160 for independent action
on averaged roots,8/160 for candidate native roots/common-first-five count, and
8/240 for the larger-request space check. Reuse the published residual/full-
action algorithms; the maintained helper's AST matches the published original.
No80-space or5/240 repeat is needed for those declared comparisons. This is
bounded root coverage/agreement, not an absolute missing-state/electronic-error
certificate. Failure stops; exact checkpoint/source/owner freeze and a new review
are required before any execution.

An eventual valid DZ→TZ B1 pole rejection would rule out treating DZ as
basis-insensitive at CAS12. A pass would only support that local basis contrast;
the resolved CAS11/CAS12 sensitivity, independent active-space saturation,
domain-wide channel/identity/neutral/embedding uncertainty remain. Any future
pole uses50meV position /5% reference FULL width with1meV floor and unchanged
operational bands, declared before that pole is evaluated. No electronic model
is chosen because it is larger, cheaper or closer to a literature candidate.

All114 holdout S matrices, equilibrium pole exclusion, nine working reference,
2048 baseline and parked new50 remain protected. No potential fit is authorized.
