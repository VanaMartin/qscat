# Prospective CO phase extraction — 17 October 2026

**The frozen analytic calibration passes its finite requirements, but the new
three-anchor extraction gate rejects.** Saved stretched nine-A2 phases pass
the provisional representation checks in both Pi sectors; compressed and
equilibrium fail required-window phase prediction. No three-anchor width curve
is qualified. All previous all-order/window rejections remain unchanged.

## What was frozen before the outcomes

The [mathematical design](../../docs/physics/co-phase-extraction.md) and
[contract](phase-extraction-contract.json) were committed at `902b76c` before
the full 30-case calibration. The native reader was committed at `4e089e3`
before opening the saved native phases under this estimator.

The unchanged variable-projection fitter uses a positive constant-width phase
step and background degrees 0–3 (`nback=1..4`). New selection reserves every
fifth full-grid point before windowing. Four inner folds select the lowest
predictively admissible order in each required window; all admissible orders,
development starts and folds enter the parameter envelope. Three seeds depend
only on the grid, not native candidates or the analytic answer. Audit phases
never initialize, select or construct alternatives. The audit can reject, not
choose a replacement. Failed/bound fits leave verdicts unset.

Model adequacy requires **0.001 rad modulo π** maximum held-phase error. This
prospective representation criterion is 50 times below the existing 0.05-rad
model-comparison budget; it is not a native-noise estimate. Position and width
envelopes retain **0.05 eV** and **5% with a 0.001-eV floor**. A fixed-width
counterexample search tries widths displaced by 1.1 times the width budget,
refitting center/background only on development points. Finding an interior
compatible alternative vetoes width identification. Not finding one is not
a global uniqueness proof.

## Analytic calibration: all 30 cases retained

The independent generator constructs unitary real-axis S factors directly,
with known pole `E_r − i Γ/2`; it never calls the phase predictor. Each of ten
families uses 31/61/121 points. Linear, cubic, exponential-background,
deterministically perturbed and three native-like scales are accompanied by
narrow-window, no-pole and neighboring-pole controls.

| Family | 31 points | 61 points | 121 points |
|---|---|---|---|
| Resolved linear background | pass | pass | pass |
| Resolved cubic background | pass | pass | pass |
| Resolved smooth exponential addition | reject | reject | reject |
| Resolved perturbed phases | pass | pass | pass |
| Equilibrium-like scale | pass | pass | pass |
| Compressed-like scale | reject | pass | pass |
| Stretched-like scale | pass | pass | pass |
| Narrow underidentified window | reject | reject | reject |
| No pole | unset | unset | unset |
| Neighboring pole, outside isolated-pole premise | unset | reject | reject |

**17 pass / 9 reject / 4 unset.** Every case satisfies its declared calibration
requirement: all required linear controls recover center/width within 0.1 meV;
every accepted isolated-pole case satisfies known-truth scientific budgets;
designated negative controls are never accepted. Unitarity errors are below
3e-15. For perturbed phases, largest accepted center/width truth errors are
**0.00243 / 0.01051 meV**. The compressed-like accepted cases retain nonzero
truth errors up to **11.011 / 21.615 meV**, within their scientific budgets.

This is deliberately a **finite calibration pass**, not universal smooth-
background qualification. Every exponential-background case rejects. Conditional
truth requirements permit those rejections; the neighboring-pole rows are stress
diagnostics, not success requirements for a single-pole model. Local replay
reproduces every case verdict. Full starts/folds/alternatives remain in the
[public payloads](qualification-evidence/phase-extraction/README.md).

## Saved nine-A2 anchor application

Before fitting, the readers reconstruct **270 original manifest entries**
(compressed/equilibrium/saved-model sets: 73/73/124), original stages/resources,
complete spectral/boundary exports, maps, grids, Pi symmetry and nonempty native
fits without MAXFIT. No scattering solve is repeated. Every candidate/start/fold
fit in the new native application succeeds away from parameter bounds.

| Anchor / R (bohr) | Required windows, requested-eV | Selected `nback` by window | Estimator verdict, both sectors |
|---|---|---|---|
| Compressed / 1.900 | 2–5; 2.2–4.8; 2.5–4.5 | none; none; 4 | **reject** |
| Equilibrium / 2.1323 | 1.8–3.3; 1.5–3.5; 1.2–3.8 | 4; 4; none | **reject** |
| Stretched / 2.500 | 0.4–1.6; 0.5–1.5; 0.6–1.4 | 3; 3; 3 | **provisional pass** |

For compressed, even the best inner-held maximum errors in the first two windows
are **0.00225936 / 0.00130390 rad**, above 0.001 rad. For equilibrium, the
broadest window's best is **0.00108757 rad**. Narrower-window fits do not rescue
a missing required window; these are model-adequacy rejections, not optimizer
failures. Their small surviving-fit parameter spans are not accepted width errors.

Stretched admits `nback=3,4` in all windows. Across every admissible development
start and inner fold:

- Position: **0.97407765–0.97415015 eV**, span **0.07251 meV**.
- Full width: **0.28902773–0.28940756 eV**, span **0.37984 meV**.
- Untouched-audit maximum: **0.000298345 rad**.
- No tested beyond-budget width alternative is compatible at 0.001 rad.

These values describe the declared isolated-pole phase representation, not a
certified native S-matrix pole or statistical confidence interval. The original
stretched nine-A2 **6.991% all-order/window rejection** and **0.0679178-rad
channel rejection** after refinement remain in force. The new rule is separately
calibrated evidence, not a replacement of those records.

## Reproducibility, costs and stopping decision

The [quoted verification](verification-phase-extraction-20261017.json) retains
all 30 calibration decisions, all native order/window scores, source hashes and
original exits. The [public verification](verification-phase-extraction-public-20261017.json)
and [replay instructions](qualification-evidence/phase-extraction/README.md)
cover 13 new payloads, 64 reader sources and the three reused raw archives.
Linux/macOS fitted-parameter drift over all 1,080 native development/inner-fold
fits is at most **7.50e-9 eV**. Platform replay shares the same fitting algorithm;
independence of the known-pole oracle comes from its direct S-factor construction.

Both bounded Docker runs used reserved physical cores 0–3, one BLAS thread,
4-GiB memory limit and the pinned tools image. Measured results on Sadaharu:

| Run | Wall time | Sampled kernel memory peak | Original command/container exit |
|---|---:|---:|---:|
| 30-case analytic calibration | 13.883 s | 0.109 GiB | 0 / 0 |
| Qualified saved-native reading | 3.969 s | 0.157 GiB | **1 / 1** |

Memory peaks are kernel high-water readings sampled every 0.2 s while the
container ran; they are not a post-exit peak measurement. Persistent sources,
outputs, logs, resource records and original exits remain under
`prepared/phase-extraction-{calibration,native}-source-20261017/` on Sadaharu.
The native rejection stops this declared estimator study. No tolerance change,
window removal, background-order ladder or expensive successor is authorized.

**Next design decision:** the missing-window prediction discrepancies are now
named and inexpensive to study. A new physically motivated background or
energy-dependent resonance representation needs a separate prospective design,
mandatory nonpolynomial known-truth successes at native-like scales, and explicit
prediction-versus-pole-identifiability checks before another saved-data reading.
No current result distinguishes all possible background/resonance decompositions.
Electronic-model, continuity and neutral-reference controls remain independent
requirements; the production sweep and potential fitting stay gated.
