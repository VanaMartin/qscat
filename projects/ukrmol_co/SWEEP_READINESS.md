# CO: validity gates for data, fits and physical claims

Apply each gate to the claim it protects. A numerical failure is not a missing
physical observable; a failed physical-insensitivity diagnostic does not erase
an accepted finite-model row. The [two-goal plan](RESEARCH_PLAN.md) recommends
local neutral and resonance datasets, **not a production physical CO model**.
No job or fit is released by a checklist.

## Data acceptance and local sweep gates

| Stage | Required evidence | Failure action |
|---|---|---|
| Every neutral/target row | Immutable source/image/input hashes; actual R/units, state, basis/core/active/ensemble, total-energy/nuclear-repulsion and reference conventions; applicable final convergence, residual, spin, density, stability/orbital/CI checks | Preserve failure/partial data; no energy from a failed iterate. Reuse successful checks where applicable, not every historical control. |
| Every scattering/outgoing row | Actual native dimensions, thresholds and channel labels; full finite closed/exterior coupling; applicable independent import/root/dipole and real-axis normalization/unitarity/reciprocity/flux/Wronskian guards | Stop that release on its first failed numerical gate; phase agreement does not excuse a failed flux guard. |
| Outgoing label | Declared continuous sheet and admissible contour/search domain; simple isolated visible pole, residual/null separation, nested moments, nonzero derivative, rank-one residue and local character witnesses; numerical error resolved | Determinant zero alone is insufficient. Unexposed/ambiguous/unresolved result remains unset or rejected, not Gamma=0. |
| Conventional estimator label | Frozen collision-energy/window/background/unwrapping rule, nearby thresholds, residuals and grid/window/background adequacy/sensitivity | Keep estimator identity; failed windows remain failed. BW/time-delay/peaks are not automatically outgoing poles. |
| One coherent sweep/branch | Fixed protocol and energy conversion; target/subspace and pole/residue correspondence in a common frame; paired B1/B2 consistency before calling a Pi curve | Split/end segment on threshold collision, lost isolation or ambiguous identity. Independent valid local rows may still be stored; no smoothing or protocol splice rescues a branch. |

Generic finite-step target-overlap triggers are diagnostics, not proof of a pole
branch jump. They motivate discriminating tests rather than automatically ending
the campaign. A frozen interpolation candidate may be prospectively tested while
common N+1 identity remains unresolved; such a test does not waive that uncertainty
or certify a unique branch. See the current disposition in the [handoff](CONTINUATION.md).

## Fit gates

- Declare domain, units, energy zero, numerical uncertainty and supported
  parameter count before fitting; test residuals and off-node representation.
  Existing neutral midpoint checks are representation evidence, not independent
  physical validation; historical neutral data are not globally blind.
- Preserve 5Z/QZ molecular-reference zeros. Current `fit_neutral` uses a
  zero-asymptote EMO with no additive offset: use an explicitly local/offset
  representation or a matching computed fragment zero, not MR atoms for CC data.
  A local table does not determine physical `De`.
- Proposed outgoing resonance fitting needs at least **four actual compatible
  training geometries**, bracketing excluded R=2.1323. B1/B2 duplication,
  interpolation and default `n_nodes=24` do not create acquisition requirements
  or independent nodes. Four is a floor, not an identifiability certificate.
- Store tables before resolving local/asymptotic API coupling. Do not invent
  `ea`, `R_inf`, tails or an asymptote to satisfy `ResonanceTarget`/fit defaults.
- Freeze model/assignment/numerics before separately releasing equilibrium
  outgoing validation. All 114 S holdouts remain untouched. A reduced pole
  prediction is not full-S or absolute electronic validation.
- Every off-node test retains its own frozen coefficients, development membership
  and prediction budgets: no validation-driven refitting. A later candidate may
  explicitly reclassify old test points as development, never as unseen validation.
  Representation subbudgets do not replace numerical or physical gates, and
  selected-point passes are not a uniform interpolation-error bound.

## Physical qualification: not entrance gates to local data collection

| Claim | Still required |
|---|---|
| Electronic-model convergence / accuracy | Quantified basis, active-space, ensemble, channel and omitted-spectrum sensitivity; matched physical comparisons and honest residual uncertainty. Compressed CAS11/CAS12 and stretched DZ/TZ insensitivity tests remain failed. |
| Global neutral / dissociation | Stable spin/core-consistent correlated treatment, overlap/size-consistency/BSSE and finite-R/tail evidence, matched fragment zero; separate `De`, `D0` and ZPE. Accepted MR R3/R4 are not a qualified dissociation limit. |
| Global resonance-to-bound branch | Threshold/state identity through crossing, confirmed bound continuation, charge/spin-consistent fragments and justified asymptote. An empty fit or unresolved positive-energy width is not binding. |
| Full scattering / VE | Independently specified reciprocal channel/background embedding, held-out predictive test and qualified QSCAT discretization. Pole curves alone omit residue/background information. |
| Attachment / DA | Applicable global branch/fragment thresholds (`D0−EA`) plus a nuclear-observable check. Lowest-Pi channel-specific dynamics does not establish total thermal attachment near 10 eV. |

External curves are optional benchmarks, not mandatory replacements for our
computed neutral data. A separate NRM/doorway Hamiltonian, exact affinity,
global bound joining and full-S validation are not prerequisites for the first
local E0/Er/Gamma rows.

## Unchanged quantitative boundaries

- Same-model numerical pole controls: **0.1 meV position/full width** and
  **1e−6 rad** in specified tests. Per-study residual/K/S/flux gates remain distinct.
- Physical comparisons: **50 meV position**, **max(5% Gamma, 1 meV) full width**,
  **0.05 rad in the specified phase metric**. Their failures reject those claims,
  not all finite-data collection. No lower-accuracy tier is adopted.
- Neutral pilot 20-meV relative-basis / 1-meV interpolation checks do not qualify
  absolute alignment. Atomic 10-meV EA / 20-meV fixed-R-offset diagnostics remain
  failed overall; fixed-R offsets are not `De`, `D0` or initial-v thresholds.

Keep nine working/new50 parked and 2048 roots retained. No silent threshold-grid
truncation, basis/root/channel escalation or repeat of accepted work. The
[results and reopening ledger](../../docs/physics/co-results-and-rejected-routes.md)
and [production decisions](PRODUCTION_DECISIONS.md) govern disposition.
