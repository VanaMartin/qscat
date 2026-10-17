# CO phase extraction: prediction, model selection and width ambiguity

## Problem and scope

CO's existing all-order/window controls reject full-width stability. A constant
background dominates several failures, but selecting orders 2–4 after seeing
those outcomes would not validate another estimator. This note defines a new,
prospective toy-stage estimator and a known-pole calibration. Original gates
and failed fits remain evidence; a new result cannot replace them.

All internal energies and full widths are Hartree, and phases are radians.
The coefficient count `nback=1..4` denotes background degree 0..3. The existing
`phase_fit.fit_phase` variable-projection fitter is retained unchanged:

\[
\eta(E)=\operatorname{atan2}(\Gamma/2,E_r-E)
 +\sum_{j=0}^{nback-1}b_j x^j,
\qquad x=(E-E_m)/E_s.
\]

Its positive-width resonance factor has pole `z=E_r-i Γ/2`. The independent
oracle constructs the unitary real-axis S factor directly,

\[
S(E)=e^{2i b(E)}\frac{E_r-E+i\Gamma/2}{E_r-E-i\Gamma/2},
\qquad \eta(E)=\arg S(E)/2.
\]

The denominator vanishes at the stated pole; the analytic exponential background
has no finite poles. A neighboring-pole stress case multiplies by a second
factor. That case is outside the isolated-pole premise. These constructed
oracles do not establish the analytic continuation of native CO eigenphases.

## Prospective estimator

`phase_extraction.estimate_phase(energy, phase, windows, settings)` consumes a
complete increasing grid, wrapped or unwrapped phases, Hartree windows and a
frozen numerical contract. No native fit or known answer supplies its seeds.

1. Reserve original grid indices `i % 5 == 0` as audit data **before windowing**.
   They are unused in order selection, initialization and ambiguity construction.
2. In every declared window, fit all four orders to development points with
   three fixed grid-only seeds `(center fraction, width/span)`:
   `(0.4,0.5)`, `(0.5,1.0)`, `(0.6,1.5)`.
3. Inner validation divides the development indices into four interleaved folds.
   Each fit starts again from the grid-only seeds, not a full-data fit. Select
   the best residual among starts within a fit; preserve every start and fold.
4. An order/window is predictively admissible only when every fit succeeds away
   from optimizer bounds and the largest inner-held error modulo π is at most
   **0.001 rad**. This stricter model-adequacy threshold is fixed before native
   use, 50 times below the 0.05-rad scientific comparison budget. It is not a
   measured native-noise confidence level.
5. Require at least one admissible order in **every** window. Choose the lowest
   admissible order separately per window; retain **all** admissible orders for
   parameter-envelope checks. Any failed/bound fit leaves the estimator verdict
   unset, even when other orders predict well. All nonadmissible rows stay saved.
6. Across all admissible development starts and inner folds, position span must
   be at most 0.05 eV and width span at most `max(0.001 eV,5% median width)`.
   Compare every admissible development fit with its untouched audit points;
   their maximum phase error must also be at most 0.001 rad. The audit may reject
   the frozen choice, but never selects a replacement order/window.

## Constructive width-ambiguity veto

Prediction alone need not identify a width. For every admissible development
best fit, fix two alternative widths displaced by **1.1 times** its width budget
and refit the center and linear background to development points. Negative
alternative widths are excluded as unphysical. Center minimization uses three
fixed grid-only starts and the same window bounds. If any successful interior
alternative predicts the development phases within 0.001 rad everywhere, it is
a constructive counterexample to width identification at that tolerance: the
estimator must reject. Save the actual alternative, residual and predictions.

Failure to find such an alternative is only a named adversarial check, **not a
global uniqueness proof**. A passing estimator is a provisional isolated-pole
phase representation under its declared tests, not a certified native S-matrix
pole or statistical uncertainty interval.

## Calibration and stopping

The executable matrix is frozen in
`projects/ukrmol_co/phase-extraction-contract.json` before its full run. It uses
31/61/121-point analytic grids, multiple background degrees, a smooth
nonpolynomial background, quantization-scale phase perturbations, a neighboring
pole, no-pole phases, and an intentionally underidentified narrow window.
The generator never calls the phase predictor under test.

Calibration requires all three resolved linear-background controls to be
accepted and recover center/full width to 0.1 meV; designated narrow-window and
no-pole controls must remain rejected or unset. Every accepted one-pole case
must recover its known center and width within the unchanged scientific budgets.
Neighboring-pole behavior is recorded as a stress diagnostic without interpreting
one of its poles as a unique physical answer. No failed case may be removed.

Only a passing calibration permits the separately frozen three-anchor native
reading. Raw source hashes, original exits, maps, grids and native fits must
qualify first. Native windows and audit splits are fixed before that reading.
Any false/unset calibration stops native application; there is no tolerance
adjustment, optimizer retry ladder or native calculation in this study.

This remains a Python toy/validation capability in `projects/ukrmol_co`.
Optimization or promotion needs separate profiling and independent review.
