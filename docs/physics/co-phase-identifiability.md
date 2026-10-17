# Local phase fidelity does not by itself identify a resonance width

## Decision to resolve

The first prospective CO estimator rejects required-window phase prediction at
compressed and equilibrium geometries, and all its smooth-exponential analytic
cases reject. More flexible backgrounds could improve prediction while making
pole parameters less identifiable. Before choosing another estimator, quantify
that tradeoff with a bounded **analytic counterexample audit**. This is not an
order-selection procedure or a new extraction estimator.

Energies and full widths are Hartree; phase is radians. On a compact positive
real-energy interval define the explicit local unitary representation

\[
S(E;z,b)=\exp[2i b(E)]\frac{E_r-E+i\Gamma/2}{E_r-E-i\Gamma/2},
\qquad z=E_r-i\Gamma/2,\quad\Gamma>0.
\]

A real polynomial background, or the entire exponential background used in
the oracle, makes the prefactor entire and nonzero. Thus it adds no finite
poles or zeros: the representation has its exact pole at `z` and zero at `z*`.
It is unitary on the real axis. These are local analytic properties; they do
**not** establish a finite-range causal Hamiltonian, correct threshold law,
high-energy behavior, or the polar multichannel asymptotics of CO. Exponentials
of polynomials generally have unsuitable behavior at infinity. The distinction
is essential: this audit tests what local fidelity/unitarity/analyticity alone
can certify, not which backgrounds are physically admissible for CO.

For a fixed alternative pole with width `Γ'`, the needed background on the
real interval is the true unwrapped phase minus the alternative resonance phase.
It is smooth there. Approximate that function with a fixed-degree Chebyshev
polynomial by linear least squares on development points. The resulting `S'`
has the alternative pole exactly, irrespective of the approximation error.
Real-axis predictive agreement can therefore be tested against a **known wrong
width**, without a nonlinear optimizer or an inverse-scattering assumption.

## Frozen experiment

`phase_identifiability.audit(contract)` returns every case, fitted coefficient
and held-grid maximum; its command writes a report before returning its gate.
The contract is `projects/ukrmol_co/phase-identifiability-contract.json`.

- Five analytic scales: resolved smooth-background, compressed-like,
  equilibrium-like, stretched-like, and intentionally narrow underidentified.
  Each has three declared windows and a nonconstant exponential addition.
- Original full grids contain 31/61/121 points. Reserve indices `i % 5 == 0`
  before windowing; only the other points enter linear background fitting.
- Fit fixed Chebyshev degrees **3, 7, 11** at each window, not a retry ladder.
  Compare exact pole widths **0.9, 1.0, 1.1** times the generating width, keeping
  center fixed. Every cell is retained: 405 total. Known poles are inputs to
  this adversarial audit, never seeds for a native estimator.
- Evaluate on 10,001 equally spaced midpoint nodes across each window, removing
  any node within `1e-12 × window span` of a development point. Also audit the
  reserved original-grid points and window endpoints. All audit phases remain
  unused in coefficient fitting or degree choice. Record the actual grids and
  coefficient vectors so predictions are reproducible.
- A width counterexample requires **all three** audit sets to agree modulo π
  within **0.001 rad**, real-axis unitarity within `3e-15`, full-rank background
  solve, and width displacement strictly larger than `max(0.001 eV, 5% Γ)`.
  This finite audit is not a rigorous supremum bound between nodes.
- Nonvacuous representation qualification requires the **true-pole degree-11
  control in every scale/window/grid** to pass all audit sets. This includes
  mandatory nonpolynomial controls at native-like scales. Other degrees stay
  recorded regardless of success. Known wrong poles are diagnostic outcomes,
  not required positive results; absence is not a uniqueness proof.

The audit must record model conditioning and rank, errors on development and
each audit set, maximum-error locations, original command/resource exits, exact
source hashes and runtime. Rank failure or nonfinite data leaves its cell
unset; it cannot become a counterexample. Any missing/failed required control
stops the study with false/unset qualification. No native data is read under
this contract, even when the representation controls pass.

## Consequence for the next physical design

If a known wrong width passes, using a broader background family and smaller
phase residual alone cannot justify a tighter width claim. A future extraction
method needs separately justified restrictions: threshold/asymptotic behavior,
causal or finite-range information, an independent background calculation, or
direct complex-energy/matrix evidence. They require their own applicability
proof for polar multichannel CO. The present scalar S factor is a local
counterexample class, not a replacement scattering calculation.

All prior extraction rejections remain preserved. A successful counterexample
is useful validation evidence for a **veto**, not release of a pole curve. A new
estimator, its mandatory known-truth controls and native application must be
declared separately. This audit stays at Python toy/validation stage.
