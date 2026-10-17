# Width identification needs more than local phase fidelity — 17 October 2026

**A flexible background can hide a 10% known width error while predicting
unused phases within 0.001 rad.** The prospectively frozen 405-cell analytic
endpoint-isolated audit finds **94 counterexamples**. All 45 mandatory true-pole controls pass,
including smooth nonpolynomial backgrounds at the three native-like scales.
Better smooth-background prediction alone is therefore insufficient motivation
for another native polynomial-fit order expansion.

## What was tested

The [design](../../docs/physics/co-phase-identifiability.md),
[contract](phase-identifiability-contract.json), implementation and two focused
tests were committed at **`c273a71` before the original complete run**. Review
then found development/endpoint overlap in 14/45 grids, invalidating the
original fully unused-data declaration. Original exit 0, producer qualification
and 98 reported matches remain [preserved](qualification-evidence/phase-identifiability/README.md).
The endpoint-exclusion repair and regression test were frozen at **`69f13b6`
before repeating** the identical matrix. Independent/public membership checks
now verify all 45 development sets disjoint from all audit sets. This is an
analytic identifiability audit, not a new extraction estimator. Each scalar
S factor has an exact known pole `E_r − i Γ/2`, exact real-axis unitarity and
an entire nonzero background. The true background includes a nonconstant
exponential term. A candidate uses a Chebyshev polynomial background and an
exact pole width 0.9, 1.0 or 1.1 times the true width, with center fixed.

All fixed degrees **3/7/11**, three windows and **31/61/121** original-grid sizes
are retained. Linear background fitting uses only development points, reserving
every fifth original index before windowing. Validation uses reserved points,
both window endpoints and a separate 10,001-midpoint grid with development
coincidences removed. No audit phase selects a degree or fits a coefficient.
The finite dense grid is not an interval-supremum certificate.

| Analytic scale | Counterexamples / 54 tested wrong-width cells | Degree 3 / 7 / 11 | True-pole degree-11 controls |
|---|---:|---|---:|
| Resolved smooth background | 0 | 0 / 0 / 0 | 9 / 9 pass |
| Compressed-like smooth background | **27** | 0 / 9 / 18 | 9 / 9 pass |
| Equilibrium-like smooth background | **13** | 0 / 5 / 8 | 9 / 9 pass |
| Stretched-like smooth background | 0 | 0 / 0 / 0 | 9 / 9 pass |
| Narrow underidentified smooth background | **54** | 18 / 18 / 18 | 9 / 9 pass |
| **Total** | **94 / 270** | **18 / 32 / 44** | **45 / 45 pass** |

The other 135 cells have the true width. Only their 45 degree-11 controls are
mandatory success requirements; every other result remains recorded. Counter-
examples require every unused audit set to pass, a full-rank background solve,
unitarity within 3e-15, and a known width error beyond the unchanged 5%-with-
0.001-eV-floor budget. Their exact finite pole positions follow analytically
from the S factor; no fit estimates those truth values.

## Representative high-density examples

Both examples below use **121 original points, degree 7, and a +10% width**.
The largest unused errors include reserved points, dense points and endpoints:

| Scale / window (Hartree) | True → candidate full width (Hartree) | Largest unused phase error |
|---|---|---:|
| Compressed-like / 0.0735–0.18375 | 0.085 → 0.0935 | **0.000801870 rad** |
| Equilibrium-like / 0.06615–0.121275 | 0.042 → 0.0462 | **0.000685684 rad** |

These examples exceed the width budget by a factor of two yet clear the
0.001-rad phase criterion. Across grid sizes, counterexample counts are
compressed **8/8/11**, equilibrium **4/5/4**, and narrow **18/18/18**. Counts
need not be monotone: they are discrete threshold diagnostics with different
development grids, not a convergence-rate claim. Resolved/stretched zero counts
mean only that this fixed-center finite family found none; uniqueness is open.

## Physical interpretation and decision

The polynomial exponential background has no finite poles or zeros, so this
failure is not caused by hiding an extra finite resonance. However, these local
analytic/unitary models are **not qualified finite-range causal CO models**:
their threshold/high-energy behavior and complex-plane growth are unconstrained.
The audit establishes a limit of local fidelity, analyticity and unitarity alone,
not a physically admissible 10% uncertainty interval for native CO.

The next extraction design needs an independent physical restriction or direct
pole evidence before fitted widths can rank expensive electronic models. A
bounded next design should first assess whether saved eigenpair/boundary exports
and outer-channel data support reconstructing the relevant energy-dependent
R/K/S matrix, including its long-range polar boundary conditions. If they do,
declare independent real-axis reconstruction and complex-pole validation gates
on a solvable scattering model before CO use. If key channel information is
missing, record that gap and cost a named outer-data/export control. An
eigenphase sum alone does not supply the full multichannel matrix. No assumed
single-p-wave threshold law or new native solve follows from this audit.

The first [prospective estimator](PHASE_EXTRACTION.md) retains its compressed/
equilibrium rejection and stretched provisional pass. All earlier all-order/
window and channel rejections remain unchanged. No native data was used in this
new study; no alternative estimator or native application is released.

## Reproducibility and cost

[Corrected quoted verification](verification-phase-identifiability-endpoint-20261017.json) retains
all case-level counts, illustrative complete fits, original exits and source
hashes. [Seven corrected public payloads](qualification-evidence/phase-identifiability-endpoint/README.md)
contain all 405 cells, exact audit grids, coefficients, host/local reports,
six source files and resource provenance. Linux/macOS agree on every cell
verdict; maximum phase-error metric drift is **3.02e-14 rad**. This is platform
replay of the same algorithm; independence of the oracle is its separate direct
S-factor construction with known pole and nonpolynomial background.

Sadaharu cost: **0.965 s wall / 419.7 cells per wall-second**, original command/
container exits **0/0**. The kernel memory high-water value sampled every 0.2 s
while running reaches **0.060 GiB**; it is not a post-exit peak. Reserved physical
cores 0–3, one BLAS thread, 4-GiB limit and pinned tools image were used. Persistent
sources/outputs/logs live in `prepared/phase-identifiability-endpoint-source-20261017/`.
The original run remains in `prepared/phase-identifiability-source-20261017/`.
**22 focused extraction/phase tests pass.** Production remains gated by model,
continuity, extraction, neutral-reference and restartable-pilot qualification.
