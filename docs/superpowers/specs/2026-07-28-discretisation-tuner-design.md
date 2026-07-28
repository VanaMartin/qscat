# Automatic discretisation tuner — Design Spec

**Date:** 2026-07-28
**Author:** Martin (martin@qscat.com) with Claude
**Status:** Approved design — spec for review
**Lifecycle:** `qm-method-lifecycle` — a NEW capability (numerical tooling): `qscat.tuning`
(deterministic primitives) + a `discretisation-tuner` skill (the supervisor). Realizes the
future-work note in the `per-molecule-discretisation` memory / `docs/physics/diatomic-ve-cross-sections.md`.

## Context

Every FEM-DVR-ECS calculation in this repo hinges on a hand-tuned grid, and **discretisation
errors have been the single most expensive class of bug**: the shared N₂-style nuclear grid
silently under-resolved the K≈58 dissociative-attachment wave (σ_DA off by ~36 orders); the LCP DA
needed the fine per-molecule grid *and* a value-vs-coefficient fix; the H₂⁺ Coulomb tail forced a
1300-bohr electronic grid. Each was found the hard way, and each grid was chosen by a human's
"good eye" for element lengths, constrained by the ease of writing round-number segments
(0.1/0.3/1.0…). Now that the grids are Python (`ElementSpec` takes an arbitrary length), that
constraint is gone: we can compute the **optimal, potential-adaptive discretisation** — minimal
memory footprint and factor-time at a fixed precision — instead of guessing it.

## Goal

Given a `ResonanceModel`, a coordinate (electronic r / nuclear R), and a target energy range,
**automatically produce the grid with the fewest DVR points that holds a target precision** — the
element mesh, quadrature order, real-region cutoff, and ECS angle/length — with the evidence that
it is converged, and validate the tool against the known-good eMoScat decks.

## Approach (the decisions from brainstorming)

- **Hybrid**: a physics **prior** lays out a good grid a-priori (no solves), then a few targeted
  **convergence probes** validate/refine it to tolerance.
- **Decoupled 1-D probes + a final 2-D spot-check**: exploit the tensor-product structure — tune
  each coordinate on cheap 1-D sub-problems; one full 2-D solve confirms.
- **Fully-adaptive equidistribution mesh**: an explicit per-element length list (each element
  sized to its local de Broglie phase), the minimal-point layout — the "no good eye" win.
- **Tune every knob**: element lengths, quadrature order (the h/p tradeoff), real cutoff, ECS
  angle + tail length/growth.
- **Precision target**: 1-D spectral proxies (nuclear vibrational eigenvalues; electronic
  resonance-pole / bound-state energies; channel-function representation) **plus a final 2-D cross
  section spot-check**, all to `rtol=1e-3` by default (user-overridable).
- **General**: works for any `ResonanceModel` (neutral + ionic) and either coordinate.

## Architecture

- **`qscat.tuning`** — deterministic, pure, unit-tested primitives (no LLM judgment):
  - `analyze` — potential analysis: the local wavenumber profile `k(x)`, turning points,
    singularities, asymptotics, channel wavenumbers over the energy range.
  - `mesh` — the equidistribution mesh generator + the quadrature (h/p) sweep.
  - `ecs` — the ECS-tail tuner: max-valid-angle scan, tail-length sizing, cutoff.
  - `probes` — the 1-D convergence probes (nuclear vibrational; electronic resonance/bound +
    channel representation) at the energy extremes, returning `(observable, converged?, cost)`.
  - `metrics` — cost model (DVR point counts; estimated 2-D matrix size / memory / factor-time)
    and the precision aggregation.
- **A `discretisation-tuner` skill** — the supervisor (the LLM-driven loop): analyze → propose →
  probe at the extremes → interpret → coarsen the over-resolved / refine the under-resolved →
  re-probe → converge → emit the config + a report. The judgment (which knob, when converged, the
  cost/precision tradeoff, edge cases) lives here; the mechanical work is in `qscat.tuning`.

## The physics (what the primitives compute)

### Equidistribution mesh
A GLL-DVR element of order `q` resolves oscillations up to ~`(q−1)/2` wavelengths (spectral
accuracy), so place element boundaries `xₖ` so the **de Broglie phase per element is constant**:

  `∫_{xₖ}^{xₖ₊₁} k(x) dx ≈ Φ_target(q) = C·(q−1)`,   `k(x) = √(2m·max(E_max − V(x), 0))`,

with `C` a safety constant **calibrated once** so the tool reproduces the eMoScat decks (below).
In classically-forbidden regions (`V > E_max`) use the decay rate `κ = √(2m(V−E))` and resolve the
decay length instead (element ~ a fixed fraction of `1/κ`), plus refinement at turning points and
singularities (the `−1/r` origin). Sweep `q ∈ {6,8,10,14}` and pick the `(mesh, q)` minimizing the
**total DVR points** — the h/p optimum (higher q ⇒ fewer, denser elements; the cost model decides).

### ECS tail
- **Max angle θ**: increase θ until (a) the potential's analytic continuation misbehaves — a
  Gaussian `e^{−αr²}` diverges for θ>45°, `−1/r` is fine but its origin sits in the real region —
  detected by the potential growing on the rotated contour; or (b) the tail eigenvalues stop being
  θ-stable (an angle scan on the 1-D problem). Take a safe fraction of the max.
- **Tail length**: the fastest outgoing wave decays `e^{−K(x−R₀)sinθ}`; extend until it reaches
  ~machine-ε, with tail elements equidistributed on the rotated contour (still oscillating +
  decaying). `K` is the largest channel wavenumber over the energy range.
- **Cutoff R₀** (real→complex pivot): where the interaction has died and the asymptotic channel
  form holds (from the potential analysis).

### The 1-D probes (validation, at both energy extremes)
- **Nuclear**: `vibrational_states` eigenvalues `eps_v` stable to `rtol` under one refinement step;
  the outgoing `riccati_bessel_en_mass(K_max)` represented to `rtol` (quadrature/projection error).
- **Electronic**: `find_resonance_pole` position + width (or the bound/Rydberg energies via
  `anion_electronic_states`) stable to `rtol`; the incident channel function
  (`riccati_bessel_en`/`coulomb_f_en`) represented to `rtol`; the ECS absorption of the fastest
  wave.
- **Extremes**: `E_max` (finest resolution) AND near-threshold `E_min` (longest wavelength, largest
  classically-allowed extent, near-threshold sensitivity). "Structures" the de Broglie prior misses
  (a narrow resonance needing local refinement in R) are caught by the resonance-pole probe and
  fed back as a local refinement.
- **Final 2-D spot-check**: one full `ve_cross_section`/`da_cross_section`/… solve at the hardest
  energy, confirming the tensor-product grid delivers `rtol` (this is the only 2-D solve; for a
  huge deck like H₂⁺ it runs on the proxy or under Docker/MUMPS and is reported as such).

## Output

- A validated **grid config** per coordinate: the `ElementSpec` length list (real + ECS tail),
  quadrature order, ECS angle, and cutoffs — directly consumable by `FemDvrEcsGrid`/the builders
  (and expressible as a committed deck).
- A **report**: achieved precision (per-probe convergence numbers), cost (DVR points, estimated 2-D
  matrix size / memory / factor-time vs the old hand grid), and the tuning decisions.

## Validating the tuner itself

The tool is only trustworthy if validated:
1. **Reproduce-or-beat the eMoScat decks** — for N₂/NO/F₂/H₂⁺, the tuner's grid must reach the
   same-or-better probe precision at same-or-fewer DVR points than the committed hand-tuned deck.
   This is also where the equidistribution constant `C` (and the ECS safety fractions) are
   **calibrated**.
2. **Flag the known failures** — the tuner must diagnose the coarse shared N₂ grid as
   *under-resolved* for F₂ DA's K≈58 wave (the bug that cost ~36 orders), and the coarse H₂⁺
   electronic grid as under-resolved for the Coulomb incident. These become regression gates.

## Sub-project decomposition (tasks for the plan)

1. **Potential analysis** (`qscat.tuning.analyze`) — `k(x)` profile, turning points, singularities,
   asymptotics, channel wavenumbers; pure, unit-tested on analytic potentials.
2. **Equidistribution mesh + h/p sweep** (`qscat.tuning.mesh`) — element-length list from a `k(x)`
   profile + order; the point-minimizing quadrature pick; forbidden-region/turning-point handling.
3. **ECS-tail tuner** (`qscat.tuning.ecs`) — max-angle scan, tail-length sizing, cutoff.
4. **1-D convergence probes** (`qscat.tuning.probes`) — nuclear + electronic, at the extremes,
   with the channel-representation + ECS-absorption checks; `(observable, converged?, cost)`.
5. **Cost/precision metrics + a `propose_grid` a-priori assembler** (`qscat.tuning.metrics`,
   `qscat.tuning.propose`) — the cost model and the one-shot a-priori grid.
6. **The `discretisation-tuner` skill** — the supervised loop procedure, calling the primitives,
   emitting the config + report.
7. **Calibration + validation** — calibrate `C`/ECS fractions and gate against the eMoScat decks +
   the known coarse-grid failures.

## Out of scope (this sub-project)

- **Coupled-channel / non-adiabatic grid coupling** — the decoupled 1-D + 2-D-spot-check model
  assumes the tensor-product structure suffices (it does for the current models).
- **Auto-running huge decks** — the 2-D spot-check for a non-laptop model (H₂⁺) runs on the proxy /
  is Docker-deferred, as elsewhere.
- **Rewiring existing models onto tuner-generated grids** — the tuner EMITS configs; adopting them
  for N₂/NO/F₂/H₂⁺ production is a separate, opt-in follow-on (the committed decks stay until then).
- **Optimizing beyond the FEM-DVR-ECS family** (spectral/other bases).

## Verification

- `uv run pytest -q -m "not slow"` pass; `uv run mypy libs/qscat/qscat` 0; `uv run ruff check .` clean.
- `qscat.tuning` primitives unit-tested on analytic potentials (known turning points, known
  equidistribution meshes, ECS-angle limits).
- The tuner reproduces-or-beats the N₂/NO/F₂/H₂⁺ eMoScat decks (same/better precision, same/fewer
  points) and flags the coarse-grid failures — gated in tests.
- The `discretisation-tuner` skill runs end-to-end on at least one molecule, emitting a validated
  config + report; `docs/physics/` gains a discretisation-tuning note; `CLAUDE.md` + the skills
  table updated.
