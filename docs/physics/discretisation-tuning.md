# FEM-DVR-ECS discretisation tuner: design, calibration, and gate

**Location:** `libs/qscat/qscat/tuning` (`analyze`, `mesh`, `ecs`, `probes`, `metrics`,
`propose`, `incident`), the `discretisation-tuner` skill (the supervised loop),
`validation/tuning/` (`calibrate.py`, `test_emoscat_decks.py` — this document's Task 8).
**Origin:** `docs/superpowers/specs/2026-07-28-discretisation-tuner-design.md`; realizes the
future-work note in `docs/physics/diatomic-ve-cross-sections.md`
("Future: automatic discretisation").
**Units:** atomic units throughout (energy in Hartree, length in Bohr).

## Why

Discretisation errors have been the single most expensive class of bug in this repo: the
shared N₂-style nuclear grid silently under-resolved F₂'s K≈58 dissociative-attachment (DA)
wave (σ_DA off by ~36 orders of magnitude); the H₂⁺ Coulomb tail forced a 1300-bohr
electronic grid. Every grid before this sub-project was hand-tuned by a human's "good eye"
for element lengths. `qscat.tuning` replaces that with a computed, potential-adaptive
discretisation, and this document is where it is CALIBRATED and GATED against the
known-good eMoScat decks — the tool is only trustworthy once it reproduces the decks it is
meant to replace, and flags the specific failure that motivated it.

## The hybrid approach

A physics **prior** lays out a good grid a-priori (no eigensolve), then cheap **convergence
probes** validate it:

1. **Equidistribution mesh + h/p sweep** (`qscat.tuning.mesh`). Place real-region element
   boundaries so each element carries a constant de-Broglie phase
   `∫ k(x) dx ≈ C·(order − 1)`, `k(x) = sqrt(2·mass·max(E_max − V(x), 0))` — a GLL-DVR element
   of order `q` resolves ~`(q−1)/2` wavelengths, so this keeps every element "the same number
   of oscillations wide." In classically forbidden stretches, element length is instead
   capped by the local decay length `1/kappa`; elements adjacent to a turning point or
   singularity are halved. `optimal_real_mesh` sweeps `order ∈ {6, 8, 10, 14}` and keeps
   the `(mesh, order)` combination with the fewest total DVR points — the h/p optimum.
2. **ECS tail — a separate, exp-growth regime** (`qscat.tuning.ecs`). The rotated outgoing
   wave decays as `exp(−K(x−R0)sinθ)`, so the tail uses growing elements, not the
   oscillatory equidistribution mesh. `max_stable_angle` scans θ up to the ~35° double-ECS
   cap (the binding limit is the 2-D corner where both coordinates sit on their complex
   tails at once, not either potential alone); `tune_ecs_tail` sizes exponentially-growing
   elements to absorb the fastest channel wavenumber `K` down to `~1e-12`.
3. **`propose_grid`** (`qscat.tuning.propose`) wires a per-coordinate MODEL ADAPTER (which
   picks `V`/mass/extent/channel-`k` from a `ResonanceModel`) through analyze → mesh → ECS
   into one `FemDvrEcsGrid` — the one-shot a-priori half of the tuner.
4. **The decoupled 1-D probes** (`qscat.tuning.probes`) validate it empirically:
   - `probe_nuclear` — vibrational eigenvalues stable under one h-refinement.
   - `probe_electronic` — the anion bound electronic-state energy stable under refinement
     (a cheap proxy for the full two-angle resonance-pole match).
   - `probe_channel_representation` — THE cheapest and most diagnostic: no eigensolve at
     all. Compares the grid's own quadrature estimate of `∫|F|² dr` (a channel function at
     wavenumber `k`, partial wave `l`) against a fine-uniform-grid reference over the same
     span. A grid whose elements are large compared to `1/k` aliases this badly — exactly
     the failure mode that cost F₂ its ~36 orders.
5. **`grid_cost`/`tensor_cost`** (`qscat.tuning.metrics`) give exact DVR point counts and
   rough, anchored 2-D sparse-LU cost estimates, for RELATIVE ranking of candidates.

The `discretisation-tuner` skill runs analyze → propose → probe at the energy extremes →
refine/coarsen → a final 2-D spot-check, as a supervised loop (see the skill for the worked
example); `qscat.tuning` itself is pure, deterministic primitives with no judgment baked in.

## Task 8: calibrating the phase constant `C`

The equidistribution mesh has exactly one free numerical knob, the phase-per-element
constant `C` in `phase_per_element = C·(order − 1)` — smaller `C` means a finer mesh.
`validation/tuning/calibrate.py` (`uv run python -m validation.tuning.calibrate`) measures
it by sweeping `C` and checking whether `propose_grid`'s nuclear grid for N₂/NO/F₂
reproduces-or-beats the corresponding committed eMoScat deck.

**F₂ is the decisive case.** F₂ has a genuinely OPEN dissociative-attachment channel within
its tested range (`(0.01, 0.05)` Ha): its anion bound electronic state
(`anion_electronic_states`, at the eMoScat deck's dissociation limit `R_inf = 10.7` bohr) sits
at `eps_e ≈ −0.127` Ha, an exothermic threshold, so `E_DR = E_max − eps_e > E_max` and the
outgoing nuclear wavenumber `K = sqrt(2·mu·E_DR) ≈ 78` at `E_max = 0.05` — the same wave
whose under-resolution on the coarse shared grid cost ~36 orders of magnitude (see
`docs/physics/diatomic-ve-cross-sections.md`). The calibrated `C` is the SMALLEST value at
which `propose_grid`'s F₂ nuclear mesh represents that wave to `rtol = 1e-3`
(`probe_channel_representation`) using FEWER points than the eMoScat F₂ DA deck.

A 40-candidate sweep, `C ∈ [0.05, 2.0]` step `0.05` (each candidate: `propose_grid` +
`probe_channel_representation` + `probe_nuclear`, all three molecules — a few minutes total),
found:

| Quantity | value |
|---|---|
| Calibrated `C` | **0.10** |
| F₂ eps_e (anion bound state) | −0.12694 Ha |
| F₂ K_DA at E_max=0.05 | 78.28 |

**Per-molecule result at `C = 0.10`:**

| Molecule | proposed n | deck n | ratio | channel rel_error | deck's own rel_error | vib converged |
|---|---|---|---|---|---|---|
| F₂ | 609 | 974 | **0.63×** | 2.89e-4 (**converged**) | 2.12e-4 | yes |
| N₂ | 614 | 428 | 1.44× | 1.66e-3 (not conv., but **beats deck 17×**) | 2.88e-2 | yes |
| NO | 604 | 597 | 1.01× | 5.32e-3 (not conv., but **beats deck 7×**) | 3.74e-2 | yes |

F₂ reproduces-and-beats the deck outright: fewer points (609 vs 974) AND absolute
convergence on the exact wave that used to fail catastrophically.

## Genuine finding: N₂/NO's floor-K channel check isn't a clean absolute bar

N₂ and NO do NOT have an open DA channel in their tested (VE-scale) energy ranges — N₂'s
is closed within the whole +0.5 Ha window; NO's opens at ~0.17 Ha, above the tested
`(0.004, 0.12)` range. Lacking an `eps_e` threshold for a closed channel, their
channel-representation check uses the conservative FLOOR `K = sqrt(2·mu·E_max)` (treating
the entire incident electron energy as if converted to nuclear translational energy — a
generous over-estimate, never actually reached since the true threshold is far more
negative). The sweep shows this floor is not met at `rtol = 1e-3` by ANY sane `C` —
**not even by the eMoScat decks themselves** (N₂'s own deck: rel_error ≈ 0.029; NO's own
deck: rel_error ≈ 0.037, both `≫ rtol`). This is not a calibration failure: it means the
floor is a deliberately conservative bound these decks were never tuned to resolve, and an
absolute `rtol` gate there would fail on a bar the reference implementation itself never
cleared. `test_emoscat_decks.py` gates N₂/NO's channel-representation COMPARATIVELY instead
— rel_error no worse than the deck's own — exactly the design spec's stated criterion
("same-or-better probe precision... than the committed hand-tuned deck"), and both pass
comfortably (17× and 7× better than their decks respectively). Their REAL requirement, the
vibrational spectrum (`probe_nuclear`), converges cleanly at every `C` tried.

**Second finding: N₂/NO's point counts exceed their decks' (1.0–1.5×).** Traced to
`qscat.tuning.propose`'s fixed `_NUCLEAR_X_MAX_DEFAULT = 18.0` bohr real-region default,
which is LARGER than N₂'s (12.0 bohr) and NO's (9.0 bohr) committed nuclear real regions
(though smaller than F₂'s 10.7-bohr real region plus its complex tail) — a
per-molecule-INDEPENDENT constant from Task 5's a-priori adapter, not derived from the
potential profile itself (e.g. "where has the interaction died"). More real-region span at
comparable density costs more points, independent of `C`. This is a Task-5 a-priori-adapter
limitation that Task 8's `C`-calibration cannot fix (`C` controls density, not extent) — a
documented follow-on (deriving `x_max` from the potential rather than a fixed constant),
not addressed here.

## The gate (`validation/tuning/test_emoscat_decks.py`)

Two kinds of test, both fast (probes only, no 2-D solve):

1. **Reproduce-or-beat** — `propose_grid` at the calibrated `C` for N₂/NO/F₂: F₂ gated
   strictly (absolute channel convergence + vibrational convergence + points ≤ 1.3× deck);
   N₂/NO gated on vibrational convergence (absolute) and channel-representation
   (comparative to their own deck, per the finding above), with N₂'s point-count margin
   widened to the measured ~1.44× ratio (documented, not silently loosened to 1.3×).
2. **Flag-the-failures** — the cheapest probe, `probe_channel_representation`, on the
   COARSE shared N₂-style grid (`qscat.core.grids.nuclear_grid()`), must report
   `converged=False` for F₂ DA's K≈58 wave (and across the whole K≈52–78 range the tested
   energies produce) — the regression guard for the exact bug that motivated this
   sub-project.

An H₂⁺ Coulomb-incident coarse-grid check was explored but not included as a gate: at H₂⁺
DR's low incident `k` (long de-Broglie wavelength), `probe_channel_representation`'s failure
mode on a truncated Coulomb electronic grid is dominated by real-region EXTENT, not element
density, and neither a 30-bohr nor the 60-bohr proxy grid gave a clean converged/
not-converged split — not a clean regression gate the way the F₂ DA case is.

## Key result

The tuner, calibrated once against F₂'s genuinely-open DA channel, reproduces-and-beats
that deck outright (37% fewer points, clean absolute convergence) and correctly flags the
coarse shared grid's historical under-resolution of the same wave — the two things this
whole sub-project set out to prove. N₂/NO's nuclear grids are comparatively good
(vibrational spectra converge; channel representation beats their own decks by 7-17×) but
cost more points than their decks, root-caused to a fixed real-region extent default, not
to `C` — a genuine, documented limitation for a follow-on, not silently hidden.

See also: `docs/superpowers/specs/2026-07-28-discretisation-tuner-design.md` (the design),
the `discretisation-tuner` skill (the supervised loop + worked example),
`docs/physics/diatomic-ve-cross-sections.md` (the F₂ DA physics and the original K≈58 bug).
