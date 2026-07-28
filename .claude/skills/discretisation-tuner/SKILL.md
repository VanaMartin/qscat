---
name: discretisation-tuner
description: Use when a qModeling calculation needs a FEM-DVR-ECS grid — supervises the qscat.tuning loop (analyze the potential → propose an adaptive grid → probe convergence at the energy extremes → refine/coarsen → 2-D spot-check → emit the minimal-cost grid at a target precision) instead of hand-picking element lengths by eye.
---

# discretisation-tuner

## Overview

Hand-tuned FEM-DVR-ECS grids have been the single most expensive class of bug in this
repo — a coarse shared nuclear grid under-resolved the K≈58 dissociative-attachment wave
(σ off by ~36 orders); the H₂⁺ Coulomb tail needed 1300 bohr. This skill replaces the
human "good eye" for element lengths with the `qscat.tuning` primitives: it computes the
**minimal-DVR-point grid that holds a target precision** for a given model, coordinate, and
energy range — an adaptive equidistribution mesh (each element carrying a ~constant de
Broglie phase), the h/p-optimal quadrature, and a double-ECS-safe absorbing tail.

You (the supervisor) run the deterministic primitives, read the convergence probes, and make
the judgment calls (which knob to move, when to stop, when to defer the 2-D check to Docker).
The physics/numerics live in `qscat.tuning`; the loop + judgment live here.

## When to Use

- Setting up (or distrusting) a grid for any qscat calculation — a new molecule, a new energy
  range, a new observable (VE / DA / DR), or a suspected under-resolution.
- Before committing a per-molecule grid deck to `validation/`.
- Diagnosing "the answer changes when I refine the grid" — the probes localize which coordinate
  and which knob is under-resolved.

Do NOT use it to re-grid a validated production deck without an explicit decision — the tuner
EMITS a config; adopting it is a separate, opt-in choice.

## What `qscat.tuning` gives you (the primitives you call)

```python
from qscat.tuning import (
    analyze_potential, PotentialProfile,        # V(x) -> local-wavenumber profile
    equidistribution_elements, optimal_real_mesh,  # adaptive mesh + h/p sweep
    max_stable_angle, tune_ecs_tail,            # ECS-tail (double-ECS-capped angle + exp absorption)
    refine, probe_nuclear, probe_electronic, probe_channel_representation,  # convergence probes
    grid_cost, tensor_cost, propose_grid,       # cost model + one-shot a-priori grid
    IncidentSpec, required_extent, tw_analysis,  # incident/test-function placement
    refine_to_2d_convergence,                   # general iterative 2-D-convergence fallback
)
```

`propose_grid(model, coordinate, energy_range, *, rtol=1e-3, incident=None)` is the a-priori
half of the hybrid — it already runs analyze → mesh → ECS and returns a `FemDvrEcsGrid`. Your
job is to VALIDATE and MINIMISE it with the probes.

## The tuning loop (the procedure)

Create a todo per step.

1. **Frame the problem.** Fix the `model`, the target `energy_range = (E_min, E_max)`, the target
   `rtol` (default 1e-3), and the observable (VE/DA/DR). Note whether it's the TI route (incident =
   channel function) or the TD route (incident = a Gaussian wavepacket).

2. **Incident / test-function placement (TW analysis) — if TD.** Either take a caller-supplied
   `IncidentSpec`, or `incident = tw_analysis(model, energy_range)` to auto-place the wavepacket
   (position/impulse/σ) so its spectrum spans the range. The incident drives BOTH the real-region
   EXTENT (`required_extent`) AND the RESOLUTION (its energy `impulse²/2` raises the mesh's effective
   `E_max`). For TI, `incident=None`.

3. **Propose the a-priori grid, per coordinate.** For each of `"nuclear"` and `"electronic"`:
   `g = propose_grid(model, coordinate, energy_range, rtol=rtol, incident=incident)`.

4. **Probe convergence at the EXTREMES.** The finest requirement is at `E_max`; the longest wave /
   largest extent is near-threshold `E_min`. At each extreme:
   - Nuclear: `probe_nuclear(model, g_R, n_vib, rtol=rtol)` (vibrational eigenvalues stable under
     one `refine`).
   - Electronic: `probe_electronic(model, g_r, R_eq, window=..., rtol=rtol)` (bound/resonance energy
     stable).
   - **Channel representation — the cheap, decisive one:** `probe_channel_representation(g, k, l,
     charge=model.charge, mass=..., rtol=rtol)` where `k` is the largest channel wavenumber the
     observable needs — the incident `k=√(2E_max)`, and for a dissociation channel the OUTGOING
     `K=√(2μ·E_DR_max)` (heavy → large; use `E_DR = E_max − ε_threshold` for an exothermic channel,
     NOT just `E_max`). This is the probe that catches the K≈58-under-resolution failures.
   Read each `ProbeResult.converged` and `.cost`.

5. **Refine / coarsen to the minimum.** For any probe with `converged == False`, the grid is
   under-resolved there — `refine` that coordinate (h) or accept the higher order the h/p sweep
   picks (p), re-propose/re-probe. For any knob COMFORTABLY over-converged (Δ ≪ rtol), try a coarser
   variant (fewer elements / lower order / smaller real extent / shorter tail) and re-probe; keep
   the coarsest variant that still holds every probe at `rtol`. **Watch the cost with `grid_cost` /
   `tensor_cost` (n_unknowns = n_r·n_R, and the est. factor memory/time) — those are relative-ranking
   estimates; pick the fewest-unknowns grid that passes.** Note: DVR point count is NON-MONOTONE
   under the h/p sweep (a higher order can win with fewer, denser elements), so rank by the probe +
   `tensor_cost`, not by element count alone.

6. **Final 2-D spot-check.** Build `TensorGrid([g_r, g_R])` and run the ACTUAL observable
   (`ve_cross_section` / `da_cross_section` / `dr_cross_section`) at the HARDEST energy, and confirm
   it agrees with a once-refined grid to `rtol`. For a non-laptop deck (H₂⁺-scale), run this on a
   reduced proxy or under Docker/MUMPS and SAY SO — do not silently skip it. Instead of a single
   manual once-refined comparison, you may call the general fallback directly: define `observable`
   as a closure over the real cross-section at that hardest energy (`observable(g_r, g_R) ->
   float`) and call `refine_to_2d_convergence(observable, g_r, g_R, rtol=..., max_iter=...)` — it
   iterates, adopting whichever coordinate (nuclear or electronic) moves the observable more, until
   both relative moves are under `rtol` or `max_iter` is hit. If `detail["iterations"]` is
   non-empty, the a-priori grid from step 3 was under-resolved in 2-D even though it may have
   passed every 1-D probe (exactly the F2 DA finding above) — report the converged grid pair
   in place of the original, plus the cost delta (`grid_cost`/`tensor_cost` before vs after).

7. **Emit the config + report.** Output the per-coordinate grid (the `ElementSpec` lists / the built
   `FemDvrEcsGrid`, expressible as a committed deck) and a report:
   - the achieved precision (each probe's convergence number),
   - the cost (`grid_cost`/`tensor_cost`) vs the previous/hand grid,
   - the tuning decisions (which knobs moved and why),
   - any deferrals (the 2-D check on Docker; an uncalibrated constant).

## Stop criteria

STOP and emit when **every** probe (both extremes, both coordinates, the channel wavenumber) holds
`rtol` under one refinement AND no single knob can be coarsened without breaking a probe. STOP and
ESCALATE if a probe cannot reach `rtol` at any feasible grid (a genuine finding — report the numbers,
don't loosen `rtol`), or if the ECS angle the potential allows can't absorb the fastest wave (report
it — the model may need a different tail representation).

## Worked example — the N₂ nuclear grid

```python
from qscat.model import N2
from qscat.tuning import propose_grid, probe_nuclear, probe_channel_representation, refine, grid_cost
g = propose_grid(N2, "nuclear", (0.04, 0.18))          # a-priori adaptive grid, calibrated C
pn = probe_nuclear(N2, g, n_vib=3)                       # vibrational eps stable?  -> converged
K = (2 * N2.mu * 0.18) ** 0.5                            # the fastest nuclear wave in-range
pc = probe_channel_representation(g, K, 0, mass=N2.mu)   # is that wave resolved?
# pn.converged is True. pc.converged is NOT (rel_error ~1.7e-3, just over rtol=1e-3) -- but this K
# is a CONSERVATIVE FLOOR (N2's DA channel is closed in-window, so there is no genuinely open
# dissociation wave to represent; sqrt(2*mu*E_max) over-estimates it), and neither does N2's own
# committed deck at this K (rel_error ~2.9e-2) -- see `validation.tuning.calibrate`'s Task 8 report.
# The comparative bar (rel_error no worse than the deck's own) IS met, by >10x -- see
# `validation/tuning/test_emoscat_decks.py`. For a molecule with a genuinely OPEN dissociation
# channel in-range (F2, K derived from its anion bound-state threshold -- `eps_e`), pc.converged
# IS the right absolute bar, and it holds.
```
Then the 2-D spot-check: `ve_cross_section(TensorGrid([elec_grid, g]), N2, ...)` at E=0.10 vs a
refined grid, agreeing to 1e-3.

## Notes

- The de Broglie phase-per-element constant `C` is calibrated (`qscat.tuning.mesh._PHASE_COEFF =
  0.10`) so the tuner reproduces-or-beats the eMoScat F2 dissociative-attachment deck -- the
  molecule with a genuinely open channel in its tested range, and the decisive calibration case;
  see `validation/tuning/calibrate.py` and `docs/physics/discretisation-tuning.md`. Trust it, don't
  re-tune it by eye. N2/NO's proposed nuclear grids cost MORE points than their committed decks
  (~1.0-1.5x) -- a `_NUCLEAR_X_MAX_DEFAULT` (fixed real-region extent) mismatch with their
  decks' own extent, not a C problem; a documented follow-on, not (yet) fixed. H2+'s proxy
  nuclear deck is a clean reproduce-and-beat (unlike N2/NO): lighter mu keeps its floor K
  modest and its deck's real-region extent is close to the fixed default.
- **Step 6 (the final 2-D spot-check) is NOT optional/a rubber stamp -- concrete proof, not just
  a warning.** `validation/tuning/test_emoscat_decks.py::test_f2_2d_da_cross_section_spot_check`
  found that F2's reproduce-and-beat NUCLEAR grid (both 1-D probes pass, fewer points than the
  deck) gives an sigma_DA that is NOT 2-D converged -- one nuclear h-refinement moves it ~5x
  (toward the eMoScat deck's own value), traced to a narrow R~2.5-2.7 bohr interaction feature
  (in `v_int`/`lambda(R)`, not `v0`) the a-priori mesh cannot see since it is built only from
  `v0`'s classical k(x) profile. The 1-D probes are NECESSARY but NOT SUFFICIENT here -- always
  run step 6 for a new grid before trusting it, even when every 1-D probe passed.
- `qscat.tuning` primitives are pure/deterministic and unit-tested on analytic potentials; the
  judgment (this procedure) is the only non-deterministic part.
- See `docs/superpowers/specs/2026-07-28-discretisation-tuner-design.md` and
  `docs/physics/discretisation-tuning.md`.
