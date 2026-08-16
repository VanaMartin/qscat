"""The general iterative 2-D-convergence fallback: given ANY scalar observable
that depends on a pair of `(electronic, nuclear)` FEM-DVR-ECS grids, refine
whichever coordinate is under-resolved until the observable stops moving.

`propose_grid`'s a-priori adiabatic mesh is built from `v0`'s classical
wavenumber profile alone, so it can miss structure that lives elsewhere in the
model (e.g. a narrow feature in `v_int`/`lambda(R)` -- see
`validation/tuning/test_emoscat_decks.py`'s F2 dissociative-attachment 2-D
spot-check). `refine_to_2d_convergence` is the general, model-agnostic
fallback for exactly that failure mode: it does not know or care what the
observable computes (`da_cross_section`, `ve_cross_section`, `dr_cross_section`,
or a synthetic test function) -- it just probes which of the two grids moves
the observable more under one `refine` step, adopts that refinement, and
repeats until the larger of the two relative moves is below `rtol` or
`max_iter` adopted steps have been taken.

- `refine_to_2d_convergence` -- the loop. Each iteration evaluates the
  observable on a once-`refine`d nuclear variant and a once-`refine`d
  electronic variant (holding the other grid fixed), adopts whichever moves
  the observable more (the larger `|Delta value| / |value|`), and records
  that step. Stops -- `converged=True` -- when the larger of the two relative
  moves is `< rtol` (no step is recorded for that final, non-adopted check).
  Stops -- `converged=False` -- after `max_iter` adopted steps without ever
  meeting `rtol`: a real, reported signal that the observable has not settled
  in the refinement budget, not a silent cap.

This module imports only `qscat.tuning.refine` and `qscat.dvr.FemDvrEcsGrid`
-- no model, no specific observable -- so it stays usable for any 2-D
cross-section (or any other coordinate-pair-dependent scalar) the caller
closes over.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

from qscat.dvr import FemDvrEcsGrid

from .probes import refine

__all__ = ["refine_to_2d_convergence"]

# Floor for the relative-move denominator, avoiding a division by ~0 when the
# observable's current value is itself (numerically) zero.
_REL_FLOOR = 1e-30


def refine_to_2d_convergence(
    observable: Callable[[FemDvrEcsGrid, FemDvrEcsGrid], float],
    g_r: FemDvrEcsGrid,
    g_R: FemDvrEcsGrid,
    *,
    rtol: float = 1e-2,
    max_iter: int = 4,
) -> tuple[FemDvrEcsGrid, FemDvrEcsGrid, dict[str, Any]]:
    """Iteratively refine whichever of `g_r` (electronic) / `g_R` (nuclear)
    grid is under-resolved for `observable`, until it stops moving.

    `observable(g_r, g_R) -> float` is an arbitrary caller-supplied closure
    (e.g. a real cross-section at a fixed energy, or a synthetic test
    function) -- this loop does not inspect what it computes.

    Each of up to `max_iter` iterations:
      1. Evaluates `observable` on a once-`refine`d nuclear variant
         (`g_r, refine(g_R)`) and a once-`refine`d electronic variant
         (`refine(g_r), g_R`), each against the CURRENT (adopted) grids.
      2. Compares each candidate's relative move `|v - current| /
         max(|current|, tiny)` off the current adopted value.
      3. If the larger of the two relative moves is `< rtol`: STOP,
         `converged=True`, and do NOT record a step for this check.
      4. Otherwise ADOPT the coordinate with the larger relative move
         (replace that grid with its refinement, update the current value),
         and record `{"coordinate": ..., "value": ..., "rel_move": ...}`.

    If the loop exhausts `max_iter` adopted steps without ever satisfying
    step 3, `converged=False` -- a genuine "did not settle" signal, not a
    silently accepted best-effort result.

    Returns `(g_r, g_R, detail)`: the (possibly refined) grid pair, and
    `detail` with exactly `"converged"` (bool), `"iterations"` (the list of
    adopted-step dicts, in order), and `"final_value"` (the last adopted
    value, or the initial `observable(g_r, g_R)` if no step was adopted).
    """
    current = observable(g_r, g_R)
    final_value = current
    iterations: list[dict[str, Any]] = []
    converged = False

    for _ in range(max_iter):
        v_nuc = observable(g_r, refine(g_R))
        v_elec = observable(refine(g_r), g_R)

        denom = max(abs(current), _REL_FLOOR)
        rel_nuc = abs(v_nuc - current) / denom
        rel_elec = abs(v_elec - current) / denom

        if max(rel_nuc, rel_elec) < rtol:
            converged = True
            break

        if rel_nuc >= rel_elec:
            g_R = refine(g_R)
            current = v_nuc
            iterations.append({"coordinate": "nuclear", "value": v_nuc, "rel_move": rel_nuc})
        else:
            g_r = refine(g_r)
            current = v_elec
            iterations.append({"coordinate": "electronic", "value": v_elec, "rel_move": rel_elec})
        final_value = current

    detail: dict[str, Any] = {
        "converged": converged,
        "iterations": iterations,
        "final_value": final_value,
    }
    return g_r, g_R, detail
