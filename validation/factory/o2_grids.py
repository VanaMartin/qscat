"""O2's FEM-DVR-ECS grids from the discretisation tuner, with the 1-D probes.

`python -m validation.factory.o2_grids [--e-max 0.10] [--spot-check]`

Builds `propose_grid(O2, ...)` for the electronic and nuclear coordinates
over the VE energy window (Alt & Houfek's Fig. 5 spans 0-2.5 eV), runs the
three 1-D probes at the window's top, and prints each grid as the
`(n_elements, endpoint)` deck `qscat.core.grids.segmented_grid` rebuilds --
the form `apps/qscat-run`'s O2 preset carries. `--spot-check` adds the 2-D
test: the exact TI `sigma_{0->1}` at one energy on the proposed pair and on
each once-refined variant (`refine_to_2d_convergence`), the step the 1-D
probes cannot replace (F2's tuner grids passed both probes and were still
5x off in sigma_DA).
"""

from __future__ import annotations

import argparse
import math
import time

import numpy as np
from qscat.core.driven import ve_cross_section
from qscat.core.grids import segmented_grid
from qscat.core.vibrational import vibrational_states
from qscat.dvr import FemDvrEcsGrid, TensorGrid
from qscat.model import O2
from qscat.tuning import (
    probe_channel_representation,
    probe_electronic,
    probe_nuclear,
    propose_grid,
    refine,
    refine_to_2d_convergence,
)

__all__ = ["deck", "o2_grids", "o2_decks", "truncate_real", "NUCLEAR_R_CUT", "main"]

E_WINDOW = (0.002, 0.10)  # Ha: 0.05-2.7 eV, the paper's Fig. 5 window
N_VIB = 12


def deck(grid: FemDvrEcsGrid) -> dict[str, object]:
    """`grid` as `segmented_grid` arguments: runs of equal-length elements
    become `(n, endpoint)` segments, real and complex separately."""
    spec = grid.spec
    real: list[tuple[int, float]] = []
    cplx: list[tuple[int, float]] = []
    x = float(spec.x_min)
    angle = 0.0
    run_len, run_n, run_cplx = None, 0, False
    for e in spec.elements:
        length, is_cplx = float(e.length), e.angle_deg != 0.0
        if is_cplx:
            angle = float(e.angle_deg)
        if run_len is not None and abs(length - run_len) < 1e-9 and is_cplx == run_cplx:
            run_n += 1
        else:
            if run_len is not None:
                (cplx if run_cplx else real).append((run_n, x))
            run_len, run_n, run_cplx = length, 1, is_cplx
        x = x + length
    if run_len is not None:
        (cplx if run_cplx else real).append((run_n, x))
    return {
        "real_segments": tuple(real),
        "complex_segments": tuple(cplx),
        "angle_deg": angle,
        "quadrature": int(spec.quadrature),
        "x_min": float(spec.x_min),
        "n": int(grid.n),
    }


def o2_grids(e_window: tuple[float, float] = E_WINDOW) -> tuple[FemDvrEcsGrid, FemDvrEcsGrid]:
    """The tuner's a-priori pair, untouched."""
    g_e = propose_grid(O2, "electronic", e_window)
    g_n = propose_grid(O2, "nuclear", e_window, channel="ve")
    return g_e, g_n


def truncate_real(grid: FemDvrEcsGrid, R_cut: float) -> FemDvrEcsGrid:
    """`grid` with its real region cut at the first element boundary at or
    beyond `R_cut`, the same ECS tail (element lengths, angle) re-attached
    there. The tuner's VE nuclear grid tiles the real region out to a fixed
    18-bohr default with 0.16-bohr elements; for O2's VE window nothing
    lives past ~4 bohr (DA closed until 3.7 eV, the anion's outer turning
    point at 2.3 eV is 4.0 bohr), so 8 bohr keeps every level the probes
    converged and drops 62 elements of empty space."""
    spec = grid.spec
    real = [e for e in spec.elements if e.angle_deg == 0.0]
    tail = [e for e in spec.elements if e.angle_deg != 0.0]
    kept: list[tuple[int, float]] = []
    x = float(spec.x_min)
    for e in real:
        x += float(e.length)
        kept.append((1, x))
        if x >= R_cut - 1e-12:
            break
    cplx: list[tuple[int, float]] = []
    for e in tail:
        x += float(e.length)
        cplx.append((1, x))
    angle = float(tail[0].angle_deg) if tail else 0.0
    return segmented_grid(
        kept, cplx, angle_deg=angle, quadrature=int(spec.quadrature), x_min=float(spec.x_min)
    )


NUCLEAR_R_CUT = 8.0
# The 2-D spot check's verdict (`--spot-check`, sigma(0->1) at 0.05 Ha):
# ONE nuclear h-refinement of the tuner's grid moved the cross section by
# 69 %, after which the pair is converged (both further refinements < 2 %).
# The 1-D probe passed that grid at rtol 1e-3 -- which for a comb of 1-8 meV
# peaks is no gate at all: a 3 meV shift of a level is most of a width, and
# the cross section AT an energy follows the peak. So the deck carried is
# the tuner's nuclear mesh refined once (every element halved).
NUCLEAR_REFINEMENTS = 1


def o2_decks(e_window: tuple[float, float] = E_WINDOW) -> tuple[FemDvrEcsGrid, FemDvrEcsGrid]:
    """The pair the O2 preset carries: the tuner's electronic grid and its
    nuclear grid truncated at `NUCLEAR_R_CUT`, then h-refined
    `NUCLEAR_REFINEMENTS` times (see above)."""
    g_e, g_n = o2_grids(e_window)
    g_n = truncate_real(g_n, NUCLEAR_R_CUT)
    for _ in range(NUCLEAR_REFINEMENTS):
        g_n = refine(g_n)
    return g_e, g_n


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--e-max", type=float, default=E_WINDOW[1])
    ap.add_argument("--spot-check", action="store_true")
    ap.add_argument("--spot-energy", type=float, default=0.05)
    a = ap.parse_args()
    window = (E_WINDOW[0], a.e_max)

    t0 = time.perf_counter()
    g_e, g_n_full = o2_grids(window)
    g_n_cut = truncate_real(g_n_full, NUCLEAR_R_CUT)
    g_e, g_n = o2_decks(window)
    print(f"[O2] electronic: {deck(g_e)}")
    print(f"[O2] nuclear (tuner, full extent, n={g_n_full.n}): {deck(g_n_full)}")
    print(f"[O2] nuclear (truncated at {NUCLEAR_R_CUT} bohr, n={g_n_cut.n}): {deck(g_n_cut)}")
    print(f"[O2] nuclear (carried: refined x{NUCLEAR_REFINEMENTS}, n={g_n.n}): {deck(g_n)}")
    print(f"[O2] 2-D unknowns: {g_e.n * g_n.n} ({time.perf_counter() - t0:.1f} s)")

    k_e = math.sqrt(2.0 * a.e_max)
    ch = probe_channel_representation(g_e, k_e, O2.ell)
    el = probe_electronic(O2, g_e, 2.6, window=None)
    nu = probe_nuclear(O2, g_n, N_VIB)
    print(f"[O2] probe channel (k={k_e:.3f}, l={O2.ell}): {ch}")
    print(f"[O2] probe electronic (anion bound state at R=2.6): {el}")
    print(f"[O2] probe nuclear ({N_VIB} levels): {nu}")

    if a.spot_check:
        E = np.array([a.spot_energy])

        def sigma_01(gr: FemDvrEcsGrid, gR: FemDvrEcsGrid) -> float:
            eps, chi = vibrational_states(gR, O2.mu, 2, O2.v0)
            tg = TensorGrid([gr, gR])
            s = ve_cross_section(tg, O2, eps, chi, 0, [1], E)
            return float(np.asarray(s)[0, 0])

        t1 = time.perf_counter()
        gr, gR, info = refine_to_2d_convergence(sigma_01, g_e, g_n, rtol=2e-2, max_iter=2)
        dt = time.perf_counter() - t1
        print(f"[O2] spot check sigma(0->1, E={a.spot_energy}): {info} ({dt:.0f} s)")
        print(f"[O2] adopted: electronic n={gr.n}, nuclear n={gR.n}")


if __name__ == "__main__":
    main()
