"""The full sigma(E) "boomerang" curve from ONE propagation, + the usable
energy window (sub-project #7, Task 5).

Task 4 (`.superpowers/sdd/task-4-report.md`) already did the convergence
work: the T-scan (E=0.10 settles to 0.93-0.95 for T>=1200; E=0.15 floors at
[1.10, 1.20] across the whole T range -- the usable-window edge, not an
unconverged transient), the wavepacket tuning (`p0` centers the incident
spectral peak between the two TI anchors), and the `F_out = hankel/2` physics
fact (the outgoing test function MUST project onto the outgoing Hankel half,
not the regular free function -- a five-order-of-magnitude effect, see the
Task 4 report). This module does not repeat that sweep; it names Task 4's
converged configuration as `TD_WORKING_GRID` and adds the one genuinely new
capability Task 5 asks for: the WHOLE sigma(E) curve from a single stored
`c(t)` trajectory (the "boomerang"), plus an honest usable-energy window
built from `|eta_incident(E)|`.

`c_{v'}(t)` does not depend on E, so `sigma_curve` is a thin wrapper around
`td_cross_section.td_ve_cross_section_2d` called with the whole `E_grid` as
its array-`E` argument: that public function already propagates ONCE and
then transforms the SAME stored trajectory at every requested energy --
the expensive part (the Crank-Nicolson propagation, ~250s at
`TD_WORKING_GRID`) happens exactly once regardless of how dense `E_grid` is.
This is TD's structural advantage over the TI solver, which needs one sparse
LU factorization per energy. `sigma_curve` does not reimplement any of this
-- it exists only to supply `TD_WORKING_GRID`'s defaults for the keyword-only
propagation parameters when the caller leaves them `None`.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from qscat.dvr import FemDvrEcsGrid, TensorGrid

from projects.n2_2d_cross_section.electronic_grid import n2_electronic_grid
from projects.n2_2d_cross_section.hamiltonian2d import ELL
from projects.n2_2d_td_cross_section import td_cross_section as td
from projects.n2_2d_td_cross_section.correlation import eta_incident
from projects.n2_ti_cross_section.nuclear_grid import n2_nuclear_grid

__all__ = ["TD_WORKING_GRID", "td_working_tgrid", "sigma_curve", "usable_window"]

# The converged configuration from Task 4 (`.superpowers/sdd/task-4-report.md`,
# "Converged configuration" + the T-scan table). Task 5 does not re-derive
# this -- it only names it, per the controller's explicit instruction not to
# redo a box/dt/T sweep. Every value is commented with its Task-4 measurement
# or justification.
TD_WORKING_GRID: dict = {
    "electronic": {
        # r_max=50 (== the ECS pivot R0) is comfortably beyond wp_in's
        # r0=25 + several sigma=5 (so the packet starts, and the resonance
        # forms, well inside the real region) and beyond wp_out's
        # r0_out=35 + sigma_out=4 (so the outgoing test function also fits
        # inside the real region). This is why the TD box is larger than
        # #6's TI working grid (r_max=16): the TD box must hold a moving
        # wavepacket, not just a static scattering solution.
        #
        # Honest caveat: r_max=50 was sized by the physical reasoning above
        # (wp_in/wp_out both fit comfortably inside the real region, and the
        # interaction V_int = -lambda(R)*exp(-alpha_c*r^2) vanishes by
        # r~5-6), and is cross-checked only INDIRECTLY, via the TD-vs-TI
        # agreement this module and `test_td_convergence.py` measure across
        # the usable window. It was NOT subjected to a direct empirical
        # r_max-convergence sweep (re-running at, say, r_max=40/50/60 and
        # checking sigma_TD is unchanged) -- that sweep is future work, out
        # of scope here.
        "r_max": 50.0,
        "order": 8,       # unchanged from #6's TI working grid.
        "n_complex": 6,   # unchanged from #6's TI working grid.
    },
    "nuclear": {
        # Matches Task 1-4's test-scale nuclear grid (same as #6's TI
        # working grid, `quadrature=10, r_max=22, n_complex=5`); the
        # vibrational bound states (R ~ 1.5-3 bohr) sit entirely inside the
        # real region, so nothing about the TD propagation needs a bigger
        # nuclear box.
        "quadrature": 10,
        "r_max": 22.0,
        "n_complex": 5,
    },
    # dt=0.5, n_steps=3000 (T=1500 a.u.): Task 4's T-scan shows
    # sigma_TD/sigma_TI at E=0.10 stops drifting with T once T>=1200
    # (0.950, 0.931, 0.945 for T=1200/1500/1800); T=1500 is the value Task 4
    # settled on. By T=1500, ||Psi|| has decayed 1.0 -> 0.024 (the resonance
    # fully depletes -- the physical "formation and decay" this method is
    # built to show).
    "dt": 0.5,
    "n_steps": 3000,
    "wp_in": {
        # r0=25: launched well inside the box (see r_max's comment above),
        # far enough out that dt=0.5/n_steps=3000 is enough time to travel
        # in, interact, and let the resonance decay.
        "r0": 25.0,
        # p0=-0.5 (inward): p0**2/2 = 0.125 Ha sits BETWEEN the two TI
        # anchor energies (0.10, 0.15), so the incident spectral weight
        # |eta_incident(E)| is large and roughly balanced at both:
        # measured |eta_in(0.10)| = 2.7034, |eta_in(0.15)| = 2.6395
        # (Task 4 report).
        "p0": -0.5,
        # sigma=5.0: narrow enough to keep the spectrum concentrated near
        # the two anchors without smearing across the whole resonance
        # structure, wide enough that both anchors sit well inside the
        # usable window (see `usable_window` below).
        "sigma": 5.0,
    },
    "wp_out": {
        # Scaled down from eMoScat's production r0=75 to fit this box;
        # still well outside the interaction range (V_int = -lambda(R)
        # exp(-alpha_c r^2) vanishes by r~5-6 since alpha_c=0.4), so the
        # outgoing test function only picks up genuine outgoing flux.
        "r0_out": 35.0,
        "p0_out": 0.5,
        "sigma_out": 4.0,
    },
}


def td_working_tgrid() -> TensorGrid:
    """The (electronic x nuclear) `TensorGrid` for `TD_WORKING_GRID`."""
    eg = TD_WORKING_GRID["electronic"]
    ng = TD_WORKING_GRID["nuclear"]
    return TensorGrid(
        [
            n2_electronic_grid(r_max=eg["r_max"], order=eg["order"], n_complex=eg["n_complex"]),
            n2_nuclear_grid(
                quadrature=ng["quadrature"], r_max=ng["r_max"], n_complex=ng["n_complex"]
            ),
        ]
    )


def sigma_curve(
    tgrid: TensorGrid,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    vprimes: list[int],
    E_grid: npt.ArrayLike,
    *,
    dt: float | None = None,
    n_steps: int | None = None,
    wp_in: dict[str, float] | None = None,
    wp_out: dict[str, float] | None = None,
) -> npt.NDArray[np.float64]:
    """The whole `sigma_{v_init->v'}(E)` curve (bohr^2) from ONE propagation.

    Thin wrapper around `td_cross_section.td_ve_cross_section_2d` -- that
    public function already does exactly this: `c_{v'}(t)` does not depend on
    E, so it runs the Crank-Nicolson propagation exactly once and then
    transforms that SAME stored trajectory once per energy in the (array)
    `E` it is given. However dense `E_grid` is, only one propagation happens
    -- the "boomerang" curve is free once that one run (~250s at
    `TD_WORKING_GRID`) is done. This function adds nothing computationally;
    it only supplies `TD_WORKING_GRID`'s defaults for the keyword-only
    propagation parameters when the caller leaves them `None`.

    Any of `dt`/`n_steps`/`wp_in`/`wp_out` left `None` default to
    `TD_WORKING_GRID`'s values.

    Returns shape `(len(E_grid), len(vprimes))`, matching
    `projects.n2_2d_cross_section.cross_section_2d.ve_cross_section_2d`'s
    array-E convention (so the two curves overlay directly) -- guaranteed
    here because `E_grid` is always passed through as an array, which is the
    branch of `td_ve_cross_section_2d` that returns that shape.
    """
    dt = TD_WORKING_GRID["dt"] if dt is None else dt
    n_steps = TD_WORKING_GRID["n_steps"] if n_steps is None else n_steps
    wp_in = TD_WORKING_GRID["wp_in"] if wp_in is None else wp_in
    wp_out = TD_WORKING_GRID["wp_out"] if wp_out is None else wp_out

    e_arr = np.atleast_1d(np.asarray(E_grid, dtype=np.float64))
    return td.td_ve_cross_section_2d(
        tgrid, eps, chi, v_init, vprimes, e_arr, dt=dt, n_steps=n_steps, wp_in=wp_in, wp_out=wp_out
    )


def usable_window(
    grid: FemDvrEcsGrid,
    E_grid: npt.ArrayLike,
    *,
    wp_in: dict[str, float] | None = None,
    l: int = ELL,
    frac: float = 0.5,
) -> tuple[tuple[float, float], npt.NDArray[np.float64]]:
    """The `(E_lo, E_hi)` sub-interval of `E_grid` where `|eta_incident(E)|`
    is at least `frac` of its peak over the grid, plus the full
    `|eta_incident(E)|` array (so a figure or test can plot/reuse it without
    recomputing).

    This is the honest usable window: outside it, the Tannor-Weeks
    deconvolution `1/eta_incident(E)` amplifies whatever residual
    truncation/discretization noise is in the stored `c(t)` -- exactly the
    effect the Task 4 report documents (E=0.15 sits farther from the
    incident wavepacket's spectral peak `p0**2/2=0.125` Ha than E=0.10 does,
    and its sigma_TD/sigma_TI ratio is correspondingly worse and does not
    shrink with propagation length T).

    `grid` is the ELECTRONIC `FemDvrEcsGrid` (i.e. `tgrid.grids[0]`),
    matching `correlation.eta_incident`'s signature. `E <= 0` entries get
    `|eta_incident| = 0` (no incident channel below threshold) rather than
    being evaluated at an imaginary `k`.

    Returns `((E_lo, E_hi), eta_abs)`. Raises `ValueError` if no energy in
    `E_grid` meets the threshold (an empty window is a configuration error,
    not a silent no-op).
    """
    wp_in = TD_WORKING_GRID["wp_in"] if wp_in is None else wp_in
    e_arr = np.atleast_1d(np.asarray(E_grid, dtype=np.float64))
    eta_abs = np.array(
        [
            abs(eta_incident(grid, float(np.sqrt(2.0 * e)), l, **wp_in)) if e > 0.0 else 0.0
            for e in e_arr
        ]
    )
    peak = eta_abs.max()
    idx = np.flatnonzero(eta_abs >= frac * peak)
    if idx.size == 0:
        raise ValueError("usable_window: no energy in E_grid meets frac*peak threshold")
    return (float(e_arr[idx[0]]), float(e_arr[idx[-1]])), eta_abs
