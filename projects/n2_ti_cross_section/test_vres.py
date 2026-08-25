"""Failing-first tests for V_d(R)/Gamma(R) recomputed per nuclear-R
(sub-project #3, Task 2).

Cross-checks `vres.vres_on_grid` against sub-project #2's own validated
pole-finder result (the development notes for
`projects/n2_resonance`): at R0 = 2.01943 bohr,
`E_pole = 0.089850 - 0.008363i Ha`, i.e. `E_res ~= 0.0898 Ha`,
`Gamma = -2*Im(E_pole) ~= 0.01673 Ha`. The nuclear grid's nearest real point
to R0 is a few thousandths of a bohr off (grid nodes don't land exactly on
R0), so a several-percent tolerance is used rather than exact equality.
"""

from __future__ import annotations

import functools

import numpy as np

from projects.n2_resonance.potential import v0
from projects.n2_ti_cross_section.nuclear_grid import n2_nuclear_grid
from projects.n2_ti_cross_section.vres import vres_on_grid

R0 = 2.01943
E_RES_R0 = 0.0898  # Ha, sub-project #2 result
GAMMA_R0 = 0.0167  # Ha, sub-project #2 result (2 * 0.00836)


@functools.cache
def _vres():
    """The nuclear grid and ONE pole walk over it, shared by all three tests.

    `vres_on_grid` re-solves the electronic eigenproblem at every nuclear
    node, which measured ~8 s; all three tests below want the same curve and
    none of them mutate it, so it is walked once rather than three times (the
    repeated-recomputation anti-pattern docs/adr/0005 names). Read-only by
    contract: copy before modifying `Vd`/`Gamma`.
    """
    grid = n2_nuclear_grid()
    Vd, Gamma = vres_on_grid(grid)
    return grid, Vd, Gamma


def test_vres_shapes_finite_and_gamma_nonnegative():
    grid, Vd, Gamma = _vres()
    assert Vd.shape == (grid.n,)
    assert Gamma.shape == (grid.n,)
    assert np.all(np.isfinite(Vd.real)) and np.all(np.isfinite(Vd.imag))
    assert np.all(np.isfinite(Gamma))
    assert np.all(Gamma >= -1e-12)


def test_vres_matches_pole_finder_at_R0():
    grid, Vd, Gamma = _vres()

    real_mask = grid.points.imag == 0.0
    real_idx = np.flatnonzero(real_mask)
    idx = real_idx[np.argmin(np.abs(grid.points[real_idx].real - R0))]
    R = grid.points[idx].real

    E_res = Vd[idx].real - v0(R)
    assert abs(E_res - E_RES_R0) / E_RES_R0 < 0.05, E_res
    assert abs(Gamma[idx] - GAMMA_R0) / GAMMA_R0 < 0.10, Gamma[idx]


def test_gamma_closes_beyond_crossing_and_on_complex_tail():
    grid, Vd, Gamma = _vres()

    real_mask = grid.points.imag == 0.0
    real_idx = np.flatnonzero(real_mask)
    real_R = grid.points[real_idx].real

    # Beyond the ~2.4 bohr resonance -> bound-state crossing, Gamma should
    # have closed to ~0.
    outer = real_idx[real_R > 4.0]
    assert np.all(Gamma[outer] < 1e-4), Gamma[outer]

    # Complex-tail points (R > 12 bohr, ECS-rotated): Gamma == 0 by the
    # documented analytic-continuation treatment (see vres.py docstring).
    tail_idx = np.flatnonzero(~real_mask)
    assert tail_idx.size > 0
    assert np.all(Gamma[tail_idx] == 0.0)

    # Vd must be finite (and equal v0(R) + a constant E_res asymptote) on
    # the tail.
    assert np.all(np.isfinite(Vd[tail_idx].real))
    assert np.all(np.isfinite(Vd[tail_idx].imag))
