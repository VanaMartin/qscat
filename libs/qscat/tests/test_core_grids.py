"""Physics/layout gates for the promoted `qscat.core.grids`/`qscat.core.vibrational`
(sub-project #A, Task 3).

These are NOT a self-comparison against the shims that now delegate to
`qscat.core` (that would be tautological, since the shims just call this
same code) -- they pin the FEM-DVR-ECS element layout by concrete invariant
(point count, ECS pivot `R0`, weight sum, known element-boundary points) and
gate `vibrational_states` against the documented N2 physics (the known
eps1-eps0 =~ 0.0124 Ha vibrational spacing, see
`projects/n2_ti_cross_section/test_vibrational.py`).
"""

from __future__ import annotations

import numpy as np
from qscat.core.grids import electronic_grid, nuclear_grid
from qscat.core.vibrational import vibrational_states
from qscat.model import N2

MU_N2 = 12766.36


def test_electronic_grid_default_layout() -> None:
    g = electronic_grid()
    assert g.n == 167
    assert g.R0 == 30.0  # ECS pivot == r_max by construction
    w = np.sum(g.weights)
    assert np.isclose(w, 78.2936269875353 + 33.81806240536799j, atol=1e-9)
    # Known inner-segment element boundaries (see `_INNER_SEGMENTS`).
    assert np.any(np.isclose(g.points.real, 1.0))
    assert np.any(np.isclose(g.points.real, 5.0))
    assert np.any(np.isclose(g.points.real, 10.0))


def test_nuclear_grid_default_layout() -> None:
    g = nuclear_grid()
    assert g.n == 428
    assert g.R0 == 12.0  # ECS pivot == 12.0 bohr by construction
    w = np.sum(g.weights)
    assert np.isclose(w, 34.92090764820161 + 16.05131596496235j, atol=1e-9)
    # Known real-segment element boundaries (see `_REAL_SEGMENTS`).
    assert np.any(np.isclose(g.points.real, 1.5))
    assert np.any(np.isclose(g.points.real, 3.0))
    assert np.any(np.isclose(g.points.real, 4.0))
    assert np.any(np.isclose(g.points.real, 12.0))


def test_vibrational_states_reproduces_the_n2_vibrational_spacing() -> None:
    grid_R = nuclear_grid()
    eps, chi = vibrational_states(grid_R, MU_N2, 4, N2.v0)

    assert eps.shape == (4,)
    assert chi.shape == (4, grid_R.n)
    assert np.all(np.diff(eps) > 0.0)

    # Documented N2 spacing (see test_vibrational.py's analytic Morse gate).
    assert np.isclose(eps[1] - eps[0], 0.0124, atol=1e-4)

    # Bound eigenvectors are real (to round-off) and Hermitian-normalized.
    assert np.max(np.abs(chi.imag)) < 1e-10
    norms = np.sum(np.abs(chi) ** 2, axis=1)
    assert np.allclose(norms, 1.0, atol=1e-10)
