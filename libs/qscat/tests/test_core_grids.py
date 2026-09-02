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
import pytest
from qscat.core.grids import electronic_grid, nuclear_grid
from qscat.core.vibrational import vibrational_states
from qscat.exceptions import GridError
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


# --- nuclear_grid public-parameter validation --------------------------------
#
# `r_max` and `n_complex` are public knobs, and the failure modes below are the
# ones a caller actually hits: `n_complex=0` used to divide by zero, a negative
# count silently dropped the ECS tail (leaving a purely real grid with no
# outgoing boundary condition), and a non-finite `r_max` built a grid of NaN
# nodes that only misbehaved much later. The gate is that each names the
# offending parameter here, at the public boundary.


@pytest.mark.parametrize("n_complex", [0, -3])
def test_nuclear_grid_rejects_non_positive_n_complex(n_complex: int) -> None:
    with pytest.raises(GridError, match="n_complex"):
        nuclear_grid(n_complex=n_complex)


def test_nuclear_grid_rejects_fractional_n_complex() -> None:
    with pytest.raises(GridError, match="n_complex"):
        nuclear_grid(n_complex=2.5)  # type: ignore[arg-type]


def test_nuclear_grid_accepts_one_complex_element() -> None:
    # The boundary of the n_complex constraint is admissible, not rejected.
    assert nuclear_grid(n_complex=1).R0 == 12.0


@pytest.mark.parametrize("r_max", [12.0, 10.0, 0.0, -5.0])
def test_nuclear_grid_rejects_r_max_inside_the_real_region(r_max: float) -> None:
    # 12.0 is the real-region endpoint itself: an empty tail, not a valid grid.
    with pytest.raises(GridError, match="r_max"):
        nuclear_grid(r_max=r_max)


@pytest.mark.parametrize("r_max", [float("nan"), float("inf"), float("-inf")])
def test_nuclear_grid_rejects_non_finite_r_max(r_max: float) -> None:
    with pytest.raises(GridError, match="r_max"):
        nuclear_grid(r_max=r_max)


def test_nuclear_grid_error_message_names_the_real_region_endpoint() -> None:
    with pytest.raises(GridError) as excinfo:
        nuclear_grid(r_max=12.0)
    assert "12.0" in str(excinfo.value)


def test_nuclear_grid_accepts_a_valid_custom_tail() -> None:
    # The small tail used by the two-angle stability families.
    g = nuclear_grid(angle_deg=18.0, r_max=14.0, n_complex=3, quadrature=8)
    assert g.R0 == 12.0
    # Nodes approach the outer edge without reaching it (it is the Dirichlet
    # endpoint), along the 18-degree ECS ray pivoted at R0.
    outer = g.points[-1] - 12.0
    assert 0.0 < abs(outer) < 2.0
    assert np.isclose(np.angle(outer), np.deg2rad(18.0))
