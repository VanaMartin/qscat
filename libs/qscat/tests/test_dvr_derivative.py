"""Tests for `qscat.dvr.derivative.dvr_first_derivative_at_node` -- the new
numerical primitive Task 3 (the `Flux` extractor) needs: the FEM-DVR
first-derivative operator's row at a node, checked against analytic
derivatives of known functions on single- and multi-element real grids.

Test-function choice matters here in a way it does not for `kinetic`: the
Dirichlet-drop basis (`FemDvrEcsGrid` docstring) structurally EXCLUDES the
two outermost grid points, so any coefficient vector built from sampling an
arbitrary smooth function implicitly assumes that function is exactly zero
at the domain edges (matching the actual propagated wavefunction, which
always is). A test function that is NOT (nearly) zero at the domain edges
sees genuine representation error there -- not a bug in this primitive, but
the SAME truncation the physical Hamiltonian eigenproblem relies on. So:
  - the `sin` case uses a single half-period `sin(pi*x/L)` (`L` = the grid's
    `R0`), which is EXACTLY zero at both domain edges by construction --
    resolved to ~1e-9 or better even at modest quadrature order.
  - the `gauss` case uses a bump placed and narrowed enough that its value
    at the domain edges is negligible (<1e-10) relative to its peak, at a
    quadrature order high enough (checked empirically, see below) for the
    GLL polynomial to resolve that width to ~1e-8 or better -- a single big
    element needs many more points to resolve a narrow bump than several
    small elements do (confirmed empirically while tuning these cases: a
    single element at length 6 needs nq~64-80 for a width-0.4 bump vs. nq~24-30
    for four 1-bohr elements resolving a width-0.15 bump).
"""

from __future__ import annotations

import numpy as np
import pytest
from qscat.dvr import ElementSpec, FemDvrEcsGrid, GridSpec
from qscat.dvr.derivative import dvr_first_derivative_at_node


def _bridge_node(grid: FemDvrEcsGrid, k: int) -> int:
    """The global index shared between element `k` and element `k+1`."""
    left = int(grid.element_maps[k][1][-1])
    right = int(grid.element_maps[k + 1][1][0])
    assert left == right
    return left


def _sin_case(grid: FemDvrEcsGrid) -> tuple[np.ndarray, np.ndarray]:
    """`f(x) = sin(pi x / L)`, `L = grid.R0` -- exactly zero at both domain edges."""
    k = np.pi / grid.R0
    r = grid.real_points
    f = np.sin(k * r)
    fp = k * np.cos(k * r)
    sqrt_w = np.sqrt(np.asarray(grid.weights, dtype=np.complex128))
    coeffs = (f * sqrt_w).astype(np.complex128)
    return coeffs, fp


def _gauss_case(grid: FemDvrEcsGrid, x0: float, s: float) -> tuple[np.ndarray, np.ndarray]:
    r = grid.real_points
    f = np.exp(-((r - x0) ** 2) / (2.0 * s**2))
    fp = -(r - x0) / s**2 * f
    sqrt_w = np.sqrt(np.asarray(grid.weights, dtype=np.complex128))
    coeffs = (f * sqrt_w).astype(np.complex128)
    return coeffs, fp


# --- single-element grid -----------------------------------------------------


def test_single_element_sin_interior_nodes() -> None:
    grid = FemDvrEcsGrid(GridSpec(quadrature=16, elements=[ElementSpec(4.0)], x_min=0.0))
    coeffs, fp = _sin_case(grid)
    for node in [0, grid.n // 2, grid.n - 1]:
        got = dvr_first_derivative_at_node(grid, node) @ coeffs
        np.testing.assert_allclose(got, fp[node], rtol=1e-8)


def test_single_element_gauss_interior_nodes() -> None:
    # nq=80 needed to resolve a width-0.4 bump across one length-6 element to
    # rtol=1e-8 (empirically tuned; see module docstring).
    grid = FemDvrEcsGrid(GridSpec(quadrature=80, elements=[ElementSpec(6.0)], x_min=0.0))
    coeffs, fp = _gauss_case(grid, x0=3.0, s=0.4)
    for node in [grid.n // 4, grid.n // 2, 3 * grid.n // 4]:
        got = dvr_first_derivative_at_node(grid, node) @ coeffs
        np.testing.assert_allclose(got, fp[node], rtol=1e-8)


# --- multi-element grid -------------------------------------------------------


def test_multi_element_sin_interior_and_border_nodes() -> None:
    grid = FemDvrEcsGrid(
        GridSpec(quadrature=10, elements=[ElementSpec(1.0) for _ in range(4)], x_min=0.0)
    )
    coeffs, fp = _sin_case(grid)
    interior_nodes = [int(grid.element_maps[1][1][3]), int(grid.element_maps[2][1][2])]
    # _bridge_node(grid, 1) sits at r=2.0 == the domain midpoint, exactly where
    # sin(pi*x/4)'s derivative vanishes (cos(pi/2)=0) -- excluded here since an
    # rtol check against a near-zero reference is a test-design artifact, not a
    # meaningful check; `test_row_is_zero_outside_the_chosen_element` still
    # exercises that bridge structurally.
    border_nodes = [_bridge_node(grid, 0), _bridge_node(grid, 2)]
    for node in interior_nodes + border_nodes:
        got = dvr_first_derivative_at_node(grid, node) @ coeffs
        np.testing.assert_allclose(got, fp[node], rtol=1e-8)


def test_multi_element_gauss_interior_and_border_nodes() -> None:
    # nq=30 needed to resolve a width-0.15 bump within a length-1 element to
    # rtol=1e-8 (empirically tuned; see module docstring).
    grid = FemDvrEcsGrid(
        GridSpec(quadrature=30, elements=[ElementSpec(1.0) for _ in range(4)], x_min=0.0)
    )
    coeffs, fp = _gauss_case(grid, x0=1.85, s=0.15)
    interior_nodes = [int(grid.element_maps[1][1][len(grid.element_maps[1][1]) // 2])]
    border_nodes = [_bridge_node(grid, 0), _bridge_node(grid, 1)]  # r=1.0, r=2.0 -- near x0=1.85
    for node in interior_nodes + border_nodes:
        got = dvr_first_derivative_at_node(grid, node) @ coeffs
        np.testing.assert_allclose(got, fp[node], rtol=1e-8)


# --- validation / edge cases --------------------------------------------------


def test_row_is_zero_outside_the_chosen_element() -> None:
    """The row only draws from ONE element's local nodes (module docstring)."""
    grid = FemDvrEcsGrid(
        GridSpec(quadrature=10, elements=[ElementSpec(1.0) for _ in range(4)], x_min=0.0)
    )
    node = _bridge_node(grid, 1)  # shared by elements 1 and 2; picks element 1
    d = dvr_first_derivative_at_node(grid, node)
    element_2_only_idx = grid.element_maps[2][1][1:]  # exclude the shared bridge index
    assert np.all(d[element_2_only_idx] == 0.0)


def test_out_of_range_raises() -> None:
    grid = FemDvrEcsGrid(GridSpec(quadrature=12, elements=[ElementSpec(4.0)], x_min=0.0))
    with pytest.raises(ValueError):
        dvr_first_derivative_at_node(grid, -1)
    with pytest.raises(ValueError):
        dvr_first_derivative_at_node(grid, grid.n)


def test_ecs_tail_node_rejected() -> None:
    grid = FemDvrEcsGrid(
        GridSpec(
            quadrature=8,
            elements=[ElementSpec(1.0), ElementSpec(1.0), ElementSpec(1.0, 35.0)],
            x_min=0.0,
        )
    )
    tail_node = int(grid.element_maps[-1][1][2])
    assert grid.real_points[tail_node] > grid.R0
    with pytest.raises(ValueError):
        dvr_first_derivative_at_node(grid, tail_node)
