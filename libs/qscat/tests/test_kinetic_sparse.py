"""Tests for `qscat.dvr.kinetic_sparse` (V2).

The existing DENSE `kinetic()` -- already validated in sub-project #1 against
analytic particle-in-a-box and bound-state theta-independence -- is retained
specifically as the differential oracle for this sparse implementation.
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp
from qscat.dvr import ElementSpec, FemDvrEcsGrid, GridSpec, kinetic, kinetic_sparse

CASES = {
    # name: (quadrature, elements)
    "all-real, uniform": (8, [ElementSpec(1.0) for _ in range(4)]),
    "all-real, graded": (6, [ElementSpec(0.5), ElementSpec(1.0), ElementSpec(2.0)]),
    "with ECS tail": (
        10,
        [ElementSpec(1.0), ElementSpec(1.0), ElementSpec(2.0, 35.0), ElementSpec(3.0, 35.0)],
    ),
    "single element": (7, [ElementSpec(1.5)]),
}


def _grid(name: str) -> FemDvrEcsGrid:
    q, els = CASES[name]
    return FemDvrEcsGrid(GridSpec(quadrature=q, elements=list(els), x_min=0.0))


@pytest.mark.parametrize("name", list(CASES))
@pytest.mark.parametrize("mass", [1.0, 12766.36])
def test_sparse_matches_dense_oracle(name: str, mass: float) -> None:
    grid = _grid(name)
    dense = kinetic(grid, mass)
    got = kinetic_sparse(grid, mass)
    assert got.shape == dense.shape
    scale = np.abs(dense).max()
    assert np.abs(got.toarray() - dense).max() <= 1e-12 * scale


# The formula below assumes the two dropped Dirichlet endpoints live in
# DIFFERENT elements, which requires tnel >= 2. For a single element both drops
# land in the same block and the count is instead (q-2)^2 -- i.e. the matrix is
# simply dense at n = q-2. The single-element grid is still covered by the
# differential test above, which is the check that actually matters.
@pytest.mark.parametrize("name", [n for n in CASES if n != "single element"])
def test_sparsity_matches_analytic_nnz_formula(name: str) -> None:
    """eMoScat KineticEnergy.cpp:95 -- nnz = q^2*tnel - 4q + 3 - tnel (tnel >= 2)."""
    grid = _grid(name)
    m = kinetic_sparse(grid, 1.0)
    m.eliminate_zeros()
    q = grid.nq
    tnel = len(grid.spec.elements)
    assert m.nnz == q**2 * tnel - 4 * q + 3 - tnel


def test_single_element_grid_is_dense() -> None:
    """tnel == 1: both Dirichlet drops fall in one block, giving a dense (q-2)^2."""
    grid = _grid("single element")
    m = kinetic_sparse(grid, 1.0)
    m.eliminate_zeros()
    assert grid.n == grid.nq - 2
    assert m.nnz == (grid.nq - 2) ** 2


def test_returns_csr_and_is_actually_sparse() -> None:
    grid = _grid("all-real, uniform")
    m = kinetic_sparse(grid, 1.0)
    assert isinstance(m, sp.csr_matrix)
    assert m.nnz < grid.n**2


def test_complex_symmetric_under_ecs() -> None:
    """ECS gives H = H^T but H != H^dagger; the kinetic term must already be so."""
    grid = _grid("with ECS tail")
    m = kinetic_sparse(grid, 1.0)
    assert abs(m - m.T).max() < 1e-12
    assert abs(m - m.conj().T).max() > 1e-6
