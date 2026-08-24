"""`ScatteringProblem` facade: delegates to the functional core, bit-for-bit."""

from __future__ import annotations

import numpy as np
import pytest
from qscat.core import ScatteringProblem, ve_cross_section, vibrational_states
from qscat.core.grids import electronic_grid, nuclear_grid
from qscat.dvr import TensorGrid
from qscat.model import N2


def _grid() -> TensorGrid:
    return TensorGrid(
        [
            electronic_grid(r_max=16.0, order=7, n_complex=5),
            nuclear_grid(r_max=22.0, quadrature=10, n_complex=5),
        ]
    )


@pytest.mark.slow
def test_problem_ve_matches_functional_api() -> None:
    tg = _grid()
    E = np.array([0.10, 0.15, 0.20])
    prob = ScatteringProblem(grid=tg, model=N2, n_vib=4, v_init=0)

    # The facade solves the same basis and calls the same solver.
    eps, chi = vibrational_states(tg.grids[1], N2.mu, 4, N2.v0)
    expected = ve_cross_section(tg, N2, eps, chi, 0, [0, 1, 2], E)
    got = prob.ve_cross_section(vprimes=[0, 1, 2], E=E)

    assert np.array_equal(got, expected)


def test_problem_exposes_basis_and_is_frozen() -> None:
    prob = ScatteringProblem(grid=_grid(), model=N2, n_vib=4)
    assert prob.eps.shape == (4,)
    assert prob.chi.shape[0] == 4
    # Basis round-trips through the NamedTuple accessors.
    assert np.array_equal(prob.basis.eps, prob.eps)
    # Frozen: cannot rebind fields.
    import dataclasses

    with pytest.raises(dataclasses.FrozenInstanceError):
        prob.n_vib = 5  # type: ignore[misc]


def test_vibrational_states_returns_named_tuple_backcompat() -> None:
    tg = _grid()
    basis = vibrational_states(tg.grids[1], N2.mu, 3, N2.v0)
    # Named access...
    assert basis.eps.shape == (3,)
    # ...and legacy tuple unpacking both work.
    eps, chi = basis
    assert np.array_equal(eps, basis.eps)
    assert chi.shape[0] == 3
