"""`ScatteringProblem` facade: delegates to the functional core, bit-for-bit."""

from __future__ import annotations

import inspect

import numpy as np
import pytest
from qscat.core import (
    ScatteringProblem,
    da_cross_section,
    td_da_cross_section,
    td_da_cross_sections_all,
    td_ve_cross_section,
    td_ve_cross_sections_all,
    ve_cross_section,
    vibrational_states,
)
from qscat.core.grids import electronic_grid, nuclear_grid
from qscat.dvr import TensorGrid
from qscat.model import N2


def _grid() -> TensorGrid:
    """A toy deck (59 x 129 = 7611), not the production working grid.

    Nothing in this module is a physics assertion: the facade check below is
    `array_equal` between two routes through the SAME solver, which holds on
    any grid the solver accepts, and the other two tests only read shapes.
    So per docs/adr/0005 point 7 the deck shrinks and the tests stay in the
    fast gate -- on the working grid (107 x 251) the facade check measured
    24.4 s, here 0.9 s.
    """
    return TensorGrid(
        [
            electronic_grid(r_max=12.0, order=5, n_complex=3),
            nuclear_grid(r_max=14.0, quadrature=6, n_complex=3),
        ]
    )


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


def test_facade_methods_have_real_signatures() -> None:
    """lib-C2: every public facade method exposes the functional solver's
    real parameters -- no `**kwargs`, no `*args`, and a real return type."""
    for name, fn in inspect.getmembers(ScatteringProblem, predicate=inspect.isfunction):
        if name.startswith("_"):
            continue
        sig = inspect.signature(fn)
        kinds = {p.kind for p in sig.parameters.values()}
        assert inspect.Parameter.VAR_KEYWORD not in kinds, f"{name} takes **kwargs"
        assert inspect.Parameter.VAR_POSITIONAL not in kinds, f"{name} takes *args"
        assert sig.return_annotation is not inspect.Signature.empty, f"{name} has no return type"


WP_IN = {"r0": 4.0, "p0": -0.5, "sigma": 1.2}
WP_OUT = {"r0_out": 6.0, "p0_out": 0.5, "sigma_out": 1.0}
NUCLEAR_WP_OUT = {"r0_out": 7.0, "p0_out": 5.0, "sigma_out": 1.0}
POSITION = 37
NUCLEAR_SURFACE = 90
DT = 0.2
N_STEPS = 3


def _problem_and_basis() -> tuple[ScatteringProblem, np.ndarray, np.ndarray, TensorGrid]:
    tg = _grid()
    prob = ScatteringProblem(grid=tg, model=N2, n_vib=4, v_init=0)
    eps, chi = vibrational_states(tg.grids[1], N2.mu, 4, N2.v0)
    return prob, eps, chi, tg


def test_problem_da_matches_functional_api() -> None:
    prob, eps, chi, tg = _problem_and_basis()
    E = np.array([0.10, 0.60])  # one closed, one open DA channel on N2
    expected = da_cross_section(tg, N2, eps, chi, 0, E, n_channels=1)
    assert np.array_equal(prob.da_cross_section(E, n_channels=1), expected)


def test_problem_dr_delegates_exact_arguments(monkeypatch) -> None:
    """dr's real solve is slow-tier (mpmath Coulomb channels), so delegation
    is checked by argument capture instead of a re-solve."""
    prob, _, _, tg = _problem_and_basis()
    sentinel = np.zeros((2, 2))
    seen: dict[str, object] = {}

    def fake_dr(tgrid, model, eps, chi, v_init, E, **kw):
        seen.update(tgrid=tgrid, model=model, eps=eps, chi=chi, v_init=v_init, E=E, **kw)
        return sentinel

    monkeypatch.setattr("qscat.core.problem.dr_cross_section", fake_dr)
    got = prob.dr_cross_section([0.01, 0.03], n_channels=2)
    assert got is sentinel
    assert seen["tgrid"] is tg
    assert seen["model"] is N2
    assert seen["v_init"] == 0
    assert np.array_equal(seen["eps"], prob.eps)
    assert np.array_equal(seen["chi"], prob.chi)
    assert seen["n_channels"] == 2
    assert seen["ordering"] == "COLAMD"


def test_problem_td_ve_matches_functional_api() -> None:
    prob, eps, chi, tg = _problem_and_basis()
    E = [0.10, 0.15]
    expected = td_ve_cross_section(
        tg, N2, eps, chi, 0, [0, 1], E, dt=DT, n_steps=N_STEPS, wp_in=WP_IN, wp_out=WP_OUT
    )
    got = prob.td_ve_cross_section([0, 1], E, dt=DT, n_steps=N_STEPS, wp_in=WP_IN, wp_out=WP_OUT)
    assert np.array_equal(got, expected)


def test_problem_td_ve_all_matches_functional_api() -> None:
    prob, eps, chi, tg = _problem_and_basis()
    expected = td_ve_cross_sections_all(
        tg,
        N2,
        eps,
        chi,
        0,
        [0, 1],
        0.10,
        dt=DT,
        n_steps=N_STEPS,
        wp_in=WP_IN,
        wp_out=WP_OUT,
        position=POSITION,
        surface=POSITION,
    )
    got = prob.td_ve_cross_sections_all(
        [0, 1],
        0.10,
        dt=DT,
        n_steps=N_STEPS,
        wp_in=WP_IN,
        wp_out=WP_OUT,
        position=POSITION,
        surface=POSITION,
    )
    assert set(got) == {"tw", "delta", "flow"}
    for key in expected:
        assert np.array_equal(got[key], expected[key])


def test_problem_td_da_matches_functional_api() -> None:
    prob, eps, chi, tg = _problem_and_basis()
    expected = td_da_cross_section(
        tg,
        N2,
        eps,
        chi,
        0,
        0.60,
        dt=DT,
        n_steps=N_STEPS,
        wp_in=WP_IN,
        method="flow",
        surface=NUCLEAR_SURFACE,
        n_channels=1,
    )
    got = prob.td_da_cross_section(
        0.60,
        dt=DT,
        n_steps=N_STEPS,
        wp_in=WP_IN,
        method="flow",
        surface=NUCLEAR_SURFACE,
        n_channels=1,
    )
    assert np.array_equal(got, expected)


def test_problem_td_da_all_matches_functional_api() -> None:
    prob, eps, chi, tg = _problem_and_basis()
    kw = dict(
        dt=DT,
        n_steps=N_STEPS,
        wp_in=WP_IN,
        surface=NUCLEAR_SURFACE,
        position=NUCLEAR_SURFACE,
        wp_out=NUCLEAR_WP_OUT,
        n_channels=1,
    )
    expected = td_da_cross_sections_all(tg, N2, eps, chi, 0, 0.60, **kw)
    got = prob.td_da_cross_sections_all(0.60, **kw)
    assert set(got) == {"flow", "delta", "tw"}
    for key in expected:
        assert np.array_equal(got[key], expected[key])
