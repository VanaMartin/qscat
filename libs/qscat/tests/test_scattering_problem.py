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
    from qscat.core.dissociation import DrResult

    prob, _, _, tg = _problem_and_basis()
    sentinel = np.zeros((2, 2))
    seen: dict[str, object] = {}

    def fake_dr_solve(tgrid, model, eps, chi, v_init, E, **kw):
        seen.update(tgrid=tgrid, model=model, eps=eps, chi=chi, v_init=v_init, E=E, **kw)
        return DrResult(sigma=sentinel, psi=None, amplitude=None)

    monkeypatch.setattr("qscat.core.problem.dr_solve", fake_dr_solve)
    got = prob.dr_cross_section([0.01, 0.03], n_channels=2)
    assert got is sentinel
    assert seen["tgrid"] is tg
    assert seen["model"] is N2
    assert seen["v_init"] == 0
    assert np.array_equal(seen["eps"], prob.eps)
    assert np.array_equal(seen["chi"], prob.chi)
    assert seen["n_channels"] == 2
    assert seen["ordering"] == "COLAMD"
    assert seen["store_wavefunction"] is False
    assert seen["store_amplitude"] is False


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


def test_problem_lcp_da_matches_functional_api() -> None:
    """Synthetic curve arrays: delegation equality needs a well-posed solve,
    not converged physics (docs/adr/0005 point 7)."""
    from qscat.core import lcp_da_cross_section

    prob, eps, chi, tg = _problem_and_basis()
    g_R = tg.grids[1]
    Vd = (0.2 * (1.0 - np.exp(-(g_R.points - 2.0))) ** 2 + 0.05).astype(np.complex128)
    Gamma = 0.01 * np.exp(-(np.abs(g_R.points - 2.4) ** 2)).astype(np.float64)
    E = np.array([0.02, 0.05])
    expected = lcp_da_cross_section(g_R, N2.mu, Vd, Gamma, eps, chi, 0, E)
    got = prob.lcp_da_cross_section(E, Vd=Vd, Gamma=Gamma)
    assert np.array_equal(got, expected)


def test_problem_resonance_levels_delegates_exact_arguments(monkeypatch) -> None:
    """The electronic pole walk is minutes-scale; delegation is checked by
    argument capture (the walk's own gates live in test_lcp_resonance_levels
    and the validation harness)."""
    prob, _, _, tg = _problem_and_basis()
    nuc_b = nuclear_grid(r_max=14.0, quadrature=6, n_complex=3, angle_deg=25.0)
    elec_b = electronic_grid(r_max=12.0, order=5, n_complex=3, angle_deg=25.0)
    sentinel = object()
    seen: dict[str, object] = {}

    def fake_levels(model, nuclear_grid_a, nuclear_grid_b, elec_grid_a, elec_grid_b, **kw):
        seen.update(
            model=model,
            nuc_a=nuclear_grid_a,
            nuc_b=nuclear_grid_b,
            elec_a=elec_grid_a,
            elec_b=elec_grid_b,
            **kw,
        )
        return sentinel

    monkeypatch.setattr("qscat.core.problem.resonance_levels", fake_levels)
    got = prob.resonance_levels(nuc_b, elec_b, n_levels=3)
    assert got is sentinel
    assert seen["model"] is N2
    assert seen["nuc_a"] is tg.grids[1] and seen["elec_a"] is tg.grids[0]
    assert seen["nuc_b"] is nuc_b and seen["elec_b"] is elec_b
    assert seen["n_levels"] == 3 and seen["return_curve"] is False


def test_problem_exact_resonance_states_delegates_exact_arguments(monkeypatch) -> None:
    """2-D pole searches are minutes of sparse factorizations; argument
    capture checks the wiring (the solver's own gates live in
    test_exact_resonance_states.py)."""
    prob, _, _, tg = _problem_and_basis()
    g_elec = TensorGrid(
        [electronic_grid(r_max=12.0, order=5, n_complex=3, angle_deg=40.0), tg.grids[1]]
    )
    g_nuc = TensorGrid(
        [tg.grids[0], nuclear_grid(r_max=14.0, quadrature=6, n_complex=3, angle_deg=25.0)]
    )
    sentinel = object()
    seen: dict[str, object] = {}

    def fake_exact(model, grid_base, grid_electronic, grid_nuclear, **kw):
        seen.update(model=model, base=grid_base, ge=grid_electronic, gn=grid_nuclear, **kw)
        return sentinel

    monkeypatch.setattr("qscat.core.problem.exact_resonance_states", fake_exact)
    got = prob.exact_resonance_states(
        g_elec, g_nuc, shifts=[-0.66 - 0.004j], window=(-0.75, -0.55, -0.05, 0.0)
    )
    assert got is sentinel
    assert seen["model"] is N2 and seen["base"] is tg
    assert seen["ge"] is g_elec and seen["gn"] is g_nuc
    assert seen["shifts"] == [-0.66 - 0.004j]
    assert seen["k"] == 8


def test_problem_nrm_methods_delegate_exact_arguments(monkeypatch) -> None:
    """The NRM ingredient build (fixed-R electronic eigenbases) is the
    expensive part and its physics gates are validation/diatomic's; argument
    capture checks the facade wiring, including the NUCLEAR-grid-first
    argument order nrm uses."""
    import qscat.core.nrm as nrm_pkg

    prob, _, _, tg = _problem_and_basis()
    phi_d = object()  # any DiscreteState; never touched by the fake
    sentinel = object()
    seen: dict[str, object] = {}

    def fake_ve(nuclear_grid, elec_grid, model, phi_d_got, eps, chi, v_init, vprimes, E, **kw):
        seen.update(
            nuc=nuclear_grid,
            elec=elec_grid,
            model=model,
            phi_d=phi_d_got,
            v_init=v_init,
            vprimes=vprimes,
            E=E,
            **kw,
        )
        return sentinel

    monkeypatch.setattr(nrm_pkg, "nrm_ve_cross_section", fake_ve)
    got = prob.nrm_ve_cross_section(phi_d, [0, 1], 0.05, n_states=20)
    assert got is sentinel
    assert seen["nuc"] is tg.grids[1] and seen["elec"] is tg.grids[0]
    assert seen["model"] is N2 and seen["phi_d"] is phi_d
    assert seen["vprimes"] == [0, 1] and seen["n_states"] == 20
    assert seen["include_background"] is True

    def fake_da(nuclear_grid, elec_grid, model, phi_d_got, eps, chi, v_init, E, **kw):
        seen.clear()
        seen.update(
            nuc=nuclear_grid,
            elec=elec_grid,
            model=model,
            phi_d=phi_d_got,
            eps=eps,
            chi=chi,
            v_init=v_init,
            E=E,
            **kw,
        )
        return sentinel

    monkeypatch.setattr(nrm_pkg, "nrm_da_cross_section", fake_da)
    got = prob.nrm_da_cross_section(phi_d, 0.05)
    assert got is sentinel
    assert seen["nuc"] is tg.grids[1] and seen["elec"] is tg.grids[0]
    assert seen["phi_d"] is phi_d and seen["ingredients"] is None
    assert np.array_equal(seen["eps"], prob.eps)
    assert np.array_equal(seen["chi"], prob.chi)
    assert seen["v_init"] == 0
