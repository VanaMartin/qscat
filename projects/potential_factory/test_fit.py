# projects/potential_factory/test_fit.py
from __future__ import annotations

import re
from dataclasses import replace

import numpy as np
import pytest
from qscat.model import N2

from projects.potential_factory.ansatz import SmoothR, from_diatomic, with_params
from projects.potential_factory.extract import extract_target
from projects.potential_factory.fit import (
    _check_omega_e,
    _joint_polish,
    _smooth_grad,
    _smooth_reparam,
    _threshold_exponent_mismatch,
    fit_coupling,
    fit_neutral,
    fit_resonance,
    model_gamma_tilde,
)
from projects.potential_factory.report import Tolerances, ecs_bounded
from projects.potential_factory.target import CouplingTarget, Curve, NeutralTarget
from projects.potential_factory.tracker import ElectronicPair


def test_fit_neutral_recovers_morse_constants_from_the_curve():
    R = np.linspace(1.4, 5.0, 60)
    target = NeutralTarget(Curve.from_table(R, N2.v0(R).real), {}, (1.5, 4.8))
    seed = with_params(from_diatomic(N2), {"D_e": 0.6, "R_e": 2.2, "beta0": 1.0})
    model, res = fit_neutral(target, seed, tol=Tolerances())
    assert res.status == "met"
    assert (
        abs(model.D_e - N2.D0) < 1e-8
        and abs(model.R_e - N2.R0) < 1e-8
        and abs(model.betas[0] - N2.alpha0) < 1e-8
    )


def test_fit_neutral_from_constants_uses_the_morse_relation():
    omega_e = N2.alpha0 * np.sqrt(2.0 * N2.D0 / N2.mu)
    target = NeutralTarget(None, {"R_e": N2.R0, "D_e": N2.D0, "omega_e": omega_e}, (1.5, 4.8))
    model, res = fit_neutral(target, from_diatomic(N2), tol=Tolerances())
    assert res.status == "met" and abs(model.betas[0] - N2.alpha0) < 1e-12
    # The ladder check compares the fitted curve's ANHARMONIC 0->1 spacing
    # against the Morse `omega_e - 2 omega_e x_e` these same constants imply,
    # not against the bare harmonic `omega_e` -- that offset is 0.83% on N2,
    # 1.96% on NO and 3.33% on F2, so comparing to `omega_e` itself would fail
    # the 1% default tolerance for an exactly reproduced curve.
    rel = float(re.search(r"omega_e_rel=([\d.eE+-]+)", res.detail).group(1))
    assert rel < 1e-10, res.detail


def test_fit_neutral_checks_the_vibrational_ladder_against_the_target_curve():
    """T0 reports the `v=0 -> v=1` spacing error, and gates "met" on it."""
    R = np.linspace(1.4, 5.0, 60)
    target = NeutralTarget(Curve.from_table(R, N2.v0(R).real), {}, (1.5, 4.8))
    _, res = fit_neutral(target, from_diatomic(N2), tol=Tolerances())
    assert res.status == "met"
    rel = float(re.search(r"omega_e_rel=([\d.eE+-]+)", res.detail).group(1))
    # The floor here is the 60-point table's own cubic-spline interpolation
    # error, not the fit: the fitted curve IS N2.v0 to round-off, while the
    # reference ladder is solved on the spline through the table.
    assert rel < 1e-5, res.detail

    # A table too narrow to confine v=0,1 cannot define a spacing at all; the
    # tier must SAY it did not check rather than quoting a meaningless number.
    Rn = np.linspace(2.0, 2.2, 20)
    narrow = NeutralTarget(Curve.from_table(Rn, N2.v0(Rn).real), {}, (2.0, 2.2))
    _, res_n = fit_neutral(narrow, from_diatomic(N2), tol=Tolerances())
    assert "omega_e: not checked" in res_n.detail and "does not confine" in res_n.detail


def test_omega_e_check_rejects_a_ladder_that_is_wrong():
    """The check has to be able to FAIL. A 20% stiffer `beta0` leaves `D_e`
    and `R_e` untouched -- the constants branch reports rms 0 either way -- so
    the ladder is the only thing that sees it."""
    omega_e = N2.alpha0 * np.sqrt(2.0 * N2.D0 / N2.mu)
    target = NeutralTarget(None, {"R_e": N2.R0, "D_e": N2.D0, "omega_e": omega_e}, (1.5, 4.8))
    stiff = with_params(from_diatomic(N2), {"beta0": N2.alpha0 * 1.2})
    ok, detail = _check_omega_e(target, stiff, Tolerances())
    assert not ok, detail
    rel = float(re.search(r"omega_e_rel=([\d.eE+-]+)", detail).group(1))
    assert rel > 0.05, detail


@pytest.mark.slow
def test_fit_resonance_round_trips_n2_curves():
    pair = ElectronicPair()
    R_desc = np.linspace(3.0, 1.6, 12)
    target = extract_target(N2, pair=pair, R_desc=R_desc, n_eps=3)
    seed = with_params(from_diatomic(N2), {"lam.f_inf": 5.5, "alpha.f_inf": 0.5})
    model, res = fit_resonance(target.resonance, seed, pair=pair, n_nodes=12, tol=Tolerances())
    assert res.status == "met", res.detail
    np.testing.assert_allclose(model.lam_R(R_desc).real, N2.lam(R_desc).real, rtol=2e-3)
    np.testing.assert_allclose(model.alpha_R(R_desc).real, N2.alpha_c, rtol=2e-3)


@pytest.mark.slow
def test_fit_coupling_round_trip_keeps_the_shell_negligible():
    pair = ElectronicPair()
    R_desc = np.linspace(3.0, 1.6, 8)
    target = extract_target(N2, pair=pair, R_desc=R_desc, n_eps=6)
    model, res = fit_coupling(
        target.coupling, from_diatomic(N2), pair=pair, n_eps=4, n_R=4, tol=Tolerances()
    )
    assert res.status == "met", res.detail
    assert res.rms < 0.02
    assert abs(model.shell_R(2.0).real) < 1e-3
    assert "shell not needed" in res.detail


@pytest.mark.slow
def test_fit_coupling_moves_the_shell_for_a_steeper_falloff():
    pair = ElectronicPair()
    R_desc = np.linspace(3.0, 1.6, 8)
    base = extract_target(N2, pair=pair, R_desc=R_desc, n_eps=6).coupling
    eps = np.geomspace(*base.eps_window, 6)
    R_asc = R_desc[::-1]
    g = model_gamma_tilde(from_diatomic(N2), pair, eps, R_asc) * np.exp(-3.0 * eps)[:, None]
    steeper = CouplingTarget.from_table(eps, R_asc, g, alpha=2.5)
    model, res = fit_coupling(
        steeper, from_diatomic(N2), pair=pair, n_eps=4, n_R=4, tol=Tolerances()
    )
    assert model.shell is not None and abs(model.shell_R(2.0).real) > 1e-3
    assert res.status == "met", res.detail
    seed_rms = float(re.search(r"seed_rms=([\d.eE+-]+)", res.detail).group(1))
    assert res.rms < seed_rms


def test_fit_coupling_reports_a_threshold_exponent_out_of_the_ansatz():
    """`ell` is fixed, so `Gamma~ ~ eps^(l+1/2)` is a property of the ansatz,
    not something a shell can fit. A target carrying a different near-threshold
    power -- a POLAR molecule, `alpha = sqrt(d + 1/4)` -- must be reported as
    out of scope instead of silently fitted with the wrong threshold law. No
    solve is run, so this test is fast."""
    model = from_diatomic(N2)
    polar = CouplingTarget.from_alt_houfek(
        a0=1.0, a1=0.0, a2=0.0, b0=1.0, b1=0.0, alpha=1.2, R_range=(1.6, 3.0)
    )
    out, res = fit_coupling(polar, model, pair=None, tol=Tolerances())
    assert res.status == "not met"
    assert np.isnan(res.rms) and np.isnan(res.max)
    assert "alpha=1.2" in res.detail and f"l+1/2={model.ell + 0.5:g}" in res.detail
    assert "out of scope" in res.detail
    assert out is model  # nothing fitted, nothing installed

    # The matching exponent is not flagged -- the guard is specific, not a
    # blanket refusal of `from_alt_houfek` targets.
    matching = CouplingTarget.from_alt_houfek(
        a0=1.0, a1=0.0, a2=0.0, b0=1.0, b1=0.0, alpha=model.ell + 0.5, R_range=(1.6, 3.0)
    )
    assert _threshold_exponent_mismatch(matching, model) is None
    assert _threshold_exponent_mismatch(polar, model) is not None


def test_smooth_grad_matches_finite_differences():
    """`_smooth_grad` is the analytic derivative `_joint_polish` hands
    `least_squares`; differential-test it against `_smooth_reparam`'s own
    `decode` in both branches (sigmoid, and the degenerate `f_0 == 0` one)."""
    for s in (
        SmoothR(18.849, -8.8746, 3.213, 1.832, coeffs=(0.3, -0.15), R_e=2.69, p=3),
        SmoothR(3.0, 0.0, 1.0, 0.0, coeffs=(), R_e=2.69, p=3),
    ):
        x0, _, _, decode = _smooth_reparam(s, 2.0, 3.6)
        for R in (2.0, 2.5333, 3.6):
            analytic = _smooth_grad(decode(x0), R)
            fd = np.zeros_like(x0)
            for i in range(x0.size):
                h = 1e-6 * max(1.0, abs(x0[i]))
                xp, xm = x0.copy(), x0.copy()
                xp[i] += h
                xm[i] -= h
                fd[i] = (float(decode(xp)(R).real) - float(decode(xm)(R).real)) / (2.0 * h)
            assert analytic.shape == x0.shape
            np.testing.assert_allclose(analytic, fd, rtol=1e-5, atol=1e-9)


def test_joint_polish_never_scores_a_model_with_no_poles_as_perfect():
    """A trial where no node yields a gated pole must be EXPENSIVE, not free.

    Scoring a missing node 0 made "find nothing anywhere" a global minimum of
    the sum of squares, and the optimizer took it -- one step to a well of
    negative depth and negative width, reported as `polish_rms=0.00e+00` with
    every node skipped. Here every node misses by construction (`alpha < 0`),
    so the polish must report that it resolved nothing (no rms to quote) and
    must not have moved off the seed on the strength of it.
    """
    pair = ElectronicPair()
    model = from_diatomic(N2)
    model = replace(model, alpha=replace(model.alpha, f_inf=-1.0))
    R_nodes = np.linspace(2.4, 2.0, 3)
    mask = np.zeros(R_nodes.size, dtype=bool)
    fitted, detail = _joint_polish(
        model,
        pair,
        R_nodes,
        lambda R: complex(-0.1, 0.0),
        np.zeros(R_nodes.size),
        mask,
        [None] * R_nodes.size,
        None,
    )
    assert "polish skipped 3/3 nodes" in detail, detail
    assert "polish_rms=nan" in detail, detail
    assert fitted.alpha.f_inf == model.alpha.f_inf
    # Not exact equality: `lam` round-trips through `log|f_0|` -> `exp` in
    # `_smooth_reparam`'s coordinates, which costs an ulp even for x == x0.
    assert fitted.lam.f_inf == model.lam.f_inf
    assert abs(fitted.lam.f_0 / model.lam.f_0 - 1.0) < 1e-12


def test_ecs_bounded_separates_analytic_bounds_from_probed_angles():
    """The 45/90 caps are properties of the ansatz, not of any grid; the
    probed angles are what was actually evaluated. Also covers the
    derive-`nuclear_deg`-from-`R_tail` path, which `fit()` does not take
    (it passes the angle it built the tail from)."""
    pair = ElectronicPair()
    model = from_diatomic(N2)
    R_tail = 12.0 + np.linspace(0.1, 6.0, 8) * np.exp(1j * np.deg2rad(35.0))
    out = ecs_bounded(model, pair, R_tail, nuclear_deg=35.0)
    assert out["electronic_max_deg"] == 45.0 and out["nuclear_max_deg"] == 90.0
    # The electronic cap is a property of the Gaussian well alone -- it does
    # not move when the grid it is probed on does.
    assert out["probed_electronic_deg"] == 35.0
    other = ecs_bounded(model, ElectronicPair(angles=(20.0, 30.0)), R_tail, nuclear_deg=35.0)
    assert other["electronic_max_deg"] == 45.0
    assert other["probed_electronic_deg"] == 20.0
    # Derived and passed nuclear angles agree to round-off.
    derived = ecs_bounded(model, pair, R_tail)
    assert abs(derived["probed_nuclear_deg"] - 35.0) < 1e-9
    assert derived["tail_growth"] == out["tail_growth"]
    with pytest.raises(ValueError, match="fewer than two points"):
        ecs_bounded(model, pair, R_tail[:1])
