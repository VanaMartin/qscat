"""Independent finite-interval and asymptotic checks; no native CO inputs."""

import mpmath as mp
import numpy as np
import pytest
from scipy.integrate import solve_ivp
from scipy.special import spherical_jn, spherical_yn

from projects.ukrmol_co.coupled_tail import exact_tail, tail_coefficient, tail_waves
from projects.ukrmol_co.exterior_numerov import (
    match_ratio,
    numerov_k,
    tail_coefficients,
    tail_columns,
)


def test_free_centrifugal_tail_has_native_sign_and_flux():
    thresholds, angular, energy = np.array([0.0, 0.4, 2.0]), np.array([0, 3, 2]), 1.0
    c2 = np.diag(angular * (angular + 1))
    coefficients, s, spread = tail_coefficients(thresholds, energy, c2, np.zeros((3, 3)), 8)
    assert spread == 0
    r, h = 7.0, 1e-4
    f, g, closed = tail_columns(thresholds, angular, energy, coefficients, s, r, r)
    for i in range(2):
        k = np.sqrt(energy - thresholds[i])
        rho = k * r
        np.testing.assert_allclose(
            f[i, i], rho * spherical_jn(angular[i], rho) / np.sqrt(k), atol=1e-14
        )
        np.testing.assert_allclose(
            g[i, i], -rho * spherical_yn(angular[i], rho) / np.sqrt(k), atol=1e-14
        )
    before = tail_columns(thresholds, angular, energy, coefficients, s, r - h, r)
    after = tail_columns(thresholds, angular, energy, coefficients, s, r + h, r)
    fp, gp = (after[i] - before[i] for i in [0, 1])
    wronskian = (f.T @ gp - fp.T @ g) / (2 * h)
    np.testing.assert_allclose(wronskian, -np.eye(2), atol=2e-9)
    assert np.all(after[2][2] < closed[2]) and np.all(before[2][2] > closed[2])


def direct_match(q, thresholds, c2, c3, energy, radius, endpoint, order):
    """Independent unscaled fundamental columns on a short nonstiff interval."""
    count = len(thresholds)

    def rhs(r, state):
        values, derivative = state.reshape(2, count, count)
        v = np.diag(thresholds - energy) + c2 / r**2 + c3 / r**3
        return np.stack([derivative, v @ values]).ravel()

    result = solve_ivp(
        rhs,
        (radius, endpoint),
        np.stack([q, np.eye(count)]).ravel(),
        method="DOP853",
        rtol=1e-12,
        atol=1e-13,
    )
    assert result.success
    value, derivative = result.y[:, -1].reshape(2, count, count)
    logarithmic = np.linalg.solve(value.T, derivative.T).T
    coefficients, s, _ = tail_coefficients(thresholds, energy, c2, c3, order)
    # Derivative of the tail columns is evaluated independently by central difference.
    delta = 1e-4
    refs = tail_columns(thresholds, [0] * count, energy, coefficients, s, endpoint, endpoint)
    before = tail_columns(
        thresholds, [0] * count, energy, coefficients, s, endpoint - delta, endpoint
    )
    after = tail_columns(
        thresholds, [0] * count, energy, coefficients, s, endpoint + delta, endpoint
    )
    residuals = [
        (right - left) / (2 * delta) - logarithmic @ reference
        for left, right, reference in zip(before, after, refs, strict=True)
    ]
    return np.linalg.solve(np.column_stack(residuals[1:]), -residuals[0])[: refs[0].shape[1]]


@pytest.mark.parametrize("thresholds", [np.array([0.0, 0.3]), np.array([0.0, 1.8])])
def test_coupled_smooth_interval_has_fourth_order_convergence(thresholds):
    q = np.array([[0.3, 0.12], [0.12, -0.2]])
    c2 = np.array([[2.0, -0.4], [-0.4, 6.0]])
    c3 = np.array([[0.1, 0.3], [0.3, -0.2]])
    reference = direct_match(q, thresholds, c2, c3, 1.0, 2.0, 8.0, 6)
    errors = []
    for step in [0.2, 0.1, 0.05]:
        k, diagnostics = numerov_k(q, thresholds, [0, 0], c2, c3, 1.0, 2.0, 8.0, step, 6)
        errors.append(np.max(abs(k - reference)))
        assert diagnostics["matching_residual"] < 1e-12
    assert errors[1] < errors[0] / 10
    assert errors[2] < errors[1] / 10
    assert errors[2] < 1e-5


def test_zero_potential_has_exact_boundary_phase_and_closed_elimination():
    q = np.array([[0.4, 0.1], [0.1, 0.3]])
    thresholds, energy, radius = np.array([0.0, 2.0]), 0.7, 2.0
    k, kap = np.sqrt(energy), np.sqrt(2 - energy)
    effective_q = q[0, 0] - q[0, 1] ** 2 * kap / (1 + q[1, 1] * kap)
    phase = np.arctan(k * effective_q) - k * radius
    expected = np.tan(phase)
    zero = np.zeros((2, 2))
    measured, _ = numerov_k(q, thresholds, [0, 0], zero, zero, energy, radius, 10.0, 0.025, 4)
    np.testing.assert_allclose(measured, [[expected]], atol=2e-8)
    # Omitting a closed channel changes the effective boundary map.
    assert abs(expected - np.tan(np.arctan(k * q[0, 0]) - k * radius)) > 0.01


def test_unequal_open_momenta_match_analytic_flux_normalized_matrix():
    thresholds, energy, radius = np.array([0.0, 0.4]), 1.0, 2.0
    momenta = np.sqrt(energy - thresholds)
    q = np.array([[0.3, 0.2], [0.2, -0.1]])
    f = np.diag(np.sin(momenta * radius) / np.sqrt(momenta))
    g = np.diag(np.cos(momenta * radius) / np.sqrt(momenta))
    fp = np.diag(np.cos(momenta * radius) * np.sqrt(momenta))
    gp = -np.diag(np.sin(momenta * radius) * np.sqrt(momenta))
    expected = np.linalg.solve(g - q @ gp, -(f - q @ fp))
    zero = np.zeros((2, 2))
    measured, _ = numerov_k(q, thresholds, [0, 0], zero, zero, energy, radius, 8.0, 0.025, 4)
    np.testing.assert_allclose(measured, expected, atol=2e-8)
    np.testing.assert_allclose(measured, measured.T, atol=2e-8)
    wrong = (np.sqrt(momenta)[:, None] * expected) / np.sqrt(momenta)
    assert np.max(abs(wrong - expected)) > 0.01


def test_closed_reference_rescaling_does_not_change_open_k():
    thresholds, energy = np.array([0.0, 2.0]), 0.7
    zero = np.zeros((2, 2))
    coefficients, s, _ = tail_coefficients(thresholds, energy, zero, zero, 4)
    values = tail_columns(thresholds, [0, 0], energy, coefficients, s, 2.0, 2.0)
    derivatives = tail_columns(
        thresholds, [0, 0], energy, coefficients, s, 2.0, 2.0, derivative=True
    )
    q = np.array([[0.4, 0.1], [0.1, 0.3]])
    expected, _ = match_ratio(np.eye(2), -q, values, derivatives)
    for scale in [1e-40, -3.0, 1e40]:
        scaled_values = (*values[:2], values[2] * scale)
        scaled_derivatives = (*derivatives[:2], derivatives[2] * scale)
        measured, _ = match_ratio(np.eye(2), -q, scaled_values, scaled_derivatives)
        np.testing.assert_allclose(measured, expected, atol=1e-13)


def test_coupled_degenerate_tail_against_noninteger_hankel_waves():
    q, dipole, radius = 1.0, 1.2, 20.0
    c2 = tail_coefficient(dipole)
    (hp, dp), _ = tail_waves(q, radius, dipole)
    errors = []
    for order in [4, 8, 12]:
        coefficients, s, _ = tail_coefficients(np.zeros(2), q**2, c2, np.zeros((2, 2)), order)
        f, g, _ = tail_columns(np.zeros(2), [1, 2], q**2, coefficients, s, radius, radius)
        fp, gp, _ = tail_columns(
            np.zeros(2), [1, 2], q**2, coefficients, s, radius, radius, derivative=True
        )
        errors.append(max(np.max(abs(g + 1j * f - hp)), np.max(abs(gp + 1j * fp - dp))))
    assert errors[1] < errors[0] / 100
    assert errors[2] < errors[1] / 10
    assert errors[2] < 1e-11


def test_hard_wall_polar_scattering_matches_independent_exact_oracle():
    zero = np.zeros((2, 2))
    k, _ = numerov_k(zero, [0, 0], [1, 2], tail_coefficient(1.2), zero, 1.0, 1.0, 20.0, 0.025, 12)
    identity = np.eye(2)
    s = np.linalg.solve(identity - 1j * k, identity + 1j * k)
    with mp.workdps(60):
        expected = np.array(exact_tail(1.0, 1.2).tolist(), complex)
    np.testing.assert_allclose(s, expected, atol=1e-7, rtol=0)
