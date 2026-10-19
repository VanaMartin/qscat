"""Independent analytic/coupled checks of cancellation-reduced Numerov algebra."""

import mpmath as mp
import numpy as np
import pytest

from projects.ukrmol_co.numerov_increment import increment_numerov_k
from projects.ukrmol_co.numerov_increment_study import exact_discrete_scalar, exact_free_k
from projects.ukrmol_co.test_exterior_numerov import direct_match


def test_defect_step_preserves_noncommuting_matrix_order():
    w = np.array([[0.3, 0.4], [0.4, -0.7]])
    ratio = np.array([[0.9, 0.02], [-0.03, 1.2]])
    identity, h = np.eye(2), 0.1
    b = identity - h**2 * w / 12
    original = np.linalg.inv(12 * np.linalg.inv(b) - 10 * identity - ratio)
    x = identity - ratio + h**2 * np.linalg.solve(b, w)
    defect = np.linalg.solve(identity + x, x)
    np.testing.assert_allclose(identity - defect, original, atol=2e-15, rtol=0)


@pytest.mark.parametrize("thresholds", [np.array([0.0, 0.3]), np.array([0.0, 1.8])])
def test_noncommuting_coupled_problem_has_fourth_order_convergence(thresholds):
    q = np.array([[0.3, 0.12], [0.12, -0.2]])
    c2 = np.array([[2.0, -0.4], [-0.4, 6.0]])
    c3 = np.array([[0.1, 0.3], [0.3, -0.2]])
    reference = direct_match(q, thresholds, c2, c3, 1.0, 2.0, 8.0, 6)
    errors = []
    for step in [0.2, 0.1, 0.05]:
        measured, _ = increment_numerov_k(q, thresholds, [0, 0], c2, c3, 1, 2, 8, step, 6)
        errors.append(np.max(abs(measured - reference)))
    assert errors[1] < errors[0] / 8
    assert errors[2] < errors[1] / 8
    assert errors[2] < 1e-5


def test_singular_boundary_matches_exact_free_hard_wall():
    zero = np.zeros((2, 2))
    measured, _ = increment_numerov_k(zero, [0, 2], [0, 0], zero, zero, 0.7, 2, 8, 0.0125, 4)
    np.testing.assert_allclose(measured, [[-np.tan(2 * np.sqrt(0.7))]], atol=2e-8, rtol=0)


def test_continuous_oracle_has_known_boundary_phase_and_closed_elimination():
    q = np.array([[0.4, 0.1], [0.1, 0.3]])
    exact = exact_free_k(q, np.array([0, 2]), 0.7, 2, 80)
    kap, k = np.sqrt(1.3), np.sqrt(0.7)
    effective_q = q[0, 0] - q[0, 1] ** 2 * kap / (1 + q[1, 1] * kap)
    assert float(exact[0, 0]) == pytest.approx(
        np.tan(np.arctan(k * effective_q) - 2 * k), abs=1e-13
    )


def test_closed_form_discrete_oracle_matches_direct_high_precision_recurrence():
    with mp.workdps(80):
        q, energy, radius, endpoint, step = map(mp.mpf, [0.4, 0.7, 2, 3, 0.125])
        w = -energy
        a = 12 / (1 - step**2 * w / 12) - 10
        derivatives = [q, mp.mpf(1)]
        for n in range(5):
            derivatives.append(w * derivatives[n])
        previous = q
        current = mp.fsum(step**n * value / mp.factorial(n) for n, value in enumerate(derivatives))
        for _ in range(1, 8):
            previous, current = current, a * current - previous
        tv = mp.fsum((-step) ** (2 * n) * w**n / mp.factorial(2 * n) for n in range(4))
        td = mp.fsum((-step) ** (2 * n + 1) * w**n / mp.factorial(2 * n + 1) for n in range(3))
        k = mp.sqrt(energy)
        coef = tv - previous / current
        expected = -(coef * mp.sin(k * endpoint) + td * k * mp.cos(k * endpoint))
        expected /= coef * mp.cos(k * endpoint) - td * k * mp.sin(k * endpoint)
        actual = exact_discrete_scalar(
            float(q), float(energy), float(radius), float(endpoint), float(step)
        )
        assert abs(actual - expected) < mp.mpf("1e-65")
