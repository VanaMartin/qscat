"""Ordered coupled algebra, retained arithmetic and independent value recurrence."""

import mpmath as mp
import numpy as np
import pytest

from projects.ukrmol_co.coupled_tail import tail_coefficient, tail_waves
from projects.ukrmol_co.numerov_coupled_oracle import direct_discrete_k, mp_tail
from projects.ukrmol_co.numerov_typed_coupled import refined_solve, typed_coupled_k, typed_tail


@pytest.mark.parametrize("n", [4, 12])
def test_large_solver_rank_one_reference(n):
    v = np.arange(1, n + 1, dtype=float) / n
    a = np.eye(n) + np.outer(v, v) / n
    rhs = np.column_stack([v, np.ones(n)])
    expected = rhs - np.outer(v, v @ rhs) / (n + v @ v)
    stats = dict(
        maximum_solve_backward_error=0.0, maximum_refinement_corrections=0, large_system_solves=0
    )
    np.testing.assert_allclose(refined_solve(a, rhs, stats), expected, atol=2e-15, rtol=0)
    assert stats["large_system_solves"] == 1


@pytest.mark.parametrize("thresholds", [[0.0, 0.3], [0.0, 1.8]])
def test_coupled_short_discrete_problem_independently_reconstructed(thresholds):
    q = np.array([[0.3, 0.12], [0.12, -0.2]])
    c2 = np.array([[2.0, -0.4], [-0.4, 6.0]])
    c3 = np.array([[0.1, 0.3], [0.3, -0.2]])
    expected = direct_discrete_k(q, thresholds, [0, 0], c2, c3, 1.0, 2.0, 4.0, 0.125, 8, 70)
    measured, stats = typed_coupled_k(q, thresholds, [0, 0], c2, c3, 1.0, 2.0, 4.0, 0.125, 8)
    np.testing.assert_allclose(measured, np.array(expected.tolist(), float), atol=1e-12, rtol=0)
    assert stats["significand_bits"] == 53 and stats["complex_significand_bits"] == 53


def test_general_tail_coefficients_and_derivatives_match_mp_equations():
    ts = np.array([0.0, 0.3, 1.8])
    angular = np.array([0, 1, 2])
    c2 = np.array([[2.0, -0.4, 0.1], [-0.4, 6.0, 0.2], [0.1, 0.2, 1.0]])
    c3 = np.array([[0.1, 0.3, -0.08], [0.3, -0.2, 0.12], [-0.08, 0.12, 0.05]])
    values, derivatives, _ = typed_tail(ts, angular, 1.0, c2, c3, 8, 20.0, 1e-10)
    with mp.workdps(70):
        expected, prime = mp_tail(
            [mp.mpf(float(t)) for t in ts],
            list(angular),
            mp.mpf(1),
            mp.matrix(c2.tolist()),
            mp.matrix(c3.tolist()),
            8,
            mp.mpf(20),
        )
    for left, right in zip((*values, *derivatives), (*expected, *prime), strict=True):
        np.testing.assert_allclose(
            left, np.array(right.tolist(), float).reshape(left.shape), atol=2e-14, rtol=0
        )


def test_degenerate_dipole_tail_against_noninteger_hankel_oracle():
    c2 = tail_coefficient(1.2)
    values, derivatives, _ = typed_tail(
        np.zeros(2), np.array([1, 2]), 1.0, c2, np.zeros((2, 2)), 12, 20.0, 1e-10
    )
    (out, prime), _ = tail_waves(1.0, 20.0, 1.2)
    np.testing.assert_allclose(values[1] + 1j * values[0], out, atol=1e-11, rtol=0)
    np.testing.assert_allclose(derivatives[1] + 1j * derivatives[0], prime, atol=1e-11, rtol=0)


@pytest.mark.skipif(np.finfo(np.longdouble).nmant <= 52, reason="Wider arithmetic absent")
def test_extended_coupled_discrete_error_is_not_rounded_to_binary64():
    q = np.array([[0.3, 0.12], [0.12, -0.2]])
    c2 = np.array([[2.0, -0.4], [-0.4, 6.0]])
    c3 = np.array([[0.1, 0.3], [0.3, -0.2]])
    expected = direct_discrete_k(q, [0.0, 1.8], [0, 0], c2, c3, 1.0, 2.0, 4.0, 0.125, 8, 90)
    k, stats = typed_coupled_k(
        q, [0.0, 1.8], [0, 0], c2, c3, 1.0, 2.0, 4.0, 0.125, 8, dtype=np.longdouble
    )
    with mp.workdps(80):
        retained = mp.matrix([[str(x) for x in row] for row in k])
        assert max(abs(a - b) for a, b in zip(retained, expected, strict=True)) < mp.mpf("1e-14")
    assert stats["significand_bits"] >= 64
