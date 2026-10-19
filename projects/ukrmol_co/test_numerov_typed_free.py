"""Pivoting, precision retention and independent finite-grid oracle controls."""

import mpmath as mp
import numpy as np
import pytest

from projects.ukrmol_co.numerov_discrete_free import discrete_free_k
from projects.ukrmol_co.numerov_typed_free import small_solve, typed_free_k


@pytest.mark.parametrize("dtype", [np.float64, np.longdouble])
def test_small_solver_requires_pivoting_and_preserves_input_precision(dtype):
    a = np.array([[0, 2, 3], [4, 0, 6], [7, 8, 1]], dtype=dtype)
    expected = np.array([[1, 2], [-2, 1], [3, -1]], dtype=dtype)
    answer = small_solve(a, a @ expected)
    assert answer.dtype == a.dtype
    np.testing.assert_allclose(answer, expected, rtol=0, atol=1e-14)


def test_singular_small_solver_fails_instead_of_returning_a_matrix():
    with pytest.raises(np.linalg.LinAlgError):
        small_solve(np.ones((2, 2)), np.eye(2))


@pytest.mark.parametrize("thresholds", [[0.0], [0.0, 0.2], [0.0, 0.0, 1.2]])
def test_typed_short_grid_matches_independent_matrix_discrete_oracle(thresholds):
    n = len(thresholds)
    q = np.full((n, n), 0.08)
    np.fill_diagonal(q, [0.4, 0.3, 0.2][:n])
    measured, diagnostics = typed_free_k(q, thresholds, 0.7, 2.0, 3.0, 0.125, np.float64)
    oracle = discrete_free_k(q, thresholds, 0.7, 2.0, 3.0, 0.125, 90)
    np.testing.assert_allclose(measured, np.array(oracle.tolist(), float), rtol=0, atol=2e-12)
    assert diagnostics["significand_bits"] == 53


@pytest.mark.skipif(np.finfo(np.longdouble).nmant <= 52, reason="No wider longdouble on this host")
def test_extended_solver_retains_bits_beyond_binary64():
    with mp.workdps(80):
        tiny = np.longdouble(2) ** -60
        a = np.array([[1, tiny], [tiny, 1]], dtype=np.longdouble)
        b = np.array([[1], [1]], dtype=np.longdouble)
        answer = small_solve(a, b)
        assert answer[0, 0] != np.longdouble(1)
        assert abs(answer[0, 0] - (1 - tiny)) < np.longdouble(2) ** -63
