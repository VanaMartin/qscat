"""Analytic boundary signs, threshold behavior and physical outgoing pole controls."""

import mpmath as mp
import numpy as np
import pytest

from projects.ukrmol_co.physical_pole import (
    boundary_s,
    exact_pole,
    exact_s,
    outgoing_residual,
    shell_boundary,
)


def test_free_nodes_and_bloch_normalization():
    for q in (0.3, np.pi, 2 * np.pi):
        state = shell_boundary(q, 0, 2, 1e-11)
        for bloch in (-1, 0, 2):
            np.testing.assert_allclose(boundary_s(q, state * 7, 2, bloch), 1, atol=1e-9, rtol=0)


def test_physical_shell_matches_exact_and_conserves_flux():
    q = 1.3
    state = shell_boundary(q, 10, 2, 1e-11)
    result = boundary_s(q, state, 2)
    with mp.workdps(60):
        expected = complex(exact_s(q, 10))
    np.testing.assert_allclose(result, expected, atol=1e-9, rtol=0)
    np.testing.assert_allclose(abs(result) ** 2, 1, atol=1e-12, rtol=0)


def test_lambert_pole_is_outgoing_not_cancelled_or_threshold_root():
    with mp.workdps(60):
        pole = exact_pole(10, 1)
        jost = 1 + 10 * mp.exp(1j * pole) * mp.sin(pole) / pole
        assert abs(jost) < mp.mpf("1e-55")
        q = complex(pole)
    assert q.real > 0 and q.imag < 0 and (q * q).imag < 0
    state = shell_boundary(q, 10, 2, 1e-11)
    assert abs(outgoing_residual(q, state, 2)) < 1e-8
    assert abs(np.exp(-2j * q) * (state[1] + 1j * q * state[0])) > 0.1
    assert abs(outgoing_residual(0, shell_boundary(0, 10, 1, 1e-11), 1) - 11) < 1e-10
    with pytest.raises(ValueError):
        exact_pole(0, 1)
