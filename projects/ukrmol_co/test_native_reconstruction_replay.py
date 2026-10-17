"""Analytic checks for the separate saved-decision replay and K-coordinate diagnosis."""

import numpy as np
import pytest

from projects.ukrmol_co.native_reconstruction_replay import compare, coordinate_diagnostic


def test_large_tangent_channel_amplifies_k_without_a_large_s_change():
    native = np.diag([-185.0, 0.1, -0.2])
    changed = native.copy()
    changed[0, 0] += 1e-5
    report = coordinate_diagnostic(changed, native)
    expected_s_change = 2e-5 / np.sqrt((1 + 185**2) * (1 + (-185 + 1e-5) ** 2))
    assert report["differences"]["k"] > 1e-6
    assert report["differences"]["s"] == pytest.approx(expected_s_change, rel=1e-6)
    assert report["large_channel_delta_frobenius_fraction"] == pytest.approx(1)
    assert report["cayley_difference_identity_max_residual"] < 1e-14


def test_cayley_difference_identity_preserves_matrix_order():
    # Noncommuting matrices distinguish the two inverse factors.
    native = np.array([[0.3, 0.2], [0.2, -0.5]])
    changed = np.array([[0.4, -0.1], [-0.1, -0.2]])
    report = coordinate_diagnostic(changed, native)
    assert report["cayley_difference_identity_max_residual"] < 1e-14


def test_phase_wrap_is_modulo_pi():
    delta = compare(np.diag([1e10, 0.0]), np.diag([-1e10, 0.0]))
    assert delta["k"] > 1e10
    assert delta["phase_rad_mod_pi"] == pytest.approx(2e-10, rel=1e-5)
