"""Independent saved-reader full-Slater contraction controls."""

import numpy as np

from projects.ukrmol_co.verify_stretched_target_overlap import full_overlap


def test_full_slater_reader_keeps_doubly_occupied_core_determinant():
    left = np.arange(1, 9).reshape(8, 1, 1).astype(complex) * (1 + 1j)
    right = np.arange(3, 11).reshape(8, 1, 1).astype(complex) * (2 - 1j)
    metric = np.diag([0.8, 0.9] + [1.0] * 12)
    actual = full_overlap(left, right, metric, (0, 0))
    expected = (0.8 * 0.9) ** 2 * np.outer(left.ravel().conj(), right.ravel())
    np.testing.assert_allclose(actual, expected, atol=1e-13, rtol=1e-13)


def test_full_slater_reader_preserves_active_orbital_signs():
    left = np.zeros((8, 12, 1))
    left[:, 0, 0] = np.arange(1, 9)
    metric = np.eye(14)
    metric[2, 2] = -1
    actual = full_overlap(left, left, metric, (1, 0))
    np.testing.assert_allclose(
        actual, -np.outer(np.arange(1, 9), np.arange(1, 9)), atol=1e-13, rtol=0
    )
