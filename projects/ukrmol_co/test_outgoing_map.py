"""Independent analytic, finite-difference/BVP and conservation controls."""

import unittest

import numpy as np
from scipy.linalg import solve

from projects.ukrmol_co.exterior_numerov import numerov_k
from projects.ukrmol_co.outgoing_map import outgoing_map, physical_s, sheet_momenta
from projects.ukrmol_co.saved_matrices import scattering_matrix


def test_free_roots_and_analytic_s_normalization():
    thresholds = np.array([0.0, 2.0])
    for energy in [0.7, 0.7 - 0.1j]:
        y, a, _ = outgoing_map(
            thresholds,
            [0, 0],
            np.zeros((2, 2)),
            np.zeros((2, 2)),
            energy,
            2,
            18,
            30,
            rtol=1e-12,
            amplitudes=energy == 0.7,
        )
        k, _ = sheet_momenta(energy, thresholds, 0.7)
        np.testing.assert_allclose(y, np.diag(1j * k), atol=1e-12, rtol=0)
        if energy == 0.7:
            q = np.diag([np.tan(k[0] * 2).real / k[0].real, 0.3])
            s, diagnostics = physical_s(q, y, a)
            np.testing.assert_allclose(s, np.ones((1, 1)), atol=2e-9, rtol=0)
            assert diagnostics["flux_defect"] < 2e-9
            wall, _ = physical_s(np.zeros((2, 2)), y, a)
            np.testing.assert_allclose(wall, [[np.exp(-4j * k[0])]], atol=2e-9, rtol=0)


def test_coupled_real_s_against_independent_boundary_numerov():
    thresholds = np.array([0.0, 0.0, 2.0])
    c2 = np.array([[2.0, 0.15, 0.1], [0.15, 6.0, 0.2], [0.1, 0.2, 2.0]])
    c3 = np.array([[0.2, 0.1, 0.05], [0.1, -0.3, 0.1], [0.05, 0.1, 0.4]])
    q = np.array([[0.4, 0.03, 0.02], [0.03, 0.5, 0.01], [0.02, 0.01, 0.3]])
    y, a, _ = outgoing_map(thresholds, [1, 2, 1], c2, c3, 0.7, 3, 35, 30, amplitudes=True)
    actual, diagnostics = physical_s(q, y, a)
    reference, _ = numerov_k(q, thresholds, np.array([1, 2, 1]), c2, c3, 0.7, 3, 35, 0.01, 8)
    np.testing.assert_allclose(actual, scattering_matrix(reference), atol=2e-8, rtol=0)
    for key in ["flux_defect", "unitarity_defect", "reciprocity_defect", "wronskian_residual"]:
        assert diagnostics[key] < 2e-8


def test_three_tolerances_against_exact_p_wave_hankel():
    energy = 0.7 - 0.1j
    k = np.sqrt(energy)
    radius = 2.0
    expected = 1j * k + (-1j / (k * radius**2)) / (1 + 1j / (k * radius))
    errors = []
    for tolerance in [1e-8, 1e-10, 1e-12]:
        y, _, _ = outgoing_map([0], [1], [[2]], [[0]], energy, radius, 30, 35, rtol=tolerance)
        errors.append(abs(y[0, 0] - expected))
    assert errors[-1] < errors[0]
    assert errors[-1] < 1e-10


def test_complex_coupled_map_against_independent_second_order_bvp():
    # Independent radial operator: no inverse-power recurrence or Riccati evolution.
    energy = 0.7 - 0.1j
    thresholds = np.array([0.0, 2.0])
    c2 = np.array([[2.0, 0.3], [0.3, 6.0]])
    c3 = np.array([[0.1, 0.2], [0.2, -0.1]])
    actual, _, _ = outgoing_map(thresholds, [1, 2], c2, c3, energy, 2, 40, 35, order=10)
    # Coefficient-space polynomial differentiation, independent of the toy's
    # barycentric grid helper (and its unrelated high-precision dependency).
    degree = 192
    x = np.cos(np.pi * np.arange(degree + 1) / degree)
    basis = np.polynomial.chebyshev.chebvander(x, degree)
    differentiated = np.column_stack(
        [
            np.polynomial.chebyshev.chebval(x, np.polynomial.chebyshev.chebder(column))
            for column in np.eye(degree + 1)
        ]
    )
    t = 1 + 39 * (x + 1) / 2
    d = solve(basis.T, differentiated.T).T * 2 / 39
    rotation = np.exp(1j * np.deg2rad(35))
    r = 2 + (t - 1) * rotation * (38 / 39)
    derivative = d * (39 / 38) / rotation
    second = derivative @ derivative
    count = len(t)
    matrix = np.kron(np.eye(2), second).astype(complex)
    k = np.array([np.sqrt(energy), 1j * np.sqrt(2 - energy)])
    for i in range(2):
        for j in range(2):
            value = c2[i, j] / r**2 + c3[i, j] / r**3
            if i == j:
                value = value + thresholds[i] - energy
            matrix[i * count : (i + 1) * count, j * count : (j + 1) * count] -= np.diag(value)
    rhs = np.zeros((2 * count, 2), complex)
    for i in range(2):
        first, last = i * count, (i + 1) * count - 1
        matrix[first] = 0
        matrix[first, i * count : (i + 1) * count] = derivative[0]
        matrix[first, first] -= 1j * k[i]
        matrix[last] = 0
        matrix[last, last] = 1
        rhs[last, i] = 1
    values = solve(matrix, rhs)
    expected = np.stack([derivative[-1] @ values[i * count : (i + 1) * count] for i in range(2)])
    np.testing.assert_allclose(actual, expected, atol=2e-7, rtol=0)


def test_invalid_sheet_contour_and_caps_reject():
    case = unittest.TestCase()
    with case.assertRaisesRegex(ValueError, "Threshold"):
        sheet_momenta(0, [0, 2], 0.5)
    with case.assertRaisesRegex(ValueError, "interval"):
        sheet_momenta(3 - 0.1j, [0, 2], 0.5)
    with case.assertRaisesRegex(ValueError, "Contour"):
        outgoing_map([0], [0], [[0]], [[0]], 0.1 - 0.2j, 2, 5, 1)
    with case.assertRaises(TimeoutError):
        outgoing_map([0], [0], [[0]], [[0]], 0.5, 2, 5, 30, max_nfev=1)


if __name__ == "__main__":
    suite = unittest.TestSuite(
        unittest.FunctionTestCase(fn)
        for name, fn in list(globals().items())
        if name.startswith("test_") and callable(fn)
    )
    raise SystemExit(0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1)
