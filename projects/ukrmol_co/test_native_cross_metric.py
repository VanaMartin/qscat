"""Independent rectangular Gaussian, reversal and uniform-ball controls."""

import numpy as np
import pytest
from scipy.linalg import block_diag

from projects.ukrmol_co.gaussian_ball_metric import CartesianGaussian, gaussian_ball_metric
from projects.ukrmol_co.native_cross_metric import all_space_cross_metric, cross_ball_witness
from projects.ukrmol_co.native_orbital_metric import (
    SphericalGaussian,
    _polynomials,
    all_space_metric,
    generator_scales,
    uniform_ball_witness,
)


def s(angular=0, center=(0.0, 0.0, 0.0), exponent=0.7, amplitude=1.0, continuum=False):
    return SphericalGaussian(center, angular, (exponent,), (amplitude,), continuum)


def cartesian(functions):
    basis = []
    changes = []
    for f in functions:
        powers, matrix = _polynomials(f.angular)
        basis.extend(
            CartesianGaussian(f.center, tuple(map(int, p)), f.exponents, f.amplitudes)
            for p in powers
        )
        changes.append(matrix.T)
    return basis, block_diag(*changes)


@pytest.mark.parametrize("extra", [0, 1, 2])
def test_rectangular_moments_hermite_and_independent_reverse(extra):
    left = [
        s(0, (0.2, -0.1, 0.3), amplitude=0.8 + 0.2j),
        s(1, exponent=0.2, continuum=True),
        s(2, (-0.4, 0.3, -0.1), 0.9, -0.2j),
    ]
    right = [
        s(1, (-0.3, 0.1, -0.2), 0.4, 0.7 - 0.1j),
        s(0, exponent=0.03, amplitude=0.1, continuum=True),
    ]
    expected = all_space_cross_metric(left, right)
    actual = all_space_cross_metric(left, right, method="hermite", extra=extra)
    assert actual.shape == (9, 4)
    np.testing.assert_allclose(actual, expected, rtol=1e-10, atol=1e-10)
    reverse = all_space_cross_metric(right, left)
    np.testing.assert_allclose(reverse.conj().T, expected, rtol=1e-12, atol=1e-12)


def test_identical_basis_reduces_to_prior_self_prescription_without_rescaling():
    functions = [
        s(0, (0.2, 0.0, 0.0), amplitude=0.6),
        s(1, exponent=0.04, amplitude=0.2 + 0.1j, continuum=True),
    ]
    scales = generator_scales(functions, 3.0)
    cross = all_space_cross_metric(functions, functions)
    original = all_space_metric(functions)
    np.testing.assert_allclose(cross, original, rtol=1e-12, atol=1e-12)
    a, ea = cross_ball_witness(functions, functions, cross, 3.0, scales, scales)
    b, eb = uniform_ball_witness(functions, original, 3.0, scales)
    np.testing.assert_allclose(a, b, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(ea, eb, rtol=1e-12, atol=1e-12)


def test_direct_finite_ball_cross_integral_is_enclosed_after_complex_mo_contraction():
    left = [
        s(0, (0.2, -0.1, 0.3), 0.6, 0.7),
        s(1, exponent=0.08, amplitude=0.2 + 0.1j, continuum=True),
    ]
    right = [s(1, (-0.1, 0.2, -0.2), 0.5, 0.6), s(0, exponent=0.05, amplitude=0.3, continuum=True)]
    radius = 1.6
    da, db = generator_scales(left, radius), generator_scales(right, radius)
    hybrid, error = cross_ball_witness(
        left, right, all_space_cross_metric(left, right), radius, da, db
    )
    a, ta = cartesian(left)
    b, tb = cartesian(right)
    exact = ta.T @ gaussian_ball_metric(a, b, radius, orders=(32, 32, 64)) @ tb
    exact = da[:, None] * exact * db[None, :]
    assert np.max(abs(exact - hybrid)) > 0.01
    assert np.all(abs(exact - hybrid) <= error + 1e-10)
    rng = np.random.default_rng(431)
    ca = rng.normal(size=(4, 3)) + 1j * rng.normal(size=(4, 3))
    cb = rng.normal(size=(4, 2)) + 1j * rng.normal(size=(4, 2))
    delta = ca.conj().T @ (hybrid - exact) @ cb
    bound = abs(ca).T @ error @ abs(cb)
    assert np.all(abs(delta) <= bound + 1e-9)
    reverse, reverse_error = cross_ball_witness(
        right, left, all_space_cross_metric(right, left), radius, db, da
    )
    np.testing.assert_allclose(reverse.conj().T, hybrid, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(reverse_error.T, error, rtol=1e-12, atol=1e-12)
    a_phase, b_phase = np.exp(0.3j), np.exp(-0.6j)
    np.testing.assert_allclose(
        (a_phase * ca).conj().T @ hybrid @ (b_phase * cb),
        a_phase.conjugate() * b_phase * (ca.conj().T @ hybrid @ cb),
        rtol=1e-12,
        atol=1e-12,
    )


def test_distinct_fixed_generator_scales_do_not_follow_the_comparison_radius():
    left = [s(0, exponent=0.04, amplitude=0.3, continuum=True)]
    right = [s(0, exponent=0.06, amplitude=0.4, continuum=True)]
    da, db = generator_scales(left, 2.0), generator_scales(right, 3.0)
    gram = all_space_cross_metric(left, right)
    small, _ = cross_ball_witness(left, right, gram, 2.0, da, db)
    large, error = cross_ball_witness(left, right, gram, 4.0, da, db)
    assert large[0, 0].real > small[0, 0].real
    np.testing.assert_allclose(error, 0, rtol=0, atol=0)


@pytest.mark.parametrize(
    "case",
    ["work", "method", "extra", "shape", "complex_scale", "negative_scale", "radius", "origin"],
)
def test_malformed_and_preallocation_work_guards(case):
    left = right = [s()]
    with pytest.raises(ValueError):
        if case == "work":
            all_space_cross_metric(left, right, max_work=1)
        elif case == "method":
            all_space_cross_metric(left, right, method="unknown")
        elif case == "extra":
            all_space_cross_metric(left, right, extra=3)
        elif case == "shape":
            cross_ball_witness(left, right, np.eye(2), 2.0, [1.0], [1.0])
        elif case == "complex_scale":
            cross_ball_witness(left, right, np.eye(1), 2.0, [1j], [1.0])
        elif case == "negative_scale":
            cross_ball_witness(left, right, np.eye(1), 2.0, [-1.0], [1.0])
        elif case == "radius":
            cross_ball_witness(left, right, np.eye(1), 0.0, [1.0], [1.0])
        else:
            all_space_cross_metric([s(center=(0.1, 0.0, 0.0), continuum=True)], right)
