"""Independent point, moment, quadrature and tail-bound controls; generated only."""

from math import gamma, pi

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import gammainc

from projects.ukrmol_co.gaussian_ball_metric import CartesianGaussian, gaussian_ball_metric
from projects.ukrmol_co.native_orbital_metric import (
    SphericalGaussian,
    _polynomials,
    all_space_metric,
    generator_scales,
    radial_overlap,
    shell_overlap,
    shells_from_export,
    solid_harmonics,
    uniform_ball_witness,
)


def shell(angular=0, center=(0.0, 0.0, 0.0), exponents=(0.7,), amplitudes=(1.0,), continuum=False):
    return SphericalGaussian(center, angular, exponents, amplitudes, continuum)


def primitive_norm(angular, exponent):
    return np.sqrt(
        2 * (2 * exponent) ** (angular + 1.5) / gamma(angular + 1.5) * (2 * angular + 1) / (4 * pi)
    )


def test_native_s_p_d_absolute_point_order_sign_and_origin_limits():
    points = np.array([[0.3, -0.7, 0.9], [-0.5, 0.4, 0.2], [0, 0, 0]])
    x, y, z = points.T
    np.testing.assert_allclose(solid_harmonics(0, points), np.ones((3, 1)), rtol=0, atol=1e-15)
    np.testing.assert_allclose(
        solid_harmonics(1, points), np.column_stack([y, z, x]), rtol=1e-13, atol=1e-14
    )
    expected = np.column_stack(
        [
            np.sqrt(3) * x * y,
            np.sqrt(3) * y * z,
            (3 * z * z - x * x - y * y - z * z) / 2,
            np.sqrt(3) * x * z,
            np.sqrt(3) / 2 * (x * x - y * y),
        ]
    )
    np.testing.assert_allclose(solid_harmonics(2, points), expected, rtol=1e-13, atol=1e-14)


@pytest.mark.parametrize("angular", range(11))
def test_all_admitted_degrees_against_explicit_cartesian_polynomials(angular):
    points = np.r_[np.random.default_rng(93).normal(size=(9, 3)) / 2, [[0, 0, 0]]]
    powers, coefficients = _polynomials(angular)
    monomials = np.prod(points[:, None, :] ** powers[None, :, :], axis=-1)
    expected = monomials @ coefficients.T
    np.testing.assert_allclose(solid_harmonics(angular, points), expected, rtol=1e-10, atol=1e-11)


@pytest.mark.parametrize("degrees", [(0, 0), (1, 2), (2, 3), (4, 1), (5, 5), (6, 3)])
def test_displaced_contracted_complex_moments_vs_three_exact_hermite_orders(degrees):
    left = shell(degrees[0], (0.3, -0.2, 0.1), (0.4, 1.2), (0.8 + 0.2j, -0.15j))
    right = shell(degrees[1], (-0.2, 0.4, -0.3), (0.3, 0.9), (-0.2, 0.6 - 0.1j))
    expected = shell_overlap(left, right)
    for extra in range(3):
        actual = shell_overlap(left, right, method="hermite", extra=extra)
        np.testing.assert_allclose(actual, expected, rtol=1e-10, atol=1e-10)
    reverse = shell_overlap(right, left)
    np.testing.assert_allclose(expected, reverse.conj().T, rtol=1e-12, atol=1e-11)


@pytest.mark.parametrize("angular", [0, 1, 2, 4, 6, 10])
def test_all_space_normalization_and_fixed_generator_sphere_restore(angular):
    exponent, radius = 0.07, 2.1
    n = primitive_norm(angular, exponent)
    f = shell(angular, exponents=(exponent,), amplitudes=(n,), continuum=True)
    full = all_space_metric([f])
    np.testing.assert_allclose(full, np.eye(2 * angular + 1), rtol=1e-12, atol=1e-12)
    fraction = gammainc(angular + 1.5, 2 * exponent * radius**2)
    scale = generator_scales([f], radius)
    np.testing.assert_allclose(
        scale, np.ones(2 * angular + 1) / np.sqrt(fraction), rtol=1e-12, atol=1e-12
    )
    restored, bound = uniform_ball_witness([f], full, radius, scale)
    np.testing.assert_allclose(restored, np.eye(2 * angular + 1), rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(bound, 0, rtol=0, atol=0)
    larger, _ = uniform_ball_witness([f], full, radius + 1, scale)
    assert np.min(np.diag(larger).real) > 1


def test_contraction_normalization_and_three_radial_orders_vs_independent_quad():
    f = shell(3, exponents=(0.2, 0.9), amplitudes=(0.7, -0.2), continuum=True)
    radius = 4.3
    reference, error = quad(
        lambda r: (
            4 * pi / 7 * r**8 * (0.7 * np.exp(-0.2 * r * r) - 0.2 * np.exp(-0.9 * r * r)) ** 2
        ),
        0,
        radius,
        epsabs=1e-13,
        epsrel=1e-13,
    )
    assert error < 1e-10
    np.testing.assert_allclose(radial_overlap(f, f, radius), reference, rtol=1e-12, atol=1e-12)
    for order in [48, 96, 192]:
        np.testing.assert_allclose(
            radial_overlap(f, f, radius, order=order), reference, rtol=1e-11, atol=1e-11
        )


def test_mixed_domain_is_not_physical_ball_and_tail_bound_propagates_through_complex_mos():
    functions = [
        shell(0, (0.2, -0.1, 0.3), (0.6,), (0.7,)),
        shell(1, (-0.2, 0.1, -0.1), (0.5,), (0.5,)),
        shell(0, exponents=(0.04,), amplitudes=(0.3,), continuum=True),
    ]
    radius = 1.4
    scales = generator_scales(functions, radius)
    full = all_space_metric(functions)
    hybrid, bound = uniform_ball_witness(functions, full, radius, scales)
    cartesian, transform = [], []
    # Expand each native solid harmonic into the older independent Cartesian
    # finite-ball oracle; identical Gaussian amplitudes, no native evaluator.
    for f in functions:
        powers, matrix = _polynomials(f.angular)
        cartesian.extend(
            [
                CartesianGaussian(f.center, tuple(map(int, power)), f.exponents, f.amplitudes)
                for power in powers
            ]
        )
        transform.append(matrix.T)
    from scipy.linalg import block_diag

    change = block_diag(*transform)
    physical = (
        change.T @ gaussian_ball_metric(cartesian, cartesian, radius, orders=(32, 32, 64)) @ change
    )
    physical = scales[:, None] * physical * scales[None, :]
    assert np.max(abs(hybrid - physical)) > 0.1
    assert np.all(abs(hybrid - physical) <= bound + 1e-10)
    rng = np.random.default_rng(17)
    c = rng.normal(size=(5, 3)) + 1j * rng.normal(size=(5, 3))
    delta = c.conj().T @ (hybrid - physical) @ c
    propagated = abs(c).T @ bound @ abs(c)
    assert np.all(abs(delta) <= propagated + 1e-9)
    assert np.linalg.eigvalsh(physical).min() >= -1e-12


def test_read_amplitudes_and_original_coefficient_array_are_not_modified():
    exponents = np.array([0.2, 0.8])
    contraction = np.array([0.7, -0.1])
    factors = np.array([0.3, 0.5])
    exported = dict(
        shells=[
            dict(
                center=np.zeros(3),
                angular=0,
                exponents=exponents,
                normalization=0.6,
                contractions=contraction,
                primitive_norms=factors,
                continuum=1,
            )
        ],
        shell_descriptors=np.array([[1], [1], [1], [1], [1], [0]]),
    )
    original = [a.copy() for a in [exponents, contraction, factors]]
    functions = shells_from_export(exported)
    np.testing.assert_allclose(functions[0].amplitudes, 0.6 * contraction * factors, rtol=0, atol=0)
    generator_scales(functions, 3.0)
    for before, after in zip(original, [exponents, contraction, factors], strict=True):
        np.testing.assert_array_equal(before, after)


@pytest.mark.parametrize(
    "case", ["degree", "exponent", "amplitude", "center", "continuum", "work", "scale"]
)
def test_malformed_and_work_cap_rejections(case):
    f = shell()
    if case == "degree":
        f = shell(11)
    elif case == "exponent":
        f = shell(exponents=(-1.0,))
    elif case == "amplitude":
        f = shell(amplitudes=(np.nan,))
    elif case == "center":
        f = shell(center=(0.1, 0, 0), continuum=True)
    elif case == "continuum":
        f = shell(continuum=1)
    with pytest.raises(ValueError):
        if case == "work":
            all_space_metric([f], max_work=1)
        elif case == "scale":
            uniform_ball_witness([f], np.eye(1), 1.0, np.zeros(1))
        else:
            all_space_metric([f])
