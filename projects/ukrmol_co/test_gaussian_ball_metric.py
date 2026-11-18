"""Analytic Gaussian product/CDF and polynomial moments versus ball quadrature."""

from itertools import product
from math import comb, gamma, pi

import numpy as np
import pytest
from scipy.special import gammainc
from scipy.stats import ncx2

from projects.ukrmol_co.four_block_overlap import InteriorExpansion, four_block_overlap
from projects.ukrmol_co.gaussian_ball_metric import CartesianGaussian, gaussian_ball_metric


def orbital(center=(0.0, 0.0, 0.0), powers=(0, 0, 0), exponents=(0.7,), coefficients=(1.0,)):
    return CartesianGaussian(center, powers, exponents, coefficients)


def s_ball(left, right, radius, center):
    result = 0j
    a, b, o = map(np.asarray, [left.center, right.center, center])
    for alpha, ca in zip(left.exponents, left.coefficients, strict=True):
        for beta, cb in zip(right.exponents, right.coefficients, strict=True):
            p = alpha + beta
            location = (alpha * a + beta * b) / p
            total = (pi / p) ** 1.5 * np.exp(-alpha * beta / p * np.sum((a - b) ** 2))
            fraction = ncx2.cdf(2 * p * radius**2, 3, 2 * p * np.sum((location - o) ** 2))
            result += np.conj(ca) * cb * total * fraction
    return result


def centered_ball(left, right, radius):
    powers = np.array(left.powers) + right.powers
    if np.any(powers % 2):
        return 0.0
    half = powers // 2
    degree = sum(half) + 1.5
    angular = np.prod([gamma(k + 0.5) for k in half])
    return sum(
        np.conj(ca)
        * cb
        * angular
        / (alpha + beta) ** degree
        * gammainc(degree, (alpha + beta) * radius**2)
        for alpha, ca in zip(left.exponents, left.coefficients, strict=True)
        for beta, cb in zip(right.exponents, right.coefficients, strict=True)
    )


def all_space(left, right):
    result = 0j
    for alpha, ca in zip(left.exponents, left.coefficients, strict=True):
        for beta, cb in zip(right.exponents, right.coefficients, strict=True):
            p = alpha + beta
            moments = []
            for a, b, ka, kb in zip(
                left.center, right.center, left.powers, right.powers, strict=True
            ):
                location = (alpha * a + beta * b) / p
                value = 0.0
                for i, j in product(range(ka + 1), range(kb + 1)):
                    if (i + j) % 2 == 0:
                        moment = gamma((i + j + 1) / 2) / p ** ((i + j + 1) / 2)
                        value += (
                            comb(ka, i)
                            * comb(kb, j)
                            * (location - a) ** (ka - i)
                            * (location - b) ** (kb - j)
                            * moment
                        )
                moments.append(value * np.exp(-alpha * beta / p * (a - b) ** 2))
            result += np.conj(ca) * cb * np.prod(moments)
    return result


@pytest.mark.parametrize("radius", [0.8, 2.7])
def test_displaced_complex_contracted_s_three_level_analytic_control(radius):
    left = orbital((0.6, -0.4, 0.3), exponents=(0.3, 1.1), coefficients=(1 + 0.2j, -0.4j))
    right = orbital((-0.3, 0.2, -0.7), exponents=(0.2, 0.9), coefficients=(-0.6, 0.4 + 0.1j))
    center = (0.1, -0.1, 0.2)
    expected = s_ball(left, right, radius, center)
    errors = []
    for orders in [(6, 6, 12), (12, 12, 24), (24, 24, 48)]:
        actual = gaussian_ball_metric([left], [right], radius, center=center, orders=orders)[0, 0]
        errors.append(abs(actual - expected))
    assert errors[-1] <= 1e-9
    assert errors[1] <= max(errors[0] * 0.1, 1e-13)
    assert errors[2] <= max(errors[1] * 0.1, 1e-13)


@pytest.mark.parametrize(
    "powers", [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (2, 0, 0), (1, 1, 0), (0, 0, 2)]
)
def test_centered_s_p_d_finite_ball_moments(powers):
    left = orbital(powers=powers, exponents=(0.4, 1.0), coefficients=(0.8, -0.2))
    basis = [left, orbital(), orbital(powers=(1, 0, 0))]
    actual = gaussian_ball_metric(basis, basis, 1.6, orders=(24, 20, 40))
    expected = np.array([[centered_ball(a, b, 1.6) for b in basis] for a in basis])
    np.testing.assert_allclose(actual, expected, rtol=1e-11, atol=1e-12)


def test_common_translation_reverse_positivity_and_complex_coefficient_covariance():
    basis = [
        orbital((0.1, -0.3, 0.2)),
        orbital((0.1, -0.3, 0.2), (1, 0, 0)),
        orbital((-0.2, 0.1, -0.4), (0, 1, 1)),
    ]
    other = [orbital((0.4, 0.2, -0.1)), orbital((0.4, 0.2, -0.1), (0, 1, 0))]
    center = np.array([0.1, 0.2, 0.3])
    metric = gaussian_ball_metric(basis, other, 2.0, center=center)
    reverse = gaussian_ball_metric(other, basis, 2.0, center=center)
    np.testing.assert_allclose(metric, reverse.conj().T, rtol=1e-13, atol=1e-13)
    shift = np.array([2.1, -0.7, 0.4])
    moved = [
        orbital(tuple(np.array(f.center) + shift), f.powers, f.exponents, f.coefficients)
        for f in basis + other
    ]
    translated = gaussian_ball_metric(moved[:3], moved[3:], 2.0, center=center + shift)
    np.testing.assert_allclose(translated, metric, rtol=1e-12, atol=1e-12)
    gram = gaussian_ball_metric(basis, basis, 2.0, center=center)
    assert np.linalg.eigvalsh(gram).min() >= -1e-12
    rng = np.random.default_rng(27)
    coefficients = rng.normal(size=(3, 2)) + 1j * rng.normal(size=(3, 2))
    transformed = coefficients.conj().T @ gram @ coefficients
    assert np.linalg.eigvalsh(transformed).min() >= -1e-12
    phases = np.exp(1j * np.array([0.3, -0.8, 1.2]))
    phased = [
        orbital(f.center, f.powers, f.exponents, tuple(phases[i] * np.array(f.coefficients)))
        for i, f in enumerate(basis)
    ]
    phased_gram = gaussian_ball_metric(phased, phased, 2.0, center=center)
    np.testing.assert_allclose(
        phased_gram, phases.conj()[:, None] * gram * phases, rtol=1e-12, atol=1e-12
    )


def test_displaced_polynomials_all_space_limit_and_nonidentity_domain():
    basis = [
        orbital((0.3, -0.2, 0.1)),
        orbital((-0.3, 0.4, -0.2), (1, 0, 0)),
        orbital((0.1, 0.2, 0.3), (0, 1, 1)),
    ]
    full = np.array([[all_space(a, b) for b in basis] for a in basis])
    large = gaussian_ball_metric(basis, basis, 7.0, orders=(48, 32, 64))
    small = gaussian_ball_metric(basis, basis, 0.8, orders=(24, 24, 48))
    np.testing.assert_allclose(large, full, rtol=1e-10, atol=1e-11)
    assert np.linalg.eigvalsh(full - small).min() >= -1e-10
    assert np.max(abs(full - small)) > 0.1
    # No whitening: rescaling the function really changes the returned metric.
    doubled = orbital(coefficients=(2.0,))
    normal = gaussian_ball_metric([orbital()], [orbital()], 0.8)
    scaled = gaussian_ball_metric([doubled], [doubled], 0.8)
    np.testing.assert_allclose(scaled, 4 * normal, rtol=1e-13, atol=1e-13)


def test_separate_radius_family_detects_diffuse_tail_not_quadrature_error():
    tight = orbital(exponents=(2.0,))
    diffuse = orbital(exponents=(0.04,))
    fractions = []
    for radius in [2.0, 4.0, 6.0]:
        fractions.append(
            [
                gaussian_ball_metric([f], [f], radius, orders=(48, 24, 48))[0, 0].real
                / all_space(f, f).real
                for f in [tight, diffuse]
            ]
        )
    expected = [
        [gammainc(1.5, 2 * f.exponents[0] * r**2) for f in [tight, diffuse]]
        for r in [2.0, 4.0, 6.0]
    ]
    np.testing.assert_allclose(fractions, expected, rtol=1e-11, atol=1e-12)
    assert fractions[0][0] > 0.99999 and fractions[0][1] < 0.2
    assert 1 - fractions[-1][1] > 0.1


def test_explicit_real_spherical_s_p_d_transform_without_native_phase_assumption():
    # These explicit Cartesian combinations define OUR angular convention.
    # They must not be inferred from an identity Gram in a native export.
    powers = [
        (0, 0, 0),
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1),
        (2, 0, 0),
        (0, 2, 0),
        (0, 0, 2),
        (1, 1, 0),
        (1, 0, 1),
        (0, 1, 1),
    ]
    alpha, radius = 0.7, 1.3
    basis = [orbital(powers=k, exponents=(alpha,)) for k in powers]
    transform = np.zeros((10, 9))
    transform[0, 0] = 1 / np.sqrt(4 * pi)
    transform[1:4, 1:4] = np.eye(3) * np.sqrt(3 / (4 * pi))
    transform[[4, 5, 6], 4] = np.sqrt(5 / (16 * pi)) * np.array([-1, -1, 2])
    transform[[4, 5], 5] = np.sqrt(15 / (16 * pi)) * np.array([1, -1])
    transform[7:10, 6:9] = np.eye(3) * np.sqrt(15 / (4 * pi))
    angular = np.array([0] + [1] * 3 + [2] * 5)
    normalization = np.array(
        [np.sqrt(2 * (2 * alpha) ** (degree + 1.5) / gamma(degree + 1.5)) for degree in angular]
    )
    transform *= normalization
    cartesian = gaussian_ball_metric(basis, basis, radius)
    gram = transform.T @ cartesian @ transform
    expected = np.diag(gammainc(angular + 1.5, 2 * alpha * radius**2))
    np.testing.assert_allclose(gram, expected, rtol=1e-11, atol=1e-12)
    wrong_phase = np.eye(9)
    wrong_phase[1, 1] = -1
    np.testing.assert_allclose(wrong_phase @ gram @ wrong_phase, gram, rtol=1e-11, atol=1e-12)
    assert np.max(abs(gram @ wrong_phase - gram)) > 0.1


def test_l2_gaussian_tail_cauchy_bound_is_not_zero_for_diffuse_orbitals():
    left = orbital((0.4, -0.3, 0.1), exponents=(0.15,))
    right = orbital((-0.2, 0.1, 0.3), exponents=(0.25,))
    for radius in [1.0, 2.0, 4.0]:
        cross_tail = all_space(left, right) - s_ball(left, right, radius, (0, 0, 0))
        self_tails = [
            all_space(f, f).real - s_ball(f, f, radius, (0, 0, 0)).real for f in [left, right]
        ]
        assert abs(cross_tail) <= np.sqrt(np.prod(self_tails)) + 1e-12
    assert abs(all_space(left, right) - s_ball(left, right, 1.0, (0, 0, 0))) > 1


def test_common_rotation_and_independent_radial_angular_refinement():
    left = orbital((0.8, -0.5, 0.4), exponents=(0.6,))
    right = orbital((-0.4, 0.2, -0.6), exponents=(0.9,))
    center = np.array([0.1, -0.1, 0.2])
    expected = s_ball(left, right, 2.5, center)
    sequences = []
    for orders_family in [
        [(6, 24, 48), (12, 24, 48), (24, 24, 48)],
        [(24, 4, 8), (24, 8, 16), (24, 16, 32)],
    ]:
        errors = [
            abs(
                gaussian_ball_metric([left], [right], 2.5, center=center, orders=orders)[0, 0]
                - expected
            )
            for orders in orders_family
        ]
        assert errors[-1] <= 1e-9
        assert errors[1] <= max(errors[0] * 0.1, 1e-13)
        assert errors[2] <= max(errors[1] * 0.1, 1e-13)
        sequences.append(errors)
    rotation = np.linalg.qr(np.random.default_rng(31).normal(size=(3, 3)))[0]
    moved = [
        orbital(tuple(rotation @ f.center), f.powers, f.exponents, f.coefficients)
        for f in [left, right]
    ]
    metric = gaussian_ball_metric([moved[0]], [moved[1]], 2.5, center=rotation @ center)
    np.testing.assert_allclose(metric, [[expected]], rtol=1e-11, atol=1e-12)
    # Moving only the domain is a different observable, not a gauge operation.
    wrong = gaussian_ball_metric([left], [right], 2.5, center=center + np.array([2, 0, 0]))[0, 0]
    assert abs(wrong - expected) > 0.1


def test_ball_metric_four_blocks_against_direct_two_electron_quadrature():
    # Independent explicit discrete two-coordinate wavefunctions on a small
    # physical grid; no determinant-overlap or removal-state code in the oracle.
    from numpy.polynomial.legendre import leggauss

    radius, orders = 1.4, (6, 6, 12)
    radial, wr = leggauss(orders[0])
    polar, wt = leggauss(orders[1])
    points, weights = [], []
    for i, r in enumerate((radial + 1) * radius / 2):
        for j, t in enumerate(polar):
            for k in range(orders[2]):
                phi = 2 * pi * k / orders[2]
                points.append(
                    [
                        r * np.sqrt(1 - t * t) * np.cos(phi),
                        r * np.sqrt(1 - t * t) * np.sin(phi),
                        r * t,
                    ]
                )
                weights.append(wr[i] * radius / 2 * wt[j] * 2 * pi / orders[2] * r * r)
    points, weights = np.array(points), np.array(weights)
    left = [orbital(), orbital(powers=(1, 0, 0)), orbital((0.3, 0.2, -0.1), exponents=(0.2,))]
    right = [
        orbital((-0.1, 0.2, 0.1)),
        orbital(powers=(1, 0, 0)),
        orbital((0.2, -0.3, 0.1), exponents=(0.3,)),
    ]
    a = InteriorExpansion(
        {(0, 1): 0.7 + 0.2j}, ({(0,): 1}, {(1,): 1}), np.array([[0.4 - 0.1j], [-0.3j]]), 2, 1
    )
    b = InteriorExpansion(
        {(0, 1): -0.2 + 0.1j}, ({(0,): 1}, {(1,): 1}), np.array([[0.2j], [0.5]]), 2, 1
    )

    def direct(functions, state):
        values = []
        for f in functions:
            delta = points - f.center
            polynomial = (
                delta[:, 0] ** f.powers[0] * delta[:, 1] ** f.powers[1] * delta[:, 2] ** f.powers[2]
            )
            radial = sum(
                c * np.exp(-alpha * np.sum(delta * delta, axis=1))
                for alpha, c in zip(f.exponents, f.coefficients, strict=True)
            )
            values.append(np.sqrt(weights) * polynomial * radial)

        def wedge(i, j):
            return (np.outer(values[i], values[j]) - np.outer(values[j], values[i])) / np.sqrt(2)

        return dict(
            q=state.bound[(0, 1)] * wedge(0, 1),
            p=-state.amplitudes[0, 0] * wedge(0, 2) - state.amplitudes[1, 0] * wedge(1, 2),
        )

    physical = four_block_overlap(a, b, gaussian_ball_metric(left, right, radius, orders=orders))
    bra, ket = direct(left, a), direct(right, b)
    for x, y in product("qp", repeat=2):
        np.testing.assert_allclose(physical[x + y], np.vdot(bra[x], ket[y]), rtol=1e-12, atol=1e-12)
    total = np.vdot(bra["q"] + bra["p"], ket["q"] + ket["p"])
    np.testing.assert_allclose(physical["total"], total, rtol=1e-12, atol=1e-12)
    norms = [
        four_block_overlap(s, s, gaussian_ball_metric(f, f, radius, orders=orders))["total"].real
        for s, f in [(a, left), (b, right)]
    ]
    assert abs(total) ** 2 <= np.prod(norms) + 1e-12


@pytest.mark.parametrize(
    "change", ["radius", "center", "orders", "basis", "exponent", "powers", "coefficient", "work"]
)
def test_malformed_and_work_caps_reject(change):
    args = dict(left=[orbital()], right=[orbital()], radius=1.0)
    if change == "radius":
        args["radius"] = -1
    elif change == "center":
        args["center"] = (np.nan, 0, 0)
    elif change == "orders":
        args["orders"] = (80, 60, 128)
    elif change == "basis":
        args["left"] = [orbital()] * 33
    elif change == "exponent":
        args["left"] = [orbital(exponents=(-1.0,))]
    elif change == "powers":
        args["left"] = [orbital(powers=(5, 0, 0))]
    elif change == "coefficient":
        args["right"] = [orbital(coefficients=(np.inf,))]
    elif change == "work":
        args["max_work"] = 1
    with pytest.raises(ValueError):
        gaussian_ball_metric(**args)
