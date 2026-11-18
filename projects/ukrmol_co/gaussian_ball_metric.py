"""Bounded Cartesian-Gaussian toy Gram on one physical ball, without rescaling.

No native spherical-CGTO adapter or saved-orbital integration is implicit here.
See docs/physics/co-common-domain-orbital-metric.md for the domain contract.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import islice

import numpy as np
from numpy.polynomial.legendre import leggauss


@dataclass(frozen=True)
class CartesianGaussian:
    """Explicit polynomial times Gaussian contraction; coefficients include all norms."""

    center: tuple[float, float, float]
    powers: tuple[int, int, int]
    exponents: tuple[float, ...]
    coefficients: tuple[complex, ...]


def _real_vector(values, size):
    array = np.asarray(values)
    if array.shape != (size,) or not np.isrealobj(array) or not np.isfinite(array).all():
        raise ValueError("Finite real vector with the declared dimension required")
    return array.astype(float)


def _basis(functions):
    # Bound even a generator before materializing a potentially large basis.
    functions = tuple(islice(iter(functions), 33))
    if not 1 <= len(functions) <= 32:
        raise ValueError("Toy basis requires one to32 functions per frame")
    primitives = 0
    for function in functions:
        if not isinstance(function, CartesianGaussian):
            raise ValueError("Explicit CartesianGaussian functions required")
        _real_vector(function.center, 3)
        if (
            len(function.powers) != 3
            or any(type(k) is not int or k < 0 for k in function.powers)
            or sum(function.powers) > 4
        ):
            raise ValueError("Cartesian nonnegative integer powers of total degree<=4 required")
        if not 1 <= len(function.exponents) <= 128:
            raise ValueError("Toy primitive cap exceeded")
        exponents = _real_vector(function.exponents, len(function.exponents))
        coefficients = np.asarray(function.coefficients)
        if (
            np.any(exponents <= 0)
            or coefficients.shape != exponents.shape
            or not np.isfinite(coefficients).all()
        ):
            raise ValueError("Positive exponents and matching finite coefficients required")
        primitives += len(exponents)
    if primitives > 128:
        raise ValueError("Toy basis requires at most128 primitives per frame")
    return functions, primitives


def _values(functions, points):
    result = np.empty((len(points), len(functions)), complex)
    for index, function in enumerate(functions):
        delta = points - function.center
        distance = np.einsum("ij,ij->i", delta, delta)
        polynomial = np.prod(delta**function.powers, axis=1)
        radial = sum(
            coefficient * np.exp(-exponent * distance)
            for exponent, coefficient in zip(function.exponents, function.coefficients, strict=True)
        )
        result[:, index] = polynomial * radial
    if not np.isfinite(result).all():
        raise ValueError("Nonfinite Gaussian values; unsupported scale")
    return result


def gaussian_ball_metric(
    left: list[CartesianGaussian],
    right: list[CartesianGaussian],
    radius: float,
    *,
    center=(0.0, 0.0, 0.0),
    orders=(24, 24, 48),
    max_work: int = 10_000_000,
) -> np.ndarray:
    """Integrate <left|right> on the same ball with positive spherical weights.

    Orders are radial, polar-cosine and azimuthal. Caps limit points, primitives
    and primitive-point work before grid allocation; values use2048-point chunks.
    Positivity/reversal are algebraic properties, not a convergence certificate.
    """
    left, left_primitives = _basis(left)
    right, right_primitives = _basis(right)
    origin = _real_vector(center, 3)
    if (
        not np.isscalar(radius)
        or not np.isrealobj(radius)
        or not np.isfinite(radius)
        or radius <= 0
    ):
        raise ValueError("Finite positive physical radius required")
    if (
        len(orders) != 3
        or any(type(n) is not int or not 2 <= n <= 128 for n in orders)
        or np.prod(orders) > 200_000
    ):
        raise ValueError("Three bounded quadrature orders with at most200000 points required")
    if type(max_work) is not int or not 1 <= max_work <= 10_000_000:
        raise ValueError("Positive bounded primitive-point budget required")
    points_count = int(np.prod(orders))
    if points_count * (left_primitives + right_primitives) > max_work:
        raise ValueError("Primitive-point work budget exceeded")
    radial, wr = leggauss(orders[0])
    polar, wt = leggauss(orders[1])
    radial, wr = (radial + 1) * radius / 2, wr * radius / 2
    azimuth = 2 * np.pi * np.arange(orders[2]) / orders[2]
    r, t, phi = np.meshgrid(radial, polar, azimuth, indexing="ij")
    transverse = np.sqrt(1 - t**2)
    points = (
        np.column_stack(
            [
                (r * transverse * np.cos(phi)).ravel(),
                (r * transverse * np.sin(phi)).ravel(),
                (r * t).ravel(),
            ]
        )
        + origin
    )
    weights = (
        (wr * radial**2)[:, None, None]
        * wt[None, :, None]
        * np.ones(orders[2])[None, None, :]
        * (2 * np.pi / orders[2])
    ).ravel()
    if not np.isfinite(points).all() or not np.isfinite(weights).all() or np.any(weights <= 0):
        raise ValueError("Nonfinite or nonpositive quadrature weights; unsupported scale")
    metric = np.zeros((len(left), len(right)), complex)
    for start in range(0, points_count, 2048):
        batch = points[start : start + 2048]
        bra, ket = _values(left, batch), _values(right, batch)
        metric += bra.conj().T @ (weights[start : start + 2048, None] * ket)
    if not np.isfinite(metric).all():
        raise ValueError("Nonfinite integrated Gram; unsupported scale")
    return metric
