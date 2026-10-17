"""Outgoing exterior DtN maps from Riccati evolution and spectral collocation."""

from __future__ import annotations

import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import solve

from projects.ukrmol_co.coupled_tail import momenta, tail_coefficient


def outgoing_momenta(energy, separation=16.0):
    """Continue the specified physical sheet and reject a threshold branch point."""
    if energy == 0 or energy == separation:
        raise ValueError("Exact threshold is outside this benchmark")
    return momenta(energy, separation=separation)


def riccati_map(energy, dipole, endpoint, angle, rtol, separation=16.0):
    """Return d/dz at z=1 for the outgoing exterior boundary-value problem."""
    rotation = np.exp(1j * np.deg2rad(angle))
    k = outgoing_momenta(energy, separation)
    if np.min((rotation * k).imag) <= 0:
        raise ValueError("Contour does not expose all outgoing/decaying channels")
    constant = np.diag([0, separation]) - energy * np.eye(2)
    coefficient = tail_coefficient(dipole)

    def rhs(t, state):
        y = state.reshape(2, 2)
        z = 1 + (t - 1) * rotation
        return (rotation * (constant + coefficient / z**2 - y @ y)).ravel()

    result = solve_ivp(
        rhs,
        (endpoint, 1),
        np.diag(1j * k).ravel(),
        method="DOP853",
        rtol=rtol,
        atol=rtol * 0.01,
    )
    if not result.success or not np.isfinite(result.y).all():
        raise ValueError(f"Riccati solve failed: {result.message}")
    return result.y[:, -1].reshape(2, 2), int(result.nfev)


def chebyshev_grid(degree, endpoint):
    """Lobatto nodes decreasing from endpoint to 1, with the t derivative."""
    index = np.arange(degree + 1)
    x = np.cos(np.pi * index / degree)
    weights = np.ones(degree + 1)
    weights[[0, -1]] = 2
    weights *= (-1.0) ** index
    differences = x[:, None] - x[None, :]
    derivative = (weights[:, None] / weights) / (differences + np.eye(degree + 1))
    derivative -= np.diag(derivative.sum(axis=1))
    return 1 + (endpoint - 1) * (x + 1) / 2, derivative * 2 / (endpoint - 1)


def spectral_map(energy, dipole, endpoint, angle, degree, separation=16.0):
    """Solve a second-order linear BVP independently of Riccati propagation."""
    rotation = np.exp(1j * np.deg2rad(angle))
    k = outgoing_momenta(energy, separation)
    if np.min((rotation * k).imag) <= 0:
        raise ValueError("Contour does not expose all outgoing/decaying channels")
    t, derivative = chebyshev_grid(degree, endpoint)
    second = derivative @ derivative / rotation**2
    z = 1 + (t - 1) * rotation
    coefficient = tail_coefficient(dipole)
    count = degree + 1
    matrix = np.kron(np.eye(2), second).astype(complex)
    for i in range(2):
        for j in range(2):
            block = slice(i * count, (i + 1) * count), slice(j * count, (j + 1) * count)
            potential = coefficient[i, j] / z**2
            if i == j:
                potential += [0, separation][i] - energy
            matrix[block] -= np.diag(potential)
    rhs = np.zeros((2 * count, 2), complex)
    for i in range(2):
        # First node is the distant endpoint; the last node is the boundary z=1.
        first, last = i * count, (i + 1) * count - 1
        matrix[first] = 0
        matrix[first, i * count : (i + 1) * count] = derivative[0] / rotation
        matrix[first, first] -= 1j * k[i]
        matrix[last] = 0
        matrix[last, last] = 1
        rhs[last, i] = 1
    values = solve(matrix, rhs)
    residual = float(np.max(abs(matrix @ values - rhs)))
    # Infinity-norm condition estimation uses the already small dense BVP matrix.
    condition = float(np.linalg.cond(matrix, p=np.inf).real)
    y = np.stack(
        [derivative[-1] @ values[i * count : (i + 1) * count] / rotation for i in range(2)]
    )
    return y, condition, residual


def exact_map(energy, dipole, separation=16.0):
    """High-precision outgoing Hankel limits for zero coupling or equal thresholds."""
    energy, dipole = mp.mpc(energy), mp.mpf(dipole)
    if separation != 0 and dipole != 0:
        raise ValueError("There is no analytic coupled oracle here")
    if separation == 0:
        coefficient = mp.matrix([[2, -2 * dipole / mp.sqrt(5)], [-2 * dipole / mp.sqrt(5), 6]])
        eigenvalues, rotation = mp.eigsy(coefficient)
        orders = [mp.sqrt(value + mp.mpf("0.25")) for value in eigenvalues]
        ks = [mp.sqrt(energy)] * 2
    else:
        rotation = mp.eye(2)
        orders = [mp.mpf("1.5"), mp.mpf("2.5")]
        ks = [
            mp.sqrt(energy),
            1j * mp.sqrt(separation - energy)
            if energy.real < separation
            else mp.sqrt(energy - separation),
        ]
    values = []
    for order, k in zip(orders, ks, strict=True):
        hankel = mp.hankel1(order, k)
        slope = (mp.hankel1(order - 1, k) - mp.hankel1(order + 1, k)) / 2
        values.append(mp.mpf("0.5") + k * slope / hankel)
    return rotation * mp.diag(values) * rotation.T
