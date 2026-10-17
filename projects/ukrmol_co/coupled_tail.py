"""Two-channel inverse-square tails and threshold shells with analytic oracles."""

from __future__ import annotations

import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import h1vp, h2vp, hankel1, hankel2


def propagate_pair(u, v, start, end, operator, rtol):
    """Propagate two regular solution columns in their original channel basis."""
    if end == start:
        return u.copy(), v.copy()

    def rhs(x, state):
        pair = state.reshape(2, 2, 2)
        return np.stack([pair[1], operator(x) @ pair[0]]).ravel()

    result = solve_ivp(
        rhs,
        (start, end),
        np.stack([u, v]).astype(complex).ravel(),
        method="DOP853",
        rtol=rtol,
        atol=rtol * 0.01,
    )
    if not result.success or not np.isfinite(result.y).all():
        raise ValueError(f"Matrix propagation failed: {result.message}")
    return tuple(result.y[:, -1].reshape(2, 2, 2))


def match_pair(u, v, plus, plus_derivative, minus, minus_derivative):
    """Match regular columns to outgoing/closed and incident-open wave columns."""
    matrix = np.block([[u, -plus], [v, -plus_derivative]])
    solution = np.linalg.solve(matrix, np.vstack([minus, minus_derivative]))
    return -solution[2:], float(np.linalg.cond(matrix))


def tail_coefficient(dipole):
    return np.array([[2.0, -2 * dipole / np.sqrt(5)], [-2 * dipole / np.sqrt(5), 6.0]])


def tail_waves(q, x, dipole):
    """Flux-normalized angular-basis waves for an uncut inverse-square tail."""
    eigenvalues, rotation = np.linalg.eigh(tail_coefficient(dipole))
    orders = np.sqrt(eigenvalues + 0.25)
    waves = []
    for sign, hankel, derivative in [(1, hankel1, h1vp), (-1, hankel2, h2vp)]:
        scale = np.sqrt(np.pi * x / 2) * np.exp(sign * 1j * (np.pi * orders / 2 + np.pi / 4))
        values = scale * hankel(orders, q * x)
        slopes = values / (2 * x) + scale * q * derivative(orders, q * x)
        phase = np.exp(-sign * 1j * np.array([1, 2]) * np.pi / 2)
        waves.append(
            ((rotation * values) @ rotation.T * phase, (rotation * slopes) @ rotation.T * phase)
        )
    return waves


def tail_scattering(q, dipole, radius, rtol, matching_dipole=None):
    coefficient = tail_coefficient(dipole)
    u, v = propagate_pair(
        np.zeros((2, 2)),
        np.eye(2),
        1,
        radius,
        lambda x: coefficient / x**2 - q * q * np.eye(2),
        rtol,
    )
    (hp, dp), (hm, dm) = tail_waves(
        q, radius, dipole if matching_dipole is None else matching_dipole
    )
    return match_pair(u, v, hp, dp, hm, dm)


def exact_tail(q, dipole):
    """High-precision direct hard-wall S; no propagation or SciPy Hankel calls."""
    q, d = mp.mpc(q), mp.mpf(dipole)
    eigenvalues, rotation = mp.eigsy(
        mp.matrix([[2, -2 * d / mp.sqrt(5)], [-2 * d / mp.sqrt(5), 6]])
    )
    waves = []
    for sign, hankel in [(1, mp.hankel1), (-1, mp.hankel2)]:
        values = []
        for value in eigenvalues:
            order = mp.sqrt(value + mp.mpf("0.25"))
            values.append(
                mp.sqrt(mp.pi / 2)
                * mp.exp(sign * 1j * (mp.pi * order / 2 + mp.pi / 4))
                * hankel(order, q)
            )
        phase = mp.diag([mp.exp(-sign * 1j * ell * mp.pi / 2) for ell in [1, 2]])
        waves.append(rotation * mp.diag(values) * rotation.T * phase)
    return waves[0] ** -1 * waves[1]


def momenta(epsilon, separation=16.0, closed_sign=1):
    """Physical real-axis branches; below upper threshold keep its decaying sheet."""
    first = np.sqrt(complex(epsilon))
    second = (
        1j * np.sqrt(complex(separation - epsilon))
        if np.real(epsilon) < separation
        else np.sqrt(complex(epsilon - separation))
    )
    return np.array([first, closed_sign * second])


def threshold_pair(epsilon, coupling, radius, rtol):
    operator = np.diag(np.array([0, 16]) - epsilon)
    u, v = propagate_pair(np.zeros((2, 2)), np.eye(2), 0, 1, lambda x: operator, rtol)
    v += np.array([[10, coupling], [coupling, 12]]) @ u
    return propagate_pair(u, v, 1, radius, lambda x: operator, rtol)


def threshold_scattering(epsilon, coupling, radius, rtol):
    k = momenta(epsilon)
    opened = np.flatnonzero(epsilon > np.array([0, 16]))
    plus = np.diag(
        [
            np.exp(1j * ki * radius) / np.sqrt(ki.real)
            if i in opened
            else np.exp(1j * ki * (radius - 1))
            for i, ki in enumerate(k)
        ]
    )
    minus = np.zeros((2, len(opened)), complex)
    for j, i in enumerate(opened):
        minus[i, j] = np.exp(-1j * k[i] * radius) / np.sqrt(k[i].real)
    u, v = threshold_pair(epsilon, coupling, radius, rtol)
    outgoing, condition = match_pair(
        u, v, plus, 1j * k[:, None] * plus, minus, -1j * k[:, None] * minus
    )
    return outgoing[opened], condition


def threshold_jost(q, coupling, radius, rtol):
    k = momenta(q * q)
    u, v = threshold_pair(q * q, coupling, radius, rtol)
    return complex(np.linalg.det(np.exp(1j * k * radius)[:, None] * (v - 1j * k[:, None] * u)))


def mp_momenta(epsilon, closed_sign=1):
    epsilon = mp.mpc(epsilon)
    return [
        mp.sqrt(epsilon),
        closed_sign * (1j * mp.sqrt(16 - epsilon) if epsilon.real < 16 else mp.sqrt(epsilon - 16)),
    ]


def exact_threshold(epsilon, coupling, closed_sign=1, flux=True):
    """Independent cotangent/Schur oracle, including deliberately wrong controls."""
    epsilon, c = mp.mpf(epsilon), mp.mpf(coupling)
    k = mp_momenta(epsilon, closed_sign)
    logder = mp.diag([ki / mp.tan(ki) for ki in k]) + mp.matrix([[10, c], [c, 12]])
    opened = [0] if epsilon < 16 else [0, 1]
    if len(opened) == 1:
        logder = mp.matrix(
            [[logder[0, 0] - logder[0, 1] * logder[1, 0] / (logder[1, 1] - 1j * k[1])]]
        )
    ko = mp.diag([k[i] for i in opened])
    plus = mp.diag([mp.exp(1j * k[i]) / (mp.sqrt(k[i]) if flux else 1) for i in opened])
    minus = mp.diag([mp.exp(-1j * k[i]) / (mp.sqrt(k[i]) if flux else 1) for i in opened])
    return plus**-1 * (logder - 1j * ko) ** -1 * (logder + 1j * ko) * minus


def exact_jost(q, coupling, incoming=False):
    k = mp_momenta(q * q)
    u = mp.diag([mp.sin(ki) / ki for ki in k])
    v = mp.diag([mp.cos(ki) for ki in k]) + mp.matrix([[10, coupling], [coupling, 12]]) * u
    if incoming:
        k[0] = -k[0]
    return mp.det(mp.diag([mp.exp(1j * ki) for ki in k]) * (v - 1j * mp.diag(k) * u))


def exact_threshold_pole(coupling, seeds):
    return mp.findroot(
        lambda q: exact_jost(q, coupling),
        tuple(mp.mpc(*s) for s in seeds),
        tol=mp.mpf("1e-50"),
        maxsteps=60,
    )
