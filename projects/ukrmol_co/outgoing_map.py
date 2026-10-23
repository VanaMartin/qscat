"""Complex outgoing Riccati/Jost maps; experimental, not a CO pole evaluator."""

from __future__ import annotations

import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import lstsq, solve


def sheet_momenta(energy, thresholds, reference):
    """Continue outgoing open / decaying closed roots from a fixed real interval."""
    thresholds = np.asarray(thresholds, float)
    if not np.isfinite(energy) or not np.isfinite(reference) or not np.isfinite(thresholds).all():
        raise ValueError("Finite energy, real reference and thresholds required")
    if np.any(thresholds == reference) or np.any(thresholds == energy):
        raise ValueError("Threshold branch point is outside this map")
    if np.any((energy.real - thresholds) * (reference - thresholds) <= 0):
        raise ValueError("Energy leaves the declared threshold interval")
    opened = thresholds < reference
    k = np.empty(len(thresholds), complex)
    k[opened] = np.sqrt(complex(energy) - thresholds[opened])
    k[~opened] = 1j * np.sqrt(thresholds[~opened] - complex(energy))
    return k, opened


def _tail_series(thresholds, c2, c3, s, order):
    """Inverse-power recurrence with explicitly continued column exponents."""
    gaps = thresholds[:, None] - thresholds
    grouped = abs(gaps) <= 1e-10
    coefficients = [np.eye(len(s), dtype=complex)]
    zero = np.zeros_like(coefficients[0])
    for p in range(1, order + 1):
        previous = coefficients[p - 1]
        before = coefficients[p - 2] if p >= 2 else zero
        earlier = coefficients[p - 3] if p >= 3 else zero
        rhs = -2 * s * (p - 1) * previous + (p - 2) * (p - 1) * before - c2 @ before - c3 @ earlier
        current = np.zeros_like(zero)
        np.divide(rhs, gaps, out=current, where=~grouped)
        same = (p * (p - 1) * previous - c2 @ previous - c3 @ before) / (2 * s * p)
        current[grouped] = same[grouped]
        coefficients.append(current)
    return np.asarray(coefficients)


def outgoing_map(
    thresholds,
    angular,
    c2,
    c3,
    energy,
    radius,
    endpoint,
    angle,
    order=8,
    *,
    reference=None,
    rtol=1e-10,
    timeout=90.0,
    max_nfev=50000,
    amplitudes=False,
):
    """Return d/dr DtN Y at the real boundary and optional flux-normalized open A.

    Energies/C2/C3 use native Rydberg radial-equation units. Endpoint is the
    real contour parameter t, not the complex coordinate's modulus. This is
    binary64 DOP853 with finite inverse-power tail data, no pole search.
    """
    started = time.monotonic()
    thresholds, angular = np.asarray(thresholds, float), np.asarray(angular)
    c2, c3 = np.asarray(c2, float), np.asarray(c3, float)
    n = len(thresholds)
    if (
        thresholds.ndim != 1
        or n == 0
        or angular.shape != (n,)
        or any(x.shape != (n, n) or not np.isfinite(x).all() for x in [c2, c3])
        or not all(np.isfinite(x) for x in [radius, endpoint, angle, rtol, timeout])
        or radius <= 0
        or endpoint <= radius
        or not 0 < angle < 90
        or not 0 < rtol < 1
        or timeout <= 0
        or max_nfev < 1
        or not isinstance(order, int)
        or order < 0
    ):
        raise ValueError("Finite compatible coefficients and a positive rotated interval required")
    for matrix in [c2, c3]:
        if not np.allclose(matrix, matrix.T, atol=1e-12, rtol=0):
            raise ValueError("Real symmetric exterior coefficients required")
    reference = float(energy.real) if reference is None else reference
    k, opened = sheet_momenta(energy, thresholds, reference)
    rotation = np.exp(1j * np.deg2rad(angle))
    if np.min((rotation * k).imag) <= 0:
        raise ValueError("Contour does not expose every outgoing/decaying channel")
    if amplitudes and complex(energy).imag != 0:
        raise ValueError("Physical S amplitudes are implemented only on the real axis")
    s = 1j * k
    coefficients = _tail_series(thresholds, c2, c3, s, order)
    end = radius + (endpoint - radius) * rotation
    powers = end ** -np.arange(order + 1)
    series = np.einsum("p,pij->ij", powers, coefficients)
    derivative = -np.einsum("p,pij->ij", np.arange(order + 1) * powers / end, coefficients)
    y = solve(series.T, (derivative + series * s).T).T
    indices = np.flatnonzero(opened)
    if amplitudes:
        inverse_rows = solve(series.T, np.eye(n)[:, indices]).T
        weights = np.sqrt(k[opened]) * np.exp(-1j * (k[opened] * end - angular[opened] * np.pi / 2))
        a = weights[:, None] * inverse_rows
    else:
        a = np.empty((0, n), complex)
    initial = np.concatenate([y.ravel(), a.ravel()])
    diagonal = np.diag(thresholds - complex(energy))
    calls = 0

    def rhs(t, state):
        nonlocal calls
        calls += 1
        if calls > max_nfev or time.monotonic() - started > timeout:
            raise TimeoutError(f"Outgoing-map cap reached after {calls} RHS calls")
        r = radius + (t - radius) * rotation
        y = state[: n * n].reshape(n, n)
        dy = diagonal + c2 / r**2 + c3 / r**3 - y @ y
        a = state[n * n :].reshape(-1, n)
        return rotation * np.concatenate([dy.ravel(), (-a @ y).ravel()])

    result = solve_ivp(
        rhs,
        (endpoint, radius),
        initial,
        method="DOP853",
        rtol=rtol,
        atol=rtol * 0.01,
        t_eval=[radius],
    )
    if not result.success or not np.isfinite(result.y).all():
        raise ValueError(f"Outgoing-map integration failed: {result.message}")
    final = result.y[:, -1]
    y = final[: n * n].reshape(n, n)
    a = final[n * n :].reshape(-1, n)
    return (
        y,
        a,
        dict(
            rhs_evaluations=calls,
            wall_seconds=time.monotonic() - started,
            transpose_defect=float(np.max(abs(y - y.T))),
            open_indices=indices.tolist(),
            endpoint=endpoint,
            angle=angle,
            order=order,
            rtol=rtol,
        ),
    )


def physical_s(q, y, a):
    """Real-axis physical open S from the bilinear Jost Wronskians; no repair."""
    q = np.asarray(q, float)
    if q.shape != y.shape or not np.allclose(q, q.T, atol=1e-10, rtol=0):
        raise ValueError("A compatible real symmetric zero-Bloch Q is required")
    wall, _, rank, singular = lstsq(a.conj().T, a.T)
    if rank != len(a):
        raise ValueError("Rank-deficient open Jost amplitude map")
    scattering = wall + 2j * a @ solve(np.eye(len(q)) - q @ y, q) @ a.T
    return scattering, dict(
        flux_defect=float(np.max(abs((y - y.conj().T) / (2j) - a.conj().T @ a))),
        wronskian_residual=float(np.max(abs(a.conj().T @ wall - a.T))),
        amplitude_condition=float(singular[0] / singular[-1]),
        unitarity_defect=float(np.max(abs(scattering.conj().T @ scattering - np.eye(len(a))))),
        reciprocity_defect=float(np.max(abs(scattering - scattering.T))),
    )
