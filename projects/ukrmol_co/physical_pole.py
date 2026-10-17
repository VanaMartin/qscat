"""Physical radial-shell boundary matching with an independent Lambert-W pole oracle."""

from __future__ import annotations

import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root


def shell_boundary(q: complex, strength: float, radius: float, rtol: float) -> np.ndarray:
    """Return (u, du/dx) at x=radius for a shell at x=1 and u_x(0)=1."""
    if radius < 1 or strength < 0 or rtol <= 0:
        raise ValueError("Require radius>=1, nonnegative shell strength and positive tolerance")

    def rhs(x, state):
        return np.array([state[1], -q * q * state[0]])

    def propagate(start, end, state):
        result = solve_ivp(rhs, (start, end), state, method="DOP853", rtol=rtol, atol=rtol * 0.01)
        if not result.success or not np.isfinite(result.y).all():
            raise ValueError(f"Radial integration failed: {result.message}")
        return result.y[:, -1]

    state = propagate(0, 1, np.array([0j, 1 + 0j]))
    state[1] += strength * state[0]
    return state if radius == 1 else propagate(1, radius, state)


def boundary_s(q: complex, state: np.ndarray, radius: float, bloch: float = 0) -> complex:
    """Match the homogeneous (u, B*u'-b*u) pair to free incoming/outgoing waves."""
    value, derivative = state
    d = radius * derivative - bloch * value
    local_q = q * radius
    return complex(
        np.exp(-2j * local_q)
        * (d + (bloch + 1j * local_q) * value)
        / (d + (bloch - 1j * local_q) * value)
    )


def outgoing_residual(q: complex, state: np.ndarray, radius: float) -> complex:
    """Normalized outgoing Jost residual: exp(i*q*x)*(u_x-i*q*u)."""
    return complex(np.exp(1j * q * radius) * (state[1] - 1j * q * state[0]))


def shoot_pole(strength: float, index: int, radius: float, rtol: float, contract: dict) -> dict:
    """One q-plane shooting solve from a declared seed, with no oracle input."""
    seed = [
        index * np.pi + contract["pole_seed_real_offset_from_n_pi"],
        contract["pole_seed_imaginary_part"],
    ]

    def residual(pair):
        q = complex(*pair)
        value = outgoing_residual(q, shell_boundary(q, strength, radius, rtol), radius)
        return [value.real, value.imag]

    solution = root(
        residual, seed, tol=contract["root_tolerance"], options={"maxfev": contract["root_maxfev"]}
    )
    q = complex(*solution.x)
    state = shell_boundary(q, strength, radius, rtol)
    incoming = np.exp(-1j * q * radius) * (state[1] + 1j * q * state[0])
    return {
        "q": [q.real, q.imag],
        "seed_q": seed,
        "solver_success": bool(solution.success),
        "solver_status": int(solution.status),
        "solver_message": str(solution.message),
        "evaluations": int(solution.nfev),
        "outgoing_residual_absolute": abs(outgoing_residual(q, state, radius)),
        "incoming_residual_absolute": float(abs(incoming)),
    }


def exact_s(q, strength):
    """Arbitrary-precision Jost-ratio oracle, without numerical radial propagation."""
    q, strength = mp.mpc(q), mp.mpf(strength)
    if not q:
        return mp.mpc(1)
    return (1 + strength * mp.exp(-1j * q) * mp.sin(q) / q) / (
        1 + strength * mp.exp(1j * q) * mp.sin(q) / q
    )


def exact_pole(strength: float, index: int):
    """Arbitrary-precision positive-real-q outgoing pole on Lambert branch -index."""
    if strength <= 0 or index < 1:
        raise ValueError("A nonzero repulsive shell and positive pole index are required")
    value = mp.mpf(strength)
    return 0.5j * (mp.lambertw(value * mp.exp(value), -index) - value)
