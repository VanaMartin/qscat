"""Time-dependent 2-D VE cross section via CN propagation + the Tannor-Weeks
energy transform (sub-project #7, Task 4 -- THE CRUX).

`.superpowers/sdd/n2-2d-exact-extraction.md` section 5.3
(`eMoScat TestFunction2d.cpp:298-307`):

    S_{v->v'}(E) = [2*pi*conj(eta_out_{v'}(E))*eta_in_v(E)]^{-1}
                   * sum_n w_n exp(i*E_tot*t_n) c_{v'}(t_n) * dt
    sigma_{v->v'}(E) = pi*|S - delta_{v,v'}|^2 / (2*E)     [bohr^2]

with `E_tot = E + eps[v_init]`. `w_n` are composite Simpson (trapezoid
fallback for an even number of samples) quadrature weights, reusing
`projects.n2_td_cross_section.td_cross_section`'s `_quadrature_weights`
pattern. `eta_in`/`eta_out` are `correlation.py`'s deconvolution factors,
projecting the incident/outgoing Gaussian wavepackets onto the SAME
energy-normalized `riccati_bessel_en` that
`projects.n2_2d_cross_section.cross_section_2d.channel_vector` uses to build
its exact TI channel functions -- the reason this converges to that module's
`ve_cross_section_2d` as `dt -> 0` and `n_steps -> infinity` (equivalently,
as `T = n_steps*dt` grows long enough for `||Psi(T)||` to decay -- the
resonance depletes and the finite-time truncation of the energy transform
vanishes). See `test_td_cross_section.py` for the convergence study and the
measured (dt, n_steps, T, norm-decay, sigma_TD/sigma_TI) numbers.

Elastic (`v' == v_init`) subtracts the Kronecker delta before squaring, per
the S-matrix convention `sigma = pi*|S - delta|^2 / (2E)`, equivalently
`4*pi^3*|T|^2/(2E)` via `S = 1 - 2*pi*i*T` (see the module docstring's
extraction reference and `cross_section_2d.py`'s `4*pi^3*|T|^2/(2E)` form --
these are the SAME formula in different conventions, not two competing ones;
do NOT use the `4*pi^3` prefactor together with this `S` (that double-counts
the `2*pi` already divided out by `eta`), and do NOT drop the `2*pi` in the
`eta` denominator (that under-counts it).
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from qscat.dvr import TensorGrid

from projects.n2_2d_cross_section.hamiltonian2d import ELL
from projects.n2_2d_td_cross_section.correlation import (
    eta_incident,
    eta_outgoing,
    outgoing_channel,
)
from projects.n2_2d_td_cross_section.td_propagation import PropagationResult, propagate
from projects.n2_2d_td_cross_section.wavepacket import initial_state

__all__ = ["td_ve_cross_section_2d"]

# Wavepacket parameter dict keys `initial_state`/`outgoing_channel` accept
# (r0/p0/sigma for the incident packet; r0_out/p0_out/sigma_out for the
# outgoing test function).
_WpIn = dict[str, float]
_WpOut = dict[str, float]


def _quadrature_weights(n_t: int) -> npt.NDArray[np.float64]:
    """Composite Simpson weights (unscaled by `dt`) for `n_t` samples.

    Identical pattern to `projects.n2_td_cross_section.td_cross_section`'s
    `_quadrature_weights`: requires an odd `n_t` (even `n_steps`) for the
    standard composite Simpson rule `dt/3*(f0+4f1+2f2+...+4f_{N-1}+fN)`;
    falls back to composite trapezoidal `dt/2*(f0+2f1+...+2f_{N-1}+fN)` for
    an even `n_t` (odd `n_steps`).
    """
    if n_t < 2:
        raise ValueError("need at least 2 time samples for a quadrature rule")
    if n_t % 2 == 1:
        w = np.ones(n_t, dtype=np.float64)
        w[1:-1:2] = 4.0
        w[2:-1:2] = 2.0
        w /= 3.0
    else:
        w = np.full(n_t, 2.0, dtype=np.float64)
        w[0] = 1.0
        w[-1] = 1.0
        w /= 2.0
    return w


def _propagate(
    tgrid: TensorGrid,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    vprimes: list[int],
    *,
    dt: float,
    n_steps: int,
    wp_in: _WpIn,
    wp_out: _WpOut,
) -> PropagationResult:
    psi0 = initial_state(tgrid, chi[v_init], **wp_in)
    out_channels = [outgoing_channel(tgrid, chi[vp], **wp_out) for vp in vprimes]
    return propagate(tgrid, psi0, out_channels, dt=dt, n_steps=n_steps)


def _sigma_from_correlations(
    tgrid: TensorGrid,
    result: PropagationResult,
    eps: npt.NDArray[np.float64],
    v_init: int,
    vprimes: list[int],
    E: float,
    dt: float,
    wp_in: _WpIn,
    wp_out: _WpOut,
) -> npt.NDArray[np.float64]:
    sigma = np.zeros(len(vprimes), dtype=np.float64)
    if E <= 0.0:
        return sigma

    weights = _quadrature_weights(result.t.size)
    e_tot = E + eps[v_init]
    k = float(np.sqrt(2.0 * E))
    eta_in = eta_incident(tgrid.grids[0], k, ELL, **wp_in)
    phase = np.exp(1j * e_tot * result.t)

    for j, vp in enumerate(vprimes):
        excess = e_tot - eps[vp]
        if excess <= 0.0:
            continue  # closed channel
        kp = float(np.sqrt(2.0 * excess))
        eta_out = eta_outgoing(tgrid.grids[0], kp, ELL, **wp_out)
        s_raw = np.sum(weights * phase * result.c[:, j]) * dt
        s = s_raw / (2.0 * np.pi * np.conj(eta_out) * eta_in)
        delta = 1.0 if vp == v_init else 0.0
        sigma[j] = np.pi * abs(s - delta) ** 2 / (2.0 * E)
    return sigma


def td_ve_cross_section_2d(
    tgrid: TensorGrid,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    vprimes: list[int],
    E: float | npt.ArrayLike,
    *,
    dt: float,
    n_steps: int,
    wp_in: _WpIn,
    wp_out: _WpOut,
) -> npt.NDArray[np.float64]:
    """sigma_{v_init->v'}(E) (bohr^2), 2-D CN propagation + Tannor-Weeks transform.

    `E` (collision energy, Hartree) may be scalar or array-like; scalar `E`
    returns shape `(len(vprimes),)`, array `E` returns `(len(E), len(vprimes))`
    -- matching `projects.n2_2d_cross_section.cross_section_2d.ve_cross_section_2d`,
    the exact differential oracle this converges to as `dt -> 0` and
    `n_steps -> infinity` (see module docstring).

    `wp_in = {"r0": ..., "p0": ..., "sigma": ...}` are the SAME incident
    Gaussian parameters used to build `Psi(0)` (via `initial_state`) and
    `eta_incident`; `wp_out = {"r0_out": ..., "p0_out": ..., "sigma_out": ...}`
    are the outgoing test function's parameters, used for both
    `outgoing_channel` (the `Phi_{v'}` propagated against) and `eta_outgoing`.
    The propagation (the expensive part) happens ONCE regardless of how many
    energies `E` are requested, since `c_{v'}(t)` does not depend on `E`.
    """
    result = _propagate(
        tgrid, eps, chi, v_init, vprimes, dt=dt, n_steps=n_steps, wp_in=wp_in, wp_out=wp_out
    )

    e_arr = np.atleast_1d(np.asarray(E, dtype=np.float64))
    out = np.stack(
        [
            _sigma_from_correlations(
                tgrid, result, eps, v_init, vprimes, float(e), dt, wp_in, wp_out
            )
            for e in e_arr
        ]
    )
    scalar = np.isscalar(E) or (isinstance(E, np.ndarray) and np.ndim(E) == 0)
    if scalar:
        return np.asarray(out[0], dtype=np.float64)
    return out
