"""Time-dependent VE cross section via Crank-Nicolson propagation + energy
transform (sub-project #4, Task 2 -- THE CRUX).

`docs/superpowers/specs/2026-07-22-n2-td-cross-section-design.md` ("Method")
and `.superpowers/sdd/n2-lcp-model-extraction.md`:

- Doorway `d_v(R) = sqrt(Gamma(R)/(2*pi)) * chi_v(R)` -- identical to
  `projects.n2_ti_cross_section.cross_section`'s doorway.
- Initial wavepacket `psi(0) = d_{v_init}`, propagated under the
  time-independent, non-Hermitian
  `H_res = T_nuc(mu) + diag(V_d(R) - i*Gamma(R)/2)` via the Crank-Nicolson
  stepper (`projects.n2_td_cross_section.propagator.make_cn_stepper`). The
  `-i*Gamma/2` term makes ``||psi||`` decay -- the resonance depletes.
- Correlation function `c_{v'}(t_n) = <d_{v'}|psi(t_n)>` -- the c-product
  (plain coefficient dot, NO conjugate), matching the TI oracle's S-matrix
  convention: the DVR basis is already 1/sqrt(weight)-normalized, and
  `psi(t_n)` is a genuinely complex ECS-driven state, not an eigenvector
  needing a Hermitian norm.
- Energy transform
  `S_{v'}(E) = (1/i) * sum_n w_n * exp(i*(E + eps[v_init])*t_n) * c_{v'}(t_n) * dt`
  (Simpson weights `w_n`, trapezoidal fallback when `n_steps` is odd), then
  `sigma_{v_init->v'}(E) = 4*pi**3*|S|**2/(2*E)`, zero if `E<=0` or the final
  channel is energetically closed (`E_tot - eps[v'] <= 0`).

Because `S_TD(E) = (1/i) * integral_0^inf exp(i*E_tot*t) * <d_v'|exp(-i*H_res*t)|d_v> dt
= <d_v'|(E_tot - H_res)^-1|d_v> = S_TI(E)` in the long-time limit, TD sigma
converges to the TI oracle's sigma
(`projects.n2_ti_cross_section.cross_section.ve_cross_section`) as
`dt -> 0` and `n_steps -> infinity` (equivalently, as the propagation time
`T = n_steps*dt` grows long enough for the correlation function to decay --
see `test_td_cross_section.py`'s V1/V2).
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from qscat.dvr import FemDvrEcsGrid, kinetic

from projects.n2_td_cross_section.propagator import make_cn_stepper

__all__ = ["td_ve_cross_section", "td_norm_ratio"]


def _quadrature_weights(n_t: int) -> npt.NDArray[np.float64]:
    """Composite Simpson weights (unscaled by `dt`) for `n_t` samples.

    Requires an odd number of samples (even number of intervals,
    i.e. `n_steps = n_t - 1` even) for the standard composite rule
    `dt/3 * (f0 + 4f1 + 2f2 + ... + 4f_{N-1} + fN)`; falls back to the
    composite trapezoidal rule `dt/2 * (f0 + 2f1 + ... + 2f_{N-1} + fN)`
    when `n_t` is even (`n_steps` odd), per the task brief's documented
    fallback.
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


def _propagate_and_correlate(
    grid: FemDvrEcsGrid,
    mu: float,
    Vd: npt.NDArray[np.complex128],
    Gamma: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    vprimes: list[int],
    *,
    dt: float,
    n_steps: int,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.complex128], float]:
    """Propagate `psi(0) = d_{v_init}` under CN and record correlations.

    Returns `(t, c, norm_ratio)`: `t` is the `(n_steps+1,)` sample-time
    array `t_n = n*dt`; `c` is the `(len(vprimes), n_steps+1)` array of
    `c_{v'}(t_n) = <d_{v'}|psi(t_n)>` (c-product, no conjugate), recorded
    BEFORE each step so `t_0 = 0` is included; `norm_ratio` is
    `||psi(T)|| / ||psi(0)||` with `T = n_steps*dt`, exposed for the
    convergence/depletion check (V2).
    """
    doorway = np.sqrt(Gamma / (2.0 * np.pi))[None, :] * chi  # (n_vib, n)

    H_res = kinetic(grid, mu) + np.diag(Vd - 0.5j * Gamma)
    step = make_cn_stepper(H_res, dt)

    psi = doorway[v_init].astype(np.complex128).copy()
    norm0 = np.linalg.norm(psi)

    n_t = n_steps + 1
    t = np.arange(n_t, dtype=np.float64) * dt
    c = np.empty((len(vprimes), n_t), dtype=np.complex128)

    for n in range(n_t):
        for k, vp in enumerate(vprimes):
            c[k, n] = np.dot(doorway[vp], psi)  # c-product: no conjugate
        if n < n_steps:
            psi = step(psi)

    norm_ratio = float(np.linalg.norm(psi) / norm0)
    return t, c, norm_ratio


def _sigma_from_correlations(
    t: npt.NDArray[np.float64],
    c: npt.NDArray[np.complex128],
    eps: npt.NDArray[np.float64],
    v_init: int,
    vprimes: list[int],
    E: float,
    dt: float,
) -> npt.NDArray[np.float64]:
    sigma = np.zeros(len(vprimes), dtype=np.float64)
    if E <= 0.0:
        return sigma

    weights = _quadrature_weights(t.size)
    E_tot = E + eps[v_init]
    phase = np.exp(1j * E_tot * t)

    for k, vp in enumerate(vprimes):
        if E_tot - eps[vp] <= 0.0:
            continue  # closed channel
        S = (1.0 / 1j) * np.sum(weights * phase * c[k]) * dt
        sigma[k] = 4.0 * np.pi**3 * np.abs(S) ** 2 / (2.0 * E)
    return sigma


def td_ve_cross_section(
    grid: FemDvrEcsGrid,
    mu: float,
    Vd: npt.NDArray[np.complex128],
    Gamma: npt.NDArray[np.float64],
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    vprimes: list[int],
    E: float | npt.ArrayLike,
    *,
    dt: float,
    n_steps: int,
    return_norm_ratio: bool = False,
) -> npt.NDArray[np.float64] | tuple[npt.NDArray[np.float64], float]:
    """sigma_{v_init->v'}(E) (bohr^2) via CN propagation + energy transform.

    `E` (collision energy, Hartree) may be a scalar or an array; scalar `E`
    returns shape `(len(vprimes),)`, array `E` returns shape
    `(len(E), len(vprimes))` -- matching
    `projects.n2_ti_cross_section.cross_section.ve_cross_section`, the
    exact differential oracle this converges to as `dt -> 0` and
    `n_steps -> infinity` (see module docstring).

    `dt` and `n_steps` control the propagation: `T = n_steps*dt` must be
    long enough for `||psi(t)||` to decay (the resonance depletes) for the
    energy transform to be well-resolved -- see `test_td_cross_section.py`.

    If `return_norm_ratio` is True, returns `(sigma, norm_ratio)` where
    `norm_ratio = ||psi(T)|| / ||psi(0)||` (for the V2 depletion check);
    otherwise returns `sigma` alone.
    """
    t, c, norm_ratio = _propagate_and_correlate(
        grid, mu, Vd, Gamma, chi, v_init, vprimes, dt=dt, n_steps=n_steps
    )

    E_arr = np.atleast_1d(np.asarray(E, dtype=np.float64))
    out = np.stack(
        [_sigma_from_correlations(t, c, eps, v_init, vprimes, float(e), dt) for e in E_arr]
    )
    if np.isscalar(E) or (isinstance(E, np.ndarray) and E.ndim == 0):
        sigma = np.asarray(out[0], dtype=np.float64)
    else:
        sigma = out

    if return_norm_ratio:
        return sigma, norm_ratio
    return sigma


def td_norm_ratio(
    grid: FemDvrEcsGrid,
    mu: float,
    Vd: npt.NDArray[np.complex128],
    Gamma: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    *,
    dt: float,
    n_steps: int,
) -> float:
    """`||psi(T)|| / ||psi(0)||` with `T = n_steps*dt`, `psi(0) = d_{v_init}`.

    A thin convenience wrapper over `_propagate_and_correlate` for callers
    (e.g. the V2 convergence test) that only need the depletion ratio, not
    a full cross section.
    """
    _, _, norm_ratio = _propagate_and_correlate(
        grid, mu, Vd, Gamma, chi, v_init, [v_init], dt=dt, n_steps=n_steps
    )
    return norm_ratio
