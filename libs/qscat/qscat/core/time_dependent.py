"""Time-dependent VE cross section: Pade propagation + the Tannor-Weeks energy
transform, for any `model: qscat.model.ResonanceModel`.

Merges `projects/n2_2d_td_cross_section/td_propagation.py` (the propagation
engine: propagate `Psi(0)` under `H` and sample) and `td_cross_section.py`
(the Tannor-Weeks transform) into one model-taking module (sub-project #A,
Task 5). The only change from those two N2 project modules is where the
Hamiltonian and partial wave `l` come from: `model.hamiltonian(tgrid)` /
`model.interaction_diag(tgrid)` / `model.ell` instead of a hardcoded N2
`build_h2d`/`interaction_diag`/`ELL` import -- everything else (the Pade
stepping, sampling cadences, the energy transform, the elastic free-reference
fix) is unchanged.

## Propagation (two sampling cadences)

  * `c_{v'}(t_n)` -- recorded at EVERY step; the raw material of the
    Tannor-Weeks energy transform.
  * density/norm snapshots -- recorded on a COARSE schedule (`sample_period`
    steps, or explicit `snapshot_times`), so the wavefunction is observed at
    static points without storing every step.

`H` is time-independent, so the sparse Pade factorization is built once and
reused; under the ECS contour `||Psi||` decays as outgoing flux is absorbed
(the resonance depletes). Norm here is the Hermitian L2 norm
`np.linalg.norm(psi)` -- the physical remaining-probability diagnostic. The
c-product is a different object, reserved for the correlations `c_{v'}(t)`.

`order` is the diagonal-Pade order of the evolution operator
(`qscat.evolution.make_pade_stepper`): `O(dt^(2*order+1))` per step. The
default 3 matches eMoScat (order-1 Crank-Nicolson under-converges badly over
a multi-thousand-step run); see `docs/physics/n2-2d-td-cross-section.md`.

## Energy transform (Tannor-Weeks)

`.superpowers/sdd/n2-2d-exact-extraction.md` section 5.3
(`eMoScat TestFunction2d.cpp:298-307`):

    S_{v->v'}(E) = [2*pi*conj(eta_out_{v'}(E))*eta_in_v(E)]^{-1}
                   * sum_n w_n exp(i*E_tot*t_n) c_{v'}(t_n) * dt
    sigma_{v->v'}(E) = pi*|S - delta_{v,v'}|^2 / (2*E)     [bohr^2]

with `E_tot = E + eps[v_init]`. `w_n` are composite Simpson (trapezoid
fallback for an even number of samples) quadrature weights. `eta_in`/`eta_out`
are `correlation.py`'s deconvolution factors, projecting the incident/
outgoing Gaussian wavepackets onto the SAME energy-normalized
`riccati_bessel_en` that `qscat.core.channels.channel_vector` uses to build
its exact TI channel functions -- the reason this converges to
`qscat.core.driven.ve_cross_section` as `dt -> 0` and `n_steps -> infinity`.

Elastic (`v' == v_init`) subtracts the unscattered reference before squaring:
`sigma = pi*|S - S_ref|^2 / (2E)`. In the standard convention `S_ref` is the
Kronecker delta (1), but that presumes the transform normalizes the free/
unscattered S-matrix to exactly 1. THIS transform does not: the outgoing
normalization factor C(E) multiplies every channel's S, so a free-particle
(`V_int=0`) propagation gives `S_free(E) = C(E) ~ 2*pi^2`, not 1. So the
elastic channel subtracts the S-matrix of a `V_int=0` reference propagation
(`_propagate(..., free=True)`), supplied via `free_result`; the literal-1
fallback (`free_result=None`) leaves a large spurious elastic background. See
`_sigma_one_energy` and the `td-elastic-wavepacket-normalization` note.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp

from qscat.dvr import TensorGrid
from qscat.evolution import make_pade_stepper
from qscat.linalg import c_product

from .correlation import eta_incident, eta_outgoing, outgoing_channel
from .wavepacket import initial_state

if TYPE_CHECKING:
    from qscat.model import ResonanceModel

__all__ = [
    "Snapshot",
    "PropagationResult",
    "propagate",
    "sigma_from_correlations",
    "td_ve_cross_section",
]

# Wavepacket parameter dict keys `initial_state`/`outgoing_channel` accept
# (r0/p0/sigma for the incident packet; r0_out/p0_out/sigma_out for the
# outgoing test function).
_WpIn = dict[str, float]
_WpOut = dict[str, float]


@dataclass(frozen=True)
class Snapshot:
    time: float
    rho_R: npt.NDArray[np.float64]  # nuclear density, sum_r |Psi|^2 (unscaled)
    rho_r: npt.NDArray[np.float64]  # electronic density, sum_R |Psi|^2 (unscaled)
    psi: npt.NDArray[np.complex128] | None  # full state, only if requested


@dataclass(frozen=True)
class PropagationResult:
    t: npt.NDArray[np.float64]  # (n_t,)  sample times n*dt
    c: npt.NDArray[np.complex128]  # (n_t, n_channels)  c_{v'}(t_n)
    norm: npt.NDArray[np.float64]  # (n_t,)  np.linalg.norm(psi) -- Hermitian L2
    snapshots: list[Snapshot]


def _densities(
    tgrid: TensorGrid, psi: npt.NDArray[np.complex128]
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    block = psi.reshape(tgrid.shape)
    dens = np.abs(block) ** 2
    r_real = tgrid.grids[0].real_points <= tgrid.grids[0].R0
    R_real = tgrid.grids[1].real_points <= tgrid.grids[1].R0
    rho_R = dens[r_real, :].sum(axis=0)
    rho_r = dens[:, R_real].sum(axis=1)
    return rho_R.astype(np.float64), rho_r.astype(np.float64)


def propagate(
    tgrid: TensorGrid,
    psi0: npt.NDArray[np.complex128],
    out_channels: list[npt.NDArray[np.complex128]],
    *,
    dt: float,
    n_steps: int,
    sample_period: int = 0,
    snapshot_times: list[float] | None = None,
    keep_psi_at: list[float] | None = None,
    hamiltonian: sp.spmatrix,
    order: int = 3,
) -> PropagationResult:
    """Propagate and sample. See module docstring for the two cadences.

    `hamiltonian` is the propagation Hamiltonian (REQUIRED here -- the
    model-agnostic engine has no default to fall back to; the caller supplies
    `model.hamiltonian(tgrid)`, or `model.hamiltonian(tgrid) -
    diag(model.interaction_diag(tgrid))` for the elastic free-reference path,
    see `_propagate`'s `free` argument).

    `order` is the diagonal-Pade order of the evolution operator
    (`qscat.evolution.make_pade_stepper`): `O(dt^(2*order+1))` per step. The
    default 3 matches eMoScat (order-1 Crank-Nicolson under-converges badly
    over a multi-thousand-step run -- ~100% accumulated error at dt=0.5/1.0 --
    and is the reason an order-1 TD cross section only reached ~10-15% of the
    TI oracle; order 3 brings it to convergence). See
    `docs/physics/n2-2d-td-cross-section.md`.
    """
    step = make_pade_stepper(hamiltonian, dt, order=order)

    n_t = n_steps + 1
    t = np.arange(n_t, dtype=np.float64) * dt
    n_ch = len(out_channels)
    c = np.empty((n_t, n_ch), dtype=np.complex128)
    norm = np.empty(n_t, dtype=np.float64)

    if snapshot_times is not None:
        snap_set = {round(x / dt) for x in snapshot_times}
    elif sample_period > 0:
        snap_set = set(range(0, n_t, sample_period)) | {n_t - 1}
    else:
        snap_set = {0, n_t - 1}
    keep_set = {round(x / dt) for x in (keep_psi_at or [])}
    snap_set |= keep_set  # a requested keep_psi_at time always gets a snapshot

    psi = psi0.astype(np.complex128).copy()
    snapshots: list[Snapshot] = []
    for n in range(n_t):
        for k in range(n_ch):
            c[n, k] = c_product(out_channels[k], psi)  # correlation: c-product
        norm[n] = float(np.linalg.norm(psi))  # Hermitian L2: physical, monotone
        if n in snap_set:
            rho_R, rho_r = _densities(tgrid, psi)
            snapshots.append(
                Snapshot(float(t[n]), rho_R, rho_r, psi.copy() if n in keep_set else None)
            )
        if n < n_steps:
            psi = step(psi)

    return PropagationResult(t=t, c=c, norm=norm, snapshots=snapshots)


def _quadrature_weights(n_t: int) -> npt.NDArray[np.float64]:
    """Composite Simpson weights (unscaled by `dt`) for `n_t` samples.

    Requires an odd `n_t` (even `n_steps`) for the standard composite Simpson
    rule `dt/3*(f0+4f1+2f2+...+4f_{N-1}+fN)`; falls back to composite
    trapezoidal `dt/2*(f0+2f1+...+2f_{N-1}+fN)` for an even `n_t` (odd
    `n_steps`).
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
    model: ResonanceModel,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    vprimes: list[int],
    *,
    dt: float,
    n_steps: int,
    wp_in: _WpIn,
    wp_out: _WpOut,
    free: bool = False,
    order: int = 3,
) -> PropagationResult:
    """Propagate the incident packet and record `c_{v'}(t)` for each `v'`.

    `free=True` propagates under `model.hamiltonian(tgrid)` with the
    interaction `V_int` removed (`model.hamiltonian(tgrid) -
    diag(model.interaction_diag(tgrid))`) -- the unscattered reference whose
    S-matrix `S_free(E)` the elastic channel subtracts instead of a literal 1
    (see `_sigma_one_energy`). Everything else -- the incident packet, the
    outgoing test functions, the grid -- is identical to the full run, so the
    spurious direct/unscattered content cancels in `S_full - S_free`.

    `order` is the diagonal-Pade evolution-operator order (default 3; see
    `propagate`). Order 1 is Crank-Nicolson.
    """
    psi0 = initial_state(tgrid, chi[v_init], **wp_in)
    out_channels = [outgoing_channel(tgrid, chi[vp], **wp_out) for vp in vprimes]
    if free:
        hamiltonian = (
            model.hamiltonian(tgrid) - sp.diags(model.interaction_diag(tgrid))
        ).tocsr()
    else:
        hamiltonian = model.hamiltonian(tgrid)
    return propagate(
        tgrid, psi0, out_channels, dt=dt, n_steps=n_steps, hamiltonian=hamiltonian, order=order
    )


def _s_vector_one_energy(
    tgrid: TensorGrid,
    model: ResonanceModel,
    result: PropagationResult,
    eps: npt.NDArray[np.float64],
    v_init: int,
    vprimes: list[int],
    E: float,
    dt: float,
    wp_in: _WpIn,
    wp_out: _WpOut,
) -> npt.NDArray[np.complex128]:
    """The complex S-matrix `S_{v_init->v'}(E)` for each `v'`, shape `(len(vprimes),)`.

    `0` for closed channels (`E_tot - eps[v'] <= 0`) and for `E <= 0`. This is
    the raw Tannor-Weeks transform (module docstring) BEFORE the `|S - ref|^2`
    step, factored out so the full run and the elastic free reference share
    one code path.
    """
    S = np.zeros(len(vprimes), dtype=np.complex128)
    if E <= 0.0:
        return S
    weights = _quadrature_weights(result.t.size)
    e_tot = E + eps[v_init]
    k = float(np.sqrt(2.0 * E))
    eta_in = eta_incident(tgrid.grids[0], k, model.ell, **wp_in)
    phase = np.exp(1j * e_tot * result.t)
    for j, vp in enumerate(vprimes):
        excess = e_tot - eps[vp]
        if excess <= 0.0:
            continue  # closed channel
        kp = float(np.sqrt(2.0 * excess))
        eta_out = eta_outgoing(tgrid.grids[0], kp, model.ell, **wp_out)
        s_raw = np.sum(weights * phase * result.c[:, j]) * dt
        S[j] = s_raw / (2.0 * np.pi * np.conj(eta_out) * eta_in)
    return S


def _sigma_one_energy(
    tgrid: TensorGrid,
    model: ResonanceModel,
    result: PropagationResult,
    eps: npt.NDArray[np.float64],
    v_init: int,
    vprimes: list[int],
    E: float,
    dt: float,
    wp_in: _WpIn,
    wp_out: _WpOut,
    free_result: PropagationResult | None = None,
) -> npt.NDArray[np.float64]:
    """`sigma_{v_init->v'}(E)` (bohr^2) at a single scalar `E`, shape `(len(vprimes),)`.

    The per-energy transform of an already-computed `PropagationResult` --
    the single-energy kernel both `sigma_from_correlations` (which adds the
    scalar-or-array `E` convention) and `td_ve_cross_section` (which also
    runs the propagation) build on.

    ELASTIC reference: `sigma = pi*|S - ref|^2/(2E)` with `ref` the free-particle
    S-matrix `S_free(E)` for the diagonal channel (`v' == v_init`) when
    `free_result` is supplied, else the literal `1`. `S_free != 1` in general:
    the outgoing normalization factor C(E) multiplies EVERY channel's S (so the
    inelastic `|S|^2` already absorbs it), but the diagonal `|S - ref|^2` only
    isolates genuine scattering if `ref` is the actual unscattered value
    `S_free = C(E)`, not `1`. `free_result` (a `V_int=0` propagation with the
    SAME wavepacket/grid) provides it; the literal-`1` fallback is only correct
    when the transform happens to normalize `S_free -> 1`, which it does not
    here -- so callers wanting a correct elastic channel MUST pass
    `free_result` (`td_ve_cross_section` does by default). Off-diagonal
    channels use `ref = 0` and are unaffected. See
    `.superpowers/sdd/` and the `td-elastic-wavepacket-normalization` note.
    """
    sigma = np.zeros(len(vprimes), dtype=np.float64)
    if E <= 0.0:
        return sigma
    s_full = _s_vector_one_energy(
        tgrid, model, result, eps, v_init, vprimes, E, dt, wp_in, wp_out
    )
    s_free = None
    if free_result is not None:
        s_free = _s_vector_one_energy(
            tgrid, model, free_result, eps, v_init, vprimes, E, dt, wp_in, wp_out
        )
    e_tot = E + eps[v_init]
    for j, vp in enumerate(vprimes):
        if e_tot - eps[vp] <= 0.0:
            continue  # closed channel
        if vp == v_init:
            ref = complex(s_free[j]) if s_free is not None else 1.0 + 0.0j
        else:
            ref = 0.0 + 0.0j
        sigma[j] = np.pi * abs(s_full[j] - ref) ** 2 / (2.0 * E)
    return sigma


def sigma_from_correlations(
    tgrid: TensorGrid,
    model: ResonanceModel,
    result: PropagationResult,
    eps: npt.NDArray[np.float64],
    v_init: int,
    vprimes: list[int],
    E: float | npt.ArrayLike,
    *,
    dt: float,
    wp_in: _WpIn,
    wp_out: _WpOut,
    free_result: PropagationResult | None = None,
) -> npt.NDArray[np.float64]:
    """sigma_{v_init->v'}(E) (bohr^2) from an ALREADY-COMPUTED `PropagationResult`.

    The "cheap sigma(E) from a stored propagation" primitive: unlike
    `td_ve_cross_section`, this never runs (or re-runs) the Pade
    propagation -- it only transforms the `c_{v'}(t)` already sitting in
    `result` (e.g. loaded back from a saved `.npz`, or a truncated/reused
    trajectory from a convergence study).

    `E` (collision energy, Hartree) may be scalar or array-like; scalar `E`
    returns shape `(len(vprimes),)`, array `E` returns `(len(E), len(vprimes))`
    -- the SAME convention as `td_ve_cross_section` and
    `qscat.core.driven.ve_cross_section`.

    `dt`, `wp_in`, `wp_out` must match the values used to produce `result`
    (the quadrature step and the incident/outgoing wavepacket parameters that
    `eta_incident`/`eta_outgoing` are evaluated with) -- this function does
    not validate that consistency, it trusts the caller.

    `free_result` is the `V_int=0` reference propagation (same wavepacket/grid,
    from `_propagate(..., free=True)`); when supplied, the diagonal/elastic
    channel subtracts its `S_free(E)` instead of a literal 1 (see
    `_sigma_one_energy`). Leave `None` to reproduce the old behavior (correct
    for the inelastic channels; the elastic channel then needs `S_free -> 1`,
    which this transform does not satisfy).
    """
    e_arr = np.atleast_1d(np.asarray(E, dtype=np.float64))
    out = np.stack(
        [
            _sigma_one_energy(
                tgrid, model, result, eps, v_init, vprimes, float(e), dt, wp_in, wp_out,
                free_result,
            )
            for e in e_arr
        ]
    )
    scalar = np.isscalar(E) or (isinstance(E, np.ndarray) and np.ndim(E) == 0)
    if scalar:
        return np.asarray(out[0], dtype=np.float64)
    return out


def td_ve_cross_section(
    tgrid: TensorGrid,
    model: ResonanceModel,
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
    order: int = 3,
    subtract_free_reference: bool = True,
) -> npt.NDArray[np.float64]:
    """sigma_{v_init->v'}(E) (bohr^2), Pade propagation + Tannor-Weeks transform.

    `E` (collision energy, Hartree) may be scalar or array-like; scalar `E`
    returns shape `(len(vprimes),)`, array `E` returns `(len(E), len(vprimes))`
    -- matching `qscat.core.driven.ve_cross_section`, the exact differential
    oracle this converges to as `dt -> 0` and `n_steps -> infinity` (see
    module docstring).

    `wp_in = {"r0": ..., "p0": ..., "sigma": ...}` are the SAME incident
    Gaussian parameters used to build `Psi(0)` (via `initial_state`) and
    `eta_incident`; `wp_out = {"r0_out": ..., "p0_out": ..., "sigma_out": ...}`
    are the outgoing test function's parameters, used for both
    `outgoing_channel` (the `Phi_{v'}` propagated against) and `eta_outgoing`.
    The propagation (the expensive part) happens ONCE regardless of how many
    energies `E` are requested, since `c_{v'}(t)` does not depend on `E`.

    `subtract_free_reference` (default `True`): when the diagonal/elastic
    channel is requested (`v_init in vprimes`), a SECOND `V_int=0` propagation
    is run to supply the free-particle reference `S_free(E)` that the elastic
    channel subtracts (instead of a literal 1) -- required for a correct
    elastic cross section, see `_sigma_one_energy`. It doubles the propagation
    cost and is a no-op (skipped) when the elastic channel is not requested;
    set `False` to force the old literal-1 behavior. The inelastic channels
    are identical either way.
    """
    result = _propagate(
        tgrid, model, eps, chi, v_init, vprimes,
        dt=dt, n_steps=n_steps, wp_in=wp_in, wp_out=wp_out, order=order,
    )
    free_result = None
    if subtract_free_reference and v_init in vprimes:
        free_result = _propagate(
            tgrid, model, eps, chi, v_init, vprimes,
            dt=dt, n_steps=n_steps, wp_in=wp_in, wp_out=wp_out, free=True, order=order,
        )
    return sigma_from_correlations(
        tgrid, model, result, eps, v_init, vprimes, E,
        dt=dt, wp_in=wp_in, wp_out=wp_out, free_result=free_result,
    )
