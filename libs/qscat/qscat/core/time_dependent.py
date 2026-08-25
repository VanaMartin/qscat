"""Time-dependent VE cross section: Pade propagation + the Tannor-Weeks energy
transform, for any `model: qscat.model.ResonanceModel`.

Merges `projects/n2_2d_td_cross_section/td_propagation.py` (the propagation
engine: propagate `Psi(0)` under `H` and sample) and `td_cross_section.py`
(the Tannor-Weeks transform) into one model-taking module. The only change
from those two N2 project modules is where the
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

Tannor & Weeks, J. Chem. Phys. 98, 3884 (1993), Eq. (39); implemented as in
`eMoScat TestFunction2d.cpp:298-307`:

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
`_sigma_one_energy` and `docs/physics/n2-2d-td-cross-section.md`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Protocol

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
    "Extractor",
    "Snapshot",
    "PropagationResult",
    "propagate",
    "free_hamiltonian",
    "sigma_from_correlations",
    "td_ve_cross_section",
    "td_ve_cross_sections_all",
    "td_da_cross_section",
    "td_da_cross_sections_all",
    "quadrature_weights",
]


class Extractor(Protocol):
    """A recorder+transform pair driven by one shared `propagate` trajectory.

    `record` is called with the current `psi(t_n)` at EVERY propagation step
    (accumulating whatever per-step datum the extractor needs -- e.g.
    `TannorWeeks` in `td_extractors.py` appends `c_product(Phi_v', psi)` per
    `v'`); `sigma` transforms the accumulated series into a cross section,
    shape `(len(E), len(vprimes))`. Letting `propagate` drive a LIST of these
    means one Pade trajectory can feed several alternative energy-extraction
    routes (Tannor-Weeks, and the delta/flow methods -- `TannorWeeks`,
    `Dirac`, `Flux` in `td_extractors.py`) without re-propagating.
    """

    def record(self, psi: npt.NDArray[np.complex128]) -> None: ...

    def sigma(
        self, E: float | npt.ArrayLike, *, free: Extractor | None = None
    ) -> npt.NDArray[np.float64]: ...


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
    extractors: list[Extractor] | None = None,
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

    This is the propagate-ONCE engine: the trajectory is computed a single
    time and `ex.record(psi)` is called on every extractor in `extractors`
    (see the `Extractor` protocol) at every step, alongside the legacy
    `out_channels` correlation bookkeeping that fills `PropagationResult.c`
    (kept for the existing callers/tests that read `.c` directly -- e.g. the
    N2 project's `td_propagation`/`td_cross_section` shims and
    `observation.py`). Pass `out_channels=[]` when only `extractors` are
    wanted (e.g. `td_ve_cross_section`'s `method="tw"` route).
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
        for ex in extractors or ():
            ex.record(psi)
        norm[n] = float(np.linalg.norm(psi))  # Hermitian L2: physical, monotone
        if n in snap_set:
            rho_R, rho_r = _densities(tgrid, psi)
            snapshots.append(
                Snapshot(float(t[n]), rho_R, rho_r, psi.copy() if n in keep_set else None)
            )
        if n < n_steps:
            psi = step(psi)

    return PropagationResult(t=t, c=c, norm=norm, snapshots=snapshots)


def quadrature_weights(n_t: int) -> npt.NDArray[np.float64]:
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


def free_hamiltonian(model: ResonanceModel, tgrid: TensorGrid) -> sp.spmatrix:
    """`model.hamiltonian(tgrid)` with the interaction `V_int` removed.

    The unscattered (`V_int=0`) reference Hamiltonian used by the elastic
    free-reference propagation (`_propagate`'s `free=True` and
    `td_ve_cross_section`'s `subtract_free_reference` path) -- see
    `_sigma_one_energy` for why the elastic channel needs this reference
    instead of a literal `1`. Public so a caller assembling its own
    `propagate(..., hamiltonian=...)` free-reference run (e.g. the
    `qscat-run` CLI) reuses the exact same reference this module does,
    rather than reaching into a private helper.

    NAME COLLISION: `qscat.core.nrm.scattering.free_hamiltonian` is a
    different function. That one is the 1-D electronic `T_r + centrifugal`
    operator with NO molecular potential at all; this one is the FULL 2-D
    `model.hamiltonian` with only the electron-molecule interaction removed.
    """
    return (model.hamiltonian(tgrid) - sp.diags(model.interaction_diag(tgrid))).tocsr()


# Alias for the same function under its original private name. This module's
# own call sites (`_propagate`, the `td_*_cross_section` free-reference paths)
# and `libs/qscat/tests/test_td_extractors.py` both import it by this name;
# keeping the alias saves churning them all for no behavioral gain.
_free_hamiltonian = free_hamiltonian


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
    hamiltonian = _free_hamiltonian(model, tgrid) if free else model.hamiltonian(tgrid)
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
    weights = quadrature_weights(result.t.size)
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
    channels use `ref = 0` and are unaffected. Tannor & Weeks state this
    directly: the transform's normalization is the wavepacket's own momentum
    amplitude, not 1 (J. Chem. Phys. 98, 3884 (1993), Eqs. 29-32). See
    `docs/physics/n2-2d-td-cross-section.md`.
    """
    sigma = np.zeros(len(vprimes), dtype=np.float64)
    if E <= 0.0:
        return sigma
    s_full = _s_vector_one_energy(tgrid, model, result, eps, v_init, vprimes, E, dt, wp_in, wp_out)
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
                tgrid,
                model,
                result,
                eps,
                v_init,
                vprimes,
                float(e),
                dt,
                wp_in,
                wp_out,
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
    method: str = "tw",
    position: int | None = None,
    surface: int | None = None,
) -> npt.NDArray[np.float64]:
    """sigma_{v_init->v'}(E) (bohr^2), Pade propagation + an energy-extraction
    method's transform: `method="tw"` (default, Tannor-Weeks -- a propagated
    Gaussian test packet), `"delta"` (a fixed-point line projection, eMoScat's
    `DiracTestFunction2d`), or `"flow"` (a fixed-surface Wronskian flux,
    eMoScat's `FluxTestFunction2d`) -- see `td_extractors.py` for the three
    `Extractor` implementations and their formulas.

    `E` (collision energy, Hartree) may be scalar or array-like; scalar `E`
    returns shape `(len(vprimes),)`, array `E` returns `(len(E), len(vprimes))`
    -- matching `qscat.core.driven.ve_cross_section`, the exact differential
    oracle this converges to as `dt -> 0` and `n_steps -> infinity` (see
    module docstring).

    `wp_in = {"r0": ..., "p0": ..., "sigma": ...}` are the SAME incident
    Gaussian parameters used to build `Psi(0)` (via `initial_state`) and
    `eta_incident`; `wp_out = {"r0_out": ..., "p0_out": ..., "sigma_out": ...}`
    are the outgoing test function's parameters -- used only by `method="tw"`
    (`outgoing_channel`/`eta_outgoing`); ignored by `"delta"`/`"flow"`. The
    propagation (the expensive part) happens ONCE regardless of how many
    energies `E` are requested, since the recorded per-step series does not
    depend on `E`.

    `position` (required, keyword-only, for `method="delta"`) and `surface`
    (required, keyword-only, for `method="flow"`) are the fixed electronic
    DVR index the `Dirac`/`Flux` extractor analyzes at -- an index in the
    real (unscaled) region past the interaction, mirroring `wp_out`'s
    asymptotic standoff (see `Dirac`/`Flux`'s own docstrings for the
    validity requirement). Omitting the one the selected `method` needs
    raises `ValueError`.

    `subtract_free_reference` (default `True`): when the diagonal/elastic
    channel is requested (`v_init in vprimes`), a SECOND `V_int=0` propagation
    is run to supply the free-particle reference `S_free(E)` that the elastic
    channel subtracts (instead of a literal 1) -- required for a correct
    elastic cross section, see `_sigma_one_energy`. It doubles the propagation
    cost and is a no-op (skipped) when the elastic channel is not requested;
    set `False` to force the old literal-1 behavior. The inelastic channels
    are identical either way.

    `method="tw"` builds a `td_extractors.TannorWeeks` extractor (and a
    second one for the free reference, if applicable), runs `propagate` with
    `out_channels=[]` (only the extractor(s) record), and returns
    `extractor.sigma(E)` -- reproducing this function's pre-refactor
    implementation (direct `_propagate` + `sigma_from_correlations`) to
    machine precision; see `libs/qscat/tests/test_td_extractors.py`'s golden
    regression test. `method="delta"`/`"flow"` mirror the SAME pattern with
    `td_extractors.Dirac`/`Flux` in place of `TannorWeeks`. Any other
    `method` raises `ValueError`.
    """
    from .td_extractors import Dirac, Flux, TannorWeeks  # deferred: avoids an import cycle

    psi0 = initial_state(tgrid, chi[v_init], **wp_in)
    hamiltonian = model.hamiltonian(tgrid)

    if method == "tw":
        tw = TannorWeeks(tgrid, model, eps, chi, v_init, vprimes, wp_out, wp_in=wp_in, dt=dt)
        propagate(
            tgrid,
            psi0,
            [],
            dt=dt,
            n_steps=n_steps,
            hamiltonian=hamiltonian,
            order=order,
            extractors=[tw],
        )

        tw_free = None
        if subtract_free_reference and v_init in vprimes:
            free_ham = _free_hamiltonian(model, tgrid)
            tw_free = TannorWeeks(
                tgrid, model, eps, chi, v_init, vprimes, wp_out, wp_in=wp_in, dt=dt
            )
            propagate(
                tgrid,
                psi0,
                [],
                dt=dt,
                n_steps=n_steps,
                hamiltonian=free_ham,
                order=order,
                extractors=[tw_free],
            )

        return tw.sigma(E, free=tw_free)

    if method == "delta":
        if position is None:
            raise ValueError("td_ve_cross_section: method='delta' requires `position`")
        dirac = Dirac(tgrid, model, eps, chi, v_init, vprimes, position, wp_in=wp_in, dt=dt)
        propagate(
            tgrid,
            psi0,
            [],
            dt=dt,
            n_steps=n_steps,
            hamiltonian=hamiltonian,
            order=order,
            extractors=[dirac],
        )

        dirac_free = None
        if subtract_free_reference and v_init in vprimes:
            free_ham = _free_hamiltonian(model, tgrid)
            dirac_free = Dirac(
                tgrid, model, eps, chi, v_init, vprimes, position, wp_in=wp_in, dt=dt
            )
            propagate(
                tgrid,
                psi0,
                [],
                dt=dt,
                n_steps=n_steps,
                hamiltonian=free_ham,
                order=order,
                extractors=[dirac_free],
            )

        return dirac.sigma(E, free=dirac_free)

    if method == "flow":
        if surface is None:
            raise ValueError("td_ve_cross_section: method='flow' requires `surface`")
        flux = Flux(tgrid, model, eps, chi, v_init, vprimes, surface, wp_in=wp_in, dt=dt)
        propagate(
            tgrid,
            psi0,
            [],
            dt=dt,
            n_steps=n_steps,
            hamiltonian=hamiltonian,
            order=order,
            extractors=[flux],
        )

        flux_free = None
        if subtract_free_reference and v_init in vprimes:
            free_ham = _free_hamiltonian(model, tgrid)
            flux_free = Flux(tgrid, model, eps, chi, v_init, vprimes, surface, wp_in=wp_in, dt=dt)
            propagate(
                tgrid,
                psi0,
                [],
                dt=dt,
                n_steps=n_steps,
                hamiltonian=free_ham,
                order=order,
                extractors=[flux_free],
            )

        return flux.sigma(E, free=flux_free)

    raise ValueError(
        f"td_ve_cross_section: unknown method {method!r} (must be one of 'tw', 'delta', 'flow')"
    )


def td_ve_cross_sections_all(
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
    position: int,
    surface: int,
    order: int = 3,
    subtract_free_reference: bool = True,
) -> dict[str, npt.NDArray[np.float64]]:
    """`{"tw": sigma, "delta": sigma, "flow": sigma}` (each bohr^2, same shape
    convention as `td_ve_cross_section`) from ONE shared propagation.

    This is the HONEST three-way comparison: `TannorWeeks`, `Dirac`, and
    `Flux` are all built up front and driven by a SINGLE `propagate(...,
    extractors=[tw, dirac, flux])` call (and a single companion `V_int=0`
    propagation, when `subtract_free_reference` applies) -- identical
    dynamics `psi(t_n)` feed all three transforms, so any spread between the
    returned cross sections reflects a genuine difference between the
    energy-extraction methods (or, at an under-converged grid/propagation
    length, a shared discretization/truncation residual all three inherit
    together -- see `docs/physics/td-extractors.md`), never a difference in
    what was propagated.

    `position`/`surface` are the `Dirac`/`Flux` fixed electronic DVR indices
    (see `td_ve_cross_section`'s docstring); both are required here (unlike
    `td_ve_cross_section`, which only needs whichever one its `method`
    selects).
    """
    from .td_extractors import Dirac, Flux, TannorWeeks  # deferred: avoids an import cycle

    psi0 = initial_state(tgrid, chi[v_init], **wp_in)
    hamiltonian = model.hamiltonian(tgrid)

    tw = TannorWeeks(tgrid, model, eps, chi, v_init, vprimes, wp_out, wp_in=wp_in, dt=dt)
    dirac = Dirac(tgrid, model, eps, chi, v_init, vprimes, position, wp_in=wp_in, dt=dt)
    flux = Flux(tgrid, model, eps, chi, v_init, vprimes, surface, wp_in=wp_in, dt=dt)
    propagate(
        tgrid,
        psi0,
        [],
        dt=dt,
        n_steps=n_steps,
        hamiltonian=hamiltonian,
        order=order,
        extractors=[tw, dirac, flux],
    )

    tw_free = dirac_free = flux_free = None
    if subtract_free_reference and v_init in vprimes:
        free_ham = _free_hamiltonian(model, tgrid)
        tw_free = TannorWeeks(tgrid, model, eps, chi, v_init, vprimes, wp_out, wp_in=wp_in, dt=dt)
        dirac_free = Dirac(tgrid, model, eps, chi, v_init, vprimes, position, wp_in=wp_in, dt=dt)
        flux_free = Flux(tgrid, model, eps, chi, v_init, vprimes, surface, wp_in=wp_in, dt=dt)
        propagate(
            tgrid,
            psi0,
            [],
            dt=dt,
            n_steps=n_steps,
            hamiltonian=free_ham,
            order=order,
            extractors=[tw_free, dirac_free, flux_free],
        )

    return {
        "tw": tw.sigma(E, free=tw_free),
        "delta": dirac.sigma(E, free=dirac_free),
        "flow": flux.sigma(E, free=flux_free),
    }


def td_da_cross_section(
    tgrid: TensorGrid,
    model: ResonanceModel,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    E: float | npt.ArrayLike,
    *,
    dt: float,
    n_steps: int,
    wp_in: _WpIn,
    method: str = "flow",
    surface: int | None = None,
    position: int | None = None,
    wp_out: _WpOut | None = None,
    n_channels: int = 1,
    order: int = 3,
) -> npt.NDArray[np.float64]:
    """sigma_DA(E) (bohr^2), Pade propagation + a NUCLEAR-axis energy-
    extraction method's transform -- the dissociative-attachment (DA)
    sibling of `td_ve_cross_section`, built on the SAME propagation engine
    but the `axis="nuclear"` `Flux`/`Dirac`/`TannorWeeks` extractors
    (`td_extractors.py`) instead of the electronic-axis ones.

    `method="flow"` (default, eMoScat's `FluxTestFunction2d` `axis_=='y'`
    branch -- the natural DA extractor: a fixed-surface Wronskian flux
    needs no propagated outgoing test packet), `"delta"` (a fixed nuclear-
    node point projection), or `"tw"` (a propagated nuclear outgoing
    Gaussian test packet) -- see `td_extractors.py`'s module docstring for
    the three nuclear-axis transforms and the `_C_DA = pi` reconciliation
    with the TI `dissociation.da_cross_section` oracle's `4*pi^3` T-matrix
    convention.

    `E` (collision energy, Hartree) may be scalar or array-like; scalar `E`
    returns shape `(n_channels,)`, array `E` returns `(len(E), n_channels)`
    -- matching `dissociation.da_cross_section`'s per-anion-channel
    contract (DA has no `v'` vibrational index; `n_channels` selects how
    many anion electronic bound states, `dissociation.
    anion_electronic_states` at `R_inf = tgrid.grids[1].R0`, are tracked as
    exit channels).

    `wp_in = {"r0": ..., "p0": ..., "sigma": ...}` is the SAME incident
    Gaussian parameters used to build `Psi(0)` (`initial_state`) and
    `eta_incident` -- identical to `td_ve_cross_section`'s `wp_in`. `surface`
    (required for `method="flow"`), `position` (required for
    `method="delta"`), and `wp_out` (required for `method="tw"`, now the
    NUCLEAR outgoing test packet's `r0_out`/`p0_out`/`sigma_out` in R) are
    the method-specific extractor parameters; omitting the one `method`
    needs raises `ValueError`. Unlike `td_ve_cross_section`, there is NO
    `wp_out`/free-reference for `"flow"`/`"delta"` (no propagated test
    packet needed there) and NO elastic free-reference subtraction for any
    method -- DA is a pure rearrangement channel with no `v'==v_init`
    diagonal to subtract a reference from (see `td_extractors.py`'s
    `Flux`/`Dirac`/`TannorWeeks(axis="nuclear")`'s `sigma` docstrings, which
    all raise on a non-`None` `free`).

    The propagation (the expensive part) happens ONCE regardless of how many
    energies `E` are requested, exactly as `td_ve_cross_section`.
    """
    from .td_extractors import Dirac, Flux, TannorWeeks  # deferred: avoids an import cycle

    psi0 = initial_state(tgrid, chi[v_init], **wp_in)
    hamiltonian = model.hamiltonian(tgrid)

    ext: Extractor
    if method == "flow":
        if surface is None:
            raise ValueError("td_da_cross_section: method='flow' requires `surface`")
        ext = Flux(
            tgrid,
            model,
            eps,
            chi,
            v_init,
            [],
            surface,
            wp_in=wp_in,
            dt=dt,
            axis="nuclear",
            n_channels=n_channels,
        )
    elif method == "delta":
        if position is None:
            raise ValueError("td_da_cross_section: method='delta' requires `position`")
        ext = Dirac(
            tgrid,
            model,
            eps,
            chi,
            v_init,
            [],
            position,
            wp_in=wp_in,
            dt=dt,
            axis="nuclear",
            n_channels=n_channels,
        )
    elif method == "tw":
        if wp_out is None:
            raise ValueError("td_da_cross_section: method='tw' requires `wp_out`")
        ext = TannorWeeks(
            tgrid,
            model,
            eps,
            chi,
            v_init,
            [],
            wp_out,
            wp_in=wp_in,
            dt=dt,
            axis="nuclear",
            n_channels=n_channels,
        )
    else:
        raise ValueError(
            f"td_da_cross_section: unknown method {method!r} (must be one of 'flow', 'delta', 'tw')"
        )

    propagate(
        tgrid,
        psi0,
        [],
        dt=dt,
        n_steps=n_steps,
        hamiltonian=hamiltonian,
        order=order,
        extractors=[ext],
    )
    return ext.sigma(E)


def td_da_cross_sections_all(
    tgrid: TensorGrid,
    model: ResonanceModel,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    E: float | npt.ArrayLike,
    *,
    dt: float,
    n_steps: int,
    wp_in: _WpIn,
    surface: int,
    position: int,
    wp_out: _WpOut,
    n_channels: int = 1,
    order: int = 3,
) -> dict[str, npt.NDArray[np.float64]]:
    """`{"flow": sigma_DA, "delta": sigma_DA, "tw": sigma_DA}` (each bohr^2,
    per anion channel, same shape convention as `td_da_cross_section`) from
    ONE shared propagation -- the DA sibling of `td_ve_cross_sections_all`.

    `Flux`, `Dirac`, and `TannorWeeks` (all `axis="nuclear"`) are built up
    front and driven by a SINGLE `propagate(..., extractors=[flux, dirac,
    tw])` call -- identical dynamics `psi(t_n)` feed all three transforms,
    so any spread between the returned cross sections reflects a genuine
    difference between the energy-extraction methods (or, at an under-
    converged grid/propagation length, a shared discretization/truncation
    residual all three inherit together -- see `docs/physics/td-da.md`),
    never a difference in what was propagated. No elastic free-reference
    run (DA has no elastic diagonal, see `td_da_cross_section`).

    `surface`/`position`/`wp_out` are all REQUIRED here (unlike
    `td_da_cross_section`, which only needs whichever one its `method`
    selects).
    """
    from .td_extractors import Dirac, Flux, TannorWeeks  # deferred: avoids an import cycle

    psi0 = initial_state(tgrid, chi[v_init], **wp_in)
    hamiltonian = model.hamiltonian(tgrid)

    flux = Flux(
        tgrid,
        model,
        eps,
        chi,
        v_init,
        [],
        surface,
        wp_in=wp_in,
        dt=dt,
        axis="nuclear",
        n_channels=n_channels,
    )
    dirac = Dirac(
        tgrid,
        model,
        eps,
        chi,
        v_init,
        [],
        position,
        wp_in=wp_in,
        dt=dt,
        axis="nuclear",
        n_channels=n_channels,
    )
    tw = TannorWeeks(
        tgrid,
        model,
        eps,
        chi,
        v_init,
        [],
        wp_out,
        wp_in=wp_in,
        dt=dt,
        axis="nuclear",
        n_channels=n_channels,
    )
    propagate(
        tgrid,
        psi0,
        [],
        dt=dt,
        n_steps=n_steps,
        hamiltonian=hamiltonian,
        order=order,
        extractors=[flux, dirac, tw],
    )

    return {
        "flow": flux.sigma(E),
        "delta": dirac.sigma(E),
        "tw": tw.sigma(E),
    }
