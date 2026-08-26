"""Exact 2-D TI dissociative attachment (DA) cross section.

DA is the transient anion's second exit channel -- `e- + AB(v=0) -> AB-* ->
A + B-`, outgoing flux in the NUCLEAR coordinate. It is computed from the SAME
driven Lippmann-Schwinger solve as VE (`qscat.core.driven.ve_cross_section`,
`return_wavefunction=True`); the two cross sections here differ only in how
they read an amplitude out of that `Psi+`.

`da_cross_section` reads the **outgoing flux**: project `Psi+` onto the anion
bound electronic state `phi_e^(n)` (at the dissociation limit `R_inf`) to get
the channel's nuclear wave `psi_n(R)`, take its VALUE at the outermost real
nuclear node `X`, and

    S_n = sqrt(K_n / 2 pi mu) psi_n(X),   sigma_n = 4 pi^3 |S_n|^2 / (2E).

`dr_cross_section` instead uses the post-form **volume T-matrix** against the
REARRANGEMENT interaction

    V_DR(r, R) = V_int(r, R) + v0(R) - V_int(r, R_inf)    (= H - H_final),

with the exit channel Phi_n(r,R) = phi_e^(n)(r) F^nuc_{K_n,0}(R) (`F^nuc` the
mass-mu energy-normalized regular nuclear Bessel), T_n = <Phi_n | V_DR | Psi+>
(c-product, masked), sigma_n = 4 pi^3 |T_n|^2/(2E). That is eMoScat's
`time_independent_model.cpp` method, GENERALIZED for H2+ over the neutral's
Rydberg electronic series with a Coulomb incident
(`channel_vector(..., charge=model.charge)`).

The two routes are algebraically the same amplitude and agree to 5e-4 on F2.
**They are not interchangeable numerically**, which is why DA uses the flux
one: the volume integral delivers `sigma_DA` as the residue of a cancellation
whose depth is set by how small `sigma_DA` is, and on NO (~1e-9 bohr^2, a
~1e6-fold cancellation) it returned answers 1e4-1e7 times too large until BOTH
integration edges were pushed far out. See `da_cross_section`'s note and
docs/physics/diatomic-ve-cross-sections.md.

`qscat.core` never imports `qscat.model` at runtime: `model` is typed against
the `ResonanceModel` protocol under `TYPE_CHECKING` only, exactly like
`driven.py`.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal, overload

import numpy as np
import numpy.typing as npt

from qscat.dvr import FemDvrEcsGrid, TensorGrid, eigen, kinetic
from qscat.linalg import Ordering, c_product
from qscat.special import riccati_bessel_en_mass

from .driven import ve_cross_section

if TYPE_CHECKING:
    from qscat.model import ResonanceModel

__all__ = [
    "DrResult",
    "anion_electronic_states",
    "da_cross_section",
    "dr_cross_section",
    "dr_solve",
    "v_dr_diag",
]

# `return_wavefunction` output types, same convention as driven.py: the driven
# Psi+ per energy (`None` below threshold), one array for scalar `E`, one list
# entry per energy for an array `E`. The SAME Psi+ the cross section is built
# from -- exposed so a caller can snapshot/animate it (e.g. qscat_run).
_Sigma = npt.NDArray[np.float64]
_Psi = npt.NDArray[np.complex128] | None
_PsiOut = _Psi | list[_Psi]
# `return_amplitude` output: the complex T-matrix amplitude `t` itself, shaped
# exactly like `sigma` -- `(n_channels,)` for scalar `E`, `(len(E), n_channels)`
# for an array `E`. Zero for a closed channel, same as `sigma` there.
_Amp = npt.NDArray[np.complex128]

# Bound-state signature on an ECS grid: true bound levels have |Im(E)| ~ 1e-15,
# ECS-continuum states jump to >= 1e-7. Same tolerance as `vibrational_states`.
_IM_TOL_HA = 1e-6


def anion_electronic_states(
    g_r: FemDvrEcsGrid,
    model: ResonanceModel,
    R_inf: float,
    n_states: int = 1,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.complex128]]:
    """Anion bound electronic state(s) at the dissociation limit `R_inf`.

    Diagonalizes `-1/2 d^2/dr^2 + surface(r, R_inf)` (electron mass 1) on the
    electronic grid. `surface` is the FULL electronic potential at `R_inf`
    (`v0(R_inf) + ell(ell+1)/2r^2 + v_int(r, R_inf)`), so `eps_e` shares the
    `H_2D` energy scale (it includes `v0(R_inf)`) and the DA threshold
    `eps_e - eps[v_init]` is correct.

    Returns `(eps_e, phi_e)`: `eps_e` the `n_states` lowest-Re eigenvalues with
    `|Im(E)| < _IM_TOL_HA` AND `Re(E) < v0(R_inf)` (the genuinely bound
    states), real, ascending; `phi_e` shape `(n_states, g_r.n)`, each
    c-product-normalized over the electronic real region. Raises `ValueError`
    if fewer than `n_states` bound states exist (e.g. `n_states` reached past
    the finite bound spectrum).
    """
    H_el = kinetic(g_r, 1.0) + np.diag(model.surface(g_r.points, R_inf))
    E, V = eigen(H_el)  # ascending Re(E)
    # Genuinely bound: near-real AND below the asymptotic electronic continuum
    # edge v0(R_inf) (as r->inf, v_int and centrifugal vanish -> the electron
    # sees only v0(R_inf)). The |Im| filter alone counts finite-basis
    # "top-of-grid numerical-junk" eigenvalues (large positive Re(E), tiny
    # |Im|) as bound; the Re(E) < e_thresh cut excludes them.
    e_thresh = float(np.real(model.v0(np.asarray(R_inf))))
    bound = np.flatnonzero((np.abs(E.imag) < _IM_TOL_HA) & (E.real < e_thresh))
    if bound.size < n_states:
        raise ValueError(
            f"anion_electronic_states(n_states={n_states}) found only "
            f"{bound.size} bound electronic state(s) (|Im(E)| < {_IM_TOL_HA} Ha "
            f"and Re(E) < v0(R_inf)={e_thresh:.6g}) at R_inf={R_inf}: the well "
            "supports fewer bound states than requested. Reduce n_states."
        )
    idx = bound[:n_states]  # E is Re-ascending, so these are the lowest-Re bound states
    eps_e = E[idx].real
    phi = V[:, idx].T.astype(np.complex128)

    real = g_r.real_points <= g_r.R0
    for i in range(n_states):
        p = phi[i].copy()
        p[~real] = 0.0
        norm2 = c_product(p, p)
        phi[i] = phi[i] / np.sqrt(norm2)
    return eps_e, phi


def v_dr_diag(tgrid: TensorGrid, model: ResonanceModel) -> npt.NDArray[np.complex128]:
    """The rearrangement interaction `V_DR = V_int(r,R) + v0(R) - V_int(r, R_inf)`,
    flat (C-order), length `tgrid.size`. `R_inf = tgrid.grids[1].R0` (the nuclear
    ECS pivot / real-region endpoint, eMoScat's `nu_inf`).

    This -- not `V_int` -- is the operator in the DA/DR T-matrix: `H - H_final`,
    where `H_final` is the asymptotic channel Hamiltonian (electron bound in
    `V_int(r, R_inf)`, free nuclei on `v0`). As `R -> R_inf` the `V_int` terms
    cancel and `V_DR -> v0(R)`.
    """
    R_inf = tgrid.grids[1].R0
    pts_r, pts_R = tgrid.points()
    v0_term = np.broadcast_to(model.v0(pts_R), tgrid.shape).ravel()
    vint_inf = np.broadcast_to(model.v_int(pts_r, R_inf), tgrid.shape).ravel()
    return np.asarray(model.interaction_diag(tgrid) + v0_term - vint_inf, dtype=np.complex128)


@overload
def da_cross_section(
    tgrid: TensorGrid,
    model: ResonanceModel,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    E: float | npt.ArrayLike,
    *,
    n_channels: int = ...,
    ordering: Ordering = ...,
    return_wavefunction: Literal[False] = ...,
) -> _Sigma: ...


@overload
def da_cross_section(
    tgrid: TensorGrid,
    model: ResonanceModel,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    E: float | npt.ArrayLike,
    *,
    n_channels: int = ...,
    ordering: Ordering = ...,
    return_wavefunction: Literal[True],
) -> tuple[_Sigma, _PsiOut]: ...


# bool catch-all (open()-style): callers holding a runtime flag forward it
# directly; the union return is narrowed by the Literal overloads above when
# the flag is literal.
@overload
def da_cross_section(
    tgrid: TensorGrid,
    model: ResonanceModel,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    E: float | npt.ArrayLike,
    *,
    n_channels: int = ...,
    ordering: Ordering = ...,
    return_wavefunction: bool = ...,
) -> _Sigma | tuple[_Sigma, _PsiOut]: ...


def da_cross_section(
    tgrid: TensorGrid,
    model: ResonanceModel,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    E: float | npt.ArrayLike,
    *,
    n_channels: int = 1,
    ordering: Ordering = "COLAMD",
    return_wavefunction: bool = False,
) -> _Sigma | tuple[_Sigma, _PsiOut]:
    """sigma_DA(E) in bohr^2, exact 2-D driven-equation DA cross section.

    *Provisional API* (docs/adr/0004-public-api-stability-policy.md): this wide
    functional signature is the layer the context-object refactor targets and
    may change in a minor release; `ScatteringProblem.da_cross_section` is the stable route.

    Reuses `ve_cross_section(..., return_wavefunction=True)` for `Psi+` (one
    `SparseLU.refactor` sweep across `E`), then reads the outgoing dissociation
    flux out of it, one `n_channels` anion channel at a time. `E` may be scalar
    (returns `(n_channels,)`) or an array (returns `(len(E), n_channels)`).
    `sigma = 0` for a closed channel (`E <= 0` or `E_DR = E_tot - eps_e <= 0`).

    The extraction projects `Psi+` onto the anion electronic state to get the
    channel's nuclear wave, then takes its VALUE at the outermost real nuclear
    node `X` (the DVR coefficient divided by `sqrt(w_b)`, never the coefficient
    itself)::

        psi_n(R) = <phi_e^(n) | Psi+>_r      [c-product over the electronic axis]
        S_n      = sqrt(K_n / 2 pi mu) psi_n(X)
        sigma_n  = 4 pi^3 |S_n|^2 / 2E

    -- algebraically the same amplitude `qscat.core.lcp.lcp_da_cross_section`
    takes off its 1-D nuclear wave, so the exact solver and the LCP are compared
    through identical arithmetic.

    .. note::
       **Why not the post-form volume T-matrix.** The textbook alternative,
       `T = <phi_e^(n) F_K | V_DR | Psi+>` with the rearrangement interaction
       `V_DR = V_int(r,R) + v0(R) - V_int(r,R_inf)` (still available as
       `v_dr_diag`, and what `dr_cross_section` uses), is formally exact and
       agrees with the flux above to 5e-4 on F2. It is nonetheless the WRONG
       tool whenever `sigma_DA` is small: `V_DR` does not decay in `r`, so the
       integrand's magnitude is set by the interaction region while the answer
       is set by how completely that region cancels. On NO -- `sigma_DA ~ 1e-9`
       bohr^2 against an integrand summing to ~2.6 -- the required cancellation
       is ~1e6-fold, and whatever has not decayed at either edge of the
       integration region survives instead of the physics. Both edges did:
       measured on the shipped decks, the T-matrix answer moved by four orders
       when the electronic real region went 16 -> 48 bohr and by a further ~800x
       when the nuclear `R_inf` went 9.0 -> 15.0 bohr (where NO's Morse `v0`
       finally reaches 1e-9 Ha), converging only there onto the flux value.
       The flux extraction needs no cancellation and is invariant under both to
       4 digits. See `docs/physics/diatomic-ve-cross-sections.md` and
       `validation/diatomic/test_no_da_thesis.py`.

    If `return_wavefunction`, also returns the driven `Psi+` (the SAME solution
    the flux is read from; `None` for `E <= 0`): one array for scalar `E`,
    one list entry per energy for an array `E` -- same convention as
    `ve_cross_section`.
    """
    e_arr = np.atleast_1d(np.asarray(E, dtype=np.float64))
    mu = model.mu
    g_r, g_R = tgrid.grids[0], tgrid.grids[1]
    R_inf = g_R.R0

    eps_e, phi = anion_electronic_states(g_r=g_r, model=model, R_inf=R_inf, n_states=n_channels)
    # `phi` is c-product-normalized over the ELECTRONIC real region, so the
    # projection that reads the channel amplitude has to be masked the same way.
    phi_real = phi.copy()
    phi_real[:, g_r.points.imag != 0.0] = 0.0

    # `b` is the outermost real NUCLEAR node -- the boundary X at which the
    # outgoing flux is read, and `sqrt(w_b)` turns its DVR coefficient into the
    # wavefunction value there. Same node `lcp_da_cross_section` uses.
    real_R = np.flatnonzero(g_R.points.imag == 0.0)
    b = int(real_R[np.argmax(g_R.points[real_R].real)])
    sqrt_w_b = np.sqrt(complex(g_R.weights[b]))

    _, psis = ve_cross_section(
        tgrid,
        model,
        eps,
        chi,
        v_init,
        [v_init],
        e_arr,
        ordering=ordering,
        return_wavefunction=True,
    )
    psi_list = psis if isinstance(psis, list) else [psis]

    out = np.zeros((len(e_arr), n_channels), dtype=np.float64)
    for ie, e in enumerate(e_arr):
        psi_plus = psi_list[ie]
        if psi_plus is None:  # E <= 0
            continue
        e_tot = float(e) + eps[v_init]
        psi_2d = psi_plus.reshape(g_r.n, g_R.n)
        for n in range(n_channels):
            e_dr = e_tot - eps_e[n]
            if e_dr <= 0.0:
                continue
            k_r = float(np.sqrt(2.0 * mu * e_dr))
            # c-product over the electronic axis: both factors are already DVR
            # coefficients, so the plain (non-conjugated) contraction IS the
            # quadrature integral -- the ECS convention, as everywhere else here.
            psi_n = phi_real[n] @ psi_2d
            s_da = np.sqrt(k_r / (2.0 * np.pi * mu)) * (psi_n[b] / sqrt_w_b)
            out[ie, n] = 4.0 * np.pi**3 * abs(s_da) ** 2 / (2.0 * float(e))

    scalar = np.isscalar(E) or (isinstance(E, np.ndarray) and np.ndim(E) == 0)
    sigma = np.asarray(out[0] if scalar else out, dtype=np.float64)
    if return_wavefunction:
        return sigma, (psi_list[0] if scalar else psi_list)
    return sigma


@dataclass(frozen=True)
class DrResult:
    """Result of the exact 2-D dissociative-recombination solve (`dr_solve`).

    - `sigma`: sigma_DR(E) in bohr^2 -- `(n_channels,)` for scalar `E`,
      `(len(E), n_channels)` for array `E`.
    - `psi`: the driven `Psi+` per energy (`None` when not stored, and `None`
      per energy below threshold) -- one array for scalar `E`, one list entry
      per energy for array `E`.
    - `amplitude`: the complex T-matrix amplitude, shaped like `sigma`
      (`None` when not stored). See `dr_solve` for the S-vs-T normalization
      note.
    """

    sigma: _Sigma
    psi: _PsiOut | None
    amplitude: _Amp | None


def dr_solve(
    tgrid: TensorGrid,
    model: ResonanceModel,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    E: float | npt.ArrayLike,
    *,
    n_channels: int = 3,
    ordering: Ordering = "COLAMD",
    store_wavefunction: bool = False,
    store_amplitude: bool = False,
) -> DrResult:
    """Exact 2-D driven-equation dissociative-recombination solve for a
    CHARGED target (e.g. H2+, `charge=-1`); returns a `DrResult`.

    *Provisional API* (docs/adr/0004-public-api-stability-policy.md): this wide
    functional signature is the layer the context-object refactor targets and
    may change in a minor release; `ScatteringProblem.dr_cross_section` is the
    stable route.

    `da_cross_section` GENERALIZED two ways: (1) the incident channel is
    Coulomb rather than free -- `channel_vector(..., charge=model.charge)`;
    (2) the exit channel is a LOOP over `n_channels` Rydberg electronic
    states of the neutral, `anion_electronic_states(..., n_states=n_channels)`
    (the SAME bound-electronic-state solver DA uses for the anion state --
    the Rydberg series is likewise bound, below the `-1/r` Coulomb
    continuum). `V_DR` (the rearrangement interaction, NOT `V_int`) is
    unchanged from DA.

    The driven Lippmann-Schwinger solve for `Psi+` reuses
    `ve_cross_section(..., return_wavefunction=True)` -- the same
    analyze-once / `SparseLU.refactor`-per-energy sweep, with `model.charge`
    forwarded to `channel_vector` so the incident channel is Coulomb --
    exactly as `da_cross_section` does. Only the exit-channel read differs:
    DR projects the post-form volume T-matrix against `V_DR` (below), DA
    reads a boundary flux.

    `E` may be scalar (returns `(n_channels,)` fields) or an array (returns
    `(len(E), n_channels)` fields). `sigma_n = 0` for a closed channel
    (`E <= 0` or `E_DR = E_tot - E_ryd(n) <= 0`, `E_ryd(n) = eps_e[n]`).

    `DrResult.psi` is populated (the driven `Psi+`, `None` for `E <= 0`) only
    when `store_wavefunction` is set: one array for scalar `E`, one list
    entry per energy for an array `E` -- same convention as
    `ve_cross_section`/`da_cross_section`.

    `DrResult.amplitude` is populated only when `store_amplitude` is set: the
    complex transition amplitude `t` the T-matrix sum is built from, shaped
    exactly like `sigma` (`(n_channels,)` for scalar `E`, `(len(E),
    n_channels)` for an array `E`, zero for a closed channel):
    `sigma = 4*pi**3 * abs(t)**2 / (2*E)`. This is the amplitude the solver
    already forms, NOT a literal unitary S-matrix element -- the thesis's
    `S_DR` differs from it by the standard `S = -2*pi*i*T` factor and its own
    normalization. That factor is a fixed rotation and rescale: it changes
    neither the zeros nor the shape of `t`'s real/imaginary crossings, so
    returning `t` as computed (rather than guessing a normalization to
    synthesize an "S") is what downstream resonance-pole fitting needs.
    """
    e_arr = np.atleast_1d(np.asarray(E, dtype=np.float64))
    mu = model.mu
    g_R = tgrid.grids[1]
    R_inf = g_R.R0

    eps_ryd, phi_ryd = anion_electronic_states(
        g_r=tgrid.grids[0], model=model, R_inf=R_inf, n_states=n_channels
    )
    v_dr = v_dr_diag(tgrid, model)
    mask = tgrid.real_mask()
    sqrt_w_R = tgrid.sqrt_weights()[1].ravel()

    # The driven Psi+ sweep is ve_cross_section's (analyze-once,
    # `SparseLU.refactor` per energy; `model.charge` forwarded to the
    # incident `channel_vector`, so H2+'s Coulomb entrance is built there)
    # -- exactly the reuse `da_cross_section` already performs. Its VE
    # sigma for the [v_init] channel is discarded; the extra cost is one
    # exit `channel_vector` + c-product per energy, marginal next to the
    # factorization at production scale -- on the small H2+ proxy deck the
    # mpmath Coulomb builds inside `channel_vector` dominate instead.
    _, psis = ve_cross_section(
        tgrid,
        model,
        eps,
        chi,
        v_init,
        [v_init],
        e_arr,
        ordering=ordering,
        return_wavefunction=True,
    )
    psi_list = psis if isinstance(psis, list) else [psis]

    out = np.zeros((len(e_arr), n_channels), dtype=np.float64)
    amp = np.zeros((len(e_arr), n_channels), dtype=np.complex128)
    for ie, e in enumerate(e_arr):
        psi_plus = psi_list[ie]
        if psi_plus is None:  # E <= 0: below threshold, sigma == 0
            continue
        e_tot = float(e) + eps[v_init]
        v_psi = v_dr * psi_plus

        for n in range(n_channels):
            e_dr = e_tot - eps_ryd[n]
            if e_dr <= 0.0:
                continue  # closed Rydberg channel
            k_r = float(np.sqrt(2.0 * mu * e_dr))
            y_coeff = riccati_bessel_en_mass(g_R.real_points, k_r, 0, mu) * sqrt_w_R
            phi_f = tgrid.outer([phi_ryd[n], y_coeff])
            phi_f[~mask] = 0.0
            t = c_product(phi_f, v_psi)
            amp[ie, n] = t
            out[ie, n] = 4.0 * np.pi**3 * abs(t) ** 2 / (2.0 * float(e))

    scalar = np.isscalar(E) or (isinstance(E, np.ndarray) and np.ndim(E) == 0)
    sigma = np.asarray(out[0] if scalar else out, dtype=np.float64)
    return DrResult(
        sigma=sigma,
        psi=(psi_list[0] if scalar else psi_list) if store_wavefunction else None,
        amplitude=(
            np.asarray(amp[0] if scalar else amp, dtype=np.complex128) if store_amplitude else None
        ),
    )


@overload
def dr_cross_section(
    tgrid: TensorGrid,
    model: ResonanceModel,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    E: float | npt.ArrayLike,
    *,
    n_channels: int = ...,
    ordering: Ordering = ...,
    return_wavefunction: Literal[False] = ...,
    return_amplitude: Literal[False] = ...,
) -> _Sigma: ...


@overload
def dr_cross_section(
    tgrid: TensorGrid,
    model: ResonanceModel,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    E: float | npt.ArrayLike,
    *,
    n_channels: int = ...,
    ordering: Ordering = ...,
    return_wavefunction: Literal[True],
    return_amplitude: Literal[False] = ...,
) -> tuple[_Sigma, _PsiOut]: ...


@overload
def dr_cross_section(
    tgrid: TensorGrid,
    model: ResonanceModel,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    E: float | npt.ArrayLike,
    *,
    n_channels: int = ...,
    ordering: Ordering = ...,
    return_wavefunction: Literal[False] = ...,
    return_amplitude: Literal[True],
) -> tuple[_Sigma, _Amp]: ...


@overload
def dr_cross_section(
    tgrid: TensorGrid,
    model: ResonanceModel,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    E: float | npt.ArrayLike,
    *,
    n_channels: int = ...,
    ordering: Ordering = ...,
    return_wavefunction: Literal[True],
    return_amplitude: Literal[True],
) -> tuple[_Sigma, _PsiOut, _Amp]: ...


# bool catch-all (open()-style): callers holding a runtime flag forward it
# directly; the union return is narrowed by the Literal overloads above when
# the flag is literal.
@overload
def dr_cross_section(
    tgrid: TensorGrid,
    model: ResonanceModel,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    E: float | npt.ArrayLike,
    *,
    n_channels: int = ...,
    ordering: Ordering = ...,
    return_wavefunction: bool = ...,
    return_amplitude: bool = ...,
) -> _Sigma | tuple[_Sigma, _PsiOut] | tuple[_Sigma, _Amp] | tuple[_Sigma, _PsiOut, _Amp]: ...


def dr_cross_section(
    tgrid: TensorGrid,
    model: ResonanceModel,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    E: float | npt.ArrayLike,
    *,
    n_channels: int = 3,
    ordering: Ordering = "COLAMD",
    return_wavefunction: bool = False,
    return_amplitude: bool = False,
) -> _Sigma | tuple[_Sigma, _PsiOut] | tuple[_Sigma, _Amp] | tuple[_Sigma, _PsiOut, _Amp]:
    """sigma_DR(E) in bohr^2, exact 2-D driven-equation dissociative-
    recombination cross section for a CHARGED target (e.g. H2+, `charge=-1`).
    See `dr_solve` for the physics; this is a thin flag-shaped-tuple wrapper
    around it.

    *Provisional API* (docs/adr/0004-public-api-stability-policy.md): this wide
    functional signature is the layer the context-object refactor targets and
    may change in a minor release; `ScatteringProblem.dr_cross_section` is the stable route.

    .. deprecated::
       `return_wavefunction`/`return_amplitude`'s flag-shaped tuple returns
       are deprecated in favor of `dr_solve`, which returns one `DrResult`
       object (`store_wavefunction`/`store_amplitude`) regardless of how many
       fields are populated. The plain `sigma`-only call (both flags `False`,
       the undisputed base case) is unaffected and stays silent.
    """
    if return_wavefunction or return_amplitude:
        warnings.warn(
            "dr_cross_section's flag-shaped tuple returns are deprecated; "
            "call dr_solve(..., store_wavefunction=..., store_amplitude=...) "
            "and read the DrResult fields",
            DeprecationWarning,
            stacklevel=2,
        )
    res = dr_solve(
        tgrid,
        model,
        eps,
        chi,
        v_init,
        E,
        n_channels=n_channels,
        ordering=ordering,
        store_wavefunction=return_wavefunction,
        store_amplitude=return_amplitude,
    )
    if return_wavefunction and return_amplitude:
        # store_amplitude=True (== return_amplitude) guarantees this at runtime;
        # the assert narrows the type past DrResult.amplitude's static Optional.
        assert res.amplitude is not None
        return res.sigma, res.psi, res.amplitude
    if return_amplitude:
        return res.sigma, res.amplitude
    if return_wavefunction:
        return res.sigma, res.psi
    return res.sigma
