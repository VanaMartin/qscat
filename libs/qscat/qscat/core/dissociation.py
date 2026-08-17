"""Exact 2-D TI dissociative attachment (DA) cross section.

DA is the transient anion's second exit channel -- `e- + AB(v=0) -> AB-* ->
A + B-`, outgoing flux in the NUCLEAR coordinate. It is computed by the SAME
driven Lippmann-Schwinger solve as VE (`qscat.core.driven.ve_cross_section`,
`return_wavefunction=True`) but projected onto the dissociation channel with
the REARRANGEMENT interaction

    V_DR(r, R) = V_int(r, R) + v0(R) - V_int(r, R_inf)    (= H - H_final),

NOT V_int. The exit channel is Phi_n(r,R) = phi_e^(n)(r) F^nuc_{K_n,0}(R),
phi_e the anion bound electronic state at the dissociation limit R_inf and
F^nuc the mass-mu energy-normalized regular nuclear Bessel; the T-matrix is
T_n = <Phi_n | V_DR | Psi+> (c-product, masked), sigma_n = 4 pi^3 |T_n|^2/(2E).

This is eMoScat's `time_independent_model.cpp` method (an earlier prototype
that used V_int instead of V_DR gave a ~1e6 unitarity violation -- that was
the bug, not a structural obstacle to a TI DA). `dr_cross_section` is the
same T-matrix GENERALIZED for H2+: looped over the neutral's Rydberg
electronic series (`n_channels` states of the same bound-electronic-state
solver) with a Coulomb incident (`channel_vector(..., charge=model.charge)`)
(sub-project D). See docs/physics/diatomic-ve-cross-sections.md and
docs/superpowers/specs/2026-07-27-da-cross-sections-design.md.

`qscat.core` never imports `qscat.model` at runtime: `model` is typed against
the `ResonanceModel` protocol under `TYPE_CHECKING` only, exactly like
`driven.py`.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal, overload

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp

from qscat.dvr import FemDvrEcsGrid, TensorGrid, eigen, kinetic
from qscat.linalg import SparseLU, c_product
from qscat.special import riccati_bessel_en_mass

from .channels import channel_vector
from .driven import ve_cross_section

if TYPE_CHECKING:
    from qscat.model import ResonanceModel

__all__ = ["anion_electronic_states", "v_dr_diag", "da_cross_section", "dr_cross_section"]

# Mirrors driven.py's re-declaration of SparseLU's private ordering Literal,
# so `ordering` passes through to ve_cross_section type-clean.
_Ordering = Literal["NATURAL", "MMD_ATA", "MMD_AT_PLUS_A", "COLAMD"]

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
    ordering: _Ordering = ...,
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
    ordering: _Ordering = ...,
    return_wavefunction: Literal[True],
) -> tuple[_Sigma, _PsiOut]: ...


def da_cross_section(
    tgrid: TensorGrid,
    model: ResonanceModel,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    E: float | npt.ArrayLike,
    *,
    n_channels: int = 1,
    ordering: _Ordering = "COLAMD",
    return_wavefunction: bool = False,
) -> _Sigma | tuple[_Sigma, _PsiOut]:
    """sigma_DA(E) in bohr^2, exact 2-D driven-equation DA cross section.

    Reuses `ve_cross_section(..., return_wavefunction=True)` for `Psi+` (one
    `SparseLU.refactor` sweep across `E`), then projects onto each of
    `n_channels` anion dissociation channels with `V_DR`. `E` may be scalar
    (returns `(n_channels,)`) or an array (returns `(len(E), n_channels)`).
    `sigma = 0` for a closed channel (`E <= 0` or `E_DR = E_tot - eps_e <= 0`).

    If `return_wavefunction`, also returns the driven `Psi+` (the SAME solution
    the T-matrix is built from; `None` for `E <= 0`): one array for scalar `E`,
    one list entry per energy for an array `E` -- same convention as
    `ve_cross_section`.
    """
    e_arr = np.atleast_1d(np.asarray(E, dtype=np.float64))
    mu = model.mu
    g_R = tgrid.grids[1]
    R_inf = g_R.R0

    eps_e, phi = anion_electronic_states(
        g_r=tgrid.grids[0], model=model, R_inf=R_inf, n_states=n_channels
    )
    v_dr = v_dr_diag(tgrid, model)
    mask = tgrid.real_mask()
    sqrt_w_R = tgrid.sqrt_weights()[1].ravel()

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
        v_psi = v_dr * psi_plus
        for n in range(n_channels):
            e_dr = e_tot - eps_e[n]
            if e_dr <= 0.0:
                continue
            k_r = float(np.sqrt(2.0 * mu * e_dr))
            y_coeff = riccati_bessel_en_mass(g_R.real_points, k_r, 0, mu) * sqrt_w_R
            phi_f = tgrid.outer([phi[n], y_coeff])
            phi_f[~mask] = 0.0
            t = c_product(phi_f, v_psi)
            out[ie, n] = 4.0 * np.pi**3 * abs(t) ** 2 / (2.0 * float(e))

    scalar = np.isscalar(E) or (isinstance(E, np.ndarray) and np.ndim(E) == 0)
    sigma = np.asarray(out[0] if scalar else out, dtype=np.float64)
    if return_wavefunction:
        return sigma, (psi_list[0] if scalar else psi_list)
    return sigma


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
    ordering: _Ordering = ...,
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
    ordering: _Ordering = ...,
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
    ordering: _Ordering = ...,
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
    ordering: _Ordering = ...,
    return_wavefunction: Literal[True],
    return_amplitude: Literal[True],
) -> tuple[_Sigma, _PsiOut, _Amp]: ...


def dr_cross_section(
    tgrid: TensorGrid,
    model: ResonanceModel,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    E: float | npt.ArrayLike,
    *,
    n_channels: int = 3,
    ordering: _Ordering = "COLAMD",
    return_wavefunction: bool = False,
    return_amplitude: bool = False,
) -> _Sigma | tuple[_Sigma, _PsiOut] | tuple[_Sigma, _Amp] | tuple[_Sigma, _PsiOut, _Amp]:
    """sigma_DR(E) in bohr^2, exact 2-D driven-equation dissociative-
    recombination cross section for a CHARGED target (e.g. H2+, `charge=-1`).

    `da_cross_section` GENERALIZED two ways: (1) the incident channel is
    Coulomb rather than free -- `channel_vector(..., charge=model.charge)`;
    (2) the exit channel is a LOOP over `n_channels` Rydberg electronic
    states of the neutral, `anion_electronic_states(..., n_states=n_channels)`
    (the SAME bound-electronic-state solver DA uses for the anion state --
    the Rydberg series is likewise bound, below the `-1/r` Coulomb
    continuum). `V_DR` (the rearrangement interaction, NOT `V_int`) is
    unchanged from DA.

    The driven Lippmann-Schwinger solve for `Psi+` is replicated inline
    (rather than reusing `ve_cross_section`) because that helper cannot pass
    `charge` through to `channel_vector`; the sparse LU is still built once
    and `refactor`-ed across the energy sweep exactly as `ve_cross_section`
    does.

    `E` may be scalar (returns `(n_channels,)`) or an array (returns
    `(len(E), n_channels)`). `sigma_n = 0` for a closed channel (`E <= 0` or
    `E_DR = E_tot - E_ryd(n) <= 0`, `E_ryd(n) = eps_e[n]`).

    If `return_wavefunction`, also returns the driven `Psi+` (`None` for
    `E <= 0`): one array for scalar `E`, one list entry per energy for an array
    `E` -- same convention as `ve_cross_section`/`da_cross_section`.

    If `return_amplitude`, also returns the complex transition amplitude `t`
    the T-matrix sum is built from, shaped exactly like `sigma` (`(n_channels,)`
    for scalar `E`, `(len(E), n_channels)` for an array `E`, zero for a closed
    channel): `sigma = 4*pi**3 * abs(t)**2 / (2*E)`. This is the amplitude the
    solver already forms, NOT a literal unitary S-matrix element -- the
    thesis's `S_DR` differs from it by the standard `S = -2*pi*i*T` factor and
    its own normalization. That factor is a fixed rotation and rescale: it
    changes neither the zeros nor the shape of `t`'s real/imaginary crossings,
    so returning `t` as computed (rather than guessing a normalization to
    synthesize an "S") is what downstream resonance-pole fitting needs.
    If both `return_wavefunction` and `return_amplitude`, the return order is
    `(sigma, psi, amplitude)`.
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

    H = model.hamiltonian(tgrid)
    v_diag = model.interaction_diag(tgrid)
    ident = sp.identity(tgrid.size, format="csc", dtype=np.complex128)

    out = np.zeros((len(e_arr), n_channels), dtype=np.float64)
    amp = np.zeros((len(e_arr), n_channels), dtype=np.complex128)
    psi_list: list[_Psi] = [None] * len(e_arr)
    lu: SparseLU | None = None
    for ie, e in enumerate(e_arr):
        if float(e) <= 0.0:
            continue  # below threshold: no driven-equation solve, sigma == 0

        e_tot = float(e) + eps[v_init]
        a = (e_tot * ident - H).tocsc()
        if lu is None:
            lu = SparseLU(a, ordering=ordering)
        else:
            lu.refactor(a)

        k = float(np.sqrt(2.0 * float(e)))
        psi_i = channel_vector(tgrid, k, chi[v_init], model.ell, charge=model.charge)
        psi_plus = psi_i + lu.solve(v_diag * psi_i)
        psi_list[ie] = psi_plus
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
    amplitude = np.asarray(amp[0] if scalar else amp, dtype=np.complex128)
    psi_out = psi_list[0] if scalar else psi_list
    if return_wavefunction and return_amplitude:
        return sigma, psi_out, amplitude
    if return_amplitude:
        return sigma, amplitude
    if return_wavefunction:
        return sigma, psi_out
    return sigma
