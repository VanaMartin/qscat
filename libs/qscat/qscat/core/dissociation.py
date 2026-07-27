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
the bug, not a structural obstacle to a TI DA). H2+ DR is the same T-matrix
looped over the neutral's Rydberg electronic series + a Coulomb incident;
deferred (sub-project D). See docs/physics/diatomic-ve-cross-sections.md and
docs/superpowers/specs/2026-07-27-da-cross-sections-design.md.

`qscat.core` never imports `qscat.model` at runtime: `model` is typed against
the `ResonanceModel` protocol under `TYPE_CHECKING` only, exactly like
`driven.py`.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

from qscat.dvr import FemDvrEcsGrid, TensorGrid, eigen, kinetic
from qscat.linalg import c_product

if TYPE_CHECKING:
    from qscat.model import ResonanceModel

__all__ = ["anion_electronic_states", "v_dr_diag"]  # Task 4 appends "da_cross_section"

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
    return np.asarray(
        model.interaction_diag(tgrid) + v0_term - vint_inf, dtype=np.complex128
    )
