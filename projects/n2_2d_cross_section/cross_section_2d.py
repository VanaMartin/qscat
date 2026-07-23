"""Exact 2-D VE cross section by the driven Lippmann-Schwinger equation.

    Psi_i    = F_{E,l}(r) chi_v(R)                     [masked to unscaled region]
    Psi_sc   = (E_tot - H_2D)^{-1} V_int Psi_i         [one sparse LU per energy]
    Psi^(+)  = Psi_i + Psi_sc
    T_{v->v'} = <chi_v' F_{E',l} | V_int | Psi^(+)>    [c-product, masked]
    sigma     = 4 pi^3 |T|^2 / k^2                     [bohr^2]

Conventions that must not be gotten wrong (each has bitten this repo or the
reference implementation):

- `chi` from `vibrational_states` is ALREADY a DVR coefficient vector; `F` is
  a FUNCTION and must be converted with `c_j = F(r_j) sqrt(w_j)` using the
  bridge-summed complex weight (`TensorGrid.sqrt_weights()`). Mixing the two
  rescales every cross section.
- Everything is paired with the C-PRODUCT (no conjugate): under ECS `H = H^T`,
  not `H^dagger`. With both sides in coefficient form the c-product IS the
  quadrature integral.
- `Psi_i` and every `Phi_f` are masked to the unscaled region. eMoScat uses a
  Hermitian dot here and is saved only by doing the same masking.

Elastic and inelastic share one formula: with `S = 1 - 2 pi i T`,
`|S - 1|^2 = 4 pi^2 |T|^2`, so Houfek's `pi |S-1|^2 / k^2` and our
`4 pi^3 |T|^2 / k^2` are the same expression. Unlike the 1-D LCP model, this
elastic T-matrix DOES contain the non-resonant background scattering.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp
from qscat.dvr import TensorGrid
from qscat.linalg import SparseLU, c_product

from projects.n2_2d_cross_section.channels import riccati_bessel_en
from projects.n2_2d_cross_section.hamiltonian2d import ELL, build_h2d, interaction_diag

__all__ = ["channel_vector", "ve_cross_section_2d"]


def channel_vector(
    tgrid: TensorGrid,
    k: float,
    chi_v: npt.NDArray[np.complex128],
    *,
    l: int = ELL,
) -> npt.NDArray[np.complex128]:
    """DVR coefficients of `F_{E,l}(r) chi_v(R)`, masked to the unscaled region.

    `chi_v` is already a coefficient vector; `F` is a function and picks up
    `sqrt(w_r)`.
    """
    g_r = tgrid.grids[0]
    f_vals = riccati_bessel_en(g_r.real_points, k, l)
    # sqrt_weights() is per-axis and broadcast-shaped ((n_r, 1) at D=2); ravel
    # it to pair elementwise with the 1-D electronic function values.
    sqrt_w_r = tgrid.sqrt_weights()[0].ravel()
    f_coeff = f_vals * sqrt_w_r

    chi = np.asarray(chi_v, dtype=np.complex128)
    chi = chi / np.sqrt(c_product(chi, chi))  # c-product normalization, not Hermitian

    psi = tgrid.outer([f_coeff, chi])
    psi[~tgrid.real_mask()] = 0.0
    return psi


def _sigma_at_one_energy(
    tgrid: TensorGrid,
    lu: SparseLU,
    v_diag: npt.NDArray[np.complex128],
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    vprimes: list[int],
    E: float,
    *,
    want_psi: bool,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.complex128] | None]:
    sigma = np.zeros(len(vprimes), dtype=np.float64)
    if E <= 0.0:
        return sigma, None

    e_tot = E + eps[v_init]
    k = float(np.sqrt(2.0 * E))

    psi_i = channel_vector(tgrid, k, chi[v_init])
    psi_plus = psi_i + lu.solve(v_diag * psi_i)
    v_psi = v_diag * psi_plus

    for j, vp in enumerate(vprimes):
        excess = e_tot - eps[vp]
        if excess <= 0.0:
            continue  # closed channel
        kp = float(np.sqrt(2.0 * excess))
        phi_f = channel_vector(tgrid, kp, chi[vp])
        t = c_product(phi_f, v_psi)
        sigma[j] = 4.0 * np.pi**3 * abs(t) ** 2 / (2.0 * E)

    return sigma, (psi_plus if want_psi else None)


def ve_cross_section_2d(
    tgrid: TensorGrid,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    vprimes: list[int],
    E: float | npt.ArrayLike,
    *,
    ordering: str = "COLAMD",
    lam_scale: float = 1.0,
    return_wavefunction: bool = False,
):
    """sigma_{v_init->v'}(E) in bohr^2, exact 2-D driven-equation solution.

    `E` may be scalar or an array; scalar returns shape `(len(vprimes),)`,
    array returns `(len(E), len(vprimes))`. One sparse LU per energy is
    reused across all `vprimes`.

    `lam_scale` scales `V_int` ONLY, for the free-particle and first-Born
    validation limits. It is a test lever, never a physics knob.
    """
    e_arr = np.atleast_1d(np.asarray(E, dtype=np.float64))
    H = build_h2d(tgrid)
    v_diag = lam_scale * interaction_diag(tgrid)
    ident = sp.identity(tgrid.size, format="csc", dtype=np.complex128)

    out = []
    psis = []
    for e in e_arr:
        e_tot = float(e) + eps[v_init]
        lu = SparseLU((e_tot * ident - H).tocsc(), ordering=ordering)
        s, psi = _sigma_at_one_energy(
            tgrid, lu, v_diag, eps, chi, v_init, vprimes, float(e),
            want_psi=return_wavefunction,
        )
        out.append(s)
        psis.append(psi)

    sigma = np.stack(out)
    scalar = np.isscalar(E) or (isinstance(E, np.ndarray) and np.ndim(E) == 0)
    result = np.asarray(sigma[0], dtype=np.float64) if scalar else sigma
    if return_wavefunction:
        return result, (psis[0] if scalar else psis)
    return result
