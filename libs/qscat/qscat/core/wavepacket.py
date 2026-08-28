"""Incident Gaussian electron wavepacket and the 2-D initial state.

`g(r) = (2 pi sigma^2)^{-1/4} exp(-(r-r0)^2/(4 sigma^2)) exp(i p0 r)`
(eMoScat `input.cpp:240`), converted to FEM-DVR coefficients on the unscaled
electronic region (`c_j = g(r_j) sqrt(w_j)`, same convention as
`qscat.core.channels.channel_vector`). `p0 < 0` launches the packet inward,
toward the molecule; the ECS tail (not evaluated here) absorbs whatever
leaves during propagation.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from qscat.dvr import FemDvrEcsGrid, TensorGrid
from qscat.linalg import c_product

__all__ = ["gaussian_coeffs", "initial_state"]


def gaussian_coeffs(
    grid: FemDvrEcsGrid, *, r0: float, p0: float, sigma: float
) -> npt.NDArray[np.complex128]:
    """DVR coefficients of `g(r)` on `grid`, zero on the ECS tail."""
    r = grid.real_points
    envelope = (2.0 * np.pi * sigma**2) ** -0.25 * np.exp(-((r - r0) ** 2) / (4.0 * sigma**2))
    g_vals = envelope * np.exp(1j * p0 * r)
    sqrt_w = np.sqrt(np.asarray(grid.weights, dtype=np.complex128))
    coeffs = (g_vals * sqrt_w).astype(np.complex128)
    coeffs[r > grid.R0] = 0.0  # unscaled region only
    return np.asarray(coeffs, dtype=np.complex128)


def initial_state(
    tgrid: TensorGrid,
    chi_v: npt.NDArray[np.complex128],
    *,
    r0: float,
    p0: float,
    sigma: float,
) -> npt.NDArray[np.complex128]:
    """`Psi(0) = g(r) chi_v(R)`, flat, masked, renormalized to unit Hermitian L2 norm.

    Renormalization uses `np.linalg.norm` (the true `sqrt(sum |psi_j|^2)`
    probability norm), NOT the c-product self-pairing: for a wavepacket
    carrying a momentum phase `exp(i p0 r)`, `c_product(g, g)` is a small
    oscillatory complex number, not a norm, and using it here would silently
    produce a state far from unit probability. The c-product is reserved for
    the chi_v self-pairing below (its ECS-basis normalization convention,
    per `vibrational_states`' docstring) and, downstream, for correlation
    functions and the S-matrix -- never for this state's overall scale.
    """
    g_coeff = gaussian_coeffs(tgrid.grids[0], r0=r0, p0=p0, sigma=sigma)
    chi = np.asarray(chi_v, dtype=np.complex128)
    chi = chi / np.sqrt(c_product(chi, chi))
    psi = tgrid.outer([g_coeff, chi])
    psi[~tgrid.real_mask()] = 0.0
    # Hermitian L2 (probability) norm -- see docstring above.
    norm = float(np.linalg.norm(psi))
    return np.asarray(psi / norm, dtype=np.complex128)
