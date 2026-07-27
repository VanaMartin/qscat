"""Tannor-Weeks eta deconvolution factors and the outgoing test function.

The correlation function `c_{v'}(t) = <Phi_{v'}|Psi(t)>` (Task 3's
`PropagationResult.c`) carries the spectral content of BOTH the incident
wavepacket `g(r)` and the outgoing test-function wavepacket `g_out(r)` --
neither of which is a pure energy eigenstate. The `eta` factors here
deconvolve that spectral content, leaving the pure single-energy S-matrix
element, per `.superpowers/sdd/n2-2d-exact-extraction.md` section 5.3
(`eMoScat TestFunction2d.cpp:298-307`):

- `eta_incident(E) = c_product(g_in_coeffs, F_{E,l}_coeffs)`: the incident
  Gaussian against the energy-normalized REGULAR free function
  `riccati_bessel_en` -- the SAME function, and the SAME `sqrt(w_r)`
  coefficient conversion, that
  `projects.n2_2d_cross_section.cross_section_2d.channel_vector` uses to
  build its exact TI incident channel function.
- `eta_outgoing(E') = c_product(g_out_coeffs, F^out_{E',l}_coeffs)`: the
  outgoing test-function Gaussian against the energy-normalized OUTGOING
  HALF of the free function, `h^{(1)}_{E',l}/2` (`j_l = (h_l^{(1)} +
  h_l^{(2)})/2`) -- per the extraction doc's `TestFunction2d.cpp:207`
  (`sphHankel1En(...)/2.0`), and settled empirically here (debug order item
  7): the REGULAR function for `F_out` gave `sigma_TD` five to six orders
  of magnitude too small against the TI oracle; the outgoing Hankel half
  brought it to within ~10-25% (see `.superpowers/sdd/task-4-report.md`).

`outgoing_channel` builds the energy-INDEPENDENT test function
`Phi_{v'} = g_out(r) chi_{v'}(R)`: the k'-dependence lives entirely in
`eta_outgoing`, evaluated once per (E, v') pair in `td_cross_section.py`,
while the (expensive) propagation against `Phi_{v'}` happens only once.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from qscat.dvr import FemDvrEcsGrid, TensorGrid
from qscat.linalg import c_product
from qscat.special import riccati_bessel_en, riccati_hankel_en

from .wavepacket import gaussian_coeffs

__all__ = ["outgoing_channel", "eta_incident", "eta_outgoing"]


def outgoing_channel(
    tgrid: TensorGrid,
    chi_v: npt.NDArray[np.complex128],
    *,
    r0_out: float,
    p0_out: float,
    sigma_out: float,
) -> npt.NDArray[np.complex128]:
    """`Phi_{v'} = g_out(r) chi_{v'}(R)`, masked, energy-independent.

    `chi_v` is already a c-product-normalized DVR coefficient vector (see
    `vibrational_states`'s docstring); no rescaling is applied here.
    """
    g_out_coeff = gaussian_coeffs(tgrid.grids[0], r0=r0_out, p0=p0_out, sigma=sigma_out)
    chi = np.asarray(chi_v, dtype=np.complex128)
    psi = tgrid.outer([g_out_coeff, chi])
    psi[~tgrid.real_mask()] = 0.0
    return psi


def _regular_coeffs(grid: FemDvrEcsGrid, k: float, l: int) -> npt.NDArray[np.complex128]:
    """`riccati_bessel_en(r, k, l) * sqrt(w_r)`, masked to the unscaled region.

    The SAME conversion `channel_vector` applies -- `F` is a function, so it
    picks up `sqrt(w_r)` to become a DVR coefficient vector; `chi_v` is
    already one and is not touched here.
    """
    f_vals = riccati_bessel_en(grid.real_points, k, l)
    sqrt_w = np.sqrt(np.asarray(grid.weights, dtype=np.complex128))
    coeffs = (f_vals * sqrt_w).astype(np.complex128)
    coeffs[grid.real_points > grid.R0] = 0.0
    return coeffs


def _outgoing_coeffs(grid: FemDvrEcsGrid, k: float, l: int) -> npt.NDArray[np.complex128]:
    """`h^{(1)}_{E,l}(r)/2 * sqrt(w_r)`, masked to the unscaled region.

    `h^{(1)}_{E,l}(r) = sqrt(2k/pi) r h_l^{(1)}(kr)` (energy-normalized,
    mass 1), `h_l^{(1)} = j_l + i*y_l`, halved -- see module docstring for
    why this (not the regular function) is `F_out`.
    """
    r = grid.real_points
    riccati_h1 = riccati_hankel_en(r, k, l)
    f_vals = riccati_h1 / 2.0
    sqrt_w = np.sqrt(np.asarray(grid.weights, dtype=np.complex128))
    coeffs = (f_vals * sqrt_w).astype(np.complex128)
    coeffs[grid.real_points > grid.R0] = 0.0
    return coeffs


def eta_incident(
    grid: FemDvrEcsGrid, k: float, l: int, *, r0: float, p0: float, sigma: float
) -> complex:
    """`eta_in(E) = c_product(g_in_coeffs, F_{E,l}_coeffs)` on the electronic grid."""
    g_coeff = gaussian_coeffs(grid, r0=r0, p0=p0, sigma=sigma)
    f_coeff = _regular_coeffs(grid, k, l)
    return c_product(g_coeff, f_coeff)


def eta_outgoing(
    grid: FemDvrEcsGrid, kp: float, l: int, *, r0_out: float, p0_out: float, sigma_out: float
) -> complex:
    """`eta_out(E') = c_product(g_out_coeffs, F^out_{E',l}_coeffs)` on the electronic grid.

    `F^out` is the outgoing Hankel half, NOT the regular function -- see
    module docstring.
    """
    g_coeff = gaussian_coeffs(grid, r0=r0_out, p0=p0_out, sigma=sigma_out)
    f_coeff = _outgoing_coeffs(grid, kp, l)
    return c_product(g_coeff, f_coeff)
