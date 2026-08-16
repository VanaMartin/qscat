"""Tannor-Weeks eta deconvolution factors and the outgoing test function.

The method is Tannor & Weeks, J. Chem. Phys. 98, 3884 (1993); the `eta`
factors are its Eq. (39) normalization, which is explicitly NOT unity -- see
`docs/physics/td-extractors.md`.

The correlation function `c_{v'}(t) = <Phi_{v'}|Psi(t)>` carries the spectral
content of BOTH the incident wavepacket `g(r)` and the outgoing test-function
wavepacket `g_out(r)` -- neither of which is a pure energy eigenstate. The
`eta` factors here deconvolve that spectral content, leaving the pure
single-energy S-matrix element (`eMoScat TestFunction2d.cpp:298-307`):

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
  (`sphHankel1En(...)/2.0`), and settled empirically here: the REGULAR
  function for `F_out` gave `sigma_TD` five to six orders of magnitude too
  small against the TI oracle; the outgoing Hankel half brought it to within
  ~10-25%.

`outgoing_channel` builds the energy-INDEPENDENT test function
`Phi_{v'} = g_out(r) chi_{v'}(R)`: the k'-dependence lives entirely in
`eta_outgoing`, evaluated once per (E, v') pair in `td_cross_section.py`,
while the (expensive) propagation against `Phi_{v'}` happens only once.

`outgoing_channel_nuclear`/nuclear `eta_outgoing` (sub-project #4/SP2, Task 4)
are the NUCLEAR-axis transpose of the above, for the `TannorWeeks(axis=
"nuclear")` dissociative-attachment (DA) extractor (`qscat.core.
td_extractors`): the outgoing Gaussian test packet moves to the NUCLEAR
coordinate `R` (mass `mu_R = model.mu`, `l=0`), projected against one of the
anion ELECTRONIC bound states `phi_c(r)` instead of a nuclear vibrational
level -- `Phi_c = phi_c(r) g_out(R)`, the direct transpose of `Phi_{v'} =
g_out(r) chi_{v'}(R)`. `eta_outgoing` gains a `mass` keyword (default `1.0`,
reproducing the electronic path bit-for-bit via `riccati_hankel_en_mass(...,
1.0) == riccati_hankel_en(...)`, see `qscat.special.radial.
riccati_hankel_en_mass`'s docstring) so the SAME function serves both axes:
electronic callers are untouched, a nuclear caller passes `mass=model.mu`.
`eta_incident` (the incident ELECTRON) is NOT generalized -- the incident
side always stays on the electronic axis, even for the nuclear DA extractor
(`td_extractors.py`'s module docstring, `Flux`/`Dirac(axis="nuclear")`
sections).
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from scipy.special import spherical_jn, spherical_yn

from qscat.dvr import FemDvrEcsGrid, TensorGrid
from qscat.linalg import c_product
from qscat.special import (
    coulomb_h1_en,
    riccati_bessel_en,
    riccati_hankel_en_mass,
)

from .wavepacket import gaussian_coeffs

__all__ = [
    "outgoing_channel",
    "outgoing_channel_nuclear",
    "eta_incident",
    "eta_outgoing",
    "hankel_point_value",
    "outgoing_surface_wave",
]


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


def outgoing_channel_nuclear(
    tgrid: TensorGrid,
    phi_c: npt.NDArray[np.complex128],
    *,
    r0_out: float,
    p0_out: float,
    sigma_out: float,
) -> npt.NDArray[np.complex128]:
    """`Phi_c = phi_c(r) g_out(R)`, masked, energy-independent.

    The NUCLEAR-axis transpose of `outgoing_channel`: the outgoing test
    packet `g_out` sits in the nuclear coordinate `R` (`tgrid.grids[1]`)
    while `phi_c` -- one of the anion electronic bound states
    (`qscat.core.dissociation.anion_electronic_states`) -- sits in the
    electronic coordinate `r` (`tgrid.grids[0]`). `phi_c` is already
    c-product-normalized (that function's docstring); no rescaling is
    applied here, mirroring `outgoing_channel`'s treatment of `chi_v`.
    """
    phi = np.asarray(phi_c, dtype=np.complex128)
    g_out_coeff = gaussian_coeffs(tgrid.grids[1], r0=r0_out, p0=p0_out, sigma=sigma_out)
    psi = tgrid.outer([phi, g_out_coeff])
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


def _outgoing_coeffs(
    grid: FemDvrEcsGrid, k: float, l: int, *, mass: float = 1.0
) -> npt.NDArray[np.complex128]:
    """`h^{(1)}_{E,l}(r)/2 * sqrt(w_r)`, masked to the unscaled region.

    `h^{(1)}_{E,l}(r) = sqrt(2*mass*k/pi) r h_l^{(1)}(kr)` (energy-normalized,
    mass `mass`), `h_l^{(1)} = j_l + i*y_l`, halved -- see module docstring
    for why this (not the regular function) is `F_out`. `mass` defaults to
    `1.0` -- `riccati_hankel_en_mass(..., 1.0)` reproduces `riccati_hankel_en`
    bit-for-bit (`qscat.special.radial.riccati_hankel_en_mass`'s docstring),
    so every existing (electronic) call site is untouched; a nuclear (DA)
    caller passes `mass=model.mu`.
    """
    r = grid.real_points
    riccati_h1 = riccati_hankel_en_mass(r, k, l, mass)
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
    grid: FemDvrEcsGrid,
    kp: float,
    l: int,
    *,
    r0_out: float,
    p0_out: float,
    sigma_out: float,
    mass: float = 1.0,
) -> complex:
    """`eta_out(E') = c_product(g_out_coeffs, F^out_{E',l}_coeffs)` on `grid`.

    `F^out` is the outgoing Hankel half, NOT the regular function -- see
    module docstring. `mass` defaults to `1.0` (electronic, byte-identical to
    the pre-`mass` code -- see `_outgoing_coeffs`'s docstring); a nuclear DA
    caller (`grid = tgrid.grids[1]`) passes `mass=model.mu`, `l=0`.
    """
    g_coeff = gaussian_coeffs(grid, r0=r0_out, p0=p0_out, sigma=sigma_out)
    f_coeff = _outgoing_coeffs(grid, kp, l, mass=mass)
    return c_product(g_coeff, f_coeff)


def hankel_point_value(
    grid: FemDvrEcsGrid, z_position: float, k: float, l: int, charge: int = 0, *, mass: float = 1.0
) -> complex:
    """`H^{(1)}_{E,l}(z_position)/2` -- the outgoing-Hankel-half VALUE at a
    single physical (real, unscaled) coordinate, e.g. `z_position =
    grid.real_points[position]` for some fixed DVR index `position`
    (`Dirac`'s analysis point, `td_extractors.py`).

    The scalar sibling of `_outgoing_coeffs`: same energy-normalized outgoing
    function -- `riccati_hankel_en_mass(z_position, k, l, mass)/2` (neutral,
    `charge == 0`) or `coulomb_h1_en(z_position, k, charge, mass, l)/2`
    (charged target) -- but evaluated at ONE point rather than converted to a
    `sqrt(w_r)`-scaled, masked DVR coefficient VECTOR: a delta-distribution
    test function needs the outgoing function's VALUE, not an integral
    against it. `grid` is accepted (unused) to keep this call-compatible
    with `_regular_coeffs`/`_outgoing_coeffs` and make "a value on THIS
    grid's real axis" explicit at call sites.

    `mass` defaults to `1.0` (the electronic reduced mass, a.u.) -- every
    existing (electronic) call site is untouched: `riccati_hankel_en_mass(
    ..., 1.0)` reproduces `riccati_hankel_en(...)` bit-for-bit (`2.0*1.0 ==
    2.0` exactly), and `coulomb_h1_en(..., 1.0, l)` is the same literal `1.0`
    the pre-`mass` code passed. A nuclear (dissociation) caller passes
    `mass=model.mu`.
    """
    del grid  # unused: kept for call-site symmetry with _regular_coeffs/_outgoing_coeffs
    if charge == 0:
        val = riccati_hankel_en_mass(np.asarray(z_position, dtype=np.float64), k, l, mass) / 2.0
    else:
        val = (
            coulomb_h1_en(np.asarray(z_position, dtype=np.complex128), k, float(charge), mass, l)
            / 2.0
        )
    return complex(np.asarray(val))


def outgoing_surface_wave(
    grid: FemDvrEcsGrid,
    z_surface: float,
    k: float,
    l: int,
    charge: float = 0.0,
    *,
    mass: float = 1.0,
) -> tuple[complex, complex]:
    """`(phi_out, dphi_out) = H^{(1)}_{E,l}(z_surface)/2` and its SPATIAL
    derivative at `z_surface` -- the `Flux` extractor's per-channel outgoing
    wave + derivative (eMoScat `FluxTestFunction2d`'s `phi_out_`/`dphi_out_`,
    confirmed by port-scout reading `FluxTestFunction2d.cpp`'s constructor:
    it samples `sphHankel1En(...)/2` -- or the Coulomb `sH1_en(...)/2` for a
    charged target -- over an element's nodes and applies the SAME DVR
    derivative `GridVector::derivative` uses; that machinery is reproduced
    here directly against the ANALYTIC function instead, see below).

    `mass` (eMoScat's `reduced_mass()`, `mu_x_=1.0` electronic / `mu_y_=mass`
    nuclear) is the mass entering the energy normalization of the outgoing
    function -- `F^{(1)}_{E,l}(r) = sqrt(2 mass k/pi) r h_l^{(1)}(kr)`
    (`riccati_hankel_en_mass`'s definition; `mass` does NOT enter the
    momentum argument `kr`, matching `riccati_bessel_en_mass`'s convention).
    Defaults to `1.0` -- every existing (electronic) call site is untouched:
    at `mass=1.0` the formulas below reproduce the pre-`mass` code bit-for-
    bit (`2.0*1.0 == 2.0` exactly). A nuclear (dissociation) caller passes
    `mass=model.mu`.

    Neutral (`charge == 0`): computed ANALYTICALLY via the product rule,

        dF/dr = sqrt(2 mass k/pi) * [h_l(kr) + kr * h_l'(kr)]

    using `scipy.special.spherical_jn`/`spherical_yn`'s `derivative=True`
    option for `h_l^{(1)}{}'(x) = j_l'(x) + i y_l'(x)` -- `qscat.special.
    radial` does not itself expose a derivative primitive (only the two
    VALUE functions), but the underlying scipy pieces it is built on already
    support one, so no finite difference is needed for this branch (checked
    against a finite difference of `riccati_hankel_en`/`riccati_hankel_en_mass`
    in `test_correlation.py`, confirming the analytic formula).

    Charged (`charge != 0`, `coulomb_h1_en`): `qscat.special.coulomb` has no
    derivative primitive for the Coulomb functions (mpmath's `coulombf`/
    `coulombg` expose no `derivative=` option, and the F_l'/G_l' recurrence
    needs extra pieces this module does not carry) -- this branch falls back
    to a high-accuracy CENTRAL finite difference (4th-order, 5-point
    stencil) of `coulomb_h1_en(..., mass, l)` itself. Kept structurally for a
    charged target (e.g. H2+); N2/F2 are neutral, so only the analytic
    branch is exercised by this sub-project's gate.
    """
    del grid  # unused: kept for call-site symmetry with hankel_point_value
    r = float(z_surface)
    if charge == 0:
        pref = np.sqrt(2.0 * mass * k / np.pi)
        x = k * r
        h_l = spherical_jn(l, x) + 1j * spherical_yn(l, x)
        h_l_prime = spherical_jn(l, x, derivative=True) + 1j * spherical_yn(
            l, x, derivative=True
        )
        phi = pref * r * h_l / 2.0
        dphi = pref * (h_l + x * h_l_prime) / 2.0
        return complex(phi), complex(dphi)

    def _h1(rr: float) -> complex:
        val = (
            coulomb_h1_en(np.asarray(rr, dtype=np.complex128), k, float(charge), mass, l) / 2.0
        )
        return complex(np.asarray(val))

    phi = _h1(r)
    h = 1e-4 * max(1.0, abs(r))
    dphi = (-_h1(r + 2 * h) + 8 * _h1(r + h) - 8 * _h1(r - h) + _h1(r - 2 * h)) / (12.0 * h)
    return phi, dphi
