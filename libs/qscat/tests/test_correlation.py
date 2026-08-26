"""Unit tests for `correlation.hankel_point_value` (sub-project #C4, Task 2)
and `correlation.outgoing_surface_wave` (Task 3): the scalar outgoing-
Hankel-half VALUE (and, for the latter, its spatial derivative) at one
physical electronic coordinate -- the `Dirac`/`Flux` extractors'
deconvolution factors.
"""

from __future__ import annotations

import numpy as np
import pytest
from qscat.core.correlation import (
    eta_outgoing,
    hankel_point_value,
    outgoing_channel_nuclear,
    outgoing_surface_wave,
)
from qscat.core.dissociation import anion_electronic_states
from qscat.core.grids import electronic_grid, nuclear_grid
from qscat.core.wavepacket import gaussian_coeffs
from qscat.dvr import TensorGrid
from qscat.linalg import c_product
from qscat.model import N2
from qscat.special import coulomb_h1_en, riccati_hankel_en, riccati_hankel_en_mass

GRID = electronic_grid(r_max=12.0, order=5, n_complex=3)


def test_hankel_point_value_matches_riccati_hankel_en_neutral() -> None:
    """`charge=0`: matches `riccati_hankel_en(z_position, k, l)/2` exactly."""
    position = 43  # r = 10.0, in the real (unscaled) region, past the interaction
    z_position = float(GRID.real_points[position])
    k, l = 0.5, 2
    got = hankel_point_value(z_position, k, l, 0)
    want = complex(riccati_hankel_en(np.asarray(z_position), k, l) / 2.0)
    assert got == want


def test_hankel_point_value_matches_coulomb_h1_en_charged() -> None:
    """`charge!=0`: matches `coulomb_h1_en(z_position, k, charge, 1.0, l)/2`."""
    position = 43
    z_position = float(GRID.real_points[position])
    k, l, charge = 0.5, 1, -1
    got = hankel_point_value(z_position, k, l, charge)
    want = complex(
        coulomb_h1_en(np.asarray(z_position, dtype=complex), k, float(charge), 1.0, l) / 2.0
    )
    assert got == want


def test_hankel_point_value_is_a_plain_scalar() -> None:
    """No `sqrt(w)` factor, no masking -- just the function's value."""
    position = 43
    z_position = float(GRID.real_points[position])
    got = hankel_point_value(z_position, 0.3, 2, 0)
    assert isinstance(got, complex)
    assert got != 0.0


# --- Task 3: outgoing_surface_wave --------------------------------------------


def test_outgoing_surface_wave_phi_matches_hankel_point_value_neutral() -> None:
    """`phi_out` must be exactly `hankel_point_value`'s VALUE (same definition)."""
    position = 43  # r = 10.0, real region
    z = float(GRID.real_points[position])
    k, l = 0.5, 2
    phi, _ = outgoing_surface_wave(z, k, l, 0)
    want = hankel_point_value(z, k, l, 0)
    assert phi == want


def test_outgoing_surface_wave_dphi_matches_finite_difference_neutral() -> None:
    """The analytic `dphi_out` (product rule + scipy `spherical_jn`/`yn`'s
    `derivative=True`) checked against an INDEPENDENT central finite
    difference of `riccati_hankel_en` itself -- confirms the analytic
    formula (module docstring)."""
    z, k, l = 10.0, 0.5, 2
    _, dphi = outgoing_surface_wave(z, k, l, 0)
    h = 1e-6
    fd = (
        complex(
            riccati_hankel_en(np.asarray(z + h), k, l) - riccati_hankel_en(np.asarray(z - h), k, l)
        )
        / (2.0 * h)
        / 2.0
    )  # /2 for the "outgoing half" convention
    np.testing.assert_allclose(dphi, fd, rtol=1e-6)


def test_outgoing_surface_wave_dphi_several_l_and_k() -> None:
    """Sanity sweep: analytic dphi_out matches an independent FD at several
    (k, l, z) combinations, not just one lucky point."""
    for k, l, z in [(0.3, 0, 6.0), (0.7, 1, 8.0), (1.2, 3, 15.0), (0.5, 5, 20.0)]:
        _, dphi = outgoing_surface_wave(z, k, l, 0)
        h = 1e-6
        fd = (
            complex(
                riccati_hankel_en(np.asarray(z + h), k, l)
                - riccati_hankel_en(np.asarray(z - h), k, l)
            )
            / (2.0 * h)
            / 2.0
        )
        np.testing.assert_allclose(dphi, fd, rtol=1e-6, atol=1e-10)


def test_outgoing_surface_wave_charged_branch_is_finite_and_matches_value() -> None:
    """Charged (Coulomb) branch: kept structurally (not gated for N2), but
    check it returns finite values and that `phi_out` matches
    `hankel_point_value`'s charged branch (same underlying value)."""
    z, k, l, charge = 10.0, 0.5, 1, -1
    phi, dphi = outgoing_surface_wave(z, k, l, charge)
    assert np.isfinite(phi)
    assert np.isfinite(dphi)
    want = hankel_point_value(z, k, l, charge)
    assert phi == want


# --- Task 2 (this task): `mass` kwarg (electronic default byte-identical +
# the nuclear-mass generalization used by `Flux(axis="nuclear")`) -----------


def test_hankel_point_value_mass_default_is_byte_identical() -> None:
    """`mass=1.0` (the default) must reproduce the pre-`mass` neutral-branch
    result bit-for-bit -- `riccati_hankel_en_mass(..., 1.0)` == `riccati_
    hankel_en(...)` exactly (`2.0*1.0 == 2.0` in IEEE754, no rounding)."""
    position = 43
    z_position = float(GRID.real_points[position])
    k, l = 0.5, 2
    got = hankel_point_value(z_position, k, l, 0)
    want = complex(riccati_hankel_en(np.asarray(z_position), k, l) / 2.0)
    assert got == want


def test_hankel_point_value_nuclear_mass_matches_riccati_hankel_en_mass() -> None:
    position = 43
    z_position = float(GRID.real_points[position])
    k, l, mu = 0.5, 0, 918.25
    got = hankel_point_value(z_position, k, l, 0, mass=mu)
    want = complex(riccati_hankel_en_mass(np.asarray(z_position), k, l, mu) / 2.0)
    assert got == want


def test_outgoing_surface_wave_mass_default_is_byte_identical() -> None:
    z, k, l = 10.0, 0.5, 2
    phi_default, dphi_default = outgoing_surface_wave(z, k, l, 0)
    phi_explicit, dphi_explicit = outgoing_surface_wave(z, k, l, 0, mass=1.0)
    assert phi_default == phi_explicit
    assert dphi_default == dphi_explicit


def test_outgoing_surface_wave_nuclear_mass_phi_matches_hankel_point_value() -> None:
    z, k, l, mu = 12.0, 0.6, 0, 918.25
    phi, _ = outgoing_surface_wave(z, k, l, 0, mass=mu)
    want = hankel_point_value(z, k, l, 0, mass=mu)
    assert phi == want


def test_outgoing_surface_wave_nuclear_mass_dphi_matches_finite_difference() -> None:
    """Same analytic-vs-FD check as the electronic (`mass=1.0`) case, at a
    nuclear-scale reduced mass."""
    z, k, l, mu = 12.0, 0.6, 0, 918.25
    _, dphi = outgoing_surface_wave(z, k, l, 0, mass=mu)
    h = 1e-6
    fd = (
        complex(
            riccati_hankel_en_mass(np.asarray(z + h), k, l, mu)
            - riccati_hankel_en_mass(np.asarray(z - h), k, l, mu)
        )
        / (2.0 * h)
        / 2.0
    )
    np.testing.assert_allclose(dphi, fd, rtol=1e-6)


# --- Task 4: `outgoing_channel_nuclear` + nuclear `eta_outgoing` (sub-project
# #4/SP2) -- the nuclear-axis transpose of `outgoing_channel`/`eta_outgoing`,
# feeding `td_extractors.TannorWeeks(axis="nuclear")`. -----------------------

TG2 = TensorGrid(
    [
        electronic_grid(r_max=12.0, order=5, n_complex=3),
        nuclear_grid(quadrature=6, r_max=14.0, n_complex=3),
    ]
)
EPS_E, PHI = anion_electronic_states(TG2.grids[0], N2, R_inf=TG2.grids[1].R0, n_states=1)
NUCLEAR_WP_OUT = {"r0_out": 7.0, "p0_out": 5.0, "sigma_out": 1.0}


def test_outgoing_channel_nuclear_shape_and_mask() -> None:
    psi = outgoing_channel_nuclear(TG2, PHI[0], **NUCLEAR_WP_OUT)
    assert psi.shape == (TG2.size,)
    assert np.all(np.isfinite(psi))
    mask = TG2.real_mask()
    assert np.all(psi[~mask] == 0.0)
    assert np.any(psi[mask] != 0.0)


def test_outgoing_channel_nuclear_matches_manual_outer_product() -> None:
    """Transpose of `outgoing_channel`: `phi_c(r)` on axis 0 (electronic),
    `g_out(R)` on axis 1 (nuclear) -- `outgoing_channel`'s `g_out(r) x
    chi_v(R)` with the two factors and axes swapped."""
    g_out = gaussian_coeffs(TG2.grids[1], r0=7.0, p0=5.0, sigma=1.0)
    want = TG2.outer([np.asarray(PHI[0], dtype=np.complex128), g_out])
    want[~TG2.real_mask()] = 0.0
    got = outgoing_channel_nuclear(TG2, PHI[0], **NUCLEAR_WP_OUT)
    np.testing.assert_allclose(got, want)


def test_eta_outgoing_mass_default_is_byte_identical() -> None:
    """`mass=1.0` (the default) must reproduce the pre-`mass` electronic
    result exactly -- `riccati_hankel_en_mass(..., 1.0)` == `riccati_hankel_
    en(...)` bit-for-bit (`_outgoing_coeffs`'s docstring)."""
    wp_out = {"r0_out": 6.0, "p0_out": 0.5, "sigma_out": 1.0}
    got = eta_outgoing(GRID, 0.5, 2, **wp_out)
    want = eta_outgoing(GRID, 0.5, 2, mass=1.0, **wp_out)
    assert got == want


def test_eta_outgoing_nuclear_mass_matches_manual_construction() -> None:
    grid = TG2.grids[1]
    k, l, mu = 0.6, 0, N2.mu
    got = eta_outgoing(grid, k, l, mass=mu, **NUCLEAR_WP_OUT)

    g_coeff = gaussian_coeffs(grid, r0=7.0, p0=5.0, sigma=1.0)
    r = grid.real_points
    f_vals = riccati_hankel_en_mass(r, k, l, mu) / 2.0
    sqrt_w = np.sqrt(np.asarray(grid.weights, dtype=np.complex128))
    f_coeff = (f_vals * sqrt_w).astype(np.complex128)
    f_coeff[r > grid.R0] = 0.0
    want = c_product(g_coeff, f_coeff)
    assert got == want


def test_eta_outgoing_nuclear_mass_is_finite_and_nonzero() -> None:
    grid = TG2.grids[1]
    k, l, mu = 0.6, 0, N2.mu
    got = eta_outgoing(grid, k, l, mass=mu, **NUCLEAR_WP_OUT)
    assert np.isfinite(got)
    assert got != 0.0


# --- lib-M16: the documented-unused `grid` parameter is dropped ------------


def test_point_value_functions_drop_grid_param() -> None:
    """lib-M16: the grid argument was documented-unused; the new signature
    drops it, the old grid-first form warns for one cycle."""
    k, l = 0.7, 0
    new = hankel_point_value(3.0, k, l)
    with pytest.warns(DeprecationWarning, match="grid"):
        old = hankel_point_value(GRID, 3.0, k, l)  # legacy call form
    assert new == old

    pv, dv = outgoing_surface_wave(3.0, k, l)
    with pytest.warns(DeprecationWarning, match="grid"):
        pv2, dv2 = outgoing_surface_wave(GRID, 3.0, k, l)
    assert (pv, dv) == (pv2, dv2)
