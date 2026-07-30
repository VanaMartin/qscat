"""Unit tests for `correlation.hankel_point_value` (sub-project #C4, Task 2)
and `correlation.outgoing_surface_wave` (Task 3): the scalar outgoing-
Hankel-half VALUE (and, for the latter, its spatial derivative) at one
physical electronic coordinate -- the `Dirac`/`Flux` extractors'
deconvolution factors.
"""

from __future__ import annotations

import numpy as np
from qscat.core.correlation import hankel_point_value, outgoing_surface_wave
from qscat.core.grids import electronic_grid
from qscat.special import coulomb_h1_en, riccati_hankel_en

GRID = electronic_grid(r_max=12.0, order=5, n_complex=3)


def test_hankel_point_value_matches_riccati_hankel_en_neutral() -> None:
    """`charge=0`: matches `riccati_hankel_en(z_position, k, l)/2` exactly."""
    position = 43  # r = 10.0, in the real (unscaled) region, past the interaction
    z_position = float(GRID.real_points[position])
    k, l = 0.5, 2
    got = hankel_point_value(GRID, z_position, k, l, 0)
    want = complex(riccati_hankel_en(np.asarray(z_position), k, l) / 2.0)
    assert got == want


def test_hankel_point_value_matches_coulomb_h1_en_charged() -> None:
    """`charge!=0`: matches `coulomb_h1_en(z_position, k, charge, 1.0, l)/2`."""
    position = 43
    z_position = float(GRID.real_points[position])
    k, l, charge = 0.5, 1, -1
    got = hankel_point_value(GRID, z_position, k, l, charge)
    want = complex(
        coulomb_h1_en(np.asarray(z_position, dtype=complex), k, float(charge), 1.0, l) / 2.0
    )
    assert got == want


def test_hankel_point_value_is_a_plain_scalar() -> None:
    """No `sqrt(w)` factor, no masking -- just the function's value."""
    position = 43
    z_position = float(GRID.real_points[position])
    got = hankel_point_value(GRID, z_position, 0.3, 2, 0)
    assert isinstance(got, complex)
    assert got != 0.0


# --- Task 3: outgoing_surface_wave --------------------------------------------


def test_outgoing_surface_wave_phi_matches_hankel_point_value_neutral() -> None:
    """`phi_out` must be exactly `hankel_point_value`'s VALUE (same definition)."""
    position = 43  # r = 10.0, real region
    z = float(GRID.real_points[position])
    k, l = 0.5, 2
    phi, _ = outgoing_surface_wave(GRID, z, k, l, 0)
    want = hankel_point_value(GRID, z, k, l, 0)
    assert phi == want


def test_outgoing_surface_wave_dphi_matches_finite_difference_neutral() -> None:
    """The analytic `dphi_out` (product rule + scipy `spherical_jn`/`yn`'s
    `derivative=True`) checked against an INDEPENDENT central finite
    difference of `riccati_hankel_en` itself -- confirms the analytic
    formula (module docstring)."""
    z, k, l = 10.0, 0.5, 2
    _, dphi = outgoing_surface_wave(GRID, z, k, l, 0)
    h = 1e-6
    fd = complex(
        riccati_hankel_en(np.asarray(z + h), k, l) - riccati_hankel_en(np.asarray(z - h), k, l)
    ) / (2.0 * h) / 2.0  # /2 for the "outgoing half" convention
    np.testing.assert_allclose(dphi, fd, rtol=1e-6)


def test_outgoing_surface_wave_dphi_several_l_and_k() -> None:
    """Sanity sweep: analytic dphi_out matches an independent FD at several
    (k, l, z) combinations, not just one lucky point."""
    for k, l, z in [(0.3, 0, 6.0), (0.7, 1, 8.0), (1.2, 3, 15.0), (0.5, 5, 20.0)]:
        _, dphi = outgoing_surface_wave(GRID, z, k, l, 0)
        h = 1e-6
        fd = complex(
            riccati_hankel_en(np.asarray(z + h), k, l) - riccati_hankel_en(np.asarray(z - h), k, l)
        ) / (2.0 * h) / 2.0
        np.testing.assert_allclose(dphi, fd, rtol=1e-6, atol=1e-10)


def test_outgoing_surface_wave_charged_branch_is_finite_and_matches_value() -> None:
    """Charged (Coulomb) branch: kept structurally (not gated for N2), but
    check it returns finite values and that `phi_out` matches
    `hankel_point_value`'s charged branch (same underlying value)."""
    z, k, l, charge = 10.0, 0.5, 1, -1
    phi, dphi = outgoing_surface_wave(GRID, z, k, l, charge)
    assert np.isfinite(phi)
    assert np.isfinite(dphi)
    want = hankel_point_value(GRID, z, k, l, charge)
    assert phi == want
