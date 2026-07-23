"""Energy-normalized regular free radial function and the electronic grid.

`F_{E,l}(r) = sqrt(2k/pi) * r * j_l(k r)` is the energy-normalized regular
solution of the free radial equation at electron mass 1 (eMoScat
`sphBesselJEn`, source/bessel.cpp:50).
"""

from __future__ import annotations

import numpy as np
import pytest
from qscat.dvr import FemDvrEcsGrid

from projects.n2_2d_cross_section.channels import riccati_bessel_en
from projects.n2_2d_cross_section.electronic_grid import n2_electronic_grid


def test_l0_reduces_to_exact_normalized_sine() -> None:
    """j_0(x) = sin(x)/x, so F_{E,0}(r) == sqrt(2/(pi k)) sin(k r) EXACTLY.

    This pins the normalization CONSTANT, not just the shape -- the whole
    cross-section scale rides on it.
    """
    r = np.linspace(0.05, 40.0, 2000)
    k = 0.63
    got = riccati_bessel_en(r, k, 0)
    want = np.sqrt(2.0 / (np.pi * k)) * np.sin(k * r)
    assert np.abs(got - want).max() < 1e-12


def test_asymptotic_envelope_matches_energy_normalization() -> None:
    """F -> sqrt(2/(pi k)) sin(k r - l pi/2): the envelope fixes the constant.

    j_l's asymptotic correction relative to sin(x-l pi/2)/x is
    O(l(l+1)/(2x^2)); with l=2 that needs k r ~ a few hundred (not ~100-130)
    to be safely below a 1e-4 tolerance, so r is pushed further out than the
    brief's original range to make the comment's "deep asymptotic" claim
    actually true at this tolerance (verified against the closed-form j_2).
    """
    k = 0.5
    r = np.linspace(2000.0, 2600.0, 20000)  # k r >> l(l+1), deep asymptotic
    f = riccati_bessel_en(r, k, 2)
    assert abs(np.abs(f).max() - np.sqrt(2.0 / (np.pi * k))) < 1e-4


def test_satisfies_free_radial_equation() -> None:
    """-F'' + l(l+1)/r^2 F = k^2 F, checked by finite differences."""
    k, ell = 0.7, 2
    h = 1e-4
    r = np.linspace(2.0, 12.0, 300)
    f = riccati_bessel_en(r, k, ell)
    fpp = (riccati_bessel_en(r + h, k, ell) - 2 * f + riccati_bessel_en(r - h, k, ell)) / h**2
    residual = -fpp + ell * (ell + 1) / r**2 * f - k**2 * f
    assert np.abs(residual).max() < 1e-5 * np.abs(f).max()


def test_regular_at_origin() -> None:
    """The REGULAR solution vanishes at r -> 0 like r^{l+1}."""
    assert abs(float(riccati_bessel_en(np.array([1e-6]), 0.5, 2)[0])) < 1e-15


def test_electronic_grid_shape_and_ecs_pivot() -> None:
    g = n2_electronic_grid(r_max=30.0, angle_deg=35.0, order=8, n_complex=8)
    assert isinstance(g, FemDvrEcsGrid)
    # R0 is x_min + sum(real element lengths), accumulated in floating point,
    # so compare approximately rather than exactly.
    assert g.R0 == pytest.approx(30.0)  # pivot at the end of the real region
    assert g.real_points.min() > 0.0  # Dirichlet endpoint at r=0 dropped
    # real region genuinely unscaled; tail genuinely scaled
    inside = g.real_points <= g.R0
    assert np.abs(g.points[inside].imag).max() < 1e-12
    assert np.abs(g.points[~inside].imag).max() > 1.0


def test_electronic_grid_is_parametrized() -> None:
    """Task 4's convergence study varies every one of these."""
    a = n2_electronic_grid(r_max=20.0, order=7, n_complex=6)
    b = n2_electronic_grid(r_max=45.0, order=9, n_complex=10)
    assert a.n != b.n
    # R0 is a floating-point accumulation of element lengths (see GridSpec),
    # so -- as the other test's comment already notes -- compare
    # approximately; r_max=45.0 with 9 outer elements of length 35/9
    # accumulates to 44.99999999999999, not exactly 45.0.
    assert a.R0 == pytest.approx(20.0) and b.R0 == pytest.approx(45.0)
    assert (
        n2_electronic_grid(angle_deg=25.0).points[-1]
        != n2_electronic_grid(angle_deg=40.0).points[-1]
    )
