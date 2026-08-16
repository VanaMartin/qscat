"""Machinery-soundness tests for the nuclear-coordinate density comparison
(sub-project #6, Task 6 -- EXPLORATORY).

Per the task brief: this is NOT a pass/fail on the LCP-vs-exact physics --
there is no assertion anywhere here on the SIZE of that difference. The
tests only check that the machinery is sound: the projected density is
real, non-negative, integrates to a finite value, is restricted to the
unscaled region, and peaks where the molecule actually is (R ~ 1.5-3 bohr).
The actual comparison (centroids, widths, shapes) is reported as data, not
gated -- see the development notes.
"""

from __future__ import annotations

import numpy as np
import pytest
from qscat.dvr import TensorGrid

from projects.n2_2d_cross_section.cross_section_2d import ve_cross_section_2d
from projects.n2_2d_cross_section.electronic_grid import n2_electronic_grid
from projects.n2_2d_cross_section.hamiltonian2d import MU
from projects.n2_2d_cross_section.nuclear_density import compare_to_lcp, nuclear_density
from projects.n2_ti_cross_section.nuclear_grid import n2_nuclear_grid
from projects.n2_ti_cross_section.vibrational import vibrational_states

# Small but physically sane grid, matching `test_cross_section_2d.py`'s TG --
# fast enough for direct `nuclear_density` unit tests; convergence itself is
# Task 4's concern, not this file's.
TG = TensorGrid(
    [
        n2_electronic_grid(r_max=16.0, order=7, n_complex=5),
        n2_nuclear_grid(quadrature=10, r_max=22.0, n_complex=5),
    ]
)
EPS, CHI = vibrational_states(TG.grids[1], MU, 4)


@pytest.fixture(scope="module")
def psi_plus_e02():
    _, psi = ve_cross_section_2d(TG, EPS, CHI, 0, [0], 0.2, return_wavefunction=True)
    assert psi is not None
    return psi


def test_density_is_real_nonnegative_and_finite(psi_plus_e02):
    R, rho = nuclear_density(TG, psi_plus_e02)
    assert rho.dtype == np.float64
    assert np.all(np.isfinite(rho))
    assert np.all(rho >= 0.0)
    integral = np.trapezoid(rho, R)
    assert np.isfinite(integral)
    assert integral > 0.0


def test_density_peaks_in_physically_sensible_range(psi_plus_e02):
    """The molecule sits at R ~ 1.5-3 bohr (equilibrium R0 = 2.01943 bohr);
    the driven-solution density should be supported there, not at some
    numerical artifact elsewhere on the grid."""
    R, rho = nuclear_density(TG, psi_plus_e02)
    r_peak = R[np.argmax(rho)]
    assert 1.5 <= r_peak <= 3.0


def test_density_is_restricted_to_the_unscaled_nuclear_region(psi_plus_e02):
    R, _rho = nuclear_density(TG, psi_plus_e02)
    g_R = TG.grids[1]
    assert np.all(R <= g_R.R0)
    # Exactly the unscaled nuclear grid points are returned -- no more, no
    # fewer -- and none of the complex-tail points leak in.
    assert R.size == int(np.sum(g_R.real_points <= g_R.R0))
    assert R.size < g_R.n  # the tail is non-empty on this grid, so this is a real check


def test_electronic_complex_tail_does_not_contribute_to_density():
    """Amplitude placed ONLY on the electronic complex-scaled tail (r > R0)
    must not appear in rho(R): that region carries outgoing flux, not
    probability density, and `nuclear_density` must mask it out before
    summing over the electronic index."""
    g_r, _g_R = TG.grids
    r_tail_mask = g_r.real_points > g_r.R0
    assert np.any(r_tail_mask)  # sanity: the electronic tail is non-empty here

    psi = np.zeros(TG.size, dtype=np.complex128)
    psi2d = psi.reshape(TG.shape)
    psi2d[r_tail_mask, :] = 1e6  # huge amplitude, entirely on the masked region
    _R, rho = nuclear_density(TG, psi)

    assert np.all(rho == 0.0)


def test_nuclear_complex_tail_is_excluded_from_the_returned_grid():
    """Amplitude on the nuclear complex-scaled tail (R > R0) must not show
    up as a returned R point at all (it is masked out AFTER the r-sum, per
    `nuclear_density`'s docstring), even though the electronic index there
    is perfectly ordinary (real, unscaled r)."""
    g_r, g_R = TG.grids
    r_real_mask = g_r.real_points <= g_r.R0
    R_tail_mask = g_R.real_points > g_R.R0
    assert np.any(R_tail_mask)

    psi = np.zeros(TG.size, dtype=np.complex128)
    psi2d = psi.reshape(TG.shape)
    # Nonzero only at (real r, complex-tail R): a "density" that exists
    # only in the region this function must discard.
    psi2d[np.ix_(r_real_mask, R_tail_mask)] = 1.0
    R, rho = nuclear_density(TG, psi)

    assert R.size == int(np.sum(g_R.real_points <= g_R.R0))
    assert np.all(rho == 0.0)


# --- Step 3: the resonance-anchor comparison itself (reported as data, no
# pass/fail on the LCP-vs-exact difference) ---------------------------------


@pytest.fixture(scope="module")
def comparison_at_resonance():
    return compare_to_lcp(0.2, v_init=0)


def test_compare_to_lcp_returns_sound_comparable_densities(comparison_at_resonance):
    d = comparison_at_resonance
    for key in ("R_exact", "density_exact", "R_lcp", "density_lcp"):
        assert key in d

    for R_key, rho_key in (("R_exact", "density_exact"), ("R_lcp", "density_lcp")):
        R = np.asarray(d[R_key])
        rho = np.asarray(d[rho_key])
        assert np.all(np.isfinite(rho))
        assert np.all(rho >= 0.0)
        # Normalized to unit area, per the task brief (not the same object
        # dimensionally otherwise).
        assert np.trapezoid(rho, R) == pytest.approx(1.0, rel=1e-6)

    for key in ("centroid_exact", "width_exact", "centroid_lcp", "width_lcp"):
        assert np.isfinite(d[key])
        assert d[key] > 0.0

    # Both centroids fall in the physically sensible bound-state region.
    assert 1.5 <= d["centroid_exact"] <= 3.0
    assert 1.5 <= d["centroid_lcp"] <= 3.0
