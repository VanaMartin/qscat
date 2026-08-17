"""Tests for the NRM discrete-state choices (PRA 77 Sec. VI A and VI B)."""

from __future__ import annotations

import numpy as np
import pytest
from qscat.core.grids import electronic_grid
from qscat.core.nrm.discrete_state import (
    AsymptoticDiscreteState,
    PhysicalDiscreteState,
    truncate,
)
from qscat.exceptions import ConvergenceError
from qscat.linalg import c_product
from qscat.model import F2


@pytest.fixture(scope="module")
def grid():
    return electronic_grid(r_max=16.0, order=8, n_complex=6)


def test_truncate_suppresses_the_tail_and_keeps_the_inner_region(grid):
    """Eq. (69)'s cutoff is ~1 well inside r_d and ~0 well outside."""
    ones = np.ones(grid.n, dtype=np.complex128)
    cut = truncate(ones, grid, r_d=10.0)
    real = grid.points.imag == 0.0
    r = grid.points[real].real
    inner = cut[real][r < 5.0]
    outer = cut[real][r > 15.0]
    assert np.all(np.abs(inner) > 0.99)
    assert np.all(np.abs(outer) < 0.01)
    assert np.all(cut[~real] == 0.0)


def test_asymptotic_state_is_c_normalized_and_localized(grid):
    """Choice B: phi_b, the R->infinity bound state."""
    ds = AsymptoticDiscreteState(grid, F2, R_inf=grid.R0)
    d = ds.phi_d(2.5)
    assert abs(c_product(d, d) - 1.0) < 1e-10
    assert np.all(np.abs(d[grid.points.imag != 0.0]) < 1e-8)


def test_asymptotic_state_is_r_independent(grid):
    """Choice B is the SAME state at every R -- that is its defining property."""
    ds = AsymptoticDiscreteState(grid, F2, R_inf=grid.R0)
    assert np.array_equal(ds.phi_d(1.8), ds.phi_d(6.0))


def test_physical_state_is_c_normalized_and_localized(grid):
    """Choice A: the scattering function at Re E_res(R), truncated."""
    R = np.array([6.0, 5.0, 4.0, 3.0, 2.6, 2.2, 1.8])
    ds = PhysicalDiscreteState(grid, F2, R, elec_grid_b=_angle_b(), r_d=10.0)
    for r_val in (6.0, 2.6, 1.8):
        d = ds.phi_d(r_val)
        assert abs(c_product(d, d) - 1.0) < 1e-10
        assert np.all(np.abs(d[grid.points.imag != 0.0]) < 1e-8)


def test_physical_state_varies_with_R(grid):
    """Choice A must NOT be R-independent -- that is what distinguishes it."""
    R = np.array([6.0, 5.0, 4.0, 3.0, 2.6, 2.2, 1.8])
    ds = PhysicalDiscreteState(grid, F2, R, elec_grid_b=_angle_b(), r_d=10.0)
    assert not np.allclose(ds.phi_d(6.0), ds.phi_d(1.8), atol=1e-6)


def test_physical_state_approaches_the_bound_state_at_large_R(grid):
    """Eq. (67): phi_d(r;R) -> phi_b(r) as R -> infinity, up to a sign."""
    R = np.array([6.0, 5.0, 4.0, 3.0, 2.6, 2.2, 1.8])
    ds = PhysicalDiscreteState(grid, F2, R, elec_grid_b=_angle_b(), r_d=10.0)
    phi_b = AsymptoticDiscreteState(grid, F2, R_inf=grid.R0).phi_d(6.0)
    overlap = abs(c_product(ds.phi_d(6.0), phi_b))
    assert overlap > 0.9


def test_physical_state_rejects_ascending_R(grid):
    """The continuation walk needs DESCENDING R (it is seeded at large R)."""
    with pytest.raises(ValueError, match="descending"):
        PhysicalDiscreteState(grid, F2, np.array([1.8, 2.6, 6.0]), elec_grid_b=_angle_b())


def test_physical_state_uses_the_scattering_branch_in_the_resonance_region(grid):
    """Choice A must actually run `scattering_state`, not only the bound
    branch. F2 has a genuine resonance (Gamma > 0) below the R ~ 2.6 bohr
    crossing (see test_lcp.py's test_gamma_positive_in_resonance_region), so
    the walk should cross from bound to scattering partway through this same
    R array -- asserted on the recorded `used_scattering` flag, not inferred
    from the output state."""
    R = np.array([6.0, 5.0, 4.0, 3.0, 2.6, 2.2, 1.8])
    ds = PhysicalDiscreteState(grid, F2, R, elec_grid_b=_angle_b(), r_d=10.0)
    assert not np.any(ds.used_scattering[:5])  # R = 6.0 .. 2.6: bound
    assert np.all(ds.used_scattering[5:])  # R = 2.2, 1.8: genuine resonance


def test_physical_state_raises_when_the_walk_freezes_onto_a_bad_sign(grid):
    """A too-narrow re_half_width lets resonance_pole_walk freeze between
    R=3.0 and R=2.6 (a documented resonance_pole_walk limitation), leaving a
    stale negative shift at R=2.2/1.8 where F2's spectrum has no bound state
    at all. The bound-branch gate must raise rather than silently return
    whatever eigenvector happens to have the smallest real part."""
    R = np.array([6.0, 5.0, 4.0, 3.0, 2.6, 2.2, 1.8])
    with pytest.raises(ConvergenceError, match="no eigenvalue"):
        PhysicalDiscreteState(
            grid, F2, R, elec_grid_b=_angle_b(), re_half_width=0.05, im_half_width=0.05
        )


def _angle_b():
    """The second ECS angle the two-angle pole matcher needs."""
    return electronic_grid(r_max=16.0, order=8, n_complex=6, angle_deg=40.0)
