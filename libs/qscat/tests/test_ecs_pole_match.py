"""Unit tests for `qscat.ecs.match_angle_stable` (multi-state two-angle matcher)."""

from __future__ import annotations

import numpy as np
import pytest
from qscat.ecs import find_resonance_pole, match_angle_stable

WIDE = (-10.0, 10.0, -10.0, 10.0)


def test_keeps_states_present_in_both_spectra():
    # Three "stable" states appear in both spectra to ~1e-12; the rest are
    # "continuum" and differ between the two by O(0.1).
    stable = np.array([-1.0 - 0.01j, -0.5 - 0.02j, 0.25 - 0.05j])
    a = np.concatenate([stable, np.array([2.0 - 1.0j, 3.0 - 1.5j])])
    b = np.concatenate([stable + 1e-12, np.array([2.4 - 1.3j, 3.6 - 1.9j])])

    energies, residuals, idx = match_angle_stable(a, b, WIDE)

    assert energies.shape == (3,)
    assert np.allclose(energies, stable, atol=1e-9)
    assert np.all(residuals < 1e-9)
    # `idx` indexes into the ORIGINAL `eigs_a`, not a filtered copy.
    assert np.allclose(np.asarray(a)[idx], stable, atol=1e-9)


def test_window_excludes_out_of_range_states():
    a = np.array([-1.0 - 0.01j, 5.0 - 0.01j])
    b = a + 1e-12
    energies, _residuals, _idx = match_angle_stable(a, b, (-2.0, 0.0, -1.0, 1.0))
    assert energies.shape == (1,)
    assert np.isclose(energies[0].real, -1.0)


def test_results_are_sorted_by_ascending_real_part():
    stable = np.array([0.25 - 0.05j, -1.0 - 0.01j, -0.5 - 0.02j])
    energies, _residuals, _idx = match_angle_stable(stable, stable + 1e-12, WIDE)
    assert np.all(np.diff(energies.real) > 0)


def test_relative_tolerance_scales_with_magnitude():
    # |E| = 100, partner off by 1e-3 -> |dE|/|E| = 1e-5 < rel_tol=1e-4: accepted.
    # NOTE: WIDE = (-10, 10, -10, 10) would exclude E=100, so this test uses
    # its own wider window (a bug in the original brief -- fixed here).
    wide_100 = (-1000.0, 1000.0, -1000.0, 1000.0)
    a = np.array([100.0 + 0.0j])
    b = np.array([100.001 + 0.0j])
    assert match_angle_stable(a, b, wide_100, rel_tol=1e-4)[0].size == 1
    # Same absolute gap at |E| = 0.1 is a 1e-2 relative gap: rejected. The
    # absolute floor `atol` must not rescue it either.
    a2 = np.array([0.1 + 0.0j])
    b2 = np.array([0.101 + 0.0j])
    assert match_angle_stable(a2, b2, WIDE, rel_tol=1e-4, atol=1e-8)[0].size == 0


def test_empty_window_raises():
    a = np.array([1.0 + 0.0j])
    with pytest.raises(ValueError, match="contains no eigenvalues"):
        match_angle_stable(a, a, (5.0, 6.0, -1.0, 1.0))


def test_no_stable_pair_returns_empty_not_an_error():
    a = np.array([1.0 + 0.0j, 2.0 + 0.0j])
    b = np.array([1.5 + 0.0j, 2.5 + 0.0j])
    energies, residuals, idx = match_angle_stable(a, b, WIDE)
    assert energies.size == 0 and residuals.size == 0 and idx.size == 0


def test_find_resonance_pole_behaviour_is_unchanged():
    # The single-pole finder returns the globally closest pair WHATEVER its
    # residual -- it must NOT inherit match_angle_stable's tolerance cut.
    a = np.array([1.0 + 0.0j, 2.0 + 0.0j])
    b = np.array([1.5 + 0.0j, 2.5 + 0.0j])
    pole, residual = find_resonance_pole(a, b, WIDE)
    assert np.isclose(pole.real, 1.25) and np.isclose(residual, 0.5)
