"""Region-split magnitude scaling for domain colouring.

The defect being fixed: `complex_to_hsv` normalises the whole field by one
scalar `mag`, so in a time-dependent wavefunction -- where the incident packet
outweighs the resonant and outgoing amplitude by orders of magnitude -- every
interesting feature renders black. See docs/physics/ for the physics; this
module is pure numpy.
"""

from __future__ import annotations

import numpy as np
import pytest
from qscat.viz import complex_to_hsv, complex_to_rgb, region_magnitudes


def _two_region_field() -> np.ndarray:
    """Left half amplitude 1.0, right half amplitude 1e-6 -- the incident-packet
    vs outgoing-tail contrast that motivates the whole feature."""
    z = np.ones((8, 10), dtype=np.complex128)
    z[:, 5:] = 1e-6
    return z


def test_scalar_mag_is_unchanged() -> None:
    # Regression guard: the widening must not alter existing scalar behaviour.
    z = np.array([[1 + 1j, 0.5], [0.0, 2.0]], dtype=np.complex128)
    hsv = complex_to_hsv(z, 2.0)
    assert hsv.shape == z.shape + (3,)
    assert np.all((hsv >= 0.0) & (hsv <= 1.0))
    # |z|=2 at mag=2 saturates value to exactly 1.0
    assert hsv[1, 1, 2] == pytest.approx(1.0)


def test_scalar_mag_leaves_the_weak_region_black() -> None:
    # This is the BUG, asserted so the fix has something to improve on.
    z = _two_region_field()
    hsv = complex_to_hsv(z, 1.0)
    assert hsv[0, 0, 2] == pytest.approx(1.0)      # strong region: full brightness
    assert hsv[0, 9, 2] < 1e-5                      # weak region: invisible


def test_region_magnitudes_rescales_each_region_independently() -> None:
    z = _two_region_field()
    mag = region_magnitudes(np.abs(z), axis=1, boundaries=[5])
    assert mag.shape == z.shape
    hsv = complex_to_hsv(z, mag)
    # Both regions now reach full brightness -- the point of the feature.
    assert hsv[0, 0, 2] == pytest.approx(1.0)
    assert hsv[0, 9, 2] == pytest.approx(1.0)


def test_region_magnitudes_uses_a_percentile_not_the_max() -> None:
    # One hot pixel must not re-flatten its region. With max-scaling the bulk
    # would sit at 1/1000; with a percentile it stays visible.
    m = np.ones((1, 100))
    m[0, 0] = 1000.0
    mag = region_magnitudes(m, axis=1, boundaries=[], percentile=90.0)
    assert mag[0, 50] == pytest.approx(1.0, rel=1e-6)


def test_region_magnitudes_is_positive_even_for_an_all_zero_region() -> None:
    # An empty region must not produce a 0 scale and a divide-by-zero downstream.
    m = np.zeros((4, 6))
    mag = region_magnitudes(m, axis=1, boundaries=[3])
    assert np.all(mag > 0.0)
    assert np.all(np.isfinite(complex_to_rgb(m.astype(np.complex128), mag)))


def test_region_magnitudes_rejects_bad_boundaries() -> None:
    m = np.ones((4, 6))
    with pytest.raises(ValueError):
        region_magnitudes(m, axis=1, boundaries=[0])        # empty leading region
    with pytest.raises(ValueError):
        region_magnitudes(m, axis=1, boundaries=[6])        # empty trailing region
    with pytest.raises(ValueError):
        region_magnitudes(m, axis=1, boundaries=[4, 2])     # not increasing


def test_complex_to_hsv_rejects_a_non_broadcastable_mag() -> None:
    z = np.ones((4, 6), dtype=np.complex128)
    with pytest.raises(ValueError):
        complex_to_hsv(z, np.ones((3, 3)))
