"""Unit test for `correlation.hankel_point_value` (sub-project #C4, Task 2):
the scalar outgoing-Hankel-half VALUE at one physical electronic coordinate,
underlying the `Dirac` delta extractor's deconvolution factor.
"""

from __future__ import annotations

import numpy as np
from qscat.core.correlation import hankel_point_value
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
