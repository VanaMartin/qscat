"""Golden-value regression for the `Extractor` refactor (sub-project #C4,
Task 1): `td_ve_cross_section(..., method="tw")` must reproduce the
pre-refactor `td_ve_cross_section(...)` output to machine precision.

`_GOLDEN_TW_SCALAR`/`_GOLDEN_TW_ARRAY` were captured by running the CURRENT
(pre-refactor) `td_ve_cross_section` on the tiny/fast N2 config below --
copied from `libs/qscat/tests/test_core_td.py`'s deliberately tiny grid
(seconds, not the `TD_WORKING_GRID`-scale minutes), which that file's own
docstring already establishes is fast and physically unconverged-by-design
(a shape/contract + regression pin, not a converged cross section). Two
back-to-back same-process runs of the pre-refactor code gave a max abs diff
of exactly 0.0 (verified while capturing these values), so `atol=1e-12,
rtol=0` is the right bar here -- this is NOT the cross-process ~1e-9 BLAS-
threading drift `test_core_td.py` documents for its own regression pin.
"""

from __future__ import annotations

import numpy as np
from qscat.core.grids import electronic_grid, nuclear_grid
from qscat.core.time_dependent import td_ve_cross_section
from qscat.core.vibrational import vibrational_states
from qscat.dvr import TensorGrid
from qscat.model import N2

TG = TensorGrid(
    [
        electronic_grid(r_max=12.0, order=5, n_complex=3),
        nuclear_grid(quadrature=6, r_max=14.0, n_complex=3),
    ]
)
EPS, CHI = vibrational_states(TG.grids[1], N2.mu, 4, N2.v0)

V_INIT = 0
VPRIMES = [0, 1]  # includes the elastic (v'=v_init) channel
WP_IN = {"r0": 4.0, "p0": -0.5, "sigma": 1.2}
WP_OUT = {"r0_out": 6.0, "p0_out": 0.5, "sigma_out": 1.0}
DT = 0.2
N_STEPS = 5  # a handful of steps -- fast, not a converged run

# Captured from the pre-refactor `td_ve_cross_section` (direct `_propagate` +
# `sigma_from_correlations`) at this exact config -- see module docstring.
_GOLDEN_TW_SCALAR = np.array([1.1294075002622328e-05, 1.891186724741104e-09])
_GOLDEN_TW_ARRAY = np.array(
    [
        [1.1294075002622328e-05, 1.891186724741104e-09],
        [9.426919301560094e-06, 1.7032730710938486e-09],
    ]
)


def test_tw_method_matches_prerefactor_golden_scalar_energy() -> None:
    sigma = td_ve_cross_section(
        TG, N2, EPS, CHI, V_INIT, VPRIMES, 0.10,
        dt=DT, n_steps=N_STEPS, wp_in=WP_IN, wp_out=WP_OUT, method="tw",
    )
    assert sigma.shape == (len(VPRIMES),)
    np.testing.assert_allclose(sigma, _GOLDEN_TW_SCALAR, rtol=0, atol=1e-12)


def test_tw_method_matches_prerefactor_golden_array_energy() -> None:
    sigma = td_ve_cross_section(
        TG, N2, EPS, CHI, V_INIT, VPRIMES, [0.10, 0.15],
        dt=DT, n_steps=N_STEPS, wp_in=WP_IN, wp_out=WP_OUT, method="tw",
    )
    assert sigma.shape == (2, len(VPRIMES))
    np.testing.assert_allclose(sigma, _GOLDEN_TW_ARRAY, rtol=0, atol=1e-12)


def test_tw_method_is_the_default() -> None:
    """Omitting `method` must give the SAME result as `method="tw"` -- the
    refactor's default did not change from the caller's point of view."""
    sigma_default = td_ve_cross_section(
        TG, N2, EPS, CHI, V_INIT, VPRIMES, [0.10, 0.15],
        dt=DT, n_steps=N_STEPS, wp_in=WP_IN, wp_out=WP_OUT,
    )
    sigma_tw = td_ve_cross_section(
        TG, N2, EPS, CHI, V_INIT, VPRIMES, [0.10, 0.15],
        dt=DT, n_steps=N_STEPS, wp_in=WP_IN, wp_out=WP_OUT, method="tw",
    )
    np.testing.assert_allclose(sigma_default, sigma_tw, rtol=0, atol=0)


def test_unknown_method_raises() -> None:
    import pytest

    with pytest.raises(ValueError, match="unknown method"):
        td_ve_cross_section(
            TG, N2, EPS, CHI, V_INIT, VPRIMES, 0.10,
            dt=DT, n_steps=N_STEPS, wp_in=WP_IN, wp_out=WP_OUT, method="delta",
        )
