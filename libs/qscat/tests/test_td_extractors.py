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
import pytest
from qscat.core.driven import ve_cross_section
from qscat.core.grids import electronic_grid, nuclear_grid
from qscat.core.td_extractors import Dirac, TannorWeeks
from qscat.core.time_dependent import _free_hamiltonian, propagate, td_ve_cross_section
from qscat.core.vibrational import vibrational_states
from qscat.core.wavepacket import initial_state
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
    with pytest.raises(ValueError, match="unknown method"):
        td_ve_cross_section(
            TG, N2, EPS, CHI, V_INIT, VPRIMES, 0.10,
            dt=DT, n_steps=N_STEPS, wp_in=WP_IN, wp_out=WP_OUT, method="delta",
        )


# --- Task 2: Dirac (delta) extractor -----------------------------------------
#
# `POSITION = 37` is `TG.grids[0]`'s DVR index at r = 6.0 bohr -- inside the
# real (unscaled) electronic region (`R0 = 12.0`), an element end colocated
# with `WP_OUT`'s r0_out = 6.0 so the delta point samples the SAME region the
# TW test packet is centered on within this short/tiny-grid propagation (a
# position further out, e.g. r=10, undersamples the wavepacket entirely at
# these step counts and gives a spuriously tiny, not-yet-arrived b_v'(t) --
# see task-2-report.md).
POSITION = 37


def _propagate_pair(
    n_steps: int, dt: float
) -> tuple[TannorWeeks, Dirac, TannorWeeks, Dirac]:
    """One full (`V_int` on) + one free (`V_int=0`) propagation, each driving
    a `TannorWeeks` AND a `Dirac` extractor from the SAME trajectory (the
    differential test's whole point: one propagate() call, two independent
    analyses of identical psi(t_n))."""
    psi0 = initial_state(TG, CHI[V_INIT], **WP_IN)

    tw = TannorWeeks(TG, N2, EPS, CHI, V_INIT, VPRIMES, WP_OUT, wp_in=WP_IN, dt=dt)
    dirac = Dirac(TG, N2, EPS, CHI, V_INIT, VPRIMES, POSITION, wp_in=WP_IN, dt=dt)
    propagate(
        TG, psi0, [], dt=dt, n_steps=n_steps, hamiltonian=N2.hamiltonian(TG),
        extractors=[tw, dirac],
    )

    tw_free = TannorWeeks(TG, N2, EPS, CHI, V_INIT, VPRIMES, WP_OUT, wp_in=WP_IN, dt=dt)
    dirac_free = Dirac(TG, N2, EPS, CHI, V_INIT, VPRIMES, POSITION, wp_in=WP_IN, dt=dt)
    propagate(
        TG, psi0, [], dt=dt, n_steps=n_steps, hamiltonian=_free_hamiltonian(N2, TG),
        extractors=[tw_free, dirac_free],
    )
    return tw, dirac, tw_free, dirac_free


# Cross-method agreement band: TW and the delta extractor are two DIFFERENT
# analyses (a propagated Gaussian test packet vs. a fixed-point line
# projection) of the identical trajectory -- they need not, and empirically do
# not, agree to machine precision. `N_STEPS_DIFF = 800` (T = 160 a.u., ~16s):
# short runs (e.g. this file's own N_STEPS = 5 golden-regression config) give
# wildly different magnitudes (the delta point hasn't been "reached" by the
# packet yet -- see task-2-report.md); by n_steps=800 the two methods'
# sigma(E) ratios have settled to [0.81, 0.83, 0.90, 0.81] across
# (E, v') in {0.10, 0.15}x{0, 1} (re-checked at n_steps=1500: [0.82, 0.84,
# 0.91, 0.85] -- stable, not still drifting toward 1). `rtol=0.20` covers the
# observed ~18.6% worst-case deviation with a small margin; if a future
# change widens it, that is a finding to report, not silently re-loosen.
N_STEPS_DIFF = 800
_DELTA_TW_RTOL = 0.20


def test_delta_agrees_with_tw_same_trajectory() -> None:
    tw, dirac, tw_free, dirac_free = _propagate_pair(N_STEPS_DIFF, DT)
    e = [0.10, 0.15]
    s_tw = tw.sigma(e, free=tw_free)
    s_delta = dirac.sigma(e, free=dirac_free)
    assert s_tw.shape == s_delta.shape == (2, len(VPRIMES))
    assert np.all(np.isfinite(s_tw))
    assert np.all(np.isfinite(s_delta))
    np.testing.assert_allclose(s_delta, s_tw, rtol=_DELTA_TW_RTOL, atol=1e-14)


@pytest.mark.slow
def test_delta_agrees_with_ti_oracle_one_anchor() -> None:
    """A looser, independent check: `Dirac.sigma` vs the exact TI
    `ve_cross_section` oracle at one anchor energy -- the SAME converged
    working grid/wavepacket `docs/physics/n2-2d-td-cross-section.md`
    establishes makes TannorWeeks match the TI oracle to ~1-2% (at T=1500).
    This uses a shorter T=1000 (~240s, one propagation) and the
    INELASTIC-only channel `v'=1` (no free-reference run needed). Measured:
    `sigma_delta/sigma_ti = 0.971` at E=0.10 (`= 1.009` at E=0.15, checked but
    not gated here -- see task-2-report.md) -- `rtol=0.10` covers the
    observed ~2.9% deviation with margin.
    """
    tg_oracle = TensorGrid(
        [
            electronic_grid(r_max=50.0, order=8, n_complex=6),
            nuclear_grid(quadrature=10, r_max=22.0, n_complex=5),
        ]
    )
    eps, chi = vibrational_states(tg_oracle.grids[1], N2.mu, 4, N2.v0)
    v_init = 0
    vprimes = [1]  # inelastic only -- no elastic free-reference propagation needed
    wp_in = {"r0": 25.0, "p0": -0.5, "sigma": 5.0}
    dt = 1.0
    n_steps = 1000
    position = 128  # r = 39.58 bohr, real region (R0=50), past the interaction

    psi0 = initial_state(tg_oracle, chi[v_init], **wp_in)
    dirac = Dirac(tg_oracle, N2, eps, chi, v_init, vprimes, position, wp_in=wp_in, dt=dt)
    propagate(
        tg_oracle, psi0, [], dt=dt, n_steps=n_steps, hamiltonian=N2.hamiltonian(tg_oracle),
        extractors=[dirac],
    )
    e = 0.10
    s_delta = dirac.sigma(e)
    s_ti = ve_cross_section(tg_oracle, N2, eps, chi, v_init, vprimes, e)
    assert np.all(np.isfinite(s_delta))
    np.testing.assert_allclose(s_delta, s_ti, rtol=0.10, atol=1e-14)
