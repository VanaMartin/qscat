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
from qscat.core.td_extractors import Dirac, Flux, TannorWeeks
from qscat.core.time_dependent import (
    _free_hamiltonian,
    propagate,
    td_ve_cross_section,
    td_ve_cross_sections_all,
)
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
        TG,
        N2,
        EPS,
        CHI,
        V_INIT,
        VPRIMES,
        0.10,
        dt=DT,
        n_steps=N_STEPS,
        wp_in=WP_IN,
        wp_out=WP_OUT,
        method="tw",
    )
    assert sigma.shape == (len(VPRIMES),)
    np.testing.assert_allclose(sigma, _GOLDEN_TW_SCALAR, rtol=0, atol=1e-12)


def test_tw_method_matches_prerefactor_golden_array_energy() -> None:
    sigma = td_ve_cross_section(
        TG,
        N2,
        EPS,
        CHI,
        V_INIT,
        VPRIMES,
        [0.10, 0.15],
        dt=DT,
        n_steps=N_STEPS,
        wp_in=WP_IN,
        wp_out=WP_OUT,
        method="tw",
    )
    assert sigma.shape == (2, len(VPRIMES))
    np.testing.assert_allclose(sigma, _GOLDEN_TW_ARRAY, rtol=0, atol=1e-12)


def test_tw_method_is_the_default() -> None:
    """Omitting `method` must give the SAME result as `method="tw"` -- the
    refactor's default did not change from the caller's point of view."""
    sigma_default = td_ve_cross_section(
        TG,
        N2,
        EPS,
        CHI,
        V_INIT,
        VPRIMES,
        [0.10, 0.15],
        dt=DT,
        n_steps=N_STEPS,
        wp_in=WP_IN,
        wp_out=WP_OUT,
    )
    sigma_tw = td_ve_cross_section(
        TG,
        N2,
        EPS,
        CHI,
        V_INIT,
        VPRIMES,
        [0.10, 0.15],
        dt=DT,
        n_steps=N_STEPS,
        wp_in=WP_IN,
        wp_out=WP_OUT,
        method="tw",
    )
    np.testing.assert_allclose(sigma_default, sigma_tw, rtol=0, atol=0)


def test_unknown_method_raises() -> None:
    with pytest.raises(ValueError, match="unknown method"):
        td_ve_cross_section(
            TG,
            N2,
            EPS,
            CHI,
            V_INIT,
            VPRIMES,
            0.10,
            dt=DT,
            n_steps=N_STEPS,
            wp_in=WP_IN,
            wp_out=WP_OUT,
            method="bogus",
        )


# --- Task 4: method="delta"/"flow" wiring + required position/surface ------


def test_tw_method_requires_wp_out() -> None:
    with pytest.raises(ValueError, match="requires `wp_out`"):
        td_ve_cross_section(
            TG,
            N2,
            EPS,
            CHI,
            V_INIT,
            VPRIMES,
            0.10,
            dt=DT,
            n_steps=N_STEPS,
            wp_in=WP_IN,
            method="tw",
        )


def test_delta_method_requires_position() -> None:
    with pytest.raises(ValueError, match="requires `position`"):
        td_ve_cross_section(
            TG,
            N2,
            EPS,
            CHI,
            V_INIT,
            VPRIMES,
            0.10,
            dt=DT,
            n_steps=N_STEPS,
            wp_in=WP_IN,
            method="delta",
        )


def test_flow_method_requires_surface() -> None:
    with pytest.raises(ValueError, match="requires `surface`"):
        td_ve_cross_section(
            TG,
            N2,
            EPS,
            CHI,
            V_INIT,
            VPRIMES,
            0.10,
            dt=DT,
            n_steps=N_STEPS,
            wp_in=WP_IN,
            method="flow",
        )


# --- Task 2: Dirac (delta) extractor -----------------------------------------
#
# `POSITION = 37` is `TG.grids[0]`'s DVR index at r = 6.0 bohr -- inside the
# real (unscaled) electronic region (`R0 = 12.0`), an element end colocated
# with `WP_OUT`'s r0_out = 6.0 so the delta point samples the SAME region the
# TW test packet is centered on within this short/tiny-grid propagation (a
# position further out, e.g. r=10, undersamples the wavepacket entirely at
# these step counts and gives a spuriously tiny, not-yet-arrived b_v'(t) --
# see docs/physics/td-da.md).
POSITION = 37


def _propagate_pair(n_steps: int, dt: float) -> tuple[TannorWeeks, Dirac, TannorWeeks, Dirac]:
    """One full (`V_int` on) + one free (`V_int=0`) propagation, each driving
    a `TannorWeeks` AND a `Dirac` extractor from the SAME trajectory (the
    differential test's whole point: one propagate() call, two independent
    analyses of identical psi(t_n))."""
    psi0 = initial_state(TG, CHI[V_INIT], **WP_IN)

    tw = TannorWeeks(TG, N2, EPS, CHI, V_INIT, VPRIMES, WP_OUT, wp_in=WP_IN, dt=dt)
    dirac = Dirac(TG, N2, EPS, CHI, V_INIT, VPRIMES, POSITION, wp_in=WP_IN, dt=dt)
    propagate(
        TG,
        psi0,
        [],
        dt=dt,
        n_steps=n_steps,
        hamiltonian=N2.hamiltonian(TG),
        extractors=[tw, dirac],
    )

    tw_free = TannorWeeks(TG, N2, EPS, CHI, V_INIT, VPRIMES, WP_OUT, wp_in=WP_IN, dt=dt)
    dirac_free = Dirac(TG, N2, EPS, CHI, V_INIT, VPRIMES, POSITION, wp_in=WP_IN, dt=dt)
    propagate(
        TG,
        psi0,
        [],
        dt=dt,
        n_steps=n_steps,
        hamiltonian=_free_hamiltonian(N2, TG),
        extractors=[tw_free, dirac_free],
    )
    return tw, dirac, tw_free, dirac_free


# Cross-method agreement band: TW and the delta extractor are two DIFFERENT
# analyses (a propagated Gaussian test packet vs. a fixed-point line
# projection) of the identical trajectory -- they need not, and empirically do
# not, agree to machine precision. `N_STEPS_DIFF = 800` (T = 160 a.u., ~16s):
# short runs (e.g. this file's own N_STEPS = 5 golden-regression config) give
# wildly different magnitudes (the delta point hasn't been "reached" by the
# packet yet -- see docs/physics/td-da.md); by n_steps=800 the two methods'
# sigma(E) ratios have settled to [0.81, 0.83, 0.90, 0.81] across
# (E, v') in {0.10, 0.15}x{0, 1} (re-checked at n_steps=1500: [0.82, 0.84,
# 0.91, 0.85] -- stable, not still drifting toward 1). `rtol=0.20` covers the
# observed ~18.6% worst-case deviation with a small margin; if a future
# change widens it, that is a finding to report, not silently re-loosen.
N_STEPS_DIFF = 800
_DELTA_TW_RTOL = 0.20


@pytest.mark.slow
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
    not gated here -- see docs/physics/td-da.md) -- `rtol=0.10` covers the
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
        tg_oracle,
        psi0,
        [],
        dt=dt,
        n_steps=n_steps,
        hamiltonian=N2.hamiltonian(tg_oracle),
        extractors=[dirac],
    )
    e = 0.10
    s_delta = dirac.sigma(e)
    s_ti = ve_cross_section(tg_oracle, N2, eps, chi, v_init, vprimes, e)
    assert np.all(np.isfinite(s_delta))
    np.testing.assert_allclose(s_delta, s_ti, rtol=0.10, atol=1e-14)


# --- Task 3: Flux (flow) extractor --------------------------------------------
#
# Reuses `POSITION = 37` (r = 6.0, colocated with `WP_OUT`) as the flux
# surface -- the same asymptotic electronic node `Dirac` uses; no reason for
# the flow extractor's surface to differ from the delta extractor's fixed
# analysis point in this differential test.


def _propagate_all_three(
    n_steps: int, dt: float
) -> tuple[TannorWeeks, Dirac, Flux, TannorWeeks, Dirac, Flux]:
    """ONE `propagate()` call drives `TannorWeeks` + `Dirac` + `Flux` together
    (plus a second call for the `V_int=0` free reference) -- the brief's
    differential-test gate: all three extractors' `sigma(E)` come from the
    IDENTICAL trajectory."""
    psi0 = initial_state(TG, CHI[V_INIT], **WP_IN)

    tw = TannorWeeks(TG, N2, EPS, CHI, V_INIT, VPRIMES, WP_OUT, wp_in=WP_IN, dt=dt)
    dirac = Dirac(TG, N2, EPS, CHI, V_INIT, VPRIMES, POSITION, wp_in=WP_IN, dt=dt)
    flux = Flux(TG, N2, EPS, CHI, V_INIT, VPRIMES, POSITION, wp_in=WP_IN, dt=dt)
    propagate(
        TG,
        psi0,
        [],
        dt=dt,
        n_steps=n_steps,
        hamiltonian=N2.hamiltonian(TG),
        extractors=[tw, dirac, flux],
    )

    tw_free = TannorWeeks(TG, N2, EPS, CHI, V_INIT, VPRIMES, WP_OUT, wp_in=WP_IN, dt=dt)
    dirac_free = Dirac(TG, N2, EPS, CHI, V_INIT, VPRIMES, POSITION, wp_in=WP_IN, dt=dt)
    flux_free = Flux(TG, N2, EPS, CHI, V_INIT, VPRIMES, POSITION, wp_in=WP_IN, dt=dt)
    propagate(
        TG,
        psi0,
        [],
        dt=dt,
        n_steps=n_steps,
        hamiltonian=_free_hamiltonian(N2, TG),
        extractors=[tw_free, dirac_free, flux_free],
    )
    return tw, dirac, flux, tw_free, dirac_free, flux_free


# Cross-method agreement band: Flux is a genuinely DIFFERENT analysis (a
# surface-current Wronskian at a FIXED electronic node vs. TW's propagated-
# test-packet volume overlap) of the identical trajectory -- no reason to
# expect machine-precision agreement. Measured (standalone script, same
# `N_STEPS_DIFF = 800` config this file already uses for the Dirac gate):
#
#   s_tw   = [[141.145, 43.485], [23.694, 2.562]]
#   s_flux = [[108.077, 35.513], [21.113, 1.970]]
#   ratio  = [[0.7657, 0.8167], [0.8911, 0.7692]]   (worst-case ~23.4% low)
#
# Re-checked at n_steps=1500: ratio [[0.7714,0.8204],[0.9098,0.8185]] -- the
# same band, not still drifting toward 1 (same qualitative pattern Dirac's
# ~18-19% band showed at these two step counts). `rtol=0.25` covers the
# observed ~23.4% worst-case deviation with a small margin; a future change
# that widens this further is a finding to report, not to silently re-loosen.
_FLUX_TW_RTOL = 0.25


@pytest.mark.slow
def test_flux_agrees_with_tw_same_trajectory() -> None:
    tw, dirac, flux, tw_free, dirac_free, flux_free = _propagate_all_three(N_STEPS_DIFF, DT)
    e = [0.10, 0.15]
    s_tw = tw.sigma(e, free=tw_free)
    s_dirac = dirac.sigma(e, free=dirac_free)  # recorded from the SAME pass; not gated here
    s_flux = flux.sigma(e, free=flux_free)
    assert s_tw.shape == s_flux.shape == s_dirac.shape == (2, len(VPRIMES))
    assert np.all(np.isfinite(s_tw))
    assert np.all(np.isfinite(s_dirac))
    assert np.all(np.isfinite(s_flux))
    np.testing.assert_allclose(s_flux, s_tw, rtol=_FLUX_TW_RTOL, atol=1e-14)


@pytest.mark.slow
def test_flux_agrees_with_ti_oracle_one_anchor() -> None:
    """A looser, independent check: `Flux.sigma` vs the exact TI
    `ve_cross_section` oracle at one anchor energy -- the SAME converged
    working grid/wavepacket/surface `test_delta_agrees_with_ti_oracle_one_
    anchor` uses (T=1000, inelastic-only channel `v'=1`, no free-reference
    run needed). See that test's docstring for the shared config rationale.

    Measured: `sigma_flux/sigma_ti = 0.970` at E=0.10 (`= 1.007` at E=0.15,
    checked but not gated here) -- essentially the SAME magnitude as
    `Dirac`'s own 0.971/1.009 at this grid (both alternative extractors
    converge to the TI oracle with a similar residual, driven by the shared
    propagation/discretization truncation rather than by which extraction
    method is used). `rtol=0.10` covers the observed ~3.0% deviation with
    margin, matching `Dirac`'s gate.
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
    surface = 128  # r = 39.58 bohr, real region (R0=50), past the interaction

    psi0 = initial_state(tg_oracle, chi[v_init], **wp_in)
    flux = Flux(tg_oracle, N2, eps, chi, v_init, vprimes, surface, wp_in=wp_in, dt=dt)
    propagate(
        tg_oracle,
        psi0,
        [],
        dt=dt,
        n_steps=n_steps,
        hamiltonian=N2.hamiltonian(tg_oracle),
        extractors=[flux],
    )
    e = 0.10
    s_flux = flux.sigma(e)
    s_ti = ve_cross_section(tg_oracle, N2, eps, chi, v_init, vprimes, e)
    assert np.all(np.isfinite(s_flux))
    np.testing.assert_allclose(s_flux, s_ti, rtol=0.10, atol=1e-14)


# --- Task 4: method="delta"/"flow" wiring + td_ve_cross_sections_all --------
#
# `td_ve_cross_section(method="delta"/"flow")` must reproduce building the
# `Dirac`/`Flux` extractor directly (as `_propagate_pair`/`_propagate_all_
# three` above already do) to machine precision -- it is the exact same
# construction + `propagate` call, just wired behind the `method=` string.


def test_delta_method_matches_direct_dirac_construction() -> None:
    sigma_method = td_ve_cross_section(
        TG,
        N2,
        EPS,
        CHI,
        V_INIT,
        VPRIMES,
        [0.10, 0.15],
        dt=DT,
        n_steps=N_STEPS,
        wp_in=WP_IN,
        method="delta",
        position=POSITION,
    )

    psi0 = initial_state(TG, CHI[V_INIT], **WP_IN)
    dirac = Dirac(TG, N2, EPS, CHI, V_INIT, VPRIMES, POSITION, wp_in=WP_IN, dt=DT)
    propagate(
        TG,
        psi0,
        [],
        dt=DT,
        n_steps=N_STEPS,
        hamiltonian=N2.hamiltonian(TG),
        extractors=[dirac],
    )
    dirac_free = Dirac(TG, N2, EPS, CHI, V_INIT, VPRIMES, POSITION, wp_in=WP_IN, dt=DT)
    propagate(
        TG,
        psi0,
        [],
        dt=DT,
        n_steps=N_STEPS,
        hamiltonian=_free_hamiltonian(N2, TG),
        extractors=[dirac_free],
    )
    sigma_direct = dirac.sigma([0.10, 0.15], free=dirac_free)
    np.testing.assert_allclose(sigma_method, sigma_direct, rtol=0, atol=0)


def test_flow_method_matches_direct_flux_construction() -> None:
    sigma_method = td_ve_cross_section(
        TG,
        N2,
        EPS,
        CHI,
        V_INIT,
        VPRIMES,
        [0.10, 0.15],
        dt=DT,
        n_steps=N_STEPS,
        wp_in=WP_IN,
        method="flow",
        surface=POSITION,
    )

    psi0 = initial_state(TG, CHI[V_INIT], **WP_IN)
    flux = Flux(TG, N2, EPS, CHI, V_INIT, VPRIMES, POSITION, wp_in=WP_IN, dt=DT)
    propagate(
        TG,
        psi0,
        [],
        dt=DT,
        n_steps=N_STEPS,
        hamiltonian=N2.hamiltonian(TG),
        extractors=[flux],
    )
    flux_free = Flux(TG, N2, EPS, CHI, V_INIT, VPRIMES, POSITION, wp_in=WP_IN, dt=DT)
    propagate(
        TG,
        psi0,
        [],
        dt=DT,
        n_steps=N_STEPS,
        hamiltonian=_free_hamiltonian(N2, TG),
        extractors=[flux_free],
    )
    sigma_direct = flux.sigma([0.10, 0.15], free=flux_free)
    np.testing.assert_allclose(sigma_method, sigma_direct, rtol=0, atol=0)


def test_cross_sections_all_matches_each_method_individually() -> None:
    """`td_ve_cross_sections_all`'s ONE-propagation result must reproduce
    calling `td_ve_cross_section` once per method (each its OWN, separately
    propagated, trajectory) to machine precision -- the propagation is
    deterministic, so driving three independent `Extractor`s from a shared
    trajectory changes nothing about what each individually records."""
    e = [0.10, 0.15]
    sigma_all = td_ve_cross_sections_all(
        TG,
        N2,
        EPS,
        CHI,
        V_INIT,
        VPRIMES,
        e,
        dt=DT,
        n_steps=N_STEPS,
        wp_in=WP_IN,
        wp_out=WP_OUT,
        position=POSITION,
        surface=POSITION,
    )
    assert set(sigma_all) == {"tw", "delta", "flow"}
    for key, method in (("tw", "tw"), ("delta", "delta"), ("flow", "flow")):
        kwargs = {"position": POSITION} if key == "delta" else {}
        if key == "flow":
            kwargs = {"surface": POSITION}
        sigma_individual = td_ve_cross_section(
            TG,
            N2,
            EPS,
            CHI,
            V_INIT,
            VPRIMES,
            e,
            dt=DT,
            n_steps=N_STEPS,
            wp_in=WP_IN,
            wp_out=WP_OUT,
            method=method,
            **kwargs,
        )
        assert sigma_all[key].shape == (2, len(VPRIMES))
        np.testing.assert_allclose(sigma_all[key], sigma_individual, rtol=0, atol=0)


# --- Task 1 (SP2): `axis` scaffolding ----------------------------------------
#
# `axis="electronic"` (default) is covered by every test above, which all
# construct extractors without an `axis` kwarg and hit the byte-identical
# golden values. `axis="nuclear"` is now implemented for all three
# extractors -- `Flux` (Task 2), `Dirac` (Task 3), `TannorWeeks` (Task 4,
# see the nuclear-TannorWeeks section further down for its own coverage) --
# so there is no remaining `NotImplementedError` stub to test here. An
# invalid `axis` is a `ValueError` for all three, and the default matches an
# explicit `axis="electronic"` for all three.

_AXIS_CTOR_ARGS: dict[str, tuple[type, tuple]] = {
    "TannorWeeks": (TannorWeeks, (TG, N2, EPS, CHI, V_INIT, VPRIMES, WP_OUT)),
    "Dirac": (Dirac, (TG, N2, EPS, CHI, V_INIT, VPRIMES, POSITION)),
    "Flux": (Flux, (TG, N2, EPS, CHI, V_INIT, VPRIMES, POSITION)),
}


@pytest.mark.parametrize("cls_name", ["TannorWeeks", "Dirac", "Flux"])
def test_axis_invalid_raises_value_error(cls_name: str) -> None:
    cls, args = _AXIS_CTOR_ARGS[cls_name]
    with pytest.raises(ValueError, match="axis must be one of"):
        cls(*args, wp_in=WP_IN, dt=DT, axis="bogus")


@pytest.mark.parametrize("cls_name", ["TannorWeeks", "Dirac", "Flux"])
def test_axis_default_matches_explicit_electronic(cls_name: str) -> None:
    cls, args = _AXIS_CTOR_ARGS[cls_name]
    default = cls(*args, wp_in=WP_IN, dt=DT)
    explicit = cls(*args, wp_in=WP_IN, dt=DT, axis="electronic")
    assert default._axis == explicit._axis == "electronic"


# --- Task 2 (SP2): nuclear-axis `Flux` (dissociative attachment) ------------
#
# Fast structural tests: builds/records/computes sigma on this file's tiny
# N2 config (`n_channels=1` -- the only bound anion electronic state that
# config's small `R_inf=grids[1].R0=12.0` supports, see
# `anion_electronic_states`'s docstring). NOT a converged DA cross section
# (same caveat the VE golden-regression config carries) -- this just checks
# "builds, records, sigma runs, threshold gating correct" (shape/finiteness/
# zero-below-threshold/zero-on-a-closed-channel/no-free-reference). The
# LOAD-BEARING convergence gate (`@slow`, pins `_C_DA`) is further below,
# against the TI `da_cross_section` oracle on F2's real (eMoScat) DA grid.

NUCLEAR_SURFACE = 90  # R=7.12 bohr on TG.grids[1] (nuclear_grid r_max=14, real region R0=12.0)


def _nuclear_flux_fixture(n_steps: int = N_STEPS) -> Flux:
    flux = Flux(
        TG,
        N2,
        EPS,
        CHI,
        V_INIT,
        [],
        NUCLEAR_SURFACE,
        wp_in=WP_IN,
        dt=DT,
        axis="nuclear",
        n_channels=1,
    )
    psi0 = initial_state(TG, CHI[V_INIT], **WP_IN)
    propagate(
        TG,
        psi0,
        [],
        dt=DT,
        n_steps=n_steps,
        hamiltonian=N2.hamiltonian(TG),
        extractors=[flux],
    )
    return flux


def test_nuclear_flux_builds_and_records() -> None:
    flux = _nuclear_flux_fixture()
    t, b, d = flux.series
    n_recorded = t.shape[0]  # N_STEPS+1: propagate() records the t=0 state too
    assert n_recorded == N_STEPS + 1
    assert b.shape == d.shape == (n_recorded, 1)
    assert np.all(np.isfinite(b))
    assert np.all(np.isfinite(d))


def test_nuclear_flux_sigma_shape_and_finite() -> None:
    flux = _nuclear_flux_fixture()
    s_scalar = flux.sigma(0.6)
    assert s_scalar.shape == (1,)
    assert np.all(np.isfinite(s_scalar))
    s_array = flux.sigma([0.10, 0.6])
    assert s_array.shape == (2, 1)
    assert np.all(np.isfinite(s_array))


def test_nuclear_flux_zero_at_or_below_threshold() -> None:
    flux = _nuclear_flux_fixture()
    np.testing.assert_allclose(flux.sigma(-0.1), [0.0])
    np.testing.assert_allclose(flux.sigma(0.0), [0.0])


def test_nuclear_flux_closed_dissociation_channel_gives_zero() -> None:
    """At this tiny grid's (too-small-`R_inf`, hence too-high) `eps_e`, E=0.10
    leaves the dissociation channel closed (`e_tot - eps_e[0] <= 0`) -- sigma
    must be exactly zero, distinct from the `E<=0` branch above."""
    flux = _nuclear_flux_fixture()
    e_tot = 0.10 + EPS[V_INIT]
    assert e_tot - flux._eps_e[0] <= 0.0  # sanity: confirms the closed-channel premise
    np.testing.assert_allclose(flux.sigma(0.10), [0.0])


def test_nuclear_flux_open_channel_path_is_finite() -> None:
    """At an energy that DOES clear this tiny grid's (unconverged) `eps_e`
    threshold, the open-channel branch executes (not skipped by `continue`)
    and returns a finite (though not physically converged) value."""
    flux = _nuclear_flux_fixture()
    e_tot = 0.6 + EPS[V_INIT]
    assert e_tot - flux._eps_e[0] > 0.0  # sanity: confirms the open-channel premise
    s = flux.sigma(0.6)
    assert np.all(np.isfinite(s))
    assert np.all(s >= 0.0)


def test_nuclear_flux_rejects_free_reference() -> None:
    flux = _nuclear_flux_fixture()
    with pytest.raises(ValueError, match="no elastic free-reference"):
        flux.sigma(0.10, free=flux)


def test_nuclear_flux_n_channels_defaults_to_one() -> None:
    flux = Flux(TG, N2, EPS, CHI, V_INIT, [], NUCLEAR_SURFACE, wp_in=WP_IN, dt=DT, axis="nuclear")
    assert flux._n_channels == 1


@pytest.mark.slow
def test_nuclear_flux_da_converges_to_ti_oracle() -> None:
    """LOAD-BEARING (pins `_C_DA`): the nuclear-`Flux` sigma_DA converges to the
    TI `da_cross_section` on the SAME grid (a differential test -- both must
    agree even where the grid is not physically converged).

    Config = a CLEAN launch-box electronic grid (r_max=25, the incident r0=12
    sits well INSIDE it, so the packet launches without ECS-tail garbage) x the
    FINE eMoScat F2 nuclear deck (which resolves the fast K_R~72 dissociation
    flux wave at the surface -- a coarse nuclear grid gives sigma_flux~0, the
    surface never sees the outgoing wave). Both failure modes were real bugs
    found while validating this: an off-box incident diverges ~1e6x; a coarse
    nuclear reads ~0.

    Measured (controller, 2026-07-31): sigma_flux/sigma_ti reaches a STABLE
    plateau ~0.86-0.97 by n>=1350 (|psi| flat at 0.0556) -- the ~3-14% gap is
    the TD-vs-TI cross-method band (as in the VE extractors), NOT a
    normalization error: a wrong `_C_DA` (e.g. 4*pi^2 off the `S=1-2*pi*i*T`
    identity) would plateau at a wildly different constant, not ~1. Heavy
    (~86k unknowns x 1500 steps, ~10 min) -- @slow. The FULL eMoScat F2 grid
    (electronic real->90 bohr, ~402k unknowns) convergence is a
    Docker/overnight-deferred production run, not this laptop gate.
    """
    from qscat.core.dissociation import da_cross_section
    from qscat.core.grids import segmented_grid
    from qscat.model import F2

    elec = electronic_grid(r_max=25.0, order=6, n_complex=3, angle_deg=40.0)
    nuc = segmented_grid(
        [(9, 1.8), (1, 2.0), (5, 2.5), (4, 2.596908), (4, 2.7), (40, 10.7)],
        [(1, 10.8), (1, 11.0), (1, 11.5), (1, 12.5), (1, 14.0), (1, 18.0), (4, 30.0), (2, 101.0)],
        angle_deg=35.0,
        quadrature=14,
    )
    tg = TensorGrid([elec, nuc])
    eps, chi = vibrational_states(nuc, F2.mu, 4, F2.v0)
    e_probe = np.array([0.03, 0.04])
    sigma_ti = np.ravel(da_cross_section(tg, F2, eps, chi, 0, e_probe))

    real = nuc.real_points
    surface = int(np.argmin(np.abs(np.where(real <= nuc.R0, real, 1e9) - 6.0)))
    wp_in = {"r0": 12.0, "p0": -0.5, "sigma": 3.0}
    psi0 = initial_state(tg, chi[0], **wp_in)
    flux = Flux(tg, F2, eps, chi, 0, [], surface, wp_in=wp_in, dt=1.0, axis="nuclear", n_channels=1)
    # `propagate` records the t=0 state itself, then every step -> 1501 samples.
    propagate(tg, psi0, [], dt=1.0, n_steps=1500, hamiltonian=F2.hamiltonian(tg), extractors=[flux])
    ratio = np.ravel(flux.sigma(e_probe)) / sigma_ti
    assert np.all(ratio > 0.7) and np.all(ratio < 1.25), (ratio, sigma_ti)


# --- Task 3 (SP2): nuclear-axis `Dirac` (delta) DA extractor ----------------
#
# Fast structural tests: builds/records/computes sigma on this file's tiny
# N2 config, mirroring the nuclear-`Flux` section above exactly (same
# `NUCLEAR_SURFACE`, same tiny-grid caveats -- not a converged DA cross
# section, just "builds, records, sigma runs, threshold gating correct").
# The LOAD-BEARING convergence gate (`@slow`) is further below, against the
# TI `da_cross_section` oracle on F2's real (eMoScat) DA grid -- a mirror of
# `test_nuclear_flux_da_converges_to_ti_oracle`.


def _nuclear_dirac_fixture(n_steps: int = N_STEPS) -> Dirac:
    dirac = Dirac(
        TG,
        N2,
        EPS,
        CHI,
        V_INIT,
        [],
        NUCLEAR_SURFACE,
        wp_in=WP_IN,
        dt=DT,
        axis="nuclear",
        n_channels=1,
    )
    psi0 = initial_state(TG, CHI[V_INIT], **WP_IN)
    propagate(
        TG,
        psi0,
        [],
        dt=DT,
        n_steps=n_steps,
        hamiltonian=N2.hamiltonian(TG),
        extractors=[dirac],
    )
    return dirac


def test_nuclear_dirac_builds_and_records() -> None:
    dirac = _nuclear_dirac_fixture()
    result = dirac.result
    n_recorded = result.t.shape[0]  # N_STEPS+1: propagate() records the t=0 state too
    assert n_recorded == N_STEPS + 1
    assert result.c.shape == (n_recorded, 1)
    assert np.all(np.isfinite(result.c))


def test_nuclear_dirac_sigma_shape_and_finite() -> None:
    dirac = _nuclear_dirac_fixture()
    s_scalar = dirac.sigma(0.6)
    assert s_scalar.shape == (1,)
    assert np.all(np.isfinite(s_scalar))
    s_array = dirac.sigma([0.10, 0.6])
    assert s_array.shape == (2, 1)
    assert np.all(np.isfinite(s_array))


def test_nuclear_dirac_zero_at_or_below_threshold() -> None:
    dirac = _nuclear_dirac_fixture()
    np.testing.assert_allclose(dirac.sigma(-0.1), [0.0])
    np.testing.assert_allclose(dirac.sigma(0.0), [0.0])


def test_nuclear_dirac_closed_dissociation_channel_gives_zero() -> None:
    """At this tiny grid's (too-small-`R_inf`, hence too-high) `eps_e`, E=0.10
    leaves the dissociation channel closed (`e_tot - eps_e[0] <= 0`) -- sigma
    must be exactly zero, distinct from the `E<=0` branch above."""
    dirac = _nuclear_dirac_fixture()
    e_tot = 0.10 + EPS[V_INIT]
    assert e_tot - dirac._eps_e[0] <= 0.0  # sanity: confirms the closed-channel premise
    np.testing.assert_allclose(dirac.sigma(0.10), [0.0])


def test_nuclear_dirac_open_channel_path_is_finite() -> None:
    """At an energy that DOES clear this tiny grid's (unconverged) `eps_e`
    threshold, the open-channel branch executes (not skipped by `continue`)
    and returns a finite (though not physically converged) value."""
    dirac = _nuclear_dirac_fixture()
    e_tot = 0.6 + EPS[V_INIT]
    assert e_tot - dirac._eps_e[0] > 0.0  # sanity: confirms the open-channel premise
    s = dirac.sigma(0.6)
    assert np.all(np.isfinite(s))
    assert np.all(s >= 0.0)


def test_nuclear_dirac_rejects_free_reference() -> None:
    dirac = _nuclear_dirac_fixture()
    with pytest.raises(ValueError, match="no elastic free-reference"):
        dirac.sigma(0.10, free=dirac)


def test_nuclear_dirac_n_channels_defaults_to_one() -> None:
    dirac = Dirac(TG, N2, EPS, CHI, V_INIT, [], NUCLEAR_SURFACE, wp_in=WP_IN, dt=DT, axis="nuclear")
    assert dirac._n_channels == 1


# Cross-method agreement band: nuclear `Dirac` and nuclear `Flux` are two
# DIFFERENT analyses (a fixed-point projection vs. a surface-current
# Wronskian) of the identical trajectory -- no reason to expect
# machine-precision agreement, mirroring the electronic Dirac-vs-Flux
# difference documented above. Measured (standalone script, this file's tiny
# N2 config, `NUCLEAR_SURFACE=90`, `n_steps=800`): `sigma_dirac/sigma_flux =
# 0.740` at E=0.6 (the only open dissociation channel on this tiny grid --
# `E=0.10` is closed for both, see `test_nuclear_dirac_closed_dissociation_
# channel_gives_zero`); re-checked at `n_steps=1500`: `0.672` (same
# qualitative band, still drifting on this deliberately tiny/toy grid -- NOT
# the convergence gate, see the `@slow` test below). `rtol=0.35` covers the
# observed ~33% worst-case deviation with a small margin; a future change
# that widens this further is a finding to report, not to silently re-loosen.
_NUCLEAR_DIRAC_FLUX_RTOL = 0.35
_NUCLEAR_DIRAC_FLUX_N_STEPS = 800


def test_nuclear_dirac_agrees_with_nuclear_flux_same_trajectory() -> None:
    """ONE `propagate()` call drives nuclear `Dirac` + nuclear `Flux`
    together -- the differential test's whole point: two independent
    analyses of the identical trajectory."""
    psi0 = initial_state(TG, CHI[V_INIT], **WP_IN)
    dirac = Dirac(
        TG,
        N2,
        EPS,
        CHI,
        V_INIT,
        [],
        NUCLEAR_SURFACE,
        wp_in=WP_IN,
        dt=DT,
        axis="nuclear",
        n_channels=1,
    )
    flux = Flux(
        TG,
        N2,
        EPS,
        CHI,
        V_INIT,
        [],
        NUCLEAR_SURFACE,
        wp_in=WP_IN,
        dt=DT,
        axis="nuclear",
        n_channels=1,
    )
    propagate(
        TG,
        psi0,
        [],
        dt=DT,
        n_steps=_NUCLEAR_DIRAC_FLUX_N_STEPS,
        hamiltonian=N2.hamiltonian(TG),
        extractors=[dirac, flux],
    )
    e = [0.10, 0.6]
    s_dirac = dirac.sigma(e)
    s_flux = flux.sigma(e)
    assert s_dirac.shape == s_flux.shape == (2, 1)
    assert np.all(np.isfinite(s_dirac))
    assert np.all(np.isfinite(s_flux))
    np.testing.assert_allclose(s_dirac[0], s_flux[0], rtol=0, atol=0)  # both closed -> exactly 0
    np.testing.assert_allclose(s_dirac[1], s_flux[1], rtol=_NUCLEAR_DIRAC_FLUX_RTOL, atol=1e-60)


@pytest.mark.slow
def test_nuclear_dirac_da_converges_to_ti_oracle() -> None:
    """LOAD-BEARING: the nuclear-`Dirac` sigma_DA converges to the TI
    `da_cross_section` on the SAME grid (a differential test -- both must
    agree even where the grid is not physically converged) -- a mirror of
    `test_nuclear_flux_da_converges_to_ti_oracle` (see that test's docstring
    for the full config rationale: the launch-box electronic grid x the
    FINE eMoScat F2 nuclear deck, the off-box-incident and coarse-nuclear
    failure modes it guards against).

    The sibling nuclear `Flux` gate plateaus at sigma_flux/sigma_ti ~
    0.86-0.97 by n>=1350 on this same config (`docs/physics/td-da.md`/the test
    above) -- `Dirac` is a different (point-value, not Wronskian-flux)
    transform of the identical propagation, so it need not land at the exact
    same ratio, but should land in the same general TD-vs-TI cross-method
    band, not at a wildly different constant (which would indicate a wrong
    `_C_DA` or a sign/prefactor error in the point-Hankel transform). Heavy
    (~86k unknowns x 1500 steps, ~10 min) -- @slow; mirrors the
    already-validated `Flux` sibling.
    """
    from qscat.core.dissociation import da_cross_section
    from qscat.core.grids import segmented_grid
    from qscat.model import F2

    elec = electronic_grid(r_max=25.0, order=6, n_complex=3, angle_deg=40.0)
    nuc = segmented_grid(
        [(9, 1.8), (1, 2.0), (5, 2.5), (4, 2.596908), (4, 2.7), (40, 10.7)],
        [(1, 10.8), (1, 11.0), (1, 11.5), (1, 12.5), (1, 14.0), (1, 18.0), (4, 30.0), (2, 101.0)],
        angle_deg=35.0,
        quadrature=14,
    )
    tg = TensorGrid([elec, nuc])
    eps, chi = vibrational_states(nuc, F2.mu, 4, F2.v0)
    e_probe = np.array([0.03, 0.04])
    sigma_ti = np.ravel(da_cross_section(tg, F2, eps, chi, 0, e_probe))

    real = nuc.real_points
    surface = int(np.argmin(np.abs(np.where(real <= nuc.R0, real, 1e9) - 6.0)))
    wp_in = {"r0": 12.0, "p0": -0.5, "sigma": 3.0}
    psi0 = initial_state(tg, chi[0], **wp_in)
    dirac = Dirac(
        tg, F2, eps, chi, 0, [], surface, wp_in=wp_in, dt=1.0, axis="nuclear", n_channels=1
    )
    # `propagate` records the t=0 state itself, then every step -> 1501 samples.
    propagate(
        tg, psi0, [], dt=1.0, n_steps=1500, hamiltonian=F2.hamiltonian(tg), extractors=[dirac]
    )
    ratio = np.ravel(dirac.sigma(e_probe)) / sigma_ti
    assert np.all(ratio > 0.7) and np.all(ratio < 1.25), (ratio, sigma_ti)


# --- Task 4 (SP2): nuclear-axis `TannorWeeks` (dissociative attachment) -----
#
# Fast structural tests: builds/records/computes sigma on this file's tiny
# N2 config, mirroring the nuclear-`Flux`/nuclear-`Dirac` sections above
# exactly (same tiny-grid caveats -- not a converged DA cross section, just
# "builds, records, sigma runs, threshold gating correct"). Unlike `Flux`/
# `Dirac`, `TannorWeeks` needs its own outgoing test-packet parameters
# (`wp_out`, now in the NUCLEAR coordinate `R`) rather than a fixed DVR
# index -- `NUCLEAR_WP_OUT` is centered near `NUCLEAR_SURFACE`'s R=7.12,
# with an outward (positive) momentum, well inside the tiny grid's real
# region (`R0=12.0`). The LOAD-BEARING convergence gate (`@slow`) is further
# below, against the TI `da_cross_section` oracle on F2's real (eMoScat) DA
# grid -- a mirror of `test_nuclear_flux_da_converges_to_ti_oracle`/
# `test_nuclear_dirac_da_converges_to_ti_oracle`.

NUCLEAR_WP_OUT = {"r0_out": 7.0, "p0_out": 5.0, "sigma_out": 1.0}


def _nuclear_tw_fixture(n_steps: int = N_STEPS) -> TannorWeeks:
    tw = TannorWeeks(
        TG,
        N2,
        EPS,
        CHI,
        V_INIT,
        [],
        NUCLEAR_WP_OUT,
        wp_in=WP_IN,
        dt=DT,
        axis="nuclear",
        n_channels=1,
    )
    psi0 = initial_state(TG, CHI[V_INIT], **WP_IN)
    propagate(
        TG,
        psi0,
        [],
        dt=DT,
        n_steps=n_steps,
        hamiltonian=N2.hamiltonian(TG),
        extractors=[tw],
    )
    return tw


def test_nuclear_tw_builds_and_records() -> None:
    tw = _nuclear_tw_fixture()
    result = tw.result
    n_recorded = result.t.shape[0]  # N_STEPS+1: propagate() records the t=0 state too
    assert n_recorded == N_STEPS + 1
    assert result.c.shape == (n_recorded, 1)
    assert np.all(np.isfinite(result.c))


def test_nuclear_tw_sigma_shape_and_finite() -> None:
    tw = _nuclear_tw_fixture()
    s_scalar = tw.sigma(0.6)
    assert s_scalar.shape == (1,)
    assert np.all(np.isfinite(s_scalar))
    s_array = tw.sigma([0.10, 0.6])
    assert s_array.shape == (2, 1)
    assert np.all(np.isfinite(s_array))


def test_nuclear_tw_zero_at_or_below_threshold() -> None:
    tw = _nuclear_tw_fixture()
    np.testing.assert_allclose(tw.sigma(-0.1), [0.0])
    np.testing.assert_allclose(tw.sigma(0.0), [0.0])


def test_nuclear_tw_closed_dissociation_channel_gives_zero() -> None:
    """At this tiny grid's (too-small-`R_inf`, hence too-high) `eps_e`, E=0.10
    leaves the dissociation channel closed (`e_tot - eps_e[0] <= 0`) -- sigma
    must be exactly zero, distinct from the `E<=0` branch above."""
    tw = _nuclear_tw_fixture()
    e_tot = 0.10 + EPS[V_INIT]
    assert e_tot - tw._eps_e[0] <= 0.0  # sanity: confirms the closed-channel premise
    np.testing.assert_allclose(tw.sigma(0.10), [0.0])


def test_nuclear_tw_open_channel_path_is_finite() -> None:
    """At an energy that DOES clear this tiny grid's (unconverged) `eps_e`
    threshold, the open-channel branch executes (not skipped by `continue`)
    and returns a finite (though not physically converged) value."""
    tw = _nuclear_tw_fixture()
    e_tot = 0.6 + EPS[V_INIT]
    assert e_tot - tw._eps_e[0] > 0.0  # sanity: confirms the open-channel premise
    s = tw.sigma(0.6)
    assert np.all(np.isfinite(s))
    assert np.all(s >= 0.0)


def test_nuclear_tw_rejects_free_reference() -> None:
    tw = _nuclear_tw_fixture()
    with pytest.raises(ValueError, match="no elastic free-reference"):
        tw.sigma(0.10, free=tw)


def test_nuclear_tw_n_channels_defaults_to_one() -> None:
    tw = TannorWeeks(
        TG, N2, EPS, CHI, V_INIT, [], NUCLEAR_WP_OUT, wp_in=WP_IN, dt=DT, axis="nuclear"
    )
    assert tw._n_channels == 1


@pytest.mark.slow
def test_nuclear_tw_da_converges_to_ti_oracle() -> None:
    """LOAD-BEARING: the nuclear-`TannorWeeks` sigma_DA converges to the TI
    `da_cross_section` on the SAME grid (a differential test -- both must
    agree even where the grid is not physically converged) -- a mirror of
    `test_nuclear_flux_da_converges_to_ti_oracle`/`test_nuclear_dirac_da_
    converges_to_ti_oracle` (see the former's docstring for the full config
    rationale: the launch-box electronic grid x the FINE eMoScat F2 nuclear
    deck, the off-box-incident and coarse-nuclear failure modes it guards
    against).

    `wp_out` is a NUCLEAR outgoing test packet: `r0_out=8.0` (position, placed
    INWARD of eMoScat's own R=9.7 so the dissociation wave reaches it well
    before this reduced grid's ECS edge at R0=10.7), `p0_out=72.0` (impulse ~
    K_R for the dissociation wave), `sigma_out=0.07` (thickness -- narrow in R
    means WIDE in K, so the packet's momentum distribution spans the K_R range
    of the probe energies; a wider-in-R packet is narrow-in-K and gives an
    ill-conditioned deconvolution -> spurious blow-up at the higher-E channels).
    Unlike `Flux`/`Dirac`'s fixed surface index, `TannorWeeks` PROPAGATES this
    Gaussian test packet against the trajectory.

    CONTROLLER-VALIDATED (2026-07-31): unlike `Flux`/`Dirac` (which plateau
    cleanly at sigma/sigma_ti ~0.86-0.97), `TannorWeeks` converges to the RIGHT
    MAGNITUDE (~1) but OSCILLATES (n=1750 ~[1.21,1.26,1.42], n=2000
    ~[0.55,1.39,1.41]) -- it is the noisiest, most test-packet-sensitive of the
    three (a propagated-Gaussian deconvolution, not a point-value/Wronskian read
    of the SAME trajectory). Hence a WIDER band `(0.4, 1.7)` at n=1750: the
    check is "converges to order ~1 (not ~1e-3 or ~1e2, which would flag a wrong
    `_C_DA` or a sign/prefactor error in the nuclear `eta_outgoing`)", not the
    tight plateau the flux/delta gates assert. Heavy (~86k unknowns x 1750
    steps, ~12 min) -- @slow. The full eMoScat 90-bohr-electronic grid is
    Docker/overnight-deferred.
    """
    from qscat.core.dissociation import da_cross_section
    from qscat.core.grids import segmented_grid
    from qscat.model import F2

    elec = electronic_grid(r_max=25.0, order=6, n_complex=3, angle_deg=40.0)
    nuc = segmented_grid(
        [(9, 1.8), (1, 2.0), (5, 2.5), (4, 2.596908), (4, 2.7), (40, 10.7)],
        [(1, 10.8), (1, 11.0), (1, 11.5), (1, 12.5), (1, 14.0), (1, 18.0), (4, 30.0), (2, 101.0)],
        angle_deg=35.0,
        quadrature=14,
    )
    tg = TensorGrid([elec, nuc])
    eps, chi = vibrational_states(nuc, F2.mu, 4, F2.v0)
    e_probe = np.array([0.03, 0.04])
    sigma_ti = np.ravel(da_cross_section(tg, F2, eps, chi, 0, e_probe))

    wp_in = {"r0": 12.0, "p0": -0.5, "sigma": 3.0}
    wp_out = {"r0_out": 8.0, "p0_out": 72.0, "sigma_out": 0.07}
    psi0 = initial_state(tg, chi[0], **wp_in)
    tw = TannorWeeks(
        tg, F2, eps, chi, 0, [], wp_out, wp_in=wp_in, dt=1.0, axis="nuclear", n_channels=1
    )
    # `propagate` records the t=0 state itself, then every step -> 1751 samples.
    propagate(tg, psi0, [], dt=1.0, n_steps=1750, hamiltonian=F2.hamiltonian(tg), extractors=[tw])
    ratio = np.ravel(tw.sigma(e_probe)) / sigma_ti
    # TW is the noisiest/most test-packet-sensitive method -> wider order-~1 band
    # (controller-measured ~[1.21,1.26,1.42] at n=1750); see the docstring.
    assert np.all(ratio > 0.4) and np.all(ratio < 1.7), (ratio, sigma_ti)


# --- `n_steps=` truncation (qscat-run Task 3) --------------------------------
#
# `apps/qscat-run`'s moment-resolved `cross_section_vs_time` artifact needs to
# read sigma(E) as of an EARLIER time from an already-completed propagation,
# without re-propagating -- these keyword-only `n_steps=` additions to
# `sigma`/`.result`/`.series` are that hook. `n_steps=None` (the default,
# used everywhere else in this file and by every pre-existing caller) MUST be
# byte-identical to the pre-addition behavior -- this file's golden/
# differential tests above (re-run unchanged, still passing) are the load-
# bearing confirmation of that; the tests below are a direct, minimal check
# of the new keyword itself: `n_steps=None` == the full series exactly,
# `n_steps=k < full` differs.


def _fresh_extractors() -> tuple[TannorWeeks, Dirac, Flux]:
    tw = TannorWeeks(TG, N2, EPS, CHI, V_INIT, VPRIMES, WP_OUT, wp_in=WP_IN, dt=DT)
    dirac = Dirac(TG, N2, EPS, CHI, V_INIT, VPRIMES, POSITION, wp_in=WP_IN, dt=DT)
    flux = Flux(TG, N2, EPS, CHI, V_INIT, VPRIMES, POSITION, wp_in=WP_IN, dt=DT)
    psi0 = initial_state(TG, CHI[V_INIT], **WP_IN)
    propagate(
        TG,
        psi0,
        [],
        dt=DT,
        n_steps=N_STEPS,
        hamiltonian=N2.hamiltonian(TG),
        extractors=[tw, dirac, flux],
    )
    return tw, dirac, flux


@pytest.mark.parametrize("extractor_index", [0, 1, 2], ids=["tw", "dirac", "flux"])
def test_sigma_n_steps_none_is_byte_identical_to_full(extractor_index: int) -> None:
    ext = _fresh_extractors()[extractor_index]
    full = ext.sigma([0.10, 0.15])
    explicit_none = ext.sigma([0.10, 0.15], n_steps=None)
    np.testing.assert_array_equal(full, explicit_none)
    # n_steps equal to the full recorded sample count is the SAME truncation
    # as n_steps=None (propagate records N_STEPS+1 samples, t=0 included).
    explicit_full_count = ext.sigma([0.10, 0.15], n_steps=N_STEPS + 1)
    np.testing.assert_array_equal(full, explicit_full_count)


@pytest.mark.parametrize("extractor_index", [0, 1, 2], ids=["tw", "dirac", "flux"])
def test_sigma_n_steps_truncated_differs_from_full(extractor_index: int) -> None:
    ext = _fresh_extractors()[extractor_index]
    full = ext.sigma([0.10, 0.15])
    truncated = ext.sigma([0.10, 0.15], n_steps=3)
    assert not np.allclose(full, truncated)
    assert np.all(np.isfinite(truncated))
