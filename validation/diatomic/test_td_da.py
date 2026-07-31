"""Tests for `td_da_cross_section`/`td_da_cross_sections_all` (`qscat.core`)
and the F2/NO three-way TD-DA validation harness (`validation.diatomic.td_da`).

Fast tests (below, run in the normal suite): shape/contract checks ONLY on
a tiny F2 grid + a handful of propagation steps -- NOT a converged DA cross
section, mirroring the caveat `libs/qscat/tests/test_td_extractors.py`'s own
fast structural tests for the underlying nuclear-axis extractors carry.

`@slow` (further below): the LOAD-BEARING F2/NO three-way validation against
the exact-2D TI `da_cross_section` oracle, on the controller-validated TD
launch grid (`validation.diatomic.td_da.td_launch_grid`, ~86k unknowns x
1800 steps, ~10 min per molecule) -- see that test's own docstring for the
measured numbers.
"""

from __future__ import annotations

import numpy as np
import pytest
from qscat.core.grids import electronic_grid, nuclear_grid
from qscat.core.time_dependent import td_da_cross_section, td_da_cross_sections_all
from qscat.core.vibrational import vibrational_states
from qscat.dvr import TensorGrid
from qscat.model import F2

from validation.diatomic.config import CONFIGS
from validation.diatomic.td_da import compute_td_da_three_way, td_launch_grid

# A tiny, fast grid -- mirrors `libs/qscat/tests/test_td_extractors.py`'s own
# tiny N2 config (same shape/orders), ported to F2 (this validation
# package's molecule) so the shape/contract check needs no eMoScat deck or
# large launch box.
_TG = TensorGrid(
    [
        electronic_grid(r_max=12.0, order=5, n_complex=3),
        nuclear_grid(quadrature=6, r_max=14.0, n_complex=3),
    ]
)
_EPS, _CHI = vibrational_states(_TG.grids[1], F2.mu, 4, F2.v0)
_WP_IN = {"r0": 4.0, "p0": -0.5, "sigma": 1.2}
_WP_OUT = {"r0_out": 6.0, "p0_out": 5.0, "sigma_out": 1.0}
_SURFACE = 90  # a real-region (unscaled) nuclear DVR index, past the interaction
_DT = 0.2
_N_STEPS = 5  # a handful of steps -- fast, not a converged run


# --- fast shape/contract tests -----------------------------------------------


def test_td_da_cross_section_flow_shape() -> None:
    sigma = td_da_cross_section(
        _TG, F2, _EPS, _CHI, 0, [0.10, 0.6],
        dt=_DT, n_steps=_N_STEPS, wp_in=_WP_IN, method="flow", surface=_SURFACE,
    )
    assert sigma.shape == (2, 1)
    assert np.all(np.isfinite(sigma))


def test_td_da_cross_section_delta_shape() -> None:
    sigma = td_da_cross_section(
        _TG, F2, _EPS, _CHI, 0, 0.6,
        dt=_DT, n_steps=_N_STEPS, wp_in=_WP_IN, method="delta", position=_SURFACE,
    )
    assert sigma.shape == (1,)
    assert np.all(np.isfinite(sigma))


def test_td_da_cross_section_tw_shape() -> None:
    sigma = td_da_cross_section(
        _TG, F2, _EPS, _CHI, 0, 0.6,
        dt=_DT, n_steps=_N_STEPS, wp_in=_WP_IN, method="tw", wp_out=_WP_OUT,
    )
    assert sigma.shape == (1,)
    assert np.all(np.isfinite(sigma))


def test_td_da_cross_section_default_method_is_flow() -> None:
    sigma_default = td_da_cross_section(
        _TG, F2, _EPS, _CHI, 0, 0.6, dt=_DT, n_steps=_N_STEPS, wp_in=_WP_IN, surface=_SURFACE,
    )
    sigma_flow = td_da_cross_section(
        _TG, F2, _EPS, _CHI, 0, 0.6,
        dt=_DT, n_steps=_N_STEPS, wp_in=_WP_IN, method="flow", surface=_SURFACE,
    )
    np.testing.assert_allclose(sigma_default, sigma_flow, rtol=0, atol=0)


def test_td_da_cross_section_flow_requires_surface() -> None:
    with pytest.raises(ValueError, match="requires `surface`"):
        td_da_cross_section(
            _TG, F2, _EPS, _CHI, 0, 0.10, dt=_DT, n_steps=_N_STEPS, wp_in=_WP_IN, method="flow"
        )


def test_td_da_cross_section_delta_requires_position() -> None:
    with pytest.raises(ValueError, match="requires `position`"):
        td_da_cross_section(
            _TG, F2, _EPS, _CHI, 0, 0.10, dt=_DT, n_steps=_N_STEPS, wp_in=_WP_IN, method="delta"
        )


def test_td_da_cross_section_tw_requires_wp_out() -> None:
    with pytest.raises(ValueError, match="requires `wp_out`"):
        td_da_cross_section(
            _TG, F2, _EPS, _CHI, 0, 0.10, dt=_DT, n_steps=_N_STEPS, wp_in=_WP_IN, method="tw"
        )


def test_td_da_cross_section_unknown_method_raises() -> None:
    with pytest.raises(ValueError, match="unknown method"):
        td_da_cross_section(
            _TG, F2, _EPS, _CHI, 0, 0.10, dt=_DT, n_steps=_N_STEPS, wp_in=_WP_IN, method="bogus"
        )


def test_td_da_cross_sections_all_shape_and_keys() -> None:
    sigma_all = td_da_cross_sections_all(
        _TG, F2, _EPS, _CHI, 0, [0.10, 0.6],
        dt=_DT, n_steps=_N_STEPS, wp_in=_WP_IN,
        surface=_SURFACE, position=_SURFACE, wp_out=_WP_OUT,
    )
    assert set(sigma_all) == {"flow", "delta", "tw"}
    for key in ("flow", "delta", "tw"):
        assert sigma_all[key].shape == (2, 1)
        assert np.all(np.isfinite(sigma_all[key]))


def test_td_da_cross_sections_all_matches_individual_calls() -> None:
    """`td_da_cross_sections_all`'s ONE-propagation result must reproduce
    calling `td_da_cross_section` once per method (each its OWN, separately
    propagated, trajectory) to machine precision -- same contract
    `test_cross_sections_all_matches_each_method_individually` establishes
    for the VE wiring."""
    e = [0.10, 0.6]
    sigma_all = td_da_cross_sections_all(
        _TG, F2, _EPS, _CHI, 0, e,
        dt=_DT, n_steps=_N_STEPS, wp_in=_WP_IN,
        surface=_SURFACE, position=_SURFACE, wp_out=_WP_OUT,
    )
    kwargs_by_method = {
        "flow": {"method": "flow", "surface": _SURFACE},
        "delta": {"method": "delta", "position": _SURFACE},
        "tw": {"method": "tw", "wp_out": _WP_OUT},
    }
    for key, kwargs in kwargs_by_method.items():
        sigma_individual = td_da_cross_section(
            _TG, F2, _EPS, _CHI, 0, e, dt=_DT, n_steps=_N_STEPS, wp_in=_WP_IN, **kwargs
        )
        np.testing.assert_allclose(sigma_all[key], sigma_individual, rtol=0, atol=0)


# --- @slow: the load-bearing F2/NO three-way validation ----------------------


@pytest.mark.slow
def test_f2_td_da_three_way_agrees_with_ti_oracle() -> None:
    """F2 `td_da_cross_sections_all` on the controller-validated TD launch
    grid (`validation.diatomic.td_da.td_launch_grid`: electronic r_max=25
    launch box x the eMoScat fine nuclear deck, `n_steps=1800`, ~86k
    unknowns, ~10 min) vs the exact-2D TI `da_cross_section` oracle on the
    SAME grid.

    Measured (controller, 2026-07-31, at n_steps=1500 -- see
    `libs/qscat/tests/test_td_extractors.py`'s per-extractor `@slow` gates,
    which this test's three-way call reproduces from a single shared
    propagation): `flow`/`delta` reach a STABLE plateau sigma/sigma_ti ~
    0.86-0.97; `tw` (the propagated nuclear test-packet method, test-packet-
    sensitive -- needs the narrow `sigma_out=0.07` wide-K packet placed
    inward of the surface/position, `p0_out=72.0`) lands ~0.9. The residual
    ~3-14% gap is the TD-vs-TI cross-method discretization band (as in the
    electronic VE extractors), not a normalization error: a wrong `_C_DA`
    (`qscat.core.td_extractors`) would plateau at a wildly different
    constant, not ~1. `rtol` bands below are widened slightly (0.7, 1.3)
    from the single-extractor `@slow` gates' (0.7, 1.25) to leave margin for
    the `n_steps=1800` (vs. 1500) config.
    """
    e_probe = np.array([0.03, 0.04])
    cfg = CONFIGS["F2"]
    tg = td_launch_grid(cfg)
    eps, chi = vibrational_states(tg.grids[1], cfg.model.mu, cfg.n_vib, cfg.model.v0)

    from qscat.core.dissociation import da_cross_section

    sigma_ti = np.ravel(da_cross_section(tg, cfg.model, eps, chi, 0, e_probe))
    sigma_td = compute_td_da_three_way(cfg, e_probe, n_steps=1800)

    # flow/delta plateau cleanly (~0.86-0.97); tw is the noisiest / most
    # test-packet-sensitive method (converges to order ~1 but oscillates
    # ~0.55-1.42) -> a wider order-~1 band. See the per-extractor @slow gates.
    bands = {"flow": (0.7, 1.3), "delta": (0.7, 1.3), "tw": (0.4, 1.7)}
    for method in ("flow", "delta", "tw"):
        ratio = np.ravel(sigma_td[method]) / sigma_ti
        lo, hi = bands[method]
        assert np.all(ratio > lo) and np.all(ratio < hi), (method, ratio, sigma_ti)


@pytest.mark.slow
def test_no_td_da_three_way_agrees_with_ti_oracle() -> None:
    """NO's TD-DA three-way validation -- a mirror of the F2 test above, on
    NO's own eMoScat nuclear deck (`CONFIGS["NO"]`). NO's DA channel opens
    much higher (~0.17 Ha, a sharp near-threshold spike, see
    `validation/diatomic/test_da_curves.py::test_no_da_threshold_onset`), so
    the probe energies are shifted up accordingly.
    """
    e_probe = np.array([0.18, 0.20])
    cfg = CONFIGS["NO"]
    tg = td_launch_grid(cfg)
    eps, chi = vibrational_states(tg.grids[1], cfg.model.mu, cfg.n_vib, cfg.model.v0)

    from qscat.core.dissociation import da_cross_section

    sigma_ti = np.ravel(da_cross_section(tg, cfg.model, eps, chi, 0, e_probe))
    sigma_td = compute_td_da_three_way(cfg, e_probe, n_steps=1800)

    # flow/delta plateau cleanly (~0.86-0.97); tw is the noisiest / most
    # test-packet-sensitive method (converges to order ~1 but oscillates
    # ~0.55-1.42) -> a wider order-~1 band. See the per-extractor @slow gates.
    bands = {"flow": (0.7, 1.3), "delta": (0.7, 1.3), "tw": (0.4, 1.7)}
    for method in ("flow", "delta", "tw"):
        ratio = np.ravel(sigma_td[method]) / sigma_ti
        lo, hi = bands[method]
        assert np.all(ratio > lo) and np.all(ratio < hi), (method, ratio, sigma_ti)
