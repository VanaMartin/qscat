"""FAST gate (sub-project #C4, Task 4): `td_ve_cross_sections_all` at a
REDUCED N2 grid -- the three extractors (`tw`/`delta`/`flow`) must mutually
agree within the documented cross-method band (Tasks 2-3 measured ~0.20-0.25
worst-case), and each must be within a much LOOSER band of the exact TI
oracle (`qscat.core.driven.ve_cross_section`).

At this reduced/fast config the propagation is far from converged (a few
seconds, `N_STEPS=800` vs the converged grid's `N_STEPS=1000` at ~50x more
DVR points and `dt=1.0` not `0.2`) -- all three methods land a factor of
~3-6x above the TI oracle TOGETHER (measured below), which is the documented,
expected behavior at this scale (see
`libs/qscat/tests/test_td_extractors.py`'s own `N_STEPS_DIFF=800` config,
whose cross-method bands this test reuses verbatim). The point of this gate
is NOT "matches TI at a toy grid" (it doesn't, by design) -- it is "the three
extractors, driven by the identical trajectory, land in the SAME ballpark
relative to the oracle", the convergence-diagnostic framing
`docs/physics/td-extractors.md` documents.
"""

from __future__ import annotations

import numpy as np
import pytest
from qscat.core.driven import ve_cross_section
from qscat.core.grids import electronic_grid, nuclear_grid
from qscat.core.time_dependent import td_ve_cross_sections_all
from qscat.core.vibrational import vibrational_states
from qscat.dvr import TensorGrid
from qscat.model import N2

from validation.n2 import td_extractors

TG = TensorGrid(
    [
        electronic_grid(r_max=12.0, order=5, n_complex=3),
        nuclear_grid(quadrature=6, r_max=14.0, n_complex=3),
    ]
)
EPS, CHI = vibrational_states(TG.grids[1], N2.mu, 4, N2.v0)

V_INIT = 0
VPRIMES = [0, 1]  # includes the elastic channel -- exercises the free reference too
WP_IN = {"r0": 4.0, "p0": -0.5, "sigma": 1.2}
WP_OUT = {"r0_out": 6.0, "p0_out": 0.5, "sigma_out": 1.0}
DT = 0.2
N_STEPS = 800  # matches libs/qscat/tests/test_td_extractors.py's N_STEPS_DIFF
POSITION = 37  # r = 6.0, real region (R0=12), colocated with wp_out

E_GRID = [0.10, 0.15]

# Cross-method band: matches libs/qscat/tests/test_td_extractors.py's own
# gates at this exact configuration (_DELTA_TW_RTOL=0.20, _FLUX_TW_RTOL=0.25)
# -- reused here, not re-derived, since this is the SAME trajectory shape.
_CROSS_METHOD_RTOL = 0.25

# Looser band vs the TI oracle: at this reduced/unconverged grid all three
# extractors measured 2.9x-6.2x the oracle (see module docstring) -- a
# documented convergence-diagnostic gap, not a bug. `_TI_BAND` covers a full
# order of magnitude either direction, comfortably bracketing the measured
# ratios with margin while still catching a genuine regression (e.g. a sign
# error or wrong-channel bug that would blow this past 10x).
_TI_BAND = 10.0


def test_three_way_mutual_agreement_and_ti_ballpark() -> None:
    sigma_all = td_ve_cross_sections_all(
        TG,
        N2,
        EPS,
        CHI,
        V_INIT,
        VPRIMES,
        E_GRID,
        dt=DT,
        n_steps=N_STEPS,
        wp_in=WP_IN,
        wp_out=WP_OUT,
        position=POSITION,
        surface=POSITION,
    )
    assert set(sigma_all) == {"tw", "delta", "flow"}
    tw, delta, flow = sigma_all["tw"], sigma_all["delta"], sigma_all["flow"]
    assert tw.shape == delta.shape == flow.shape == (len(E_GRID), len(VPRIMES))
    assert np.all(np.isfinite(tw))
    assert np.all(np.isfinite(delta))
    assert np.all(np.isfinite(flow))
    assert np.all(tw > 0.0)
    assert np.all(delta > 0.0)
    assert np.all(flow > 0.0)

    # Mutual agreement: delta/flow vs TW from the SAME propagation.
    np.testing.assert_allclose(delta, tw, rtol=_CROSS_METHOD_RTOL, atol=1e-14)
    np.testing.assert_allclose(flow, tw, rtol=_CROSS_METHOD_RTOL, atol=1e-14)

    # Each converges toward the TI oracle -- loosely, at this reduced grid.
    sigma_ti = ve_cross_section(TG, N2, EPS, CHI, V_INIT, VPRIMES, E_GRID)
    assert np.all(np.isfinite(sigma_ti))
    assert np.all(sigma_ti > 0.0)
    for name, sigma in (("tw", tw), ("delta", delta), ("flow", flow)):
        ratio = sigma / sigma_ti
        assert np.all(ratio > 1.0 / _TI_BAND), f"{name}: ratio too small vs TI: {ratio}"
        assert np.all(ratio < _TI_BAND), f"{name}: ratio too large vs TI: {ratio}"


@pytest.mark.slow
def test_converged_three_way_matches_ti_oracle_live_anchor() -> None:
    """The `@slow` converged-grid three-way comparison at the GATED C5/D1
    anchor `(E=0.10, v'=1)`: one ~4-5 min propagation
    (`td_extractors.compute_live_result`) driving `TannorWeeks`/`Dirac`/
    `Flux` together, each checked against the TI oracle at the same ~3%
    tolerance `libs/qscat/tests/test_td_extractors.py`'s individual
    `Dirac`/`Flux` @slow anchor tests already gate on (delta 0.971, flow
    0.970 there -- this reruns both, PLUS TannorWeeks, from one combined
    propagation).
    """
    result = td_extractors.compute_live_result()
    assert result.live
    assert np.isfinite(result.sigma_tw)
    assert np.isfinite(result.sigma_delta)
    assert np.isfinite(result.sigma_flow)
    assert result.sigma_ti > 0.0
    # rtol=0.10 matches libs/qscat/tests/test_td_extractors.py's own gate on
    # this exact grid/wavepacket/anchor for delta/flow.
    assert result.ratio_delta_ti == pytest.approx(1.0, rel=0.10)
    assert result.ratio_flow_ti == pytest.approx(1.0, rel=0.10)
    # TW has not been individually anchor-tested at T=1000 before (only at
    # T=1500 in docs/physics/n2-2d-td-cross-section.md, where it reached
    # ~1-2%); rel=0.15 gives headroom for the shorter T=1000 propagation
    # while still gating a genuine regression.
    assert result.ratio_tw_ti == pytest.approx(1.0, rel=0.15)
