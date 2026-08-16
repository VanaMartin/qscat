from __future__ import annotations

import numpy as np
from qscat.model import F2, N2
from qscat.tuning import interaction_region


def test_f2_interaction_region_brackets_the_coupling() -> None:
    R_lo, R_hi = interaction_region(F2)
    assert 0.5 < R_lo < 2.5 and 3.0 < R_hi < 8.0  # ~[1.5,4]-ish, where lambda(R) is significant
    assert R_lo < R_hi


def test_region_is_where_vint_is_still_transitioning() -> None:
    # N2's lambda(R) saturates to a substantial (non-zero) plateau on BOTH
    # sides -- verified numerically, see resonance.py's module docstring --
    # so "outside [R_lo, R_hi] is small" (the F2-style check) does not hold
    # here. What DOES hold: past R_hi, s(R) has essentially stopped
    # changing (it's within a small fraction of its peak of its own
    # asymptote) -- the property the downstream sampler actually relies on
    # to freeze at a single far value.
    R_lo, R_hi = interaction_region(N2, frac=0.05)
    r = np.linspace(0.1, 15, 200)
    s = lambda R: np.max(np.abs(np.real(N2.v_int(r[:, None], np.array([[R]])))))
    peak = max(s(R) for R in np.linspace(0.5, 6, 60))
    assert abs(s(R_hi + 1.0) - s(R_hi)) < 0.05 * peak * 1.5  # already saturated
    assert R_lo < R_hi


def _elec_grids():
    from qscat.core.grids import electronic_grid

    return (
        electronic_grid(r_max=16.0, order=7, n_complex=6, angle_deg=35.0),
        electronic_grid(r_max=16.0, order=7, n_complex=6, angle_deg=44.0),
    )


def test_resonance_curve_dense_interaction_sparse_far():
    from qscat.model import F2
    from qscat.tuning import interaction_region, resonance_curve

    ga, gb = _elec_grids()
    R_lo, R_hi = interaction_region(F2)
    R, Vd, G = resonance_curve(F2, ga, gb, R_max=22.0, n_dense=20)
    # most samples land inside the interaction region; the far region is sparse (~1 pt near R_max)
    inside = (R >= R_lo) & (R <= R_hi)
    far = R > R_hi + 1.0
    assert inside.sum() >= 15  # dense inside
    assert far.sum() <= 3  # sparse far
    assert np.all(np.isfinite(Vd)) and np.all(G >= 0.0)
    # Gamma peaks inside the interaction region (the resonance), ~0 far
    far_peak = G[far].max() if far.any() else 0.0
    assert G[inside].max() > 10 * far_peak + 1e-12
