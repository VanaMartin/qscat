from __future__ import annotations

import numpy as np
from qscat.tuning import (
    analyze_potential,
    equidistribution_elements,
    optimal_real_mesh,
    order_for_wavenumber,
    refine_elements_in_window,
)
from qscat.tuning.mesh import _clamp_lengths_span_preserving


def _flat_profile(k_const, x_max, m):
    E = k_const**2 / (2 * m)
    return analyze_potential(lambda x: np.zeros_like(np.asarray(x, float)), 0.0, x_max, m, E)


def test_uniform_k_gives_uniform_elements():
    # constant k -> equidistribution => (nearly) uniform element lengths
    p = _flat_profile(2.0, 20.0, 1.0)
    els = equidistribution_elements(p, 8, phase_per_element=5.0, min_len=1e-3, max_len=10.0)
    assert np.std(els) / np.mean(els) < 0.05  # ~uniform
    assert abs(2.0 * np.mean(els) - 5.0) < 0.3  # k*h ~ phase_per_element


def test_finer_where_k_larger():
    # step: k large for x<5, small for x>5 -> smaller elements in the fast region
    def V(x):
        x = np.asarray(x, float)
        return np.where(x < 5.0, -2.0, 0.0)

    p = analyze_potential(V, 0.0, 20.0, 1.0, 0.5)
    els = equidistribution_elements(p, 8, phase_per_element=5.0, min_len=1e-3, max_len=10.0)
    # reconstruct boundaries; mean element length in [0,5] < mean in [5,20]
    b = np.concatenate([[0.0], np.cumsum(els)])
    lo = np.diff(b)[b[:-1] < 5.0].mean()
    hi = np.diff(b)[b[:-1] >= 5.0].mean()
    assert lo < hi


def test_optimal_mesh_minimizes_points():
    p = _flat_profile(2.0, 20.0, 1.0)
    els, order = optimal_real_mesh(p, phase_coeff=1.0, min_len=1e-3, max_len=10.0)
    assert order in (6, 8, 10, 14)
    assert len(els) * (order - 1) <= 400  # a sane point budget for this easy case


def test_uniform_case_preserves_domain_span():
    # Locks the general invariant: elements must tile the domain exactly,
    # regardless of whether min/max_len clamping fires.
    p = _flat_profile(2.0, 20.0, 1.0)
    els = equidistribution_elements(p, 8, phase_per_element=5.0, min_len=1e-3, max_len=10.0)
    assert abs(sum(els) - (p.x[-1] - p.x[0])) < 1e-9


def test_oversized_element_is_subdivided_not_dropped():
    # k tiny + phase_per_element large -> a single natural element spans the
    # whole domain (20), far past max_len=10: must be SUBDIVIDED, not
    # clipped (clipping to 10 would silently drop half the domain).
    p = _flat_profile(0.1, 20.0, 1.0)
    span = p.x[-1] - p.x[0]
    els = equidistribution_elements(p, 8, phase_per_element=5.0, min_len=1e-3, max_len=10.0)
    assert abs(sum(els) - span) < 1e-9
    assert all(e <= 10.0 + 1e-9 for e in els)
    assert all(e >= 1e-3 - 1e-9 for e in els)


def test_undersized_elements_are_merged_not_inflated():
    # k large + phase_per_element small -> many natural elements far below
    # min_len=1.0: must be MERGED forward, not clipped up (clipping each of
    # ~400 elements up to 1.0 would inflate the total span ~20x).
    p = _flat_profile(10.0, 20.0, 1.0)
    span = p.x[-1] - p.x[0]
    els = equidistribution_elements(p, 8, phase_per_element=0.5, min_len=1.0, max_len=10.0)
    assert abs(sum(els) - span) < 1e-9
    assert all(e >= 1.0 - 1e-9 for e in els)
    assert all(e <= 10.0 + 1e-9 for e in els)


def test_order_for_wavenumber_resolves_fast_wave():
    # K_exit=78, lambda~0.08; at 0.15-bohr elements order 6 gives 3.2 ppw
    # (too few), 14 gives 7.5 -- the F2 exit-wave anchor from the design
    # brief.
    assert order_for_wavenumber(78.0, 0.15, target_ppw=6.0) == 14
    assert order_for_wavenumber(5.0, 0.15, target_ppw=6.0) <= 8  # slow wave, low order ok


def test_order_for_wavenumber_nonpositive_k_returns_cheapest_order():
    assert order_for_wavenumber(0.0, 0.15) == 6
    assert order_for_wavenumber(-1.0, 0.15) == 6


def test_refine_elements_in_window_only_local_and_span_preserving():
    lengths = [0.5] * 10  # span 5.0, from x_min=0
    out = refine_elements_in_window(lengths, 0.0, 2.0, 3.0, 0.1)
    assert abs(sum(out) - 5.0) < 1e-12  # span preserved
    assert max(out) <= 0.5 + 1e-12  # nothing coarsened
    # elements straddling [2,3] are now <=0.1; those outside stay 0.5
    assert min(out) <= 0.1 + 1e-12


def test_refine_elements_in_window_elements_outside_untouched():
    lengths = [0.5] * 10
    out = refine_elements_in_window(lengths, 0.0, 2.0, 3.0, 0.1)
    # the elements before/after the refined window are exactly the
    # original 0.5-length ones, unchanged.
    n_refined = sum(1 for h in out if h < 0.5 - 1e-9)
    assert n_refined > 0
    assert all(abs(h - 0.5) < 1e-12 or h <= 0.1 + 1e-12 for h in out)


def test_refine_elements_in_window_noop_when_already_fine():
    lengths = [0.05] * 10  # already below target_len everywhere
    out = refine_elements_in_window(lengths, 0.0, 2.0, 3.0, 0.1)
    assert out == lengths


def test_merge_never_breaches_max_len_hard_cap():
    # A small element (0.5) sitting right next to a near-max_len element
    # (9.9): gluing them (0.5 + 9.9 = 10.4) would breach the HARD max_len
    # cap of 10.0. The merge must instead emit them separately -- max_len
    # is inviolable, min_len is only a best-effort floor, so we deliberately
    # do NOT assert every element >= min_len here (0.5 stays below it).
    out = _clamp_lengths_span_preserving(np.array([0.5, 9.9]), min_len=1.0, max_len=10.0)
    assert all(e <= 10.0 + 1e-9 for e in out)
    assert abs(sum(out) - 10.4) < 1e-9
