"""BO quasi-bound levels in the H2+ Rydberg curves."""

from __future__ import annotations

import numpy as np
import pytest
from qscat.exceptions import GridError
from qscat.model import H2P

from validation.h2plus.config import full_grid, proxy_grid
from validation.h2plus.reference_levels import EPS0, LEVELS
from validation.h2plus.rydberg_levels import rydberg_levels

# n_vib=3 on the proxy grid: curve 0 (n=0) is a
# MUCH deeper well than curves 1, 2 -- it plunges from ~-0.26 Ha near R~2 (the
# ion equilibrium) down to a shallow, ~0.006 Ha-deep minimum near R~6.9 before
# flattening at ~-1.385 Ha for large R (`v_int` does not vanish as R -> inf;
# `v_dr_diag`'s explicit `- V_int(r, R_inf)` subtraction elsewhere in the
# codebase already documents this). That flat, shallow well genuinely supports
# only 3 numerically clean (|Im(E)| < 1e-6 Ha) bound vibrational levels on this
# proxy grid before the 4th picks up ~1e-5 Ha of contamination -- a real
# capacity limit of this curve+grid, not a bug. curves 1/2 (the actual
# near-threshold Rydberg series) support >=12 clean levels each on the same
# grid. Against `full_grid()`'s nuclear grid the same split holds with curve 0
# reaching 5 rather than 3 -- see `rydberg_levels`' `allow_partial`.
_N_VIB_PROXY = 3


def test_curves_are_ordered_and_below_the_ion() -> None:
    tg = proxy_grid()
    res = rydberg_levels(H2P, tg.grids[0], tg.grids[1], n_curves=3, n_vib=_N_VIB_PROXY)
    assert res.curves.shape == (3, tg.grids[1].n)
    # Rydberg curves are ordered in energy at every R.
    real = np.abs(tg.grids[1].points.imag) < 1e-12
    ordered = np.diff(res.curves[:, real].real, axis=0)
    assert np.all(ordered > 0.0)


def test_levels_are_real_and_ascending() -> None:
    tg = proxy_grid()
    res = rydberg_levels(H2P, tg.grids[0], tg.grids[1], n_curves=2, n_vib=_N_VIB_PROXY)
    assert res.energies.shape == (2, _N_VIB_PROXY)
    assert np.all(np.isreal(res.energies))
    assert np.all(np.diff(res.energies, axis=1) > 0.0)  # vibrational ladder


def test_allow_partial_pads_the_curve_that_cannot_supply_n_vib() -> None:
    """A count curve 0 cannot meet must raise by default, and must come back
    NaN-padded under `allow_partial` -- with the curves that CAN meet it
    still full. Without the flag the only way to include curves 1-2's upper
    levels would be a uniform `n_vib` that curve 0 rejects, and the only way
    to satisfy curve 0 would be to truncate curves 1-2; the DR windows are
    built out of exactly those upper levels, so both are wrong answers.
    """
    tg = proxy_grid()
    g_r, g_R = tg.grids[0], tg.grids[1]
    too_many = _N_VIB_PROXY + 3

    with pytest.raises(GridError):
        rydberg_levels(H2P, g_r, g_R, n_curves=3, n_vib=too_many)

    res = rydberg_levels(H2P, g_r, g_R, n_curves=3, n_vib=too_many, allow_partial=True)
    assert res.energies.shape == (3, too_many)

    # Curve 0 is the shallow one: some levels present, the rest padded.
    finite_0 = np.isfinite(res.energies[0])
    assert finite_0.any() and not finite_0.all()
    # Padding is a suffix, never a hole in the middle of a ladder.
    assert np.all(np.diff(finite_0.astype(int)) <= 0)
    # The levels curve 0 does supply are the same ones the strict call gives.
    strict = rydberg_levels(H2P, g_r, g_R, n_curves=3, n_vib=_N_VIB_PROXY)
    np.testing.assert_allclose(res.energies[0, :_N_VIB_PROXY], strict.energies[0])

    # Curves 1-2 are deep enough to fill the whole row.
    assert np.all(np.isfinite(res.energies[1:]))


@pytest.mark.slow
def test_levels_match_the_published_omega_table() -> None:
    """Every level, against the author's own published `omega_i^j` table.

    Stronger than `test_curve_asymptotes_match_table_4_1`, which pins only
    the three `R -> infinity` plateaus: this compares the full vibrational
    ladder inside each curve -- 5 levels for Ry_0 plus 12 each for
    Ry_1..Ry_4 -- so it catches a nuclear-grid or reduced-mass error that
    leaves the electronic asymptotes untouched.

    This is the well-conditioned way to validate this model against the
    published data. The cross section itself cannot do the job: its
    resonances are ~2e-5 Ha wide and the published sweep samples at 1e-5 Ha,
    so a few-uHa position difference swings sigma by tens of percent (see
    docs/physics/h2plus-dr.md). The levels are where the same agreement is
    directly readable.

    `atol=1e-5` Ha: measured agreement is <=4e-6 across all 53 levels, and
    1e-5 leaves headroom without being loose enough to hide a regression --
    it is still a fifth of a typical resonance width.

    Uses the FULL nuclear grid (`full_grid().grids[1]`), not the proxy's,
    because the published levels are converged ones; the electronic side
    stays on the proxy grid, which `rydberg_levels`' docstring records as
    sufficient for curves 0-4. Runtime is a few minutes -- 818 nuclear
    points x one dense electronic eigensolve each.
    """
    res = rydberg_levels(
        H2P,
        proxy_grid().grids[0],
        full_grid().grids[1],
        n_curves=len(LEVELS),
        n_vib=max(len(r) for r in LEVELS),
        allow_partial=True,
    )
    worst = 0.0
    for j, published in enumerate(LEVELS):
        mine = res.energies[j, : len(published)] - EPS0  # to the electron-energy frame
        assert np.all(np.isfinite(mine)), f"curve {j} supplied fewer levels than published"
        np.testing.assert_allclose(mine, published, atol=1e-5)
        worst = max(worst, float(np.max(np.abs(mine - np.asarray(published)))))
    # Guard the headroom claim above: if agreement ever degrades toward the
    # tolerance, this fails while there is still margin to investigate.
    assert worst < 6e-6, f"level agreement degraded to {worst:.2e} Ha"


@pytest.mark.slow
def test_curve_asymptotes_match_table_4_1() -> None:
    """The large-R plateau of curves Ry_0/Ry_1/Ry_2 against the thesis's
    Table 4.1 (Vana 2017, p. 62), the published asymptotic electron energy
    spectrum for the same fixed-nuclei electronic problem. This is the
    strongest evidence the curve construction (looping `eigen` per nuclear
    point) is right -- an external, independently-published oracle, not
    just an internal self-consistency check -- so it must be a gate, not
    docstring prose.

    `proxy_grid()`'s real nuclear region reaches R=14 bohr, well past where
    Fig. 4.2 shows the curves have already flattened (by R~8-12), so the
    curve value at the grid's largest real R is a fair stand-in for the
    R->infinity asymptote. Runtime ~35-40s (510 nuclear points x the
    377-point electronic grid, one dense electronic eigensolve per point) --
    kept off the default (non-slow) run for that reason.
    """
    tg = proxy_grid()
    g_r, g_R = tg.grids[0], tg.grids[1]
    res = rydberg_levels(H2P, g_r, g_R, n_curves=3, n_vib=1)

    real = np.abs(g_R.points.imag) < 1e-12
    real_idx = np.flatnonzero(real)
    r_max_idx = real_idx[np.argmax(g_R.points[real_idx].real)]

    # Table 4.1: Ry_0, Ry_1, Ry_2 asymptotic electron energies (Ha).
    published = np.array([-1.38492776, -0.12499996, -0.05481037])
    measured = res.curves[:, r_max_idx].real
    # Measured agreement on this deck is ~5e-5; 2e-4
    # leaves headroom without being loose enough to hide a real regression.
    np.testing.assert_allclose(measured, published, atol=2e-4)
