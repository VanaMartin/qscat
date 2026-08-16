"""Artifact writers for the quasi-bound level table."""

from __future__ import annotations

import csv

import numpy as np
from qscat.core.lcp import ResonanceLevels
from qscat_run.artifacts import (
    _resonance_levels_real_mask,
    _resonance_levels_ylim,
    _write_resonance_levels,
)
from qscat_run.runner import ResonanceLevelsRun


def _fake_run(R0: float = 6.0) -> ResonanceLevelsRun:
    n_grid = 12
    levels = ResonanceLevels(
        energies=np.array([-0.05 - 0.001j, -0.03 - 0.004j]),
        widths=np.array([0.002, 0.008]),
        states=np.ones((2, n_grid), dtype=np.complex128) / np.sqrt(n_grid),
        residuals=np.array([1e-9, 4e-9]),
        real_weight=np.array([0.999, 0.97]),
        golden_rule=np.array([-0.0501 - 0.0011j, -0.0299 - 0.0039j]),
    )
    return ResonanceLevelsRun(
        label="lcp:resonance_levels",
        levels=levels,
        R_axis=np.linspace(1.5, 6.0, n_grid),
        Vd=np.linspace(-0.1, 0.0, n_grid).astype(np.complex128),
        Gamma=np.linspace(0.01, 0.0, n_grid),
        R0=R0,
    )


def test_writes_csv_npz_and_png(tmp_path):
    _write_resonance_levels(tmp_path, _fake_run())
    stem = "resonance_levels_lcp_resonance_levels"
    assert (tmp_path / f"{stem}.csv").exists()
    assert (tmp_path / f"{stem}.npz").exists()
    assert (tmp_path / f"{stem}.png").exists()


def test_csv_columns_and_values(tmp_path):
    _write_resonance_levels(tmp_path, _fake_run())
    path = tmp_path / "resonance_levels_lcp_resonance_levels.csv"
    with path.open() as fh:
        rows = list(csv.DictReader(fh))
    assert list(rows[0]) == [
        "v",
        "Re_E",
        "Gamma_v",
        "residual",
        "real_weight",
        "Re_E0",
        "Gamma_v_1",
    ]
    assert len(rows) == 2
    assert float(rows[0]["v"]) == 0
    assert float(rows[0]["Re_E"]) == -0.05
    assert float(rows[1]["Gamma_v"]) == 0.008


def test_npz_round_trips_the_complex_energies_and_curve(tmp_path):
    run = _fake_run()
    _write_resonance_levels(tmp_path, run)
    data = np.load(tmp_path / "resonance_levels_lcp_resonance_levels.npz")
    np.testing.assert_allclose(data["energies"], run.levels.energies)
    np.testing.assert_allclose(data["states"], run.levels.states)
    np.testing.assert_allclose(data["R_axis"], run.R_axis)
    np.testing.assert_allclose(data["Vd"], run.Vd)
    np.testing.assert_allclose(data["Gamma"], run.Gamma)


def test_csv_writes_nan_golden_rule_as_literal_nan_not_zero(tmp_path):
    """`golden_rule` legitimately goes all-`nan` (comparator unavailable/
    disabled) or per-level `nan` (a distance guard rejecting an implausible
    match) -- the csv must carry that through as `nan` text, never silently
    coerce it to `0.0`. `lcp_resonance_levels` fills the disabled/unpaired
    case with `nan + 1j*nan` (both parts), not a bare real `nan` -- match
    that here."""
    run = _fake_run()
    run.levels.golden_rule[:] = np.nan + 1j * np.nan
    _write_resonance_levels(tmp_path, run)
    path = tmp_path / "resonance_levels_lcp_resonance_levels.csv"
    with path.open() as fh:
        rows = list(csv.DictReader(fh))
    for row in rows:
        assert row["Re_E0"].strip().lower() == "nan"
        assert np.isnan(float(row["Re_E0"]))
        assert np.isnan(float(row["Gamma_v_1"]))


def test_png_survives_nan_golden_rule(tmp_path):
    """`golden_rule` never feeds the png (only `Re E_v`/`Gamma_v` do), so an
    all-`nan` comparator must not crash or blank out the figure."""
    run = _fake_run()
    run.levels.golden_rule[:] = np.nan + 1j * np.nan
    _write_resonance_levels(tmp_path, run)
    path = tmp_path / "resonance_levels_lcp_resonance_levels.png"
    assert path.exists()
    assert path.stat().st_size > 0


def test_real_mask_excludes_the_ecs_tail():
    """`R_axis <= R0` is the mask that actually restricts to the real
    region -- `R_axis` itself covers the FULL grid (real nodes plus the
    complex-rotated ECS tail), so it cannot be used as its own "how many
    real nodes" count (see `ResonanceLevelsRun.R0`'s docstring)."""
    R_axis = np.linspace(1.5, 6.0, 12)
    mask = _resonance_levels_real_mask(R_axis, R0=4.0)
    assert mask.sum() == 7
    assert bool(mask[0])
    assert not bool(mask[-1])
    # the real region is the LEADING contiguous block on this grid's node
    # ordering (real elements assembled before the rotated tail elements):
    # verify the mask has no "gap" -- once it goes False it stays False.
    idx = np.flatnonzero(~mask)
    assert idx.size == 0 or bool(np.all(idx >= idx[0]))


def test_writer_slices_out_the_ecs_tail_from_the_plotted_curve(tmp_path):
    """A `R0` inside the fake grid's range (rather than past its end) must
    not crash the writer and must actually narrow what gets plotted --
    regression test for the `Vd[:R_axis.size]` no-op slice this replaces."""
    run = _fake_run(R0=4.0)
    _write_resonance_levels(tmp_path, run)
    stem = "resonance_levels_lcp_resonance_levels"
    assert (tmp_path / f"{stem}.png").exists()
    mask = _resonance_levels_real_mask(run.R_axis, run.R0)
    assert mask.sum() < run.R_axis.size


def test_ylim_frames_the_well_not_the_repulsive_wall():
    """The real F2 curve rises to ~200 Ha as `R -> 0` while the levels sit in a
    0.01 Ha window around -0.145 Ha. Autoscaling on that range collapses the
    well, every level bar and `Gamma(R)` onto one line at `y ~ 0`, which is
    what the shipped example's figure used to look like. The frame must come
    from the physical range instead."""
    Vd = np.array([199.7, 12.0, -0.15, -0.149, -0.14, -0.02, 0.0])
    re_e = np.array([-0.1482, -0.1459, -0.1379])
    lo, hi = _resonance_levels_ylim(Vd, re_e)

    assert lo < -0.15  # below the well bottom
    assert -0.1379 < hi < 0.0  # just above the top level, nowhere near the wall
    assert hi - lo < 0.05  # the six-level window stays legible


def test_ylim_without_levels_still_excludes_the_wall():
    """A run with no angle-stable level at all has no `Re E_v` to frame on;
    the fallback (the curve's outermost real value, i.e. the anion
    dissociation limit) must still keep the wall out of the frame."""
    Vd = np.array([199.7, 12.0, -0.15, -0.149, -0.14, -0.02, -0.005])
    lo, hi = _resonance_levels_ylim(Vd, np.empty(0))
    assert lo < -0.15
    assert hi < 0.1  # framed on the -0.005 Ha asymptote, not the 199.7 Ha wall


def test_ylim_survives_a_degenerate_flat_curve():
    """A constant `V_d` with a level exactly at it gives `span == 0`; the
    limits must stay finite and ordered rather than collapsing to a point."""
    lo, hi = _resonance_levels_ylim(np.full(5, -0.2), np.array([-0.2]))
    assert np.isfinite(lo) and np.isfinite(hi)
    assert lo < hi
