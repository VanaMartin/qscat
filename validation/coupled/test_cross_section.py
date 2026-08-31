"""Campaign structure, on a grid small enough to run in the fast tier.

The `test_archived_*` tests below lock the numbers of the BARE `s = 0.3`
campaign, whose interpretations are withdrawn -- that model's anion is unbound
at every R, so its cross section is not a property of partial-wave coupling.
They are kept deliberately: the archived record in the physics note quotes
those numbers, and a withdrawn result still has to be reproducible from the
committed data or the archive cannot be audited. They guard the RECORD, not a
claim.
"""

from __future__ import annotations

import json

import numpy as np
import pytest
from qscat.core.grids import segmented_grid
from qscat.dvr import TensorGrid
from scipy.signal import find_peaks

from validation.coupled.cross_section import (
    KAPPA_RUN,
    N_CHANNEL_VALUES,
    RESULTS,
    S_RUN,
    VPRIMES,
    main,
)

# The note's two headline tables (docs/physics/coupled-partial-waves.md, "The
# cross section"), transcribed here so a silent drift between the committed
# data and the prose it is quoted from fails a test rather than waiting for
# the next review. Loads the already-committed `results/cross_section.json`
# -- no solve, so this belongs in the fast tier.
_PEAK_SIGMA = {
    "1": [67.032, 5.970, 0.754, 0.165, 0.040],
    "4": [70.918, 11.833, 3.223, 1.202, 0.474],
}
_PEAK_E = {
    "1": [0.07425, 0.06925, 0.07375, 0.07850, 0.08375],
    "4": [0.05800, 0.05050, 0.05400, 0.05850, 0.06400],
}
_TV_RATIO = [1.062, 1.020, 1.018, 1.009, 1.007]
# median |sigma4-sigma1|/sigma4 over E > 0.10 Ha (v'=0..4) -- what the note
# and CLAUDE.md both call the "tail does NOT converge" claim, corrected once
# already in the note (round 2) and once more in CLAUDE.md (round 4) after
# surviving a whole review cycle there. Guarded here so a reintroduced "the
# tail is where the models converge" can no longer pass silently.
_TAIL_MEDIAN = [0.251, 0.240, 0.157, 0.137, 0.144]


def _committed_report() -> dict:
    return json.loads((RESULTS / "cross_section.json").read_text())


def test_archived_bare_campaign_has_one_peak_per_curve_per_channel() -> None:
    """One local maximum per curve per channel, at both prominence floors
    the note cites (1% and 0.1% of the channel's own peak height) -- the
    structural null the ARCHIVED bare-campaign result rested on."""
    d = _committed_report()
    for n_ch in ("1", "4"):
        sigma = np.asarray(d["sigma"][n_ch]["total"], dtype=float)
        for vp in range(5):
            s = sigma[:, vp]
            for frac in (0.01, 0.001):
                peaks, _props = find_peaks(s, prominence=frac * s.max())
                assert len(peaks) == 1, (n_ch, vp, frac, len(peaks))


def test_archived_bare_campaign_peak_table_matches_the_note() -> None:
    """The peak heights, positions, ratio and shift the note's ARCHIVE quotes --
    all reproduced from the committed JSON alone, so the withdrawn result
    stays auditable."""
    d = _committed_report()
    E = np.asarray(d["sigma"]["1"]["E"], dtype=float)
    for n_ch in ("1", "4"):
        sigma = np.asarray(d["sigma"][n_ch]["total"], dtype=float)
        for vp in range(5):
            i = sigma[:, vp].argmax()
            assert abs(sigma[i, vp] - _PEAK_SIGMA[n_ch][vp]) < 5e-3
            assert abs(E[i] - _PEAK_E[n_ch][vp]) < 1e-6

    sigma1 = np.asarray(d["sigma"]["1"]["total"], dtype=float)
    sigma4 = np.asarray(d["sigma"]["4"]["total"], dtype=float)
    ratio_2sf = [round(sigma4[:, vp].max() / sigma1[:, vp].max(), 1) for vp in range(5)]
    assert ratio_2sf == [1.1, 2.0, 4.3, 7.3, 11.8]


def test_archived_bare_campaign_tv_ratios_and_channel_truncation() -> None:
    """Total variation (normalised by peak-to-peak range) and the N_l=3-vs-4
    truncation check the median/max table reports."""
    d = _committed_report()
    for vp in range(5):
        s1 = np.asarray(d["sigma"]["1"]["total"], dtype=float)[:, vp]
        s3 = np.asarray(d["sigma"]["3"]["total"], dtype=float)[:, vp]
        s4 = np.asarray(d["sigma"]["4"]["total"], dtype=float)[:, vp]
        tv1 = np.sum(np.abs(np.diff(s1))) / (s1.max() - s1.min())
        tv4 = np.sum(np.abs(np.diff(s4))) / (s4.max() - s4.min())
        assert abs(tv4 / tv1 - _TV_RATIO[vp]) < 2e-3

        open_mask = s4 > 0
        med = np.median(np.abs(s4[open_mask] - s3[open_mask]) / s4[open_mask])
        assert med < 1e-3


def test_archived_bare_campaign_tail_does_not_converge() -> None:
    """The tail (E > 0.10 Ha) settles at a 14-25 % median relative
    difference -- it does NOT shrink toward agreement. A few percent
    relative tolerance: this reads committed JSON (argmax/diff/divide, no
    BLAS path), so there is no cross-architecture concern, but a bound
    tight enough to be brittle buys nothing either."""
    d = _committed_report()
    E = np.asarray(d["sigma"]["1"]["E"], dtype=float)
    s1 = np.asarray(d["sigma"]["1"]["total"], dtype=float)
    s4 = np.asarray(d["sigma"]["4"]["total"], dtype=float)
    tail = E > 0.10
    for vp in range(5):
        rel = np.abs(s4[tail, vp] - s1[tail, vp]) / s4[tail, vp]
        med = np.median(rel)
        assert med == pytest.approx(_TAIL_MEDIAN[vp], rel=0.05), (vp, med)


def _tiny_grid() -> TensorGrid:
    el = segmented_grid(((4, 8.0),), ((2, 12.0),), angle_deg=35.0, quadrature=6)
    nu = segmented_grid(((3, 4.0),), ((2, 6.0),), angle_deg=30.0, quadrature=6, x_min=1.0)
    return TensorGrid([el, nu])


def test_the_declared_parameters_are_the_ones_measured() -> None:
    """S_RUN is the anisotropy the ARCHIVED campaign ran at. It was chosen as
    the point where all 41 R are resonant -- which is the same statement as
    "the anion is nowhere bound", the reason that campaign is superseded.
    N_l = 2 is deliberately absent: it was measured 30 % from converged
    against N_l = 4."""
    assert (S_RUN, KAPPA_RUN) == (0.3, 0.5)
    assert 2 not in N_CHANNEL_VALUES
    assert set(N_CHANNEL_VALUES) == {1, 3, 4}
    assert VPRIMES == [0, 1, 2, 3, 4]


def test_main_writes_every_model_on_one_mesh(tmp_path) -> None:
    """Two energies and a tiny grid: this checks the report's shape and that
    every model saw the SAME mesh, not the physics."""
    report = main(results=tmp_path, energies=np.array([0.02, 0.05]), grid=_tiny_grid())
    assert set(report["n_channels"]) == set(N_CHANNEL_VALUES)
    meshes = {tuple(report["sigma"][str(n)]["E"]) for n in N_CHANNEL_VALUES}
    assert len(meshes) == 1, "the branches must share one energy mesh"
    for n in N_CHANNEL_VALUES:
        s = report["sigma"][str(n)]
        assert np.shape(s["total"]) == (2, len(VPRIMES))
        assert np.shape(s["restricted"]) == (2, len(VPRIMES))
    assert (tmp_path / "cross_section.json").exists()


@pytest.mark.slow
def test_probe_reports_the_real_cost() -> None:
    """The campaign's cost estimate is an extrapolation; this is the measurement."""
    from validation.coupled.cross_section import probe_one_energy

    out = probe_one_energy(1)
    assert out["unknowns"] == 78804
    assert out["factor_s"] > 0.0
