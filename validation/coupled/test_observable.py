"""The gate: three criteria, one verdict, and a decision that is recorded
either way. A closed gate is a result, not a failure."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from validation.coupled.observable import E_SWEEP, gate_decision


def _report(dgamma: float, dsigma: float, n_poles_max: int) -> dict:
    return {
        "max_relative_gamma_shift": dgamma,
        "median_relative_sigma_shift": dsigma,
        "max_n_poles": n_poles_max,
    }


def test_a_second_genuine_pole_opens_the_gate_on_its_own() -> None:
    out = gate_decision(_report(0.0, 0.0, 2))
    assert out["open"] is True
    assert "second genuine pole" in out["reason"]


def test_the_ubiquitous_artefact_does_not_open_the_gate() -> None:
    """`n_stable` is 2 everywhere because a spurious near-threshold state is
    always present. The gate must read `n_poles`, so a campaign that found one
    genuine pole per R leaves criterion (a) shut however many stable states
    were counted."""
    out = gate_decision(_report(0.0, 0.0, 1))
    assert out["open"] is False


def test_a_large_width_shift_opens_the_gate() -> None:
    out = gate_decision(_report(0.08, 0.0, 1))
    assert out["open"] is True


def test_a_large_cross_section_shift_opens_the_gate() -> None:
    out = gate_decision(_report(0.0, 0.09, 1))
    assert out["open"] is True


def test_small_effects_leave_the_gate_shut() -> None:
    out = gate_decision(_report(0.02, 0.03, 1))
    assert out["open"] is False
    assert "not run" in out["reason"]


def test_a_zero_curve_difference_leaves_the_cross_section_unchanged() -> None:
    """The differential structure: if the coupled curve equals the fixed-l
    curve, the two cross sections must be identical, not merely close. Anything
    else means the two branches are not going through the same code."""
    from qscat.core.grids import electronic_grid

    from validation.coupled.observable import lcp_from_curve
    from validation.coupled.screen import CoupledCurve

    # `electronic_grid` stands in for a generic FemDvrEcsGrid here -- the test
    # only exercises the LCP resolvent's differential structure, not any real
    # electronic physics -- but its fixed inner segments require `r_max` to
    # exceed 10 bohr (`qscat.core.grids.electronic_grid`), so 6.0 (as given)
    # raises `GridError`; 12.0 is the smallest round value that clears it.
    grid = electronic_grid(r_max=12.0, angle_deg=30.0, order=6, n_complex=3)
    R = np.linspace(1.8, 4.0, 9)
    curve = CoupledCurve(
        R=R,
        E_res=np.full(R.size, 0.03 - 0.005j),
        residual=np.zeros(R.size),
        n_stable=np.ones(R.size, dtype=np.intp),
        n_poles=np.ones(R.size, dtype=np.intp),
    )
    vd = np.asarray(-0.05 + 0.0 * grid.real_points, dtype=np.complex128)
    gamma = np.full(grid.n, 0.01)
    out = lcp_from_curve(curve, curve, grid, vd, gamma)
    np.testing.assert_array_equal(out["full"], out["fixed"])


def test_the_energy_sweep_is_the_declared_one() -> None:
    assert E_SWEEP.size == 41
    assert E_SWEEP[0] == pytest.approx(0.020)
    assert E_SWEEP[-1] == pytest.approx(0.100)
    # NO's DA channel opens at +0.172 Ha, well above this sweep -- which is
    # why the gate observable is VE and never DA.
    assert float(np.max(E_SWEEP)) < 0.172


# `_summarize` and `main` produced every number in `gate.json` and had no
# test of their own -- the nested-depth walking (`s_curves` is {s: payload},
# `kappa_curves` is {kappa: {s: payload}}), the shared-`s` matching, and the
# perturbation descent are the most defect-prone code in this project. These
# tests build a small SYNTHETIC report with the same shape `screen.main()`
# produces -- they do not read the committed campaign file, so they exercise
# `_summarize`'s own logic rather than re-asserting one fixed dataset's
# numbers.
_R_SYNTH = [2.0, 2.5, 3.0, 3.5, 4.0]


def _curve_payload(v_d: float, gamma: float) -> dict[str, list[float]]:
    """A `CoupledCurve`-shaped payload, constant over `_R_SYNTH`."""
    n = len(_R_SYNTH)
    return {
        "v_d": [v_d] * n,
        "gamma": [gamma] * n,
        "residual": [1e-9] * n,
        "n_stable": [1] * n,
        "n_poles": [1] * n,
    }


def _synthetic_screen_report(walk: dict[str, tuple[float, float]]) -> dict:
    """A minimal `screen.main()`-shaped report for `kappa = 0.5`,
    `n_channels` 1 (fixed-l) vs. 2 (full) -- `_summarize` only ever reads
    `kappa_curves[str(max(n_channels))]` and `kappa_curves["1"]`, so a
    two-channel report exercises exactly the same code path as the real
    four-channel campaign.

    `walk` is `{s: (gamma_full, gamma_fixed)}`; `v_d` is held fixed at 0.03
    Ha (NO's resonance window) for both branches at every `s`, so the curve
    difference driving `_perturbation_fraction`/`lcp_from_curve` is in
    `Gamma` alone -- exactly what these tests need to control.
    """
    fixed_by_s = {s: _curve_payload(0.03, g_fixed) for s, (_, g_fixed) in walk.items()}
    full_by_s = {s: _curve_payload(0.03, g_full) for s, (g_full, _) in walk.items()}
    return {
        "n_channels": [1, 2],
        "R": _R_SYNTH,
        "s_curves": {"1": fixed_by_s, "2": full_by_s},
        "kappa_curves": {"1": {"0.5": fixed_by_s}, "2": {"0.5": full_by_s}},
        "n_channels_5_check": {},
    }


def test_summarize_picks_the_largest_shared_s_for_gamma_and_descends_for_sigma() -> None:
    """`gamma_s` (criterion b) needs no construction, so it must be the
    LARGEST shared `s` regardless of how big the difference is. `sigma_s`
    (criterion c) rides the curve-difference-on-shipped-LCP construction, so
    it must walk DOWN from the largest shared `s` to the largest one where
    that construction is still a perturbation (see the module docstring and
    `docs/physics/coupled-partial-waves.md`).

    Gamma differs by orders of magnitude (50 Ha vs. 0.01 Ha) at s = 0.3 and
    0.2 -- enormously past any real `Gamma_shipped(R)` scale, so the
    perturbation guard must reject both -- and is identical (no difference)
    at s = 0.1, which must therefore be where criterion (c) is evaluated.
    """
    from validation.coupled.observable import _summarize

    report = _synthetic_screen_report(
        {
            "0.3": (50.0, 0.01),
            "0.2": (50.0, 0.01),
            "0.1": (0.01, 0.01),
        }
    )
    summary = _summarize(report)
    assert summary["gamma_s"] == pytest.approx(0.3)
    assert summary["sigma_s"] == pytest.approx(0.1)
    # criterion (b) reads the raw payload difference directly: (50-0.01)/0.01.
    assert summary["max_relative_gamma_shift"] == pytest.approx((50.0 - 0.01) / 0.01)
    # criterion (c) is evaluated where the two curves are IDENTICAL, so the
    # cross section difference must be exactly zero.
    assert summary["median_relative_sigma_shift"] == pytest.approx(0.0, abs=1e-12)


def test_summarize_raises_when_no_shared_s_passes_the_perturbation_limit() -> None:
    """If the curve difference exceeds the perturbation limit at EVERY
    shared `s`, criterion (c) has no `s` at which it can be evaluated, and
    `_summarize` must say so rather than silently reporting a number built
    from a construction that has already collapsed."""
    from validation.coupled.observable import _summarize

    report = _synthetic_screen_report({"0.2": (50.0, 0.01)})
    with pytest.raises(ValueError, match="every shared s"):
        _summarize(report)


def test_main_reads_the_screen_report_and_writes_the_gate_from_results(tmp_path: Path) -> None:
    """`main` end to end, entirely inside `tmp_path`: it must read
    `screen.json` from and write `gate.json` to the `results` directory it
    is given, never the committed `validation/coupled/results/` -- the fix
    for the same tracked-file-overwrite problem `screen.main` has."""
    from validation.coupled.observable import main

    report = _synthetic_screen_report(
        {
            "0.3": (50.0, 0.01),
            "0.2": (50.0, 0.01),
            "0.1": (0.01, 0.01),
        }
    )
    (tmp_path / "screen.json").write_text(json.dumps(report))

    out = main(results=tmp_path)

    assert (tmp_path / "gate.json").exists()
    written = json.loads((tmp_path / "gate.json").read_text())
    assert written == out
    assert out["summary"]["gamma_s"] == pytest.approx(0.3)
    assert out["summary"]["sigma_s"] == pytest.approx(0.1)
    # Gamma moved by orders of magnitude at the criterion-(b) point, so the
    # gate must be open.
    assert out["verdict"]["open"] is True
