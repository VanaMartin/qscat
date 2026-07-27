"""Task 7: Group F -- the time-dependent (TD) 2-D cross section
(`projects/n2_2d_td_cross_section`, sub-project #7) reported into the
harness WITHOUT a live TD propagation, but WITH a live TI cross-check.

Why no live TD check (measured, not assumed): a full propagation at
`TD_WORKING_GRID` (N=47188, T=1500 a.u., order-3 Pade / dt=1.0) takes
several minutes wall -- and BOTH a full run and a V_int=0 free-reference run
are needed for the elastic channel, doubling that. The order-3 Pade operator
applies 3 LU back-substitutions per step (vs Crank-Nicolson's 1), so even at
dt=1.0/n_steps=1500 the propagation is well over the harness's ~60s per-group
budget. (Order-1 CN was cheaper per step but under-converged, capping
sigma_TD/sigma_TI at ~0.93/1.10; the order-3 operator brings both anchors to
~0.97/0.99 -- see docs/physics/n2-2d-td-cross-section.md -- at the cost of the
extra factors, reinforcing the don't-run-live decision.) Group F therefore
reports the ALREADY-VALIDATED sigma_TD as a literal,
cited constant: no TD propagation runs when the harness executes this
module. See `docs/physics/n2-2d-td-cross-section.md` for the full method
and `.superpowers/sdd/task-4-report.md` / `task-6-report.md` for the
underlying runs this constant is taken from verbatim.

sigma_TI, by contrast, IS obtained live: it is cheap (sparse solves, ~seconds,
not a ~250s propagation), so hardcoding it would be needless. It comes from
`validation.n2.exact2d.compute_exact2d_results()` -- Group E's exact 2-D TI
solver, already run (and `lru_cache`d) earlier in the same harness pass --
reusing that cache is free. E=0.10/v'=1 is itself one of Group E's anchors
(`reference.ANCHOR_COORDS`), so its sigma_TI is read straight out of the
cached result list. E=0.15/v'=1 is not one of Group E's anchors, so for that
one row this module makes one extra, equally cheap live call to
`qscat.core.driven.ve_cross_section` on the SAME cached working grid/vibrational
states `exact2d.build_system()` already built for Group E.

These rows are NOTE, never PASS/FAIL: they report a validated fact, not a
live gate. The genuine PASS/FAIL gate on this same comparison lives in
`projects/n2_2d_td_cross_section/test_td_cross_section.py`'s `@pytest.mark.slow`
tests (run outside the default harness; run explicitly via
`uv run pytest projects/n2_2d_td_cross_section -m slow`).
"""

from __future__ import annotations

from dataclasses import dataclass

from qscat.core.driven import ve_cross_section
from qscat.model import N2

from validation.n2 import exact2d

__all__ = ["TdExact2dResult", "compute_td_exact2d_results"]

# Recorded sigma_TD anchors (TD_WORKING_GRID: order-3 Pade, dt=1.0, T=1500).
# There is no cheap live path for these -- a live check costs ~200s, over the
# harness's ~60s/group budget (see module docstring) -- so sigma_TD stays a
# cited, literal constant. sigma_TI is NOT included here: it is obtained live
# in `compute_td_exact2d_results()` below, from Group E's already-cached
# exact 2-D TI solver. (Order-1 Crank-Nicolson previously gave 5.6973/0.6904
# here, ratios 0.93/1.10; the order-3 Pade operator brings both to ~0.97/0.99.)
_RECORDED_TD_ANCHORS = (
    # (energy_ha, channel, sigma_td, rtol, source)
    (
        0.10,
        1,
        5.9595,
        0.06,
        "test_td_cross_section.py::test_v2a_td_matches_ti_at_e010 (@slow) / "
        "order-3 Pade TD_WORKING_GRID run",
    ),
    (
        0.15,
        1,
        0.6185,
        0.06,
        "test_td_cross_section.py::test_v2a_td_matches_ti_at_e015 (@slow) / "
        "order-3 Pade TD_WORKING_GRID run",
    ),
)

_ENERGY_MATCH_TOL_HA = 1e-6


@dataclass(frozen=True)
class TdExact2dResult:
    """One TD-vs-TI anchor: sigma_TD is a RECORDED constant from a
    previously-measured, validated propagation at `TD_WORKING_GRID` (Task 4's
    T-scan and Task 6's independent reproduction of it) -- NOT recomputed
    live by the harness (see module docstring for why). sigma_TI, by
    contrast, IS computed live, every harness run, via
    `validation.n2.exact2d` -- see `compute_td_exact2d_results()`.

    `rtol` is the tolerance the corresponding `@slow` test in
    `projects/n2_2d_td_cross_section/test_td_cross_section.py` gates on
    (imported nowhere here to avoid a hard dependency on that test module;
    restated as a literal, matching it by construction -- see that module's
    own docstring/tests for the authoritative value).
    """

    energy_ha: float
    channel: int  # v' (v_init = 0 in every case here)
    sigma_td: float  # bohr^2, recorded (Task 4/6 TD_WORKING_GRID run)
    sigma_ti: float  # bohr^2, LIVE (this harness run's exact 2-D TI solve)
    ratio_td_ti: float  # sigma_td / sigma_ti (recomputed from the live sigma_ti)
    rtol: float  # the @slow test's gate on |ratio - 1|
    source: str  # where sigma_td comes from


def _sigma_ti_live(energy_ha: float, channel: int) -> float:
    """sigma_TI at (energy_ha, channel), computed live.

    Reuses Group E's `exact2d.compute_exact2d_results()` (already `lru_cache`d
    from Group E's own run earlier in this same harness pass, so this lookup
    is free) whenever the coordinate is one of its anchors. Otherwise falls
    back to one live `ve_cross_section` call (`qscat.core.driven`, against
    `qscat.model.N2` directly) on the SAME cached working grid/vibrational
    states `exact2d.build_system()` built for Group E -- still cheap (sparse
    solves, ~seconds), just not pre-tabulated.
    """
    for r in exact2d.compute_exact2d_results():
        if r.channel == channel and abs(r.energy_ha - energy_ha) < _ENERGY_MATCH_TOL_HA:
            return r.sigma_exact
    tgrid, eps, chi = exact2d.build_system()
    sigma = ve_cross_section(tgrid, N2, eps, chi, 0, [channel], energy_ha)
    return float(sigma[0])


def compute_td_exact2d_results() -> list[TdExact2dResult]:
    """The 2 TD-vs-TI anchors at `TD_WORKING_GRID`'s converged T=1500
    configuration: sigma_TD recorded (see module docstring for why a live TD
    propagation is not run in-harness), sigma_TI computed LIVE every run via
    `validation.n2.exact2d` (cheap -- see `_sigma_ti_live`).

    sigma_TD is reproduced independently across Task 4 (the T-scan that
    established T=1500 as converged), Task 5/6 (the full sigma(E) curve and
    the committed `.npz`/figures), and `test_td_cross_section.py`'s `@slow`
    V2a test (the actual gating pytest check) -- all four report the SAME
    number to the digits shown, so restating it here is not a fresh
    measurement subject to run-to-run drift, it is citing an already
    multiply-reproduced result.
    """
    results = []
    for energy_ha, channel, sigma_td, rtol, source in _RECORDED_TD_ANCHORS:
        sigma_ti = _sigma_ti_live(energy_ha, channel)
        results.append(
            TdExact2dResult(
                energy_ha=energy_ha,
                channel=channel,
                sigma_td=sigma_td,
                sigma_ti=sigma_ti,
                ratio_td_ti=sigma_td / sigma_ti,
                rtol=rtol,
                source=source,
            )
        )
    return results
