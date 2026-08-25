"""Sub-project #C4, Task 4: the honest three-way N2 comparison of the
`TannorWeeks`/`Dirac`/`Flux` TD energy extractors, all driven by the SAME
propagation via `qscat.core.time_dependent.td_ve_cross_sections_all`.

`main()` (the `@slow` path, run via `python -m validation.n2.td_extractors`)
computes sigma_TW/sigma_delta/sigma_flow at the converged N2 working grid
`libs/qscat/tests/test_td_extractors.py`'s `@slow` anchor tests already
validate individually (`Dirac`/`Flux` vs the TI oracle, `sigma_delta/sigma_ti
= 0.971`, `sigma_flow/sigma_ti = 0.970` at E=0.10) -- this module reruns that
SAME grid/wavepacket/energy but with all three extractors sharing ONE
propagation (`td_ve_cross_sections_all`), so the reported three-way spread is
the honest, identical-dynamics comparison the brief asks for, and adds the
Houfek cross-check (`validation.n2.loader`) alongside the TI oracle.

Budget: one propagation at this grid (T=1000 a.u., order-3 Pade, dt=1.0)
costs ~240-300s wall (measured in the individual `Dirac`/`Flux` `@slow`
tests); running THREE extractors from the same trajectory adds negligible
overhead (the sparse LU back-substitutions dominate, not the O(1)-per-step
c-product bookkeeping each extractor performs) -- so the three-way run costs
about the SAME as one of the individual single-extractor runs, not three
times as much. `vprimes=[1]` (inelastic only, matching the individual
anchor tests) means `v_init=0` is never requested, so no elastic
free-reference propagation is needed either.

Only the E=0.10 anchor (the GATED C5/D1 anchor `(0.1, 1)` in
`validation.n2.reference.ANCHOR_COORDS`) is run live here, within the
harness's documented patience budget; E=0.15 (used elsewhere as a validated,
non-C5-gated anchor -- see `docs/physics/n2-2d-td-cross-section.md`) is
RECORDED as a note from the same already-published measurements
(`libs/qscat/tests/test_td_extractors.py`'s docstrings) rather than re-run,
following `validation/n2` Group F's own precedent
(`validation/n2/td_exact2d.py`) for citing an already-validated, expensive
propagation instead of re-running it inside a budget that cannot afford a
second ~250s propagation. This is a recorded, cited NOTE, not a silent
fabrication -- see the module docstring of `td_exact2d.py` for the same
pattern.
"""

from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np
from qscat.core.driven import ve_cross_section
from qscat.core.grids import electronic_grid, nuclear_grid
from qscat.core.time_dependent import td_ve_cross_sections_all
from qscat.core.vibrational import vibrational_states
from qscat.dvr import TensorGrid
from qscat.model import N2

from validation.figures import FIGURE_DIR
from validation.n2 import loader

__all__ = [
    "DT",
    "E_LIVE",
    "FIGURE_PATH",
    "N_STEPS",
    "POSITION",
    "SURFACE",
    "VPRIMES",
    "V_INIT",
    "WP_IN",
    "WP_OUT",
    "ThreeWayResult",
    "compute_live_result",
    "converged_tgrid",
    "main",
    "recorded_note_result",
]

FIGURE_PATH = FIGURE_DIR / "n2-td-extractors-comparison.png"

# Converged N2 working grid -- SAME as `libs/qscat/tests/test_td_extractors.py`'s
# `test_delta_agrees_with_ti_oracle_one_anchor`/`test_flux_agrees_with_ti_oracle_
# one_anchor` (@slow), so this module's live E=0.10 measurement is directly
# comparable to (and reproduces) those already-recorded numbers.
V_INIT = 0
VPRIMES = [1]  # inelastic only -- no elastic free-reference propagation needed
WP_IN = {"r0": 25.0, "p0": -0.5, "sigma": 5.0}
WP_OUT = {"r0_out": 35.0, "p0_out": 0.5, "sigma_out": 4.0}
DT = 1.0
N_STEPS = 1000  # T = 1000 a.u.
POSITION = 128  # r = 39.58 bohr, real region (R0=50), past the interaction
SURFACE = 128

E_LIVE = 0.10  # the GATED C5/D1 anchor `(0.1, 1)` -- run live, within budget

# Recorded E=0.15 ratios from `libs/qscat/tests/test_td_extractors.py`'s @slow
# docstrings (same grid/wavepacket, T=1000 propagation): sigma_delta/sigma_ti =
# 1.009, sigma_flow/sigma_ti = 1.007. Those tests drove ONE extractor at a time,
# not the combined three-way path -- this module's own `compute_live_result` is
# what exercises the ONE-propagation three-way path (run live at E=0.10);
# E=0.15 is cited here from those single-extractor runs, not re-run.
_RECORDED_E15_DELTA_RATIO = 1.009
_RECORDED_E15_FLOW_RATIO = 1.007
_RECORDED_SOURCE = (
    "libs/qscat/tests/test_td_extractors.py::test_delta_agrees_with_ti_oracle_one_anchor "
    "/ test_flux_agrees_with_ti_oracle_one_anchor (@slow docstrings, single-extractor "
    "propagations at the same grid/wavepacket)"
)


def converged_tgrid() -> TensorGrid:
    return TensorGrid(
        [
            electronic_grid(r_max=50.0, order=8, n_complex=6),
            nuclear_grid(quadrature=10, r_max=22.0, n_complex=5),
        ]
    )


@dataclass(frozen=True)
class ThreeWayResult:
    """One energy's sigma_TW/sigma_delta/sigma_flow (bohr^2) vs the TI
    oracle and Houfek's golden data, from ONE shared propagation (`live`)
    or cited from a previously-published, already-validated measurement
    (`live=False`)."""

    energy_ha: float
    channel: int
    sigma_tw: float
    sigma_delta: float
    sigma_flow: float
    sigma_ti: float
    sigma_houfek: float
    ratio_tw_ti: float
    ratio_delta_ti: float
    ratio_flow_ti: float
    ratio_delta_tw: float
    ratio_flow_tw: float
    live: bool
    wall_seconds: float | None  # None for a recorded (not-run-here) row
    source: str


def _houfek_sigma(energy_ha: float, channel: int) -> float:
    d = loader.load()
    i = int(np.argmin(np.abs(d.energy - energy_ha)))
    return float(d.sigma[i, channel])


def compute_live_result(energy_ha: float = E_LIVE, channel: int = VPRIMES[0]) -> ThreeWayResult:
    """The ONE live three-way measurement this module runs: a single
    `td_ve_cross_sections_all` propagation (T=1000, order-3 Pade) at the
    converged grid, transformed at `energy_ha`, compared against the TI
    oracle (`qscat.core.driven.ve_cross_section`, cheap -- one sparse solve)
    and Houfek's golden data.
    """
    tgrid = converged_tgrid()
    eps, chi = vibrational_states(tgrid.grids[1], N2.mu, 4, N2.v0)

    t0 = time.time()
    sigma_all = td_ve_cross_sections_all(
        tgrid,
        N2,
        eps,
        chi,
        V_INIT,
        [channel],
        energy_ha,
        dt=DT,
        n_steps=N_STEPS,
        wp_in=WP_IN,
        wp_out=WP_OUT,
        position=POSITION,
        surface=SURFACE,
    )
    wall = time.time() - t0

    sigma_ti = float(ve_cross_section(tgrid, N2, eps, chi, V_INIT, [channel], energy_ha)[0])
    sigma_houfek = _houfek_sigma(energy_ha, channel)

    sigma_tw = float(sigma_all["tw"][0])
    sigma_delta = float(sigma_all["delta"][0])
    sigma_flow = float(sigma_all["flow"][0])

    return ThreeWayResult(
        energy_ha=energy_ha,
        channel=channel,
        sigma_tw=sigma_tw,
        sigma_delta=sigma_delta,
        sigma_flow=sigma_flow,
        sigma_ti=sigma_ti,
        sigma_houfek=sigma_houfek,
        ratio_tw_ti=sigma_tw / sigma_ti if sigma_ti else float("inf"),
        ratio_delta_ti=sigma_delta / sigma_ti if sigma_ti else float("inf"),
        ratio_flow_ti=sigma_flow / sigma_ti if sigma_ti else float("inf"),
        ratio_delta_tw=sigma_delta / sigma_tw if sigma_tw else float("inf"),
        ratio_flow_tw=sigma_flow / sigma_tw if sigma_tw else float("inf"),
        live=True,
        wall_seconds=wall,
        source="this module, live run (td_ve_cross_sections_all, one propagation)",
    )


def recorded_note_result(energy_ha: float = 0.15, channel: int = VPRIMES[0]) -> ThreeWayResult:
    """E=0.15's row, cited from already-published single-extractor
    measurements (module docstring) rather than re-run -- a second full
    propagation does not fit this module's patience budget alongside the
    live E=0.10 run. `sigma_tw` is left as `nan` (TW was not part of those
    individually-run @slow tests); only the delta/flow ratios are recorded.
    """
    sigma_houfek = _houfek_sigma(energy_ha, channel)
    return ThreeWayResult(
        energy_ha=energy_ha,
        channel=channel,
        sigma_tw=float("nan"),
        sigma_delta=float("nan"),
        sigma_flow=float("nan"),
        sigma_ti=float("nan"),
        sigma_houfek=sigma_houfek,
        ratio_tw_ti=float("nan"),
        ratio_delta_ti=_RECORDED_E15_DELTA_RATIO,
        ratio_flow_ti=_RECORDED_E15_FLOW_RATIO,
        ratio_delta_tw=float("nan"),
        ratio_flow_tw=float("nan"),
        live=False,
        wall_seconds=None,
        source=_RECORDED_SOURCE,
    )


def _print_report(results: list[ThreeWayResult]) -> None:
    print("N2 TD extractor three-way comparison (TannorWeeks / Dirac(delta) / Flux(flow))")
    print("=" * 78)
    for r in results:
        tag = "LIVE" if r.live else "RECORDED"
        print(f"\nE = {r.energy_ha} Ha, v'={r.channel}  [{tag}]")
        if r.live:
            print(f"  wall time        : {r.wall_seconds:.1f}s")
            print(f"  sigma_TW         : {r.sigma_tw:.4f} bohr^2")
            print(f"  sigma_delta      : {r.sigma_delta:.4f} bohr^2")
            print(f"  sigma_flow       : {r.sigma_flow:.4f} bohr^2")
            print(f"  sigma_TI (oracle): {r.sigma_ti:.4f} bohr^2")
            print(f"  sigma_Houfek     : {r.sigma_houfek:.4f} bohr^2")
            print(f"  ratio TW/TI      : {r.ratio_tw_ti:.4f}")
            print(f"  ratio delta/TI   : {r.ratio_delta_ti:.4f}")
            print(f"  ratio flow/TI    : {r.ratio_flow_ti:.4f}")
            print(f"  ratio delta/TW   : {r.ratio_delta_tw:.4f}")
            print(f"  ratio flow/TW    : {r.ratio_flow_tw:.4f}")
        else:
            print(f"  sigma_Houfek     : {r.sigma_houfek:.4f} bohr^2")
            print(f"  ratio delta/TI   : {r.ratio_delta_ti:.4f}  (recorded, not re-run)")
            print(f"  ratio flow/TI    : {r.ratio_flow_ti:.4f}  (recorded, not re-run)")
        print(f"  source: {r.source}")


def _write_figure(results: list[ThreeWayResult]) -> None:
    """Accuracy/cost figure: per-method sigma-vs-TI ratio at the live
    anchor, plus the recorded E=0.15 ratios for delta/flow, and a
    qualitative cost bar (TW/delta/flow relative per-step cost)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    live = next(r for r in results if r.live)
    recorded = next((r for r in results if not r.live), None)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5))

    methods = ["TW", "delta", "flow"]
    ratios_e010 = [live.ratio_tw_ti, live.ratio_delta_ti, live.ratio_flow_ti]
    ax1.bar(methods, ratios_e010, color=["#4C72B0", "#DD8452", "#55A868"])
    ax1.axhline(1.0, color="k", linestyle="--", linewidth=1.0, label="TI oracle")
    if recorded is not None:
        ax1.scatter(
            [1, 2],
            [recorded.ratio_delta_ti, recorded.ratio_flow_ti],
            marker="x",
            color="k",
            zorder=5,
            label=f"E={recorded.energy_ha} (recorded)",
        )
    ax1.set_ylabel(r"$\sigma_{\rm method} / \sigma_{\rm TI}$")
    ax1.set_title(f"sigma-vs-TI accuracy at E={live.energy_ha} Ha, v'={live.channel}")
    ax1.legend(fontsize="small")
    ax1.grid(True, alpha=0.3)

    # Qualitative per-step cost: TW does a full outgoing-channel c-product
    # (volume overlap) each step; delta a single point projection; flow a
    # point projection PLUS its DVR first-derivative row -- roughly double
    # delta's per-step bookkeeping (both dwarfed by the shared O(N) sparse
    # LU back-substitution, so this bar is qualitative, not a measured
    # profile).
    cost_labels = ["delta\n(point)", "flow\n(point+deriv)", "TW\n(full overlap)"]
    cost_relative = [1.0, 1.5, 3.0]
    ax2.barh(cost_labels, cost_relative, color=["#DD8452", "#55A868", "#4C72B0"])
    ax2.set_xlabel("qualitative per-step extraction cost (relative)")
    ax2.set_title("Cost ranking (shared O(N) LU solve dominates all three)")

    fig.tight_layout()
    FIGURE_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE_PATH)
    plt.close(fig)
    print(f"\nWrote {FIGURE_PATH}")


def main() -> list[ThreeWayResult]:
    live = compute_live_result()
    recorded = recorded_note_result()
    results = [live, recorded]
    _print_report(results)
    _write_figure(results)
    return results


if __name__ == "__main__":
    main()
