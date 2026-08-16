"""Convergence study for the exact 2-D N2 solver, and the working grid it picks.

eMoScat asserts a 35-degree ECS angle and a 98-bohr electronic box without
documenting the study that justified them. This module redoes that study and
records the numbers, so the grid used for the benchmark is EARNED rather than
inherited.

Run: `uv run python -m projects.n2_2d_cross_section.convergence`
"""

from __future__ import annotations

import time
from typing import Any

from qscat.dvr import TensorGrid

from projects.n2_2d_cross_section.cross_section_2d import ve_cross_section_2d
from projects.n2_2d_cross_section.electronic_grid import n2_electronic_grid
from projects.n2_2d_cross_section.hamiltonian2d import MU
from projects.n2_ti_cross_section.nuclear_grid import n2_nuclear_grid
from projects.n2_ti_cross_section.vibrational import vibrational_states

__all__ = ["BASELINE", "WORKING_GRID", "build_tgrid", "working_tgrid", "convergence_table"]

# The anchor the study is run at: the resonance region, a well-behaved VE channel.
STUDY_E, STUDY_VP, N_VIB = 0.2, 1, 4

BASELINE: dict[str, Any] = {
    "r_max": 30.0,
    "angle_deg": 35.0,
    "order": 8,
    "n_complex": 8,
    "nuc_r_max": 40.0,
    "nuc_quadrature": 14,
    "nuc_n_complex": 10,
}

# Chosen from the measured table (the measured convergence table):
# every axis is converged to <=3.4e-7 relative (individually, about BASELINE)
# at N=71476 -- four to six orders of magnitude inside the ~1% criterion --
# so the cheap end of each axis was combined and VERIFIED directly (sweeps
# vary one axis at a time, so additivity is not assumed): at this combined
# grid, N=26857 (2.7x smaller than BASELINE), sigma = 1.256450927036e-01,
# a relative deviation of 2.368e-06 from BASELINE's 1.256447951966e-01 --
# still four orders of magnitude inside 1%, and a single sparse LU solve
# here costs ~3.1s wall (vs BASELINE's 38.0s).
WORKING_GRID: dict[str, Any] = {
    "r_max": 16.0,  # r_max 16->22->30->45: rel change 3.3e-10/1.1e-9/1.2e-10;
    # flat well below r_max=16 already -- the interaction lives at r<~3 bohr
    "angle_deg": 35.0,  # angle 25->30->35->40 @ BASELINE: rel change
    # 7.8e-10/6.1e-10/3.5e-10 (all tiny) -- BUT a SEPARATE measurement on
    # this same small grid (r_max=16, order=7, n_complex=5) found 25 deg
    # gives 6.75e-5 relative deviation from the converged value, ~30x worse
    # than the other three angles there: a shallow ECS contour combined
    # with few complex tail elements (n_complex=5, chosen below) under-
    # resolves the rotated continuum. 30-40 deg are all safely converged at
    # n_complex=5; 35 deg is chosen (matches eMoScat, sits mid-range).
    "order": 7,  # order 7->8: rel change 3.4e-7 (the LARGEST single-axis
    # change measured anywhere in the sweep); 8->9: 3.8e-9. Even the
    # largest measured change is ~4 orders of magnitude under the 1%
    # criterion, so order=7 (N=61204 vs 71476 @ BASELINE) is affordable.
    "n_complex": 5,  # n_complex 5->8: rel change 4.1e-8; 8->11: 2.6e-10.
    # Cheapest tested (N=62488 vs 71476 @ BASELINE), safely converged --
    # PROVIDED angle_deg stays >=30 deg (see angle_deg's comment above).
    "nuc_r_max": 20.0,  # nuc_r_max 20->30->40: rel change 6.2e-13/1.5e-13
    # (fully converged) -- and does not change nuclear grid point count at
    # all (only the ECS tail element length), so this choice costs nothing
    # either way; the smallest tested value is used.
    "nuc_quadrature": 10,  # nuc_quadrature 10->12->14: rel change
    # 3.5e-8/1.3e-10. Cheapest tested (nuclear n=296 vs 428 @ BASELINE,
    # N=49432 vs 71476), safely converged.
    "nuc_n_complex": 5,  # nuc_n_complex 5->10->13: rel change
    # 3.4e-13/3.6e-14 (fully converged, and -- like nuc_r_max -- barely
    # moves sigma at all). Cheapest tested (N=60621 vs 71476 @ BASELINE).
}

SWEEPS: dict[str, list[Any]] = {
    "r_max": [16.0, 22.0, 30.0, 45.0],
    "angle_deg": [25.0, 30.0, 35.0, 40.0],  # the sharpest check
    "order": [7, 8, 9],
    "n_complex": [5, 8, 11],
    "nuc_r_max": [20.0, 30.0, 40.0],
    "nuc_quadrature": [10, 12, 14],
    "nuc_n_complex": [5, 10, 13],
}


def build_tgrid(params: dict[str, Any]) -> TensorGrid:
    """Both grids share `angle_deg` -- one ECS contour angle for the problem."""
    return TensorGrid(
        [
            n2_electronic_grid(
                r_max=params["r_max"],
                angle_deg=params["angle_deg"],
                order=params["order"],
                n_complex=params["n_complex"],
            ),
            n2_nuclear_grid(
                r_max=params["nuc_r_max"],
                angle_deg=params["angle_deg"],
                quadrature=params["nuc_quadrature"],
                n_complex=params["nuc_n_complex"],
            ),
        ]
    )


def working_tgrid() -> TensorGrid:
    return build_tgrid(WORKING_GRID)


def _one(params: dict[str, Any], ordering: str = "COLAMD") -> dict[str, Any]:
    t0 = time.perf_counter()
    tg = build_tgrid(params)
    eps, chi = vibrational_states(tg.grids[1], MU, N_VIB)
    sigma = float(ve_cross_section_2d(tg, eps, chi, 0, [STUDY_VP], STUDY_E, ordering=ordering)[0])
    return {"N": tg.size, "sigma": sigma, "seconds": time.perf_counter() - t0}


def convergence_table() -> list[dict[str, Any]]:
    """Vary ONE axis at a time about BASELINE; report sigma and its drift.

    `rel_change` is a RAW FRACTION (not a percentage): sigma is converged to
    ~1e-6 relative on this problem (see WORKING_GRID's comments), four
    orders of magnitude tighter than the ~1% acceptance criterion, so a
    `%.2f` PERCENTAGE display rounds every row to `0.00%` and makes a
    converged sweep visually indistinguishable from a broken one that
    varies nothing at all. Report/print this in scientific notation instead
    (see `main`).
    """
    rows: list[dict[str, Any]] = []
    for key, values in SWEEPS.items():
        prev: float | None = None
        for v in values:
            params = {**BASELINE, key: v}
            row = {"axis": key, "value": v, **_one(params)}
            row["rel_change"] = None if prev is None else abs(row["sigma"] - prev) / abs(prev)
            prev = row["sigma"]
            rows.append(row)
    return rows


def main() -> None:
    rows = convergence_table()
    print("| axis | value | N | sigma (bohr^2) | rel. change | s |")
    print("|---|---|---|---|---|---|")
    for r in rows:
        rc = "-" if r["rel_change"] is None else f"{r['rel_change']:.3e}"
        print(
            f"| {r['axis']} | {r['value']} | {r['N']} | "
            f"{r['sigma']:.12e} | {rc} | {r['seconds']:.1f} |"
        )
    # Ordering comparison on the REAL Hamiltonian (a small random matrix
    # suggested MMD_AT_PLUS_A roughly halves fill; confirm or refute here).
    # Measured: COLAMD 38.0s vs MMD_AT_PLUS_A 586.6s at BASELINE (N=71476,
    # same sigma to displayed precision) -- MMD is ~15x SLOWER here, the
    # opposite of the small-random-matrix result. COLAMD wins; this is why
    # the spec said to measure on the real problem rather than assume.
    for ordering in ("COLAMD", "MMD_AT_PLUS_A"):
        r = _one(BASELINE, ordering=ordering)
        print(f"ordering={ordering:<14} N={r['N']} sigma={r['sigma']:.12e} {r['seconds']:.1f}s")


if __name__ == "__main__":
    main()
