"""Profile the qscat hot paths to confirm where time actually goes.

Roadmap Part 5 ("profile before cutting"): the optimize-in-Rust stage should be
driven by measurement, not intuition. This runs cProfile over the two dominant
code paths on a representative N2 problem and prints the top functions by
cumulative time:

  * TI  — `ScatteringProblem.ve_cross_section` over an energy sweep. The cost is
          dominated by the sparse LU (`qscat.linalg.SparseLU`) factor+solve per
          energy; this is optimization target #1.
  * TD  — a short `td_ve_cross_section` propagation. The cost is the per-step
          cached-factor solve plus the extractor `record` projections
          (optimization target #2).

Run:  uv run python -m benchmarks.profile_hotpaths [--td] [--top N]

This is a measurement tool, not a gate; it imports only qscat (the library),
never validation/projects, so it profiles the shipped code paths.
"""

from __future__ import annotations

import argparse
import cProfile
import pstats
from io import StringIO

import numpy as np
from qscat.core import ScatteringProblem
from qscat.core.grids import electronic_grid, nuclear_grid
from qscat.dvr import TensorGrid
from qscat.model import N2


def _problem(order: int, n_complex: int) -> ScatteringProblem:
    grid = TensorGrid(
        [
            electronic_grid(r_max=16.0, order=order, n_complex=n_complex),
            nuclear_grid(r_max=22.0, quadrature=10, n_complex=n_complex),
        ]
    )
    return ScatteringProblem(grid=grid, model=N2, n_vib=4, v_init=0)


def _run_ti(prob: ScatteringProblem) -> None:
    E = np.linspace(0.04, 0.20, 9)
    prob.ve_cross_section(vprimes=[0, 1, 2], E=E)


def _run_td(prob: ScatteringProblem) -> None:
    prob.td_ve_cross_section(
        vprimes=[0, 1],
        E=np.array([0.10, 0.15]),
        dt=1.0,
        n_steps=400,
        wp_in={"r0": 8.0, "p0": -0.35, "sigma": 4.0},
        wp_out={"r0_out": 8.0, "p0_out": 0.35, "sigma_out": 4.0},
        order=3,
    )


def _profile(fn, prob: ScatteringProblem, top: int, label: str) -> None:
    pr = cProfile.Profile()
    pr.enable()
    fn(prob)
    pr.disable()
    buf = StringIO()
    stats = pstats.Stats(pr, stream=buf).sort_stats("cumulative")
    stats.print_stats(top)
    print(f"\n===== {label}: top {top} by cumulative time =====")
    print(buf.getvalue())


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--td", action="store_true", help="also profile a TD propagation")
    ap.add_argument("--top", type=int, default=25, help="number of functions to show")
    ap.add_argument("--order", type=int, default=7, help="electronic DVR order")
    ap.add_argument("--n-complex", type=int, default=5, help="ECS tail elements")
    args = ap.parse_args()

    prob = _problem(args.order, args.n_complex)
    _profile(_run_ti, prob, args.top, "TI ve_cross_section (SparseLU factor+solve)")
    if args.td:
        _profile(_run_td, prob, args.top, "TD propagation (per-step solve + record)")


if __name__ == "__main__":
    main()
