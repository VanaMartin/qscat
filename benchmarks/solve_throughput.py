"""Repeated-solve (propagation-pattern) throughput: SuperLU vs MUMPS.

`profile_hotpaths.py --td` showed the TD propagation is dominated by the
per-step sparse SOLVE (~82%), not the one-shot factorization (~17%) — the
opposite balance from a TI energy sweep, where `refactor` reuse makes the factor
the cost. MUMPS's documented win is factorization (72×); its solve/back-
substitution is a separate question. A time propagation factors once per Padé
pole and then solves thousands of times against that factor, so the SOLVE cost is
what sets TD wall-clock.

This benchmark isolates that: build a representative complex-symmetric ECS
propagation matrix `A = I - (dt/2)i·H` (the Crank-Nicolson / Padé shift of the
N2 2-D Hamiltonian), factor it once with each backend, then time K repeated
solves. It reports factor time, per-solve time, and the break-even solve count
(when a slower-factoring / faster-solving backend wins for a K-solve run).

Run (SuperLU only on a MUMPS-less box; both in the Docker `test` image):

    uv run python -m benchmarks.solve_throughput [--order N] [--solves K]

Depends only on qscat (matrix built from qscat.model.N2 on a qscat grid).
"""

from __future__ import annotations

import argparse
import time

import numpy as np
import scipy.sparse as sp
from qscat.core.grids import electronic_grid, nuclear_grid
from qscat.dvr import TensorGrid
from qscat.linalg import SparseLU
from qscat.linalg.sparse_lu import mumps_available
from qscat.model import N2


def _propagation_matrix(order: int, n_complex: int, dt: float) -> sp.csr_matrix:
    tg = TensorGrid([
        electronic_grid(r_max=16.0, order=order, n_complex=n_complex),
        nuclear_grid(r_max=22.0, quadrature=10, n_complex=n_complex),
    ])
    H = N2.hamiltonian(tg).tocsr()
    n = H.shape[0]
    # One Crank-Nicolson / Padé-pole shift: complex-symmetric like H itself.
    A = (sp.identity(n, dtype=np.complex128, format="csr") - (0.5j * dt) * H).tocsr()
    return A


def _measure(A: sp.csr_matrix, backend: str, n_solves: int) -> dict[str, float]:
    rng = np.random.default_rng(0)
    n = A.shape[0]
    t0 = time.perf_counter()
    lu = SparseLU(A, backend=backend)
    factor_s = time.perf_counter() - t0

    # A fresh RHS each solve, mimicking a propagation (state changes every step).
    rhs = [
        rng.standard_normal(n) + 1j * rng.standard_normal(n) for _ in range(n_solves)
    ]
    t0 = time.perf_counter()
    for b in rhs:
        lu.solve(b)
    solve_total_s = time.perf_counter() - t0
    return {
        "n": float(n),
        "factor_s": factor_s,
        "per_solve_ms": 1e3 * solve_total_s / n_solves,
        "solve_total_s": solve_total_s,
        "backend_used": lu.backend_used,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--order", type=int, default=7, help="electronic DVR order")
    ap.add_argument("--n-complex", type=int, default=5, help="ECS tail elements")
    ap.add_argument("--dt", type=float, default=1.0)
    ap.add_argument("--solves", type=int, default=200, help="repeated solves to time")
    args = ap.parse_args(argv)

    A = _propagation_matrix(args.order, args.n_complex, args.dt)
    backends = ["scipy"] + (["mumps"] if mumps_available() else [])
    print(f"matrix N={A.shape[0]}, nnz={A.nnz}, {args.solves} solves\n")

    results = {}
    for backend in backends:
        r = _measure(A, backend, args.solves)
        results[backend] = r
        print(
            f"{backend:>6}: factor {r['factor_s']:7.3f} s | "
            f"per-solve {r['per_solve_ms']:7.3f} ms | "
            f"{args.solves} solves {r['solve_total_s']:7.3f} s"
        )

    if "mumps" in results:
        s, m = results["scipy"], results["mumps"]
        d_factor = m["factor_s"] - s["factor_s"]
        d_solve = s["per_solve_ms"] - m["per_solve_ms"]  # >0 if MUMPS solve faster
        print()
        print(f"factor: MUMPS {'faster' if d_factor < 0 else 'slower'} by {abs(d_factor):.3f} s")
        if d_solve > 0:
            breakeven = d_factor / (d_solve / 1e3) if d_solve else float("inf")
            print(f"solve : MUMPS faster by {d_solve:.3f} ms/solve")
            print(f"break-even: MUMPS wins a run of > {breakeven:.0f} solves")
        else:
            print(f"solve : MUMPS SLOWER by {abs(d_solve):.3f} ms/solve (SuperLU wins the solve)")
    else:
        print("\n(mumps unavailable here — run in the Docker `test` image for the comparison)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
