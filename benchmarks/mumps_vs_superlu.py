"""V3 (Task 4) -- THE benchmark: does MUMPS actually beat SuperLU on the real
N2 2-D matrices?

Builds a real N2 2-D Hamiltonian `H_2D` via
`projects.n2_2d_cross_section.hamiltonian2d.build_h2d` on a representative
grid, forms the driven complex-symmetric matrix `A = (E_tot*I - H)`, and
factors + solves it with `SparseLU(A, backend="scipy")` (SuperLU) versus
`SparseLU(A, backend="mumps")` (MUMPS complex-symmetric SYM=2). It reports, per
(grid, backend): factor time, solve time, peak RSS, fill_factor, and the
ordering each engine used.

MEASURE the win -- this script asserts NOTHING about the speedup. If MUMPS-seq
does not beat SuperLU, that is a real, reportable finding; SuperLU remains the
safe fallback (`backend="auto"` -> SuperLU when MUMPS is absent), so nothing
regresses either way.

Memory: each (grid, backend) measurement runs in a FRESH subprocess so
`resource.getrusage(RUSAGE_SELF).ru_maxrss` (the process high-water mark) is
attributable to that one factorization and not contaminated by the other
backend's factors. This deliberately does NOT call `SparseLU.memory_bytes()`,
which materializes SuperLU's L/U factors and costs +6 GB at production scale
(see `qscat.linalg.sparse_lu`'s module docstring).

Run (in the Docker `test` image / a container with system MUMPS + qscat[mumps]):

    uv run python -m benchmarks.mumps_vs_superlu                # working + td grids
    uv run python -m benchmarks.mumps_vs_superlu --include-production
    uv run python -m benchmarks.mumps_vs_superlu --grids working

The driver writes a Markdown table to
prints it; `--out PATH` also writes it to a file.
"""

from __future__ import annotations

import argparse
import json
import platform
import resource
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

# --- grid registry -----------------------------------------------------------
# Each entry builds a real N2 TensorGrid. Kept as thunks so the driver never
# builds a grid itself (only the worker subprocess does), and so importing this
# module is cheap.

_REPO_ROOT = Path(__file__).resolve().parent.parent

# Representative collision energy for the driven matrix A = (E_tot*I - H); the
# same E=0.2 Ha the convergence study and Group E anchors use. E_tot = E +
# eps[0] (neutral ground vibrational energy), matching ve_cross_section_2d.
_STUDY_E = 0.2


def _build_grid(name: str) -> Any:
    """Build one named N2 TensorGrid (worker side only)."""
    from qscat.dvr import TensorGrid

    from projects.n2_2d_cross_section.convergence import working_tgrid
    from projects.n2_2d_cross_section.electronic_grid import n2_electronic_grid
    from projects.n2_2d_td_cross_section.convergence import td_working_tgrid
    from projects.n2_ti_cross_section.nuclear_grid import n2_nuclear_grid

    if name == "working":  # sub-project #6 TI working grid (~27k)
        return working_tgrid()
    if name == "td":  # sub-project #7 TD working grid (~47k)
        return td_working_tgrid()
    if name == "production":  # real N2 production deck (~143k)
        return TensorGrid(
            [
                n2_electronic_grid(r_max=98.0, order=8, n_complex=15),
                n2_nuclear_grid(quadrature=14, r_max=40.0, n_complex=10),
            ]
        )
    raise ValueError(f"unknown grid {name!r}")


def _maxrss_mb() -> float:
    """Process peak resident set size in MB (ru_maxrss is KB on Linux, bytes on macOS)."""
    kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    raw_bytes = kb if platform.system() == "Darwin" else kb * 1024
    return raw_bytes / (1024.0 * 1024.0)


# --- worker: one (grid, backend) measurement in an isolated process ----------


def _measure(grid_name: str, backend: str) -> dict[str, Any]:
    """Build A, factor with `backend`, solve one RHS; return timing/memory."""
    import numpy as np
    import scipy.sparse as sp
    from qscat.linalg import SparseLU

    from projects.n2_2d_cross_section.hamiltonian2d import MU, build_h2d
    from projects.n2_ti_cross_section.vibrational import vibrational_states

    tgrid = _build_grid(grid_name)
    n = int(tgrid.size)
    eps, _chi = vibrational_states(tgrid.grids[1], MU, 4)
    e_tot = _STUDY_E + float(eps[0])

    H = build_h2d(tgrid)
    ident = sp.identity(n, format="csc", dtype=np.complex128)
    A = (e_tot * ident - H).tocsc()
    nnz = int(A.nnz)

    rss_before = _maxrss_mb()

    t0 = time.perf_counter()
    lu = SparseLU(A, backend=backend)
    factor_s = time.perf_counter() - t0

    rss_after_factor = _maxrss_mb()

    rng = np.random.default_rng(0)
    b = (rng.standard_normal(n) + 1j * rng.standard_normal(n)).astype(np.complex128)
    t1 = time.perf_counter()
    x = lu.solve(b)
    solve_s = time.perf_counter() - t1

    residual = float(np.linalg.norm(A @ x - b) / np.linalg.norm(b))

    # The MUMPS matrix type actually driven: SYM=2 (complex-symmetric, single
    # upper triangle) when `symmetric` was detected/forced True, else SYM=0
    # (general unsymmetric, full matrix). SuperLU has no symmetric mode, so
    # its SYM cell is "-". This is the audit column: it makes visible which
    # mode ran, so a fill_factor is never misread as the wrong storage count.
    if lu.backend_used == "mumps":
        sym = "SYM=2" if lu.symmetric else "SYM=0"
    else:
        sym = "-"

    return {
        "grid": grid_name,
        "backend": backend,
        "backend_used": lu.backend_used,
        "symmetric": bool(lu.symmetric),
        "sym": sym,
        "N": n,
        "nnz": nnz,
        "factor_s": factor_s,
        "solve_s": solve_s,
        "fill_factor": float(lu.fill_factor),
        "ordering_used": lu.ordering_used,
        "peak_rss_mb": rss_after_factor,
        "factor_rss_delta_mb": rss_after_factor - rss_before,
        "residual": residual,
    }


# --- driver: spawn a subprocess per measurement, assemble the table ----------


def _run_worker(grid_name: str, backend: str) -> dict[str, Any] | None:
    """Run one measurement in a fresh `python -m benchmarks.mumps_vs_superlu
    --worker` subprocess; return its parsed JSON, or None on failure."""
    cmd = [sys.executable, "-m", "benchmarks.mumps_vs_superlu", "--worker", grid_name, backend]
    print(f"  [{grid_name}/{backend}] running...", flush=True)
    proc = subprocess.run(cmd, cwd=_REPO_ROOT, capture_output=True, text=True)
    if proc.returncode != 0:
        print(f"  [{grid_name}/{backend}] FAILED (rc={proc.returncode}):", flush=True)
        print(proc.stderr[-2000:], flush=True)
        return None
    line = proc.stdout.strip().splitlines()[-1]
    return json.loads(line)


def _backends() -> list[str]:
    """`scipy` always; `mumps` only where importable."""
    from qscat.linalg._mumps_backend import mumps_available

    return ["scipy", "mumps"] if mumps_available() else ["scipy"]


def _format_table(rows: list[dict[str, Any]]) -> str:
    lines: list[str] = []
    lines.append("# Task 4 V3 benchmark: MUMPS (SYM=2) vs SuperLU on real N2 2-D matrices")
    lines.append("")
    lines.append(
        f"Machine: {platform.platform()}; Python {platform.python_version()}. "
        f"Matrix: A = (E_tot*I - H_2D), complex-symmetric, E={_STUDY_E} Ha. "
        "Each row measured in an isolated subprocess (peak RSS = process "
        "ru_maxrss high-water mark; NOT SparseLU.memory_bytes())."
    )
    lines.append("")
    header = (
        "| grid | backend | SYM | N | nnz | factor (s) | solve (s) | "
        "peak RSS (MB) | factor RSS delta (MB) | fill_factor | ordering | residual |"
    )
    sep = "|---|---|---|---|---|---|---|---|---|---|---|---|"
    lines.append(header)
    lines.append(sep)
    for r in rows:
        lines.append(
            f"| {r['grid']} | {r['backend_used']} | {r.get('sym', '-')} | "
            f"{r['N']:,} | {r['nnz']:,} | "
            f"{r['factor_s']:.3f} | {r['solve_s']:.4f} | "
            f"{r['peak_rss_mb']:.0f} | {r['factor_rss_delta_mb']:.0f} | "
            f"{r['fill_factor']:.2f} | {r['ordering_used']} | {r['residual']:.2e} |"
        )
    lines.append("")

    # Per-grid speed/memory comparison (measured, no assertion).
    by_grid: dict[str, dict[str, dict[str, Any]]] = {}
    for r in rows:
        by_grid.setdefault(r["grid"], {})[r["backend_used"]] = r
    lines.append("## Measured MUMPS-vs-SuperLU comparison (per grid)")
    lines.append("")
    for grid, per in by_grid.items():
        if "mumps" in per and "scipy" in per:
            m, s = per["mumps"], per["scipy"]
            fac = s["factor_s"] / m["factor_s"] if m["factor_s"] > 0 else float("inf")
            mem = s["peak_rss_mb"] / m["peak_rss_mb"] if m["peak_rss_mb"] > 0 else float("inf")
            won = "MUMPS faster" if fac > 1.0 else "SuperLU faster"
            lines.append(
                f"- **{grid}** (N={m['N']:,}): MUMPS mode {m.get('sym', '?')}; "
                f"factor speedup SuperLU/MUMPS = "
                f"{fac:.2f}x ({won}); peak-RSS ratio SuperLU/MUMPS = {mem:.2f}x; "
                f"MUMPS ordering={m['ordering_used']}, SuperLU ordering={s['ordering_used']}."
            )
        else:
            only = next(iter(per))
            lines.append(
                f"- **{grid}** (N={per[only]['N']:,}): only `{only}` measured "
                "(MUMPS unavailable in this environment -- run in the Docker test image)."
            )
    lines.append("")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--worker",
        nargs=2,
        metavar=("GRID", "BACKEND"),
        help="internal: run one measurement and print a JSON line",
    )
    parser.add_argument(
        "--grids",
        nargs="+",
        default=["working", "td"],
        help="grids to benchmark (default: working td)",
    )
    parser.add_argument(
        "--include-production",
        action="store_true",
        help="also benchmark the ~143k production deck (needs >~14 GB RAM for SuperLU)",
    )
    parser.add_argument(
        "--out",
        default=None,
        help="also write the results table to this path (default: print only)",
    )
    args = parser.parse_args(argv)

    if args.worker is not None:
        grid_name, backend = args.worker
        result = _measure(grid_name, backend)
        print(json.dumps(result))
        return 0

    grids = list(args.grids)
    if args.include_production and "production" not in grids:
        grids.append("production")

    backends = _backends()
    print(f"Benchmarking grids={grids} backends={backends}", flush=True)
    rows: list[dict[str, Any]] = []
    for grid_name in grids:
        for backend in backends:
            r = _run_worker(grid_name, backend)
            if r is not None:
                rows.append(r)

    table = _format_table(rows)
    print("\n" + table, flush=True)
    if args.out is not None:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(table)
        print(f"\nWrote table to {out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
