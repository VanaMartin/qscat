"""V3 (sub-project #9, Task 2) -- THE sweep-reuse benchmark: how much does
`SparseLU.refactor` (reuse the symbolic analysis, re-run only the numeric
factorization) save when sweeping an energy-independent-pattern matrix
`A(E) = E*I - H` across many collision energies?

This is the sweep that `ve_cross_section_2d` now performs internally: `H` is
fixed and the identity only shifts the already-present diagonal, so every
`A(E)` shares one sparsity pattern. Analyze once at the first energy, then
`refactor` at each subsequent energy.

Two routes, per backend, over the real N2 working-grid `H_2D` (`build_h2d`):

  (a) REUSE    -- one `SparseLU(A(E0), backend=...)`, then `refactor(A(E))`
                  for every subsequent energy (symbolic analysis reused).
  (b) NO-REUSE -- a fresh `SparseLU(A(E), backend=...)` per energy (full
                  analyze + factor every time).

On the MUMPS backend the analysis is a real, bounded fraction of the
factorization, so REUSE is faster than NO-REUSE and the saving grows with the
sweep length `M`. On the scipy (SuperLU) backend there is no symbolic-reuse
hook, so `refactor` re-runs `splu`: REUSE and NO-REUSE cost the same (recorded
here to confirm the expected no-speedup, at a small `M` because each SuperLU
factorization of this deck is seconds, not milliseconds).

MEASURE -- this script asserts NOTHING about the speedup. The saving is real
but its size depends on the analysis/factor ratio of the specific matrix and
backend; the point is to report it, not to gate on it.

Run (in the Docker `test` image / a container with system MUMPS + qscat[mumps]):

    uv run python -m benchmarks.sweep_reuse
    uv run python -m benchmarks.sweep_reuse --mumps-m 50 100 --scipy-m 8

The driver writes a Markdown table to
prints it; `--out PATH` also writes it to a file.
"""

from __future__ import annotations

import argparse
import platform
import time
from pathlib import Path
from typing import Any

import numpy as np
import scipy.sparse as sp
from qscat.linalg import SparseLU
from qscat.linalg._mumps_backend import mumps_available

from projects.n2_2d_cross_section.convergence import working_tgrid
from projects.n2_2d_cross_section.hamiltonian2d import MU, build_h2d
from projects.n2_ti_cross_section.vibrational import vibrational_states

_REPO_ROOT = Path(__file__).resolve().parent.parent

# Physical open-channel energy window for the sweep (Ha): every energy is > 0
# so every `A(E)` is actually factored -- the benchmark measures factorization
# reuse, so no energy is allowed to short-circuit below threshold.
_E_LO, _E_HI = 0.01, 1.0


class _Deck:
    """The real N2 working-grid matrix family `A(E) = e_tot*I - H`, built once."""

    def __init__(self) -> None:
        tg = working_tgrid()
        self.n = int(tg.size)
        eps, _chi = vibrational_states(tg.grids[1], MU, 4)
        self._eps0 = float(eps[0])
        self._H = build_h2d(tg)
        self._I = sp.identity(self.n, format="csc", dtype=np.complex128)
        self.nnz = int(self._H.nnz)

    def a(self, e: float) -> sp.csc_matrix[np.complex128]:
        """`A(E) = (E + eps[0])*I - H`, the same driven matrix
        `ve_cross_section_2d` forms (`e_tot = E + eps[v_init=0]`)."""
        return ((e + self._eps0) * self._I - self._H).tocsc()


def _energies(m: int) -> np.ndarray:
    return np.linspace(_E_LO, _E_HI, m)


def _sweep_reuse(deck: _Deck, energies: np.ndarray, backend: str) -> float:
    """(a) One `SparseLU`, then `refactor` per subsequent energy. Wall seconds."""
    t0 = time.perf_counter()
    lu: SparseLU | None = None
    for e in energies:
        a = deck.a(float(e))
        if lu is None:
            lu = SparseLU(a, backend=backend)
        else:
            lu.refactor(a)
    return time.perf_counter() - t0


def _sweep_no_reuse(deck: _Deck, energies: np.ndarray, backend: str) -> float:
    """(b) A fresh `SparseLU` (full analyze+factor) per energy. Wall seconds."""
    t0 = time.perf_counter()
    for e in energies:
        SparseLU(deck.a(float(e)), backend=backend)
    return time.perf_counter() - t0


def _measure(deck: _Deck, backend: str, m: int) -> dict[str, Any]:
    energies = _energies(m)
    reuse_s = _sweep_reuse(deck, energies, backend)
    no_reuse_s = _sweep_no_reuse(deck, energies, backend)
    saved_s = no_reuse_s - reuse_s
    frac = saved_s / no_reuse_s if no_reuse_s > 0 else 0.0
    # Per-energy analysis saving: the reuse route pays one full analyze+factor
    # (the first energy) and M-1 cheap refactors; the saving is spread over the
    # M-1 reused energies.
    per_energy_saved = saved_s / (m - 1) if m > 1 else 0.0
    return {
        "backend": backend,
        "M": m,
        "reuse_s": reuse_s,
        "no_reuse_s": no_reuse_s,
        "saved_s": saved_s,
        "frac_saved": frac,
        "per_energy_saved_s": per_energy_saved,
    }


def _format_table(deck: _Deck, rows: list[dict[str, Any]]) -> str:
    lines: list[str] = []
    lines.append("# Task 2 V3 benchmark: energy-sweep symbolic reuse (`SparseLU.refactor`)")
    lines.append("")
    lines.append(
        f"Machine: {platform.platform()}; Python {platform.python_version()}. "
        f"Matrix: A(E) = (E + eps[0])*I - H_2D on the N2 working grid "
        f"(N={deck.n:,}, H nnz={deck.nnz:,}), complex-symmetric, "
        f"E in [{_E_LO}, {_E_HI}] Ha. REUSE = one SparseLU + refactor per "
        "energy; NO-REUSE = fresh SparseLU per energy. Measured, not asserted."
    )
    lines.append("")
    header = (
        "| backend | M | reuse total (s) | no-reuse total (s) | "
        "saved (s) | fraction saved | per-energy analysis saved (s) |"
    )
    sep = "|---|---|---|---|---|---|---|"
    lines.append(header)
    lines.append(sep)
    for r in rows:
        lines.append(
            f"| {r['backend']} | {r['M']} | {r['reuse_s']:.3f} | "
            f"{r['no_reuse_s']:.3f} | {r['saved_s']:.3f} | "
            f"{r['frac_saved'] * 100:.1f}% | {r['per_energy_saved_s']:.4f} |"
        )
    lines.append("")
    lines.append(
        "The MUMPS rows show the reuse saving (symbolic analysis skipped on "
        "every refactor) growing with M; the scipy rows show ~no saving "
        "(SuperLU has no symbolic-reuse hook, so refactor re-runs `splu`) -- "
        "the expected control."
    )
    lines.append("")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mumps-m",
        nargs="+",
        type=int,
        default=[50, 100],
        help="sweep lengths M for the MUMPS backend (default: 50 100)",
    )
    parser.add_argument(
        "--scipy-m",
        nargs="+",
        type=int,
        default=[8],
        help=(
            "sweep lengths M for the scipy backend (default: 8 -- kept small "
            "because each SuperLU factorization of this deck is seconds)"
        ),
    )
    parser.add_argument(
        "--out",
        default=None,
        help="also write the results table to this path (default: print only)",
    )
    args = parser.parse_args(argv)

    print("Building N2 working-grid deck...", flush=True)
    deck = _Deck()
    print(f"  N={deck.n:,}, H nnz={deck.nnz:,}", flush=True)

    rows: list[dict[str, Any]] = []

    if mumps_available():
        # Warm up the MUMPS library (first factorization pays one-time init
        # that must not land inside a timed sweep).
        SparseLU(deck.a(0.2), backend="mumps")
        for m in args.mumps_m:
            print(f"  [mumps/M={m}] measuring reuse vs no-reuse...", flush=True)
            rows.append(_measure(deck, "mumps", m))
    else:
        print("  MUMPS unavailable -- skipping (run in the Docker test image).", flush=True)

    # scipy warm-up + small-M control (no symbolic reuse -> no speedup).
    SparseLU(deck.a(0.2), backend="scipy")
    for m in args.scipy_m:
        print(f"  [scipy/M={m}] measuring reuse vs no-reuse...", flush=True)
        rows.append(_measure(deck, "scipy", m))

    table = _format_table(deck, rows)
    print("\n" + table, flush=True)
    if args.out is not None:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(table)
        print(f"\nWrote table to {out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
