"""Cost model for the FEM-DVR-ECS discretisation tuner.

`grid_cost` reports a single 1-D grid's point count. `tensor_cost` reports
ROUGH, ANCHORED-BUT-NOT-PREDICTIVE estimates of the tensor-product (2-D)
problem size and sparse-LU factorization cost -- their job is RELATIVE
RANKING (is candidate grid A cheaper than candidate grid B), never an
absolute runtime/memory prediction. See the docstrings below for exactly
what they are anchored to and how rough that anchor is.

Both reports are fixed-shape: `GridCost` and `TensorCost` name every key and
its type, so a caller can annotate a report it stores or passes on, and a
type checker catches a mistyped key rather than a `KeyError` at run time.
They are plain `TypedDict`s -- what the two functions return IS a `dict`,
subscripted and iterated exactly as before.
"""

from __future__ import annotations

from typing import TypedDict

from qscat.dvr import FemDvrEcsGrid

__all__ = ["GridCost", "TensorCost", "grid_cost", "tensor_cost"]


class GridCost(TypedDict):
    """`grid_cost`'s report: one 1-D grid's discretisation cost.

    - `n_points`: the grid's DVR point count (`FemDvrEcsGrid.n`). EXACT --
      the number of unknowns a 1-D problem on this grid carries, and the
      quantity a grid-vs-deck size comparison is made on.
    """

    n_points: int


class TensorCost(TypedDict):
    """`tensor_cost`'s report: the tensor-product (2-D) problem's size, plus
    ROUGH anchored estimates of its sparse-LU cost. Only `n_unknowns` is
    exact; the three `est_*` fields are for RELATIVE ranking of candidate
    grids, never absolute prediction -- see `tensor_cost`.

    - `n_unknowns`: `g_r.n * g_R.n`, the tensor-grid unknown count. EXACT.
    - `est_nnz`: estimated nonzeros in the assembled 2-D matrix, at a
      constant `~22.9` nonzeros per row.
    - `est_factor_gib`: estimated peak resident memory of a SuperLU
      factorization, in GiB.
    - `est_factor_seconds`: estimated SuperLU factorization wall-clock time,
      in seconds.
    """

    n_unknowns: int
    est_nnz: int
    est_factor_gib: float
    est_factor_seconds: float


# --- anchor: the N2 2-D production deck, docs/physics/mumps-sparse-backend.md
# ("N = 335*428 = 143,380, nnz = 3,276,450, ~22.9 nonzeros/row") and CLAUDE.md
# ("the 143k production deck (3.6 s / 0.8 GB vs 260 s / 7.4 GB)" for
# MUMPS/SuperLU factor time and peak RSS). We anchor to the SuperLU numbers
# (CLAUDE.md's rounded 260 s / 7.4 GiB) since scipy/SuperLU is the always-
# available fallback backend -- a MUMPS-backed factorization of the same
# matrix is MUCH cheaper (per CLAUDE.md, ~72x faster / ~9x less peak memory),
# so treat `est_factor_seconds`/`est_factor_gib` below as a conservative
# (worst-case-backend) estimate, not what an actual `backend="auto"` run
# would cost if MUMPS is installed.
_ANCHOR_N = 143_000.0
_ANCHOR_SECONDS = 260.0
_ANCHOR_GIB = 7.4

# nnz/row measured on that same anchor deck (3,276,450 / 143,380 ~= 22.85);
# used as a simple constant-bandwidth-ish factor, not a function of grid
# structure (a rougher estimate still, but the ROUGH point stands).
_NNZ_PER_ROW = 22.9

# Empirical fit of SuperLU factor time/peak-RSS vs N across the three
# measured decks in docs/physics/mumps-sparse-backend.md (N = 26,857 /
# 47,188 / 143,380): factor time grows close to N**2.3, peak RSS closer to
# N**1.5 (plausible for a 2-D sparse LU with COLAMD fill-in growing faster
# than O(N), between the O(N log N) nested-dissection ideal and O(N**2)
# dense-like blowup). These exponents are chosen to roughly track that
# empirical trend, not derived from first principles -- see the module
# docstring: RELATIVE RANKING is the job, not an absolute prediction.
_SECONDS_EXPONENT = 2.3
_GIB_EXPONENT = 1.5

_SECONDS_COEFF = _ANCHOR_SECONDS / _ANCHOR_N**_SECONDS_EXPONENT
_GIB_COEFF = _ANCHOR_GIB / _ANCHOR_N**_GIB_EXPONENT


def grid_cost(grid: FemDvrEcsGrid) -> GridCost:
    """A single grid's discretisation cost: just its DVR point count."""
    return {"n_points": grid.n}


def tensor_cost(g_r: FemDvrEcsGrid, g_R: FemDvrEcsGrid) -> TensorCost:
    """ROUGH, ANCHORED cost ESTIMATES for the tensor-product 2-D problem.

    `n_unknowns = g_r.n * g_R.n` is exact (the tensor-grid unknown count).
    `est_nnz`, `est_factor_gib`, and `est_factor_seconds` are simple monotone
    scalings of `n_unknowns`, calibrated so the CLAUDE.md/`docs/physics/
    mumps-sparse-backend.md` ~143k-unknown N2 production-deck SuperLU
    numbers (260 s factor time, 7.4 GiB peak RSS, ~22.9 nonzeros/row) are
    ROUGHLY reproduced at that anchor size. They are NOT a validated
    performance model: no assembly-time-dependence on quadrature order,
    ECS-tail bandwidth, or ordering algorithm is captured, and a MUMPS-
    backed factorization of the same matrix would be far cheaper (see
    `grid_cost`'s module docstring). Use these ONLY to rank candidate grids
    against each other (bigger `n_unknowns` always costs more here by
    construction), never to predict a real wall-clock time or memory
    footprint.
    """
    n_unknowns = g_r.n * g_R.n
    est_nnz = round(n_unknowns * _NNZ_PER_ROW)
    est_factor_gib = _GIB_COEFF * float(n_unknowns) ** _GIB_EXPONENT
    est_factor_seconds = _SECONDS_COEFF * float(n_unknowns) ** _SECONDS_EXPONENT
    return {
        "n_unknowns": n_unknowns,
        "est_nnz": est_nnz,
        "est_factor_gib": est_factor_gib,
        "est_factor_seconds": est_factor_seconds,
    }
