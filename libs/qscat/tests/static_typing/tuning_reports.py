"""Static-typing fixture: the tuner's fixed-shape reports keep their fields.

This module is NEVER executed. `test_static_types.py` runs a type checker over
it, and every `assert_type` below fails there if an inferred type drifts. A
runtime test cannot make these assertions -- key types and narrowing exist only
at check time, and at run time each report is an ordinary `dict` that would
answer `isinstance(..., dict)` no matter how it is annotated.

Each function takes its inputs as parameters rather than building a grid, so
nothing here needs a real solve to be type-checked.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Literal, assert_type

from qscat.dvr import FemDvrEcsGrid
from qscat.tuning import (
    GridCost,
    Refine2dReport,
    Refine2dStep,
    RefinementCoordinate,
    TensorCost,
    grid_cost,
    refine_to_2d_convergence,
    tensor_cost,
)


def check_grid_cost(g: FemDvrEcsGrid) -> None:
    cost = grid_cost(g)
    assert_type(cost, GridCost)
    assert_type(cost["n_points"], int)

    # A key the report does not carry is a type error now, instead of a
    # `KeyError` at run time. The suppression IS the assertion: strict mode
    # reports an UNUSED ignore -- and so fails this fixture -- the moment
    # `grid_cost` goes back to returning an open `dict[str, Any]`.
    _ = cost["n_pointz"]  # type: ignore[typeddict-item]


def check_tensor_cost(g_r: FemDvrEcsGrid, g_R: FemDvrEcsGrid) -> None:
    cost = tensor_cost(g_r, g_R)
    assert_type(cost, TensorCost)
    # The exact count is an int; the three anchored estimates are floats.
    assert_type(cost["n_unknowns"], int)
    assert_type(cost["est_nnz"], int)
    assert_type(cost["est_factor_gib"], float)
    assert_type(cost["est_factor_seconds"], float)
    _ = cost["est_factor_minutes"]  # type: ignore[typeddict-item]


def check_refine_report(
    observable: Callable[[FemDvrEcsGrid, FemDvrEcsGrid], float],
    g_r: FemDvrEcsGrid,
    g_R: FemDvrEcsGrid,
) -> None:
    g_r2, g_R2, detail = refine_to_2d_convergence(observable, g_r, g_R)
    assert_type(g_r2, FemDvrEcsGrid)
    assert_type(g_R2, FemDvrEcsGrid)
    assert_type(detail, Refine2dReport)
    assert_type(detail["converged"], bool)
    assert_type(detail["final_value"], float)
    assert_type(detail["iterations"], list[Refine2dStep])

    step = detail["iterations"][0]
    assert_type(step["value"], float)
    assert_type(step["rel_move"], float)
    # The coordinate is the two-value alias, not a bare `str`, so a caller
    # dispatching on it has both spellings checked.
    assert_type(step["coordinate"], RefinementCoordinate)
    assert_type(step["coordinate"], Literal["electronic", "nuclear"])
