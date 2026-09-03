from __future__ import annotations

from qscat.core.grids import electronic_grid, nuclear_grid
from qscat.tuning import grid_cost, tensor_cost


def test_grid_cost_reports_point_count():
    g = nuclear_grid(r_max=22.0, n_complex=8, quadrature=12)
    assert grid_cost(g) == {"n_points": g.n}


def test_tensor_cost_reports_exact_unknown_count():
    g_r = electronic_grid(r_max=16.0, order=7, n_complex=6)
    g_R = nuclear_grid(r_max=22.0, n_complex=8, quadrature=12)
    cost = tensor_cost(g_r, g_R)
    assert cost["n_unknowns"] == g_r.n * g_R.n
    assert cost["est_nnz"] > cost["n_unknowns"]  # at least diagonal + neighbors


def test_tensor_cost_estimates_are_monotone_in_grid_size():
    # A coarser and a finer nuclear grid over the same electronic grid: the
    # bigger tensor product must cost more on every estimated axis.
    g_r = electronic_grid(r_max=16.0, order=7, n_complex=6)
    coarse = nuclear_grid(r_max=22.0, n_complex=4, quadrature=8)
    fine = nuclear_grid(r_max=22.0, n_complex=8, quadrature=14)
    assert fine.n > coarse.n

    cost_coarse = tensor_cost(g_r, coarse)
    cost_fine = tensor_cost(g_r, fine)

    assert cost_fine["n_unknowns"] > cost_coarse["n_unknowns"]
    assert cost_fine["est_nnz"] > cost_coarse["est_nnz"]
    assert cost_fine["est_factor_gib"] > cost_coarse["est_factor_gib"]
    assert cost_fine["est_factor_seconds"] > cost_coarse["est_factor_seconds"]


def test_tensor_cost_roughly_reproduces_the_143k_anchor():
    # docs/physics/mumps-sparse-backend.md / CLAUDE.md: ~143k-unknown N2
    # production deck -> SuperLU ~260 s / ~7.4 GiB. The estimate is rough
    # (relative-ranking tool, not a validated performance model), so check
    # it lands within an order of magnitude, not close agreement.
    class _FakeGrid:
        def __init__(self, n: int) -> None:
            self.n = n

    cost = tensor_cost(_FakeGrid(335), _FakeGrid(428))  # 335*428 = 143,380
    assert 26.0 <= cost["est_factor_seconds"] <= 2600.0
    assert 0.74 <= cost["est_factor_gib"] <= 74.0


def test_cost_reports_are_still_plain_dicts():
    # `GridCost`/`TensorCost` name the report keys for a type checker and
    # construct nothing at run time. A caller that stores a report, merges it,
    # or writes it to JSON has to keep working, so what comes back must be an
    # ordinary `dict` with exactly the documented keys -- not a new class that
    # merely looks like one.
    g_r = electronic_grid(r_max=16.0, order=7, n_complex=6)
    g_R = nuclear_grid(r_max=22.0, n_complex=8, quadrature=12)

    one = grid_cost(g_R)
    two = tensor_cost(g_r, g_R)
    assert type(one) is dict
    assert type(two) is dict
    assert set(one) == {"n_points"}
    assert set(two) == {"n_unknowns", "est_nnz", "est_factor_gib", "est_factor_seconds"}
    assert isinstance(one["n_points"], int)
    assert isinstance(two["est_factor_gib"], float)
    assert {**one, "extra": 1}["extra"] == 1  # still merges like any dict
