"""Tests for `qscat.dvr.TensorGrid`: geometry, broadcasting, the ECS real-region
mask, and separable-state construction.
"""

from __future__ import annotations

import numpy as np
import pytest
from qscat.dvr import ElementSpec, FemDvrEcsGrid, GridSpec, TensorGrid


def _real_grid(q: int, n_el: int, length: float = 1.0) -> FemDvrEcsGrid:
    return FemDvrEcsGrid(
        GridSpec(quadrature=q, elements=[ElementSpec(length) for _ in range(n_el)])
    )


def _ecs_grid(q: int, n_real: int, n_cplx: int) -> FemDvrEcsGrid:
    els = [ElementSpec(1.0) for _ in range(n_real)]
    els += [ElementSpec(1.0, 35.0) for _ in range(n_cplx)]
    return FemDvrEcsGrid(GridSpec(quadrature=q, elements=els))


def test_shape_size_and_ndim() -> None:
    ga, gb, gc = _real_grid(6, 2), _real_grid(5, 3), _real_grid(4, 2)
    tg = TensorGrid([ga, gb, gc])
    assert tg.ndim == 3
    assert tg.shape == (ga.n, gb.n, gc.n)
    assert tg.size == ga.n * gb.n * gc.n
    assert tg.grids == (ga, gb, gc)


def test_points_are_broadcastable_not_meshgrid() -> None:
    ga, gb = _real_grid(6, 2), _real_grid(5, 3)
    tg = TensorGrid([ga, gb])
    pa, pb = tg.points()
    assert pa.shape == (ga.n, 1)
    assert pb.shape == (1, gb.n)
    # broadcasting them together reproduces the full grid
    assert np.broadcast_shapes(pa.shape, pb.shape) == (ga.n, gb.n)
    assert np.allclose(pa.ravel(), ga.points)
    assert np.allclose(pb.ravel(), gb.points)


def test_points_for_d1_is_plain_1d() -> None:
    g = _real_grid(6, 3)
    (p,) = TensorGrid([g]).points()
    assert p.shape == (g.n,)
    assert np.allclose(p, g.points)


def test_weights_are_broadcastable_and_match_points_shape() -> None:
    ga, gb = _real_grid(6, 2), _real_grid(5, 3)
    tg = TensorGrid([ga, gb])
    wa, wb = tg.weights()
    pa, pb = tg.points()
    assert wa.shape == pa.shape == (ga.n, 1)
    assert wb.shape == pb.shape == (1, gb.n)
    assert np.allclose(wa.ravel(), ga.weights)
    assert np.allclose(wb.ravel(), gb.weights)


def test_weights_are_complex_on_the_ecs_tail() -> None:
    """A weight that coerced to float would silently kill the ECS tail."""
    g = _ecs_grid(6, 3, 2)
    (w,) = TensorGrid([g]).weights()
    assert np.abs(w.imag).max() > 1e-6


def test_sqrt_weights_is_sqrt_of_weights_per_axis() -> None:
    ga, gb = _ecs_grid(6, 3, 2), _ecs_grid(5, 2, 2)
    tg = TensorGrid([ga, gb])
    for w, sw in zip(tg.weights(), tg.sqrt_weights(), strict=True):
        assert np.allclose(sw**2, w)
        assert sw.shape == w.shape


def test_weights_product_over_axes_matches_explicit_outer_product() -> None:
    """The broadcast product `w_a * w_b` (what a physics routine actually
    computes) must equal the explicit outer product of the flat 1-D weights,
    in the same C-order flattening `points()`/`real_mask()` use.
    """
    ga, gb = _real_grid(6, 2), _real_grid(5, 3)
    tg = TensorGrid([ga, gb])
    wa, wb = tg.weights()
    broadcast_flat = np.broadcast_to(wa * wb, tg.shape).ravel()
    explicit_outer = np.outer(ga.weights, gb.weights).ravel()
    assert np.allclose(broadcast_flat, explicit_outer)


def test_sqrt_weights_converts_a_separable_function_to_basis_coefficients() -> None:
    """The intended use: `c = f(points) * sqrt_weight_product` converts a
    known function to FEM-DVR basis coefficients. For a separable
    `f(x, y) = f_a(x) * f_b(y)`, the coefficient array must itself be
    separable into the two axes' own `f_d(x_d) * sqrt(w_d)` factors.
    """
    ga, gb = _ecs_grid(7, 3, 2), _ecs_grid(6, 2, 2)
    tg = TensorGrid([ga, gb])
    pa, pb = tg.points()
    swa, swb = tg.sqrt_weights()

    def f_a(x: np.ndarray) -> np.ndarray:
        return np.exp(-0.1 * x**2)

    def f_b(y: np.ndarray) -> np.ndarray:
        return np.exp(-0.2 * (y - 0.5) ** 2)

    ca_flat = (f_a(pa) * swa).ravel()
    cb_flat = (f_b(pb) * swb).ravel()
    got = np.outer(ca_flat, cb_flat).ravel()

    want = ((f_a(pa) * swa) * (f_b(pb) * swb)).ravel()
    assert np.allclose(got, want)


def test_real_mask_is_and_across_dimensions() -> None:
    ga = _ecs_grid(6, 3, 2)
    gb = _ecs_grid(5, 2, 2)
    tg = TensorGrid([ga, gb])
    mask = tg.real_mask()
    assert mask.shape == (tg.size,)
    assert mask.dtype == np.bool_
    ma = ga.real_points <= ga.R0
    mb = gb.real_points <= gb.R0
    assert np.array_equal(mask, np.outer(ma, mb).ravel())
    # a genuine mixture -- otherwise the test proves nothing
    assert 0 < int(mask.sum()) < mask.size


def test_real_mask_all_true_when_no_ecs() -> None:
    tg = TensorGrid([_real_grid(6, 2), _real_grid(5, 3)])
    assert tg.real_mask().all()


def test_outer_builds_separable_state_in_c_order() -> None:
    ga, gb = _real_grid(6, 2), _real_grid(5, 3)
    tg = TensorGrid([ga, gb])
    rng = np.random.default_rng(0)
    a = rng.standard_normal(ga.n) + 1j * rng.standard_normal(ga.n)
    b = rng.standard_normal(gb.n) + 1j * rng.standard_normal(gb.n)
    psi = tg.outer([a, b])
    assert psi.shape == (tg.size,)
    assert np.allclose(psi.reshape(tg.shape), np.outer(a, b))


def test_outer_rejects_wrong_length() -> None:
    tg = TensorGrid([_real_grid(6, 2), _real_grid(5, 3)])
    with pytest.raises(ValueError, match="expected 2"):
        tg.outer([np.ones(3)])


def test_rejects_empty_grid_list() -> None:
    with pytest.raises(ValueError, match="at least one"):
        TensorGrid([])
