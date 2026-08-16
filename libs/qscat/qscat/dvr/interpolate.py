"""Evaluate a FEM-DVR-ECS field at arbitrary points (Lagrange interpolation).

The DVR representation is not diagonal for off-node evaluation: the value of a
state at an arbitrary point ``x`` is

    f(x) = sum_i state_i * b_i(z(x)),   b_i(z) = L_i(z) / sqrt(w_i),

where ``L_i`` is the element-local Lagrange cardinal function of node ``i``,
``z(x)`` is the exterior-complex-scaling image of ``x`` (`qscat.ecs.ecs_map`),
and ``w_i`` is the bridge-summed DVR weight. Because ``L_i`` is element-local,
each evaluation point only couples to the ``nq`` nodes of the one element that
contains it -- so the operator mapping a state to its values on a set of sample
points is SPARSE. `dvr_interpolation_matrix` assembles that sparse operator once;
applying it to any state (or, tensored per-axis, to a 2-D state) is the fast,
repeatable projection used for wavefunction visualisation.

Ported from eMoScat's ``FemDvrEcsGrid::basis_function_value`` +
``EquidistantProjector`` (the bridge-factor treatment of the two dropped
Dirichlet endpoints is reproduced exactly). Sibling of
`qscat.dvr.dvr_first_derivative_at_node` (value-at-a-point vs derivative-at-a-node).
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp

from qscat.ecs import ecs_map
from qscat.exceptions import GridError

from .grid import FemDvrEcsGrid

__all__ = ["dvr_interpolation_matrix"]


def _element_boundaries(grid: FemDvrEcsGrid) -> npt.NDArray[np.float64]:
    """Cumulative real element boundaries ``ar`` (length tnel+1); ``ar[0]=x_min``."""
    lengths = np.array([el.length for el in grid.spec.elements], dtype=np.float64)
    return grid.spec.x_min + np.concatenate(([0.0], np.cumsum(lengths)))


def dvr_interpolation_matrix(grid: FemDvrEcsGrid, x: npt.ArrayLike) -> sp.csr_matrix:
    """Sparse operator ``P`` s.t. ``(P @ state)[k]`` is the field value at ``x[k]``.

    Parameters
    ----------
    grid : FemDvrEcsGrid
        The FEM-DVR-ECS radial grid the state is represented on.
    x : array_like
        Real sample coordinates (unscaled), each within ``[x_min, x_max]``.

    Returns
    -------
    scipy.sparse.csr_matrix
        Complex matrix of shape ``(len(x), grid.n)``. ``P @ state`` gives the
        state's values at ``x``; ``P`` has at most ``nq`` nonzeros per row (the
        nodes of the containing element) and is built once, reused per frame.

    Raises
    ------
    GridError
        If any sample point lies outside ``[x_min, x_max]``.
    """
    xs = np.atleast_1d(np.asarray(x, dtype=np.float64))
    spec = grid.spec
    ar = _element_boundaries(grid)
    x_min, x_max = ar[0], ar[-1]
    tol = 1e-9 * max(1.0, abs(x_max))
    if np.any(xs < x_min - tol) or np.any(xs > x_max + tol):
        raise GridError(f"sample points must lie within [{x_min}, {x_max}] (the grid extent)")
    xs = np.clip(xs, x_min, x_max)

    R0 = grid.R0
    pts = grid.points  # complex ECS node positions (retained)
    wsqrt = np.sqrt(grid.weights)  # 1/sqrt(w_i) scaling of each basis function
    first_angle = spec.elements[0].angle_deg
    last_angle = spec.elements[-1].angle_deg
    z_min = ecs_map(x_min, R0, first_angle)  # dropped Dirichlet endpoint images
    z_max = ecs_map(x_max, R0, last_angle)

    # Element index of each sample point (ar strictly increasing).
    elem = np.clip(np.searchsorted(ar, xs, side="right") - 1, 0, len(ar) - 2)

    rows: list[int] = []
    cols: list[int] = []
    vals: list[complex] = []
    for k, (xk, e) in enumerate(zip(xs, elem, strict=True)):
        gidx = grid.element_maps[e][1]  # consecutive global node indices of element e
        z = ecs_map(float(xk), R0, spec.elements[e].angle_deg)
        zl = pts[gidx]  # this element's complex node positions
        is_first = int(gidx[0]) == 0
        is_last = int(gidx[-1]) == grid.n - 1
        for pos, i in enumerate(gidx):
            # Lagrange cardinal L_i(z) over this element's nodes (exclude self).
            other = np.arange(len(gidx)) != pos
            lag = np.prod((z - zl[other]) / (zl[pos] - zl[other]))
            # Bridge factors for the two globally-dropped Dirichlet endpoints.
            if is_first:
                lag *= (z - z_min) / (pts[i] - z_min)
            if is_last:
                lag *= (z - z_max) / (pts[i] - z_max)
            rows.append(k)
            cols.append(int(i))
            vals.append(lag / wsqrt[i])

    return sp.csr_matrix(
        (np.array(vals, dtype=np.complex128), (rows, cols)),
        shape=(xs.shape[0], grid.n),
    )
