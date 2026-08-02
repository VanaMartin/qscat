"""Project a 2-D FEM-DVR-ECS state onto an equidistant grid (cached, reusable).

The projection coefficients (the sparse `dvr_interpolation_matrix` per axis) are
computed ONCE at construction and applied per state -- the design that makes
frame-by-frame wavefunction animation cheap: build the projector, then call
`project(psi(t))` for every time step.

Ported from eMoScat's ``EquidistantProjector2d``. For a 2-D tensor state
``M[i, j]`` (electronic node i x nuclear node j) and per-axis interpolation
operators ``P0`` (rows = axis-0 samples) and ``P1`` (rows = axis-1 samples), the
projected field is the separable contraction ``P0 @ M @ P1.T`` -- equivalent to
the full Kronecker operator ``P0 (x) P1`` but without ever forming it.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from qscat.dvr import TensorGrid, dvr_interpolation_matrix
from qscat.exceptions import GridError

__all__ = ["EquidistantProjector"]


class EquidistantProjector:
    """Cached projector from a 2-D `TensorGrid` state to a uniform sampling grid.

    Parameters
    ----------
    tgrid : TensorGrid
        A 2-D tensor grid (``len(tgrid.grids) == 2``); axis 0 is the first grid
        (e.g. electronic r), axis 1 the second (e.g. nuclear R).
    samples : int or tuple[int, int]
        Number of uniform samples per axis (a single int applies to both).
    extent : tuple[tuple[float, float], tuple[float, float]] or None
        ``((a0, b0), (a1, b1))`` real sampling ranges per axis; defaults to each
        grid's full real region ``[x_min, R0]`` when None.
    """

    def __init__(
        self,
        tgrid: TensorGrid,
        *,
        samples: int | tuple[int, int] = 400,
        extent: tuple[tuple[float, float], tuple[float, float]] | None = None,
    ) -> None:
        if len(tgrid.grids) != 2:
            raise GridError(
                f"EquidistantProjector needs a 2-D TensorGrid, got {len(tgrid.grids)}-D"
            )
        self.tgrid = tgrid
        n0, n1 = (samples, samples) if isinstance(samples, int) else samples
        g0g, g1g = tgrid.grids
        if extent is None:
            extent = (
                (float(g0g.spec.x_min), float(g0g.R0)),
                (float(g1g.spec.x_min), float(g1g.R0)),
            )
        (a0, b0), (a1, b1) = extent
        self.axis0 = np.linspace(a0, b0, n0)
        self.axis1 = np.linspace(a1, b1, n1)
        # The two 1-D sparse interpolation operators (built once).
        self._p0 = dvr_interpolation_matrix(tgrid.grids[0], self.axis0)
        self._p1 = dvr_interpolation_matrix(tgrid.grids[1], self.axis1)

    def project(self, state: npt.NDArray[np.complex128]) -> npt.NDArray[np.complex128]:
        """Project a flat 2-D state to the uniform grid, shape ``(n0, n1)`` complex."""
        expected = self.tgrid.grids[0].n * self.tgrid.grids[1].n
        state = np.asarray(state, dtype=np.complex128).reshape(-1)
        if state.shape[0] != expected:
            raise GridError(
                f"state has {state.shape[0]} entries, expected {expected} for this grid"
            )
        m = state.reshape(self.tgrid.grids[0].n, self.tgrid.grids[1].n)
        # P0 @ M @ P1.T  ->  (n0, n1); sparse @ dense @ sparse.T
        return np.asarray(self._p0 @ m @ self._p1.T)

    def project_values(
        self, field: npt.NDArray[np.complex128]
    ) -> npt.NDArray[np.complex128]:
        """Project a nodal-VALUE field (e.g. a potential) on the same grid.

        `project` interpolates a √w-scaled DVR state (a wavefunction). A field
        given as plain nodal values ``f(x_i)`` (e.g. a potential surface) has
        interpolant ``sum_i f(x_i) L_i(x)`` -- no 1/√w -- so it must be scaled by
        ``√w`` before going through the same operator. Returns shape ``(n0, n1)``.
        """
        g0, g1 = self.tgrid.grids
        f = np.asarray(field, dtype=np.complex128).reshape(g0.n, g1.n)
        sw = np.sqrt(np.outer(g0.weights, g1.weights))  # √(w_i w_j) node scaling
        return self.project((sw * f).reshape(-1))
