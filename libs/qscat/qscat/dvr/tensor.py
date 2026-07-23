"""Tensor products of FEM-DVR-ECS grids, and the N-dimensional Hamiltonian
assembled on them.

A separable-kinetic, diagonal-potential Hamiltonian on a tensor product of DVR
grids is a Kronecker sum plus a diagonal:

    H = sum_d  I x ... x T_d x ... x I  +  diag(V(x_0, ..., x_{D-1}))

Nothing about that is specific to two dimensions, or to molecular scattering.
This module is the dimension-general form; `qscat.linalg.kron_sum` does the
Kronecker algebra and knows nothing about grids.

Index convention: numpy-native C order, LAST axis fastest, so a state of shape
`tgrid.shape` ravels to the vector the Hamiltonian acts on. eMoScat uses the
opposite convention (first coordinate fastest); the two are physically
identical and differ only in basis ordering.

See `docs/physics/nd-tensor-hamiltonian.md`.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp

from qscat.linalg import kron_sum

from .grid import FemDvrEcsGrid
from .kinetic import kinetic_sparse

__all__ = ["TensorGrid", "kinetic_nd", "potential_nd", "hamiltonian_nd"]


class TensorGrid:
    """Tensor product of D FEM-DVR-ECS grids (C order: last axis fastest)."""

    def __init__(self, grids: Sequence[FemDvrEcsGrid]) -> None:
        tup = tuple(grids)
        if not tup:
            raise ValueError("TensorGrid requires at least one grid")
        self._grids = tup

    @property
    def grids(self) -> tuple[FemDvrEcsGrid, ...]:
        return self._grids

    @property
    def ndim(self) -> int:
        return len(self._grids)

    @property
    def shape(self) -> tuple[int, ...]:
        return tuple(g.n for g in self._grids)

    @property
    def size(self) -> int:
        return int(np.prod(self.shape))

    def _broadcast_shape(self, d: int) -> tuple[int, ...]:
        """Shape that puts axis `d`'s data on axis `d` and 1 everywhere else."""
        return (1,) * d + (-1,) + (1,) * (self.ndim - d - 1)

    def points(self) -> tuple[npt.NDArray[np.complex128], ...]:
        """The D complex (ECS-scaled) coordinate arrays, shaped to broadcast.

        For D = 2 these are `(n_0, 1)` and `(1, n_1)`, so a potential can be
        written naturally as `V(r, R)` without materializing a full meshgrid.
        For D = 1 the single array is plain 1-D of shape `(n_0,)`.
        """
        return tuple(
            np.asarray(g.points, dtype=np.complex128).reshape(self._broadcast_shape(d))
            for d, g in enumerate(self._grids)
        )

    def real_mask(self) -> npt.NDArray[np.bool_]:
        """Flat boolean mask, True where EVERY coordinate is in the unscaled region.

        Under exterior complex scaling a driving term or channel projection is
        only meaningful on the unscaled region, so anything of that kind must be
        masked with this before use. Making it a property of the grid is
        deliberate: the physics layer should not have to remember.
        """
        mask: npt.NDArray[np.bool_] | None = None
        for d, g in enumerate(self._grids):
            md = np.asarray(g.real_points <= g.R0, dtype=bool).reshape(
                self._broadcast_shape(d)
            )
            mask = md if mask is None else (mask & md)
        assert mask is not None  # ndim >= 1 guaranteed by __init__
        return np.asarray(np.broadcast_to(mask, self.shape).ravel(), dtype=bool)

    def outer(self, vectors: Sequence[npt.ArrayLike]) -> npt.NDArray[np.complex128]:
        """Separable state `⊗_d vectors[d]`, flattened to length `size`."""
        vecs = list(vectors)
        if len(vecs) != self.ndim:
            raise ValueError(f"expected {self.ndim} vectors, got {len(vecs)}")
        for d, v in enumerate(vecs):
            got = np.asarray(v).shape
            if got != (self.shape[d],):
                raise ValueError(f"vector {d} has shape {got}, expected {(self.shape[d],)}")
        out = np.asarray(vecs[0], dtype=np.complex128)
        for v in vecs[1:]:
            out = np.multiply.outer(out, np.asarray(v, dtype=np.complex128))
        return np.asarray(out.ravel(), dtype=np.complex128)


def kinetic_nd(tgrid: TensorGrid, masses: Sequence[float]) -> sp.csr_matrix:
    """`sum_d I x ... x T_d x ... x I`, with `T_d` built at mass `masses[d]`."""
    ms = list(masses)
    if len(ms) != tgrid.ndim:
        raise ValueError(f"expected {tgrid.ndim} masses, got {len(ms)}")
    return kron_sum([kinetic_sparse(g, m) for g, m in zip(tgrid.grids, ms, strict=True)])


def potential_nd(
    tgrid: TensorGrid, V: Callable[..., npt.ArrayLike]
) -> npt.NDArray[np.complex128]:
    """Evaluate `V` at the D-dimensional COMPLEX points, flattened.

    `V` is called as `V(x_0, ..., x_{D-1})` with the broadcastable arrays from
    `TensorGrid.points()`. It MUST NOT coerce its arguments to a real dtype:
    the points are complex on the ECS tail, and discarding the imaginary part
    silently destroys the analytic continuation the method depends on.
    """
    vals = np.asarray(V(*tgrid.points()), dtype=np.complex128)
    return np.asarray(np.broadcast_to(vals, tgrid.shape).ravel(), dtype=np.complex128)


def hamiltonian_nd(
    tgrid: TensorGrid, masses: Sequence[float], V: Callable[..., npt.ArrayLike]
) -> sp.csr_matrix:
    """`H = kinetic_nd(tgrid, masses) + diag(potential_nd(tgrid, V))` as CSR.

    Complex symmetric (`H = H^T`), NOT Hermitian, whenever any grid has an ECS
    tail. Use general algorithms only.
    """
    T = kinetic_nd(tgrid, masses)
    V_diag = potential_nd(tgrid, V)
    return sp.csr_matrix(T + sp.diags(V_diag, format="csr"))
