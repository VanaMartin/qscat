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

    def weights(self) -> tuple[npt.NDArray[np.complex128], ...]:
        """The D bridge-summed, COMPLEX (ECS-scaled) quadrature weight arrays,
        shaped to broadcast exactly like `points()` (same `_broadcast_shape`
        convention: axis `d`'s weight varies along axis `d`, size 1 elsewhere).

        Converting a *function* `f` to FEM-DVR basis coefficients (not
        evaluating an already-DVR-diagonal potential, which is what
        `potential_nd` does) requires `c_j = f(x_j) * sqrt(w_j)`, using the
        GLOBAL, bridge-summed, COMPLEX weight `FemDvrEcsGrid.weights` --
        `qscat.dvr.kinetic`'s module docstring calls getting this wrong (using
        a real or per-element weight instead) "the classic assembly trap".
        `TensorGrid` already spares a caller from remembering `real_mask()`;
        the same principle applies here, since a driving term or an
        asymptotic channel function must go through this exact factor at
        every coordinate to land correctly on the pre-normalized DVR basis.
        See also `sqrt_weights()`, which supplies that factor directly.
        """
        return tuple(
            np.asarray(g.weights, dtype=np.complex128).reshape(self._broadcast_shape(d))
            for d, g in enumerate(self._grids)
        )

    def sqrt_weights(self) -> tuple[npt.NDArray[np.complex128], ...]:
        """`sqrt(weights())` per axis -- the exact factor a basis-coefficient
        conversion `c_j = f(x_j) * sqrt(w_j)` needs, so a caller never has to
        write `np.sqrt` on the bridge-summed complex weight itself (a detail
        easy to get subtly wrong under ECS, since the weight is complex and
        `sqrt` of a complex number requires choosing a branch -- NumPy's
        principal branch, `Re(sqrt(z)) >= 0`, is the one used here and is the
        correct one for this weight, which never crosses the negative real
        axis for a valid ECS tail angle).
        """
        return tuple(np.sqrt(w) for w in self.weights())

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
    """`sum_d I x … x T_d x … x I`, with `T_d` built at mass `masses[d]`."""
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

    TWO BROADCASTING TRAPS to know about, because `np.broadcast_to` will
    happily paper over both rather than raising:

    - A RANK-DEFICIENT return value. If `V` returns an array with fewer axes
      than `tgrid.ndim` (e.g. a plain `(n,)` array at D=2, perhaps because `V`
      evaluated only one coordinate, or returned a precomputed 1-D array by
      mistake), NumPy broadcasting aligns it against the LAST axis and tiles
      it silently along every other axis -- a wrong potential with a
      plausible shape, not a crash. This function guards against exactly
      that: any non-scalar `V` result whose `ndim` is not `tgrid.ndim` raises
      `ValueError`. A genuine scalar (`ndim == 0`, a spatially constant `V`)
      is fine and is deliberately exempt, since a true scalar broadcasts
      unambiguously to every axis.
    - AXIS TRANSPOSITION when two axes have the SAME size. This one is NOT
      caught by any shape check, because the shapes are literally equal:
      `np.meshgrid` defaults to `indexing="xy"` (the first two returned
      arrays are transposed relative to the `"ij"`/broadcastable convention
      `TensorGrid.points()` uses), so a `V` that builds its own meshgrid
      internally -- or receives an externally precomputed surface of shape
      `(n1, n0)` instead of `(n0, n1)` -- silently transposes the potential
      whenever `n0 == n1`. There is no shape-based guard for this: call
      `V(*tgrid.points())` (or pass `indexing="ij"` explicitly to any
      internal `np.meshgrid` call) rather than reconstructing the grid by
      hand.
    """
    vals = np.asarray(V(*tgrid.points()), dtype=np.complex128)
    if vals.ndim not in (0, tgrid.ndim):
        raise ValueError(
            f"V(*tgrid.points()) returned an array of shape {vals.shape} "
            f"(ndim={vals.ndim}), but the tensor grid has ndim={tgrid.ndim}. "
            "A rank-deficient result broadcasts silently along the trailing "
            "axes rather than raising (numpy's default broadcasting rule) -- "
            "this is almost never the intended potential. A true scalar "
            "(ndim=0, spatially constant V) is exempt; anything else must "
            "have exactly tgrid.ndim axes."
        )
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
