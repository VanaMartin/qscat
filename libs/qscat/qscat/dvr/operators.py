"""Hamiltonian assembly + eigensolver helpers for the FEM-DVR-ECS grid.

See `docs/physics/femdvr-ecs.md`. H = T + diag(V) is complex-symmetric but non-Hermitian in
general (ECS elements), so we use the general eigensolver (np.linalg.eig,
c.f. ZGEEV) rather than a Hermitian one, and sort by ascending Re(E).
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import numpy.typing as npt

from .grid import FemDvrEcsGrid
from .kinetic import kinetic

__all__ = ["hamiltonian", "eigen"]

PotentialLike = (
    Callable[[npt.NDArray[np.complex128]], npt.ArrayLike] | npt.ArrayLike
)


def hamiltonian(
    grid: FemDvrEcsGrid, V: PotentialLike, mass: float
) -> npt.NDArray[np.complex128]:
    """Assemble H = T + diag(V) on the FEM-DVR basis (diagonal-potential approx)."""
    T = kinetic(grid, mass)
    raw_vals = V(grid.points) if callable(V) else V
    Vals = np.broadcast_to(np.asarray(raw_vals), (grid.n,)).astype(np.complex128)
    result: npt.NDArray[np.complex128] = T + np.diag(Vals)
    return result


def eigen(
    H: npt.NDArray[np.complex128],
) -> tuple[npt.NDArray[np.complex128], npt.NDArray[np.complex128]]:
    """Eigendecompose H, sorted by ascending Re(E).

    Eigenvectors carry numpy's `v†v=1` (Hermitian) normalization; for ECS
    c-product observables (residues/widths in sub-project #2), re-normalize
    to `vᵀv=1`.
    """
    E, vecs = np.linalg.eig(H)  # complex, non-Hermitian
    order = np.argsort(E.real)
    return E[order], vecs[:, order]
