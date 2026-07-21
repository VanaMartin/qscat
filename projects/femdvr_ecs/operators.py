"""Hamiltonian assembly + eigensolver helpers for the FEM-DVR-ECS grid.

See .superpowers/sdd/femdvr-ecs-extraction.md section 3 (Eigensolver) and
task-2-brief.md. H = T + diag(V) is complex-symmetric but non-Hermitian in
general (ECS elements), so we use the general eigensolver (np.linalg.eig,
c.f. ZGEEV) rather than a Hermitian one, and sort by ascending Re(E).
"""

import numpy as np

from kinetic import kinetic


def hamiltonian(grid, V, mass: float) -> np.ndarray:
    """Assemble H = T + diag(V) on the FEM-DVR basis (diagonal-potential approx)."""
    T = kinetic(grid, mass)
    Vals = V(grid.points) if callable(V) else np.asarray(V)
    Vals = np.broadcast_to(Vals, (grid.n,)).astype(complex)
    return T + np.diag(Vals)


def eigen(H: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Eigendecompose H, sorted by ascending Re(E)."""
    E, vecs = np.linalg.eig(H)           # complex, non-Hermitian
    order = np.argsort(E.real)
    return E[order], vecs[:, order]
