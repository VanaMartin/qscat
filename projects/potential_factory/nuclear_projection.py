"""Sparse toy assembly in a real orthonormal neutral nuclear DVR basis."""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp


def projected_hamiltonian(radial, thresholds, basis, interaction) -> sp.csr_matrix:
    """Return radial-major H; thresholds include V0, radial includes centrifugal.

    Basis columns are real DVR coefficients, not point values. Interaction has
    shape (radial points, nuclear DVR points); no extra weights are applied.
    See docs/physics/finite-box-nuclear-projection.md for conventions and limits.
    """
    kinetic = sp.csr_matrix(radial, dtype=complex)
    eps, B, values = np.asarray(thresholds), np.asarray(basis), np.asarray(interaction)
    nr = kinetic.shape[0]
    if (
        kinetic.shape != (nr, nr)
        or nr < 1
        or eps.ndim != 1
        or len(eps) < 1
        or not np.isrealobj(eps)
        or not np.isrealobj(B)
        or B.ndim != 2
        or B.shape[1] != len(eps)
        or values.shape != (nr, B.shape[0])
        or values.size > 2_000_000
        or not all(np.isfinite(x).all() for x in (kinetic.data, eps, B, values))
    ):
        raise ValueError(
            "Finite radial matrix, real thresholds/basis and matching interaction required"
        )
    if not np.allclose(B.T @ B, np.eye(len(eps)), rtol=0, atol=1e-10):
        raise ValueError("Nuclear DVR coefficient columns must be orthonormal")
    blocks = [sp.csr_matrix(B.T @ (row[:, None] * B)) for row in values]
    return (
        sp.kron(kinetic, sp.eye(len(eps)), format="csr")
        + sp.kron(sp.eye(nr), sp.diags(eps), format="csr")
        + sp.block_diag(blocks, format="csr")
    ).tocsr()
