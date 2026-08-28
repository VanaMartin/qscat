"""Channel-outermost assembly of a coupled-channel Hamiltonian.

The state vector is [psi_{l=1}(x), psi_{l=2}(x), ...] -- channels OUTERMOST.
That layout is what makes every off-diagonal block a plain diagonal matrix
(the coupling is a local potential), so the assembly is one `sp.bmat` and the
nonzero count grows as

    nnz = n_ch * nnz(H_block) + (n_ch^2 - n_ch) * n,

NOT as n_ch^2 times the whole matrix. It also makes the zero-anisotropy limit
exactly block-diagonal, with each block the corresponding single-channel
Hamiltonian -- which is how the embedding gate can be an identity rather than
a tolerance.

Nothing here knows about grids, dimensions, or physics: the same function
assembles the fixed-R electronic problem and the full 2-D one.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp

__all__ = ["assemble_coupled"]


def assemble_coupled(
    diagonal: Sequence[sp.spmatrix],
    coupling: Sequence[Sequence[npt.NDArray[np.complex128] | None]],
) -> sp.csr_matrix:
    """Assemble the channel-block matrix from per-channel blocks and couplings.

    `diagonal[i]` is channel `i`'s COMPLETE block (kinetic plus its own
    diagonal potential). `coupling[i][j]`, `i != j`, is the flattened
    off-diagonal potential; `None` means no coupling. `coupling[i][i]` must be
    `None`: a diagonal potential passed there would be added on top of the one
    already inside `diagonal[i]`.
    """
    n_ch = len(diagonal)
    if len(coupling) != n_ch or any(len(row) != n_ch for row in coupling):
        raise ValueError(f"coupling must be {n_ch}x{n_ch}, got {[len(r) for r in coupling]}")
    n = diagonal[0].shape[0]
    for i, blk in enumerate(diagonal):
        if blk.shape != (n, n):
            raise ValueError(f"diagonal[{i}] has size {blk.shape}, expected {(n, n)}")

    rows: list[list[sp.spmatrix | None]] = []
    for i in range(n_ch):
        row: list[sp.spmatrix | None] = []
        for j in range(n_ch):
            entry = coupling[i][j]
            if i == j:
                if entry is not None:
                    raise ValueError(
                        f"coupling[{i}][{i}] must be None -- channel {i}'s diagonal "
                        "potential belongs in diagonal[i], and passing it here would "
                        "add it twice"
                    )
                row.append(sp.csr_matrix(diagonal[i]))
            elif entry is None:
                row.append(None)
            else:
                vals = np.asarray(entry, dtype=np.complex128).ravel()
                if vals.size != n:
                    raise ValueError(f"coupling[{i}][{j}] has size {vals.size}, expected {n}")
                row.append(sp.diags(vals, format="csr"))
        rows.append(row)
    return sp.csr_matrix(sp.bmat(rows, format="csr"))
