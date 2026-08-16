"""FEM-DVR kinetic-energy operator assembly.

Ported from eMoScat's KineticEnergy.cpp:15-87; see `docs/physics/femdvr-ecs.md`.
`kinetic_sparse` is the COO/CSR sibling, with the dense assembly as its
differential oracle.

T = -(1/2*mass) d^2/dz^2, assembled element-by-element on the global,
Dirichlet-trimmed FEM-DVR basis exposed by `FemDvrEcsGrid`:

  Per element k, with reference GLL weights `wl` (shared by all elements,
  same quadrature order) and complex half-length `hz = grid.hz[k]`:
    wze[l]      = hz * wl[l]                         (scaled quadrature weight)
    dBF[a, l]   = grid.dLp[l, a] / hz                 (scaled derivative,
                  basis index a, node index l; dLp[j,i] = L_i'(x_j))
    dBF[a, :]  /= sqrt(grid.weights[global_idx(a)])   (normalize by the GLOBAL
                  bridge-summed weight at basis a's global index -- NOT the
                  local element weight; this is the classic assembly trap)
    T_local[a, b] = (1/(2*mass)) * sum_l wze[l] * dBF[a, l] * dBF[b, l]

  T_local is computed over ALL nq local basis indices, then the RETAINED
  sub-block (per `grid.element_maps[k]`) is combined into the global matrix.
  Adjacent elements share exactly one bridge global index, so accumulating
  (scatter-add for the dense path, COO duplicate-summation for the sparse
  path) reproduces the bridge-corner coupling automatically -- no special
  casing needed.

  The per-element `T_local` computation (the `wze`/`dBF`/`norm`/`einsum`
  block) is shared by both assembly paths via the private `_element_block`
  helper below. The bridge ACCUMULATION itself is deliberately left as two
  independent implementations (dense scatter-add vs. COO triplet emission +
  summation on CSR conversion): that accumulation is the bug-prone part of
  this method, so keeping it duplicated means the differential test between
  `kinetic` and `kinetic_sparse` still exercises two genuinely different code
  paths. The einsum math itself is already pinned by the analytic
  particle-in-a-box benchmarks, so sharing *that* part costs no coverage.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp

from .gll import gll_nodes_weights
from .grid import FemDvrEcsGrid

__all__ = ["kinetic", "kinetic_sparse"]


def _element_block(
    grid: FemDvrEcsGrid, k: int, mass: float, wl: npt.NDArray[np.float64]
) -> npt.NDArray[np.complex128]:
    """Retained (local, local) sub-block of element `k`'s local T_local.

    Shared by `kinetic` and `kinetic_sparse`; see module docstring for why the
    surrounding bridge accumulation is NOT shared.
    """
    local, global_idx = grid.element_maps[k]
    hz = grid.hz[k]

    wze = hz * wl  # (nq,)
    dBF = grid.dLp.T / hz  # dBF[a, l] = dLp[l, a] / hz

    # Normalize each basis row by the GLOBAL bridge-summed weight at that
    # basis function's global index. Only the retained local indices have
    # a defined global index; the dropped Dirichlet endpoint(s) never
    # participate in T_local's retained sub-block, so leave them alone.
    norm = np.ones(grid.nq, dtype=complex)
    norm[local] = 1.0 / np.sqrt(grid.weights[global_idx])
    dBF_n = dBF * norm[:, np.newaxis]  # broadcast over node axis

    T_local = (1.0 / (2.0 * mass)) * np.einsum("l,al,bl->ab", wze, dBF_n, dBF_n)
    block: npt.NDArray[np.complex128] = T_local[np.ix_(local, local)]
    return block


def kinetic(grid: FemDvrEcsGrid, mass: float) -> npt.NDArray[np.complex128]:
    """Assemble the (n, n) complex FEM-DVR kinetic-energy matrix (dense)."""
    n = grid.n
    T: npt.NDArray[np.complex128] = np.zeros((n, n), dtype=complex)

    _, wl = gll_nodes_weights(grid.nq)  # reference GLL weights on (-1, 1)

    for k, (_local, global_idx) in enumerate(grid.element_maps):
        T[np.ix_(global_idx, global_idx)] += _element_block(grid, k, mass, wl)

    return T


def kinetic_sparse(grid: FemDvrEcsGrid, mass: float) -> sp.csr_matrix:
    """Sparse (CSR) FEM-DVR kinetic-energy matrix -- the sparse sibling of `kinetic`.

    Identical mathematics to the dense `kinetic()`, which is retained as this
    function's differential oracle. The only structural difference is that
    per-element blocks are emitted as COO triplets instead of scatter-added into
    a dense array: `coo_matrix` SUMS duplicate `(row, col)` entries on
    conversion, which reproduces the dense version's `+=` bridge accumulation
    exactly. No bridge special-casing is needed or wanted.

    Nonzero count is `nq**2 * tnel - 4*nq + 3 - tnel` (eMoScat
    `KineticEnergy.cpp:95`) -- the union of the per-element `nq x nq` blocks,
    overlapping by one index at each bridge, less the two dropped Dirichlet
    endpoints.
    """
    n = grid.n
    nq = grid.nq

    _, wl = gll_nodes_weights(nq)  # reference GLL weights on (-1, 1)

    rows: list[npt.NDArray[np.intp]] = []
    cols: list[npt.NDArray[np.intp]] = []
    vals: list[npt.NDArray[np.complex128]] = []

    for k, (_local, global_idx) in enumerate(grid.element_maps):
        block = _element_block(grid, k, mass, wl)

        gi, gj = np.meshgrid(global_idx, global_idx, indexing="ij")
        rows.append(gi.ravel())
        cols.append(gj.ravel())
        vals.append(block.ravel())

    coo = sp.coo_matrix(
        (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
        shape=(n, n),
        dtype=np.complex128,
    )
    return sp.csr_matrix(coo)  # duplicate (row, col) entries are summed here
