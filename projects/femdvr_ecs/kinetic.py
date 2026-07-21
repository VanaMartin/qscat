"""FEM-DVR kinetic-energy operator assembly.

Ported from eMoScat's KineticEnergy.cpp:15-87; see
.superpowers/sdd/femdvr-ecs-extraction.md section 2 and task-2-brief.md.

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
  sub-block (per `grid.element_maps[k]`) is scatter-added (`+=`) into the
  global matrix. Adjacent elements share exactly one bridge global index, so
  the `+=` accumulates the bridge-corner coupling automatically -- no special
  casing needed.
"""

import numpy as np

from projects.femdvr_ecs import gll


def kinetic(grid, mass: float) -> np.ndarray:
    """Assemble the (n, n) complex FEM-DVR kinetic-energy matrix."""
    n = grid.n
    nq = grid.nq
    T = np.zeros((n, n), dtype=complex)

    _, wl = gll.gll_nodes_weights(nq)   # reference GLL weights on (-1, 1)

    for k, (local, global_idx) in enumerate(grid.element_maps):
        hz = grid.hz[k]

        wze = hz * wl                              # (nq,)
        dBF = grid.dLp.T / hz                       # dBF[a, l] = dLp[l, a] / hz

        # Normalize each basis row by the GLOBAL bridge-summed weight at that
        # basis function's global index. Only the retained local indices have
        # a defined global index; the dropped Dirichlet endpoint(s) never
        # participate in T_local's retained sub-block, so leave them alone.
        norm = np.ones(nq, dtype=complex)
        norm[local] = 1.0 / np.sqrt(grid.weights[global_idx])
        dBF_n = dBF * norm[:, np.newaxis]           # broadcast over node axis

        T_local = (1.0 / (2.0 * mass)) * np.einsum("l,al,bl->ab", wze, dBF_n, dBF_n)

        T[np.ix_(global_idx, global_idx)] += T_local[np.ix_(local, local)]

    return T
