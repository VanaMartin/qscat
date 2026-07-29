"""FEM-DVR first-derivative-at-a-node operator.

The genuinely new numerical primitive `qscat.core.td_extractors.Flux` (the
flow time-dependent energy extractor) needs: a ROW `d` such that
`d @ psi_coeffs ~= d/dx psi(x_node)` at a real (unscaled) grid node, built
from the SAME element-local Gauss-Lobatto Lagrange-derivative matrix
`kinetic._element_block` uses to assemble the kinetic-energy operator
(`grid.dLp`, `grid.hz`, `grid.element_maps`) -- see that module's docstring
for the shared conventions this mirrors.

Ported from eMoScat's `GridVector::derivative` (`FemDvrEcs/GridVector.cpp:240`,
confirmed by port-scout reading for this task): that function stores a
coefficient-convention vector internally (`vector_[i] = value * sqrt(w(i))`,
`GridVector::f`) and computes

    d/dx f(x_i) = (2/length) * sum_k dlp(k, i) * vector_[k] / sqrt(w(k))

over the ONE element containing the target node (`length = 2*hz`, so
`2/length == 1/hz`). Rewritten as a row `d` acting on a COEFFICIENT vector
(the same `coeff = f(x)*sqrt(w)` convention `qscat.core.wavepacket.
gaussian_coeffs`/`qscat.core.td_extractors.Dirac` use):

    d[global_idx[a]] = dLp[node_local, a_local] / (hz * sqrt(weights[global_idx[a]]))

summed over the element's local nodes `a` (only one element's worth of
entries is nonzero -- an interior node touches exactly one element; a bridge
node's row still only draws from ONE of its two adjacent elements, see
below), so that `d @ psi_coeffs ~= d/dx psi(x_node)`.

For a BRIDGE node (the shared border of two elements), this picks the FIRST
element (ascending element index, i.e. the element to the LEFT of the
border) whose `element_maps` entry contains `node_index`. A genuinely
smooth, well-resolved function gives (to the GLL quadrature's near-spectral
accuracy) the SAME derivative from either side, so this is a convention, not
a source of error -- confirmed by the border-node cases in
`test_dvr_derivative.py`.

Only meaningful in the real (unscaled) region: the ECS tail's Jacobian
(`hz`) is complex, and this primitive -- like `Dirac`'s fixed-point
projection -- is only used where the analysis surface actually lives (past
the interaction, still inside `R0`).
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from .grid import FemDvrEcsGrid

__all__ = ["dvr_first_derivative_at_node"]


def dvr_first_derivative_at_node(
    grid: FemDvrEcsGrid, node_index: int
) -> npt.NDArray[np.complex128]:
    """Row `d` (length `grid.n`) s.t. `d @ psi_coeffs ~= d/dx psi(x_node)`.

    `node_index` must be a real (unscaled) grid node
    (`grid.real_points[node_index] <= grid.R0`) -- see module docstring.
    """
    n = grid.n
    if not (0 <= node_index < n):
        raise ValueError(f"node_index {node_index} out of range for grid of size {n}")
    if grid.real_points[node_index] > grid.R0:
        raise ValueError(
            f"node_index {node_index} (r={grid.real_points[node_index]}) is not in the real "
            f"(unscaled) region (R0={grid.R0}) -- dvr_first_derivative_at_node requires a real "
            "hz/sqrt(w) surface"
        )

    for k, (local, global_idx) in enumerate(grid.element_maps):
        matches = np.flatnonzero(global_idx == node_index)
        if matches.size == 0:
            continue
        node_local = int(local[matches[0]])
        hz = grid.hz[k]
        sqrt_w = np.sqrt(np.asarray(grid.weights[global_idx], dtype=np.complex128))
        d = np.zeros(n, dtype=np.complex128)
        d[global_idx] = grid.dLp[node_local, local] / (hz * sqrt_w)
        return d

    raise AssertionError(f"node_index {node_index} not found in any element (grid bug)")
