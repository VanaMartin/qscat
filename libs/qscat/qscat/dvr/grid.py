"""FEM-DVR-ECS grid geometry: complex element boundaries, node placement,
bridge-summed weights, and the Dirichlet endpoint drop.

Ported from eMoScat's FemDvrEcsGrid.cpp (`initialize`, the real-`aa` overload)
and DvrGrid.cpp; see .superpowers/sdd/femdvr-ecs-extraction.md sections 1-3.

Construction (mirrors FemDvrEcsGrid.cpp:84-112 exactly):
  - Per element i (0-indexed), with reference GLL nodes/weights (xi, w) on
    (-1, 1) shared by all elements (same quadrature order):
      hz[i]      = 0.5 * length[i] * exp(i * angle[i])      (complex half-length)
      hz_real[i] = 0.5 * length[i]                           (real half-length)
      all_xr[ii+j] = hz_real[i]*xi[j] + hz_real[i] + ar[i]   (unscaled bookkeeping)
      all_xz[ii+j] = ecs_map(all_xr[ii+j], R0, angle[i])     (ECS-scaled point)
      all_wz[ii+j] += hz[i]*w[j]                             (bridge-summed weight)
    where ii = i*(nq-1) is the shared index of element i's first local node
    with element i-1's last local node (this is exactly how the bridge sum
    happens: all_wz is accumulated with `+=` into the same global slot).
  - ar is the real (unscaled) cumulative element boundary (cumulative sum of
    per-element lengths, starting at x_min). Point placement routes through
    `qscat.ecs.ecs_map` -- the single source of the ECS coordinate transform
    -- rather than re-deriving `R0 + (x - R0) e^{i theta}` locally; this is
    numerically identical to the original per-element complex-boundary
    (`az`) construction for the validated case of a single shared tail angle
    (see `qscat.ecs.ecs_map` and the `GridSpec` docstring caveat about
    multiple different tail angles being unverified).
  - The weight Jacobian `hz` (used only for `all_wz`, not point placement)
    still carries the per-element complex half-length directly: it is the
    local derivative dz/dxi of the ECS map's linear tail segment, which is
    independent of how the point coordinate itself is computed.
  - The two outermost global points (x_min and x_max) are dropped (Dirichlet
    psi=0): nb = tnel*(nq-1) + 1 - 2.

Local-to-global index mapping (`element_maps`):
  The retained local GLL node indices differ at the first/last element because
  of the Dirichlet drop, and are NOT simply `range(nq)` shifted by a constant
  offset -- a consumer must not assume `global = start + local`. Instead each
  element carries an EXPLICIT `(local_idx, global_idx)` pair of equal-length
  int arrays (see `element_maps` below); `global_idx[k]` is where local node
  `local_idx[k]` of that element scatters to in the length-`nb` global basis.
  Adjacent elements deliberately share their boundary global index (the last
  entry of element i's global_idx equals the first entry of element i+1's
  global_idx), which is what makes a `+=` scatter-accumulate assemble the
  bridge coupling correctly in the kinetic-operator build (Task 2).

  Worked example (nq=3, tnel=2, single-element-local nodes = [0, 1, 2]):
    element 0 (first): drop local 0            -> local_idx=[1, 2], global_idx=[0, 1]
    element 1 (last):  drop local nq-1=2        -> local_idx=[0, 1], global_idx=[1, 2]
    Shared global index 1 (the bridge node) appears in both elements' global_idx,
    exactly where a consumer should accumulate contributions from both.
    nb = tnel*(nq-1) + 1 - 2 = 2*2+1-2 = 3, matching the union {0, 1, 2}.
  With a single element (tnel=1), BOTH endpoints (local 0 and local nq-1) are
  dropped: local_idx=[1, ..., nq-2], global_idx=[0, ..., nq-3].
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from qscat.ecs import ecs_map

from . import gll
from .spec import GridSpec

__all__ = ["FemDvrEcsGrid"]


class FemDvrEcsGrid:
    """FEM-DVR-ECS grid geometry built from a validated `GridSpec`.

    See the module docstring for the `element_maps` local->global index
    convention that the kinetic-operator assembly relies on.
    """

    spec: GridSpec
    nq: int
    n: int
    R0: float
    points: npt.NDArray[np.complex128]
    weights: npt.NDArray[np.complex128]
    real_points: npt.NDArray[np.float64]
    dLp: npt.NDArray[np.float64]
    hz: npt.NDArray[np.complex128]
    element_ranges: list[tuple[int, int]]
    element_maps: list[tuple[npt.NDArray[np.intp], npt.NDArray[np.intp]]]

    def __init__(self, spec: GridSpec) -> None:
        nq = spec.quadrature
        elements = spec.elements
        tnel = len(elements)

        xi, w = gll.gll_nodes_weights(nq)  # reference nodes/weights on (-1, 1)
        dLp = gll.diff_matrix(xi)

        # real (unscaled) cumulative element boundaries: ar[0] = x_min
        ar = np.empty(tnel + 1, dtype=float)
        ar[0] = spec.x_min
        for i, el in enumerate(elements):
            ar[i + 1] = ar[i] + el.length

        # complex half-lengths (weight Jacobian only -- point placement goes
        # through ecs_map, see module docstring)
        hz = np.empty(tnel, dtype=complex)
        for i, el in enumerate(elements):
            theta = np.deg2rad(el.angle_deg)
            hz[i] = 0.5 * el.length * np.exp(1j * theta)

        nall = tnel * (nq - 1) + 1  # all grid points, including x_min and x_max
        nb = nall - 2  # basis functions after Dirichlet drop

        all_xz = np.zeros(nall, dtype=complex)
        all_wz = np.zeros(nall, dtype=complex)
        all_xr = np.zeros(nall, dtype=float)
        element_span_all: list[tuple[int, int]] = []  # inclusive (start, end) in "all" index space

        for i, el in enumerate(elements):
            ii = i * (nq - 1)
            hzi = hz[i]
            hz_real = 0.5 * el.length
            for j in range(nq):
                x_real = hz_real * xi[j] + hz_real + ar[i]
                all_xr[ii + j] = x_real
                all_xz[ii + j] = ecs_map(x_real, spec.R0, el.angle_deg)
                all_wz[ii + j] += hzi * w[j]  # bridge sum: shared index accumulates
            element_span_all.append((ii, ii + nq - 1))

        # Dirichlet drop: remove index 0 (x_min) and index nall-1 (x_max).
        # Retained global index = "all" index - 1.
        points = all_xz[1 : nall - 1]
        weights = all_wz[1 : nall - 1]
        real_points = all_xr[1 : nall - 1]

        # element_ranges: half-open [start, stop) slices into the global basis
        # spanned by each element, INCLUSIVE of shared boundary indices with
        # neighbors (i.e. adjacent ranges overlap by one index). This is a
        # coarse, position-independent summary kept for backward-compat /
        # quick range checks; it does NOT say which local index maps to which
        # global index, so a kinetic-assembly consumer should use
        # `element_maps` instead (see module docstring).
        element_ranges: list[tuple[int, int]] = []
        for s, e in element_span_all:
            gs = max(0, s - 1)
            ge = min(nb - 1, e - 1)
            element_ranges.append((gs, ge + 1))  # half-open [start, stop)

        # element_maps: explicit per-element (local_idx, global_idx) pairs.
        # local_idx = retained local GLL node indices (0..nq-1), dropping the
        # Dirichlet endpoint(s): local 0 for the first element, local nq-1 for
        # the last element (both, if tnel == 1). global_idx = the matching
        # row/column index into the length-nb global basis ("all" index - 1).
        # See module docstring for the convention and a worked example.
        element_maps: list[tuple[npt.NDArray[np.intp], npt.NDArray[np.intp]]] = []
        for i, (s, _e) in enumerate(element_span_all):
            local = np.arange(nq)
            if i == 0:
                local = local[local != 0]
            if i == tnel - 1:
                local = local[local != nq - 1]
            global_idx = (s + local) - 1
            element_maps.append((local, global_idx))

        self.spec = spec
        self.nq = nq
        self.n = nb
        self.R0 = spec.R0
        self.points = points
        self.weights = weights
        self.real_points = real_points
        self.dLp = dLp
        self.hz = hz
        self.element_ranges = element_ranges
        self.element_maps = element_maps
