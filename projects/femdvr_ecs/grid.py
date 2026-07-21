"""FEM-DVR-ECS grid geometry: complex element boundaries, node placement,
bridge-summed weights, and the Dirichlet endpoint drop.

Ported from eMoScat's FemDvrEcsGrid.cpp (`initialize`, the real-`aa` overload)
and DvrGrid.cpp; see .superpowers/sdd/femdvr-ecs-extraction.md sections 1-3.

Construction (mirrors FemDvrEcsGrid.cpp:84-112 exactly):
  - Per element i (0-indexed), with reference GLL nodes/weights (xi, w) on
    (-1, 1) shared by all elements (same quadrature order):
      hz[i]      = 0.5 * length[i] * exp(i * angle[i])      (complex half-length)
      hz_real[i] = 0.5 * length[i]                           (real half-length)
      all_xz[ii+j] = hz[i]*xi[j] + hz[i] + az[i]             (ECS-scaled point)
      all_wz[ii+j] += hz[i]*w[j]                             (bridge-summed weight)
      all_xr[ii+j] = hz_real[i]*xi[j] + hz_real[i] + ar[i]   (unscaled bookkeeping)
    where ii = i*(nq-1) is the shared index of element i's first local node
    with element i-1's last local node (this is exactly how the bridge sum
    happens: all_wz is accumulated with `+=` into the same global slot).
  - az / ar are the complex/real cumulative element boundaries (cumulative
    sums of the per-element complex/real lengths, starting at x_min).
  - The two outermost global points (x_min and x_max) are dropped (Dirichlet
    psi=0): nb = tnel*(nq-1) + 1 - 2.
"""

import gll
import numpy as np
from spec import GridSpec


class FemDvrEcsGrid:
    """FEM-DVR-ECS grid geometry built from a validated `GridSpec`."""

    def __init__(self, spec: GridSpec) -> None:
        nq = spec.quadrature
        elements = spec.elements
        tnel = len(elements)

        xi, w = gll.gll_nodes_weights(nq)   # reference nodes/weights on (-1, 1)
        dLp = gll.diff_matrix(xi)

        # real (unscaled) cumulative element boundaries: ar[0] = x_min
        ar = np.empty(tnel + 1, dtype=float)
        ar[0] = spec.x_min
        for i, el in enumerate(elements):
            ar[i + 1] = ar[i] + el.length

        # complex (ECS-scaled) cumulative element boundaries: az[0] = x_min
        az = np.empty(tnel + 1, dtype=complex)
        az[0] = spec.x_min
        hz = np.empty(tnel, dtype=complex)
        for i, el in enumerate(elements):
            theta = np.deg2rad(el.angle_deg)
            Lk = el.length * np.exp(1j * theta)
            hz[i] = 0.5 * Lk
            az[i + 1] = az[i] + Lk

        nall = tnel * (nq - 1) + 1   # all grid points, including x_min and x_max
        nb = nall - 2                # basis functions after Dirichlet drop

        all_xz = np.zeros(nall, dtype=complex)
        all_wz = np.zeros(nall, dtype=complex)
        all_xr = np.zeros(nall, dtype=float)
        element_span_all: list[tuple[int, int]] = []   # inclusive (start, end) in "all" index space

        for i, el in enumerate(elements):
            ii = i * (nq - 1)
            hzi = hz[i]
            hz_real = 0.5 * el.length
            for j in range(nq):
                all_xz[ii + j] = hzi * xi[j] + hzi + az[i]
                all_wz[ii + j] += hzi * w[j]     # bridge sum: shared index accumulates
                all_xr[ii + j] = hz_real * xi[j] + hz_real + ar[i]
            element_span_all.append((ii, ii + nq - 1))

        # Dirichlet drop: remove index 0 (x_min) and index nall-1 (x_max).
        # Retained global index = "all" index - 1.
        points = all_xz[1: nall - 1]
        weights = all_wz[1: nall - 1]
        real_points = all_xr[1: nall - 1]

        element_ranges: list[tuple[int, int]] = []
        for s, e in element_span_all:
            gs = max(0, s - 1)
            ge = min(nb - 1, e - 1)
            element_ranges.append((gs, ge + 1))   # half-open [start, stop)

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
