"""FEM-DVR-ECS: finite-element DVR radial grid with exterior complex scaling.

Public API:
  - `ElementSpec`, `GridSpec` -- grid specification dataclasses.
  - `FemDvrEcsGrid` -- grid geometry (points, bridge-summed weights, the
    local-to-global `element_maps`) built from a validated `GridSpec`.
  - `kinetic`, `kinetic_sparse` -- assemble the FEM-DVR kinetic-energy matrix
    on a grid, dense or sparse (CSR). The dense one is the sparse one's
    differential oracle.
  - `dvr_first_derivative_at_node` -- the FEM-DVR first-derivative operator's
    row at a real grid node (`d @ psi_coeffs ~= d/dx psi(x_node)`), built
    from the same element-local Lagrange-derivative matrix `kinetic` uses;
    the flow (flux) time-dependent energy extractor's numerical primitive.
  - `hamiltonian`, `eigen` -- diagonal-potential Hamiltonian assembly and the
    complex-symmetric (non-Hermitian) eigensolver.
  - `gll_nodes_weights`, `diff_matrix` -- the underlying Gauss-Lobatto-
    Legendre quadrature/differentiation building blocks (reusable outside
    FEM-DVR-ECS).
  - `TensorGrid` -- tensor product of D FEM-DVR-ECS grids (C order, last axis
    fastest), with the ECS real-region mask and separable-state construction.
  - `kinetic_nd`, `potential_nd`, `hamiltonian_nd` -- the N-dimensional
    Kronecker-sum Hamiltonian assembled on a `TensorGrid`, sparse (CSR).

See `docs/physics/femdvr-ecs.md` for the method and its validation benchmarks.
The method is Rescigno & McCurdy, Phys. Rev. A 62, 032706 (2000); this
implementation was ported from eMoScat.
"""

from __future__ import annotations

from .derivative import dvr_first_derivative_at_node
from .gll import diff_matrix, gll_nodes_weights
from .grid import FemDvrEcsGrid
from .interpolate import dvr_interpolation_matrix
from .kinetic import kinetic, kinetic_sparse
from .operators import eigen, hamiltonian
from .spec import ElementSpec, GridSpec
from .tensor import TensorGrid, hamiltonian_nd, kinetic_nd, potential_nd

__all__ = [
    "ElementSpec",
    "GridSpec",
    "FemDvrEcsGrid",
    "kinetic",
    "kinetic_sparse",
    "dvr_first_derivative_at_node",
    "dvr_interpolation_matrix",
    "hamiltonian",
    "eigen",
    "gll_nodes_weights",
    "diff_matrix",
    "TensorGrid",
    "kinetic_nd",
    "potential_nd",
    "hamiltonian_nd",
]
