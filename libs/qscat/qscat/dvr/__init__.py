"""FEM-DVR-ECS: finite-element DVR radial grid with exterior complex scaling.

Public API:
  - `ElementSpec`, `GridSpec` -- grid specification dataclasses.
  - `FemDvrEcsGrid` -- grid geometry (points, bridge-summed weights, the
    local-to-global `element_maps`) built from a validated `GridSpec`.
  - `kinetic`, `kinetic_sparse` -- assemble the FEM-DVR kinetic-energy matrix
    on a grid, dense or sparse (CSR). The dense one is the sparse one's
    differential oracle.
  - `hamiltonian`, `eigen` -- diagonal-potential Hamiltonian assembly and the
    complex-symmetric (non-Hermitian) eigensolver.
  - `gll_nodes_weights`, `diff_matrix` -- the underlying Gauss-Lobatto-
    Legendre quadrature/differentiation building blocks (reusable outside
    FEM-DVR-ECS).

See `docs/physics/femdvr-ecs.md` for the method and its validation
benchmarks, and `.superpowers/sdd/femdvr-ecs-extraction.md` for the port-scout
extraction from eMoScat this implementation is based on.
"""

from __future__ import annotations

from .gll import diff_matrix, gll_nodes_weights
from .grid import FemDvrEcsGrid
from .kinetic import kinetic, kinetic_sparse
from .operators import eigen, hamiltonian
from .spec import ElementSpec, GridSpec

__all__ = [
    "ElementSpec",
    "GridSpec",
    "FemDvrEcsGrid",
    "kinetic",
    "kinetic_sparse",
    "hamiltonian",
    "eigen",
    "gll_nodes_weights",
    "diff_matrix",
]
