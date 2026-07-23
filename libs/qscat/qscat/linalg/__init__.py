"""Linear-algebra helpers: dimension-general Kronecker sums, cached sparse
factorizations, and the exterior-complex-scaling c-product.

Pure linear algebra -- nothing here knows about grids, potentials or physics,
so it composes with any discretization.

Public API:
  - `kron_sum` -- `sum_d I x ... x A_d x ... x I` for arbitrary D.
  - `c_product` -- the bilinear (non-conjugated) ECS inner product.
  - `SparseLU` -- cached sparse LU factorization (factor once, solve many),
    with fill-in and memory diagnostics.

See `docs/physics/nd-tensor-hamiltonian.md`.
"""

from __future__ import annotations

from .inner import c_product
from .kron import kron_sum
from .sparse_lu import SparseLU

__all__ = ["kron_sum", "c_product", "SparseLU"]
