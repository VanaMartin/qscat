"""Linear-algebra helpers: dimension-general Kronecker sums, cached sparse
factorizations, and the exterior-complex-scaling c-product.

Pure linear algebra -- nothing here knows about grids, potentials or physics,
so it composes with any discretization.

Public API:
  - `kron_sum` -- `sum_d I x ... x A_d x ... x I` for arbitrary D.
  - `c_product` -- the bilinear (non-conjugated) ECS inner product.
  - `SparseLU` -- cached sparse LU factorization (factor once, solve many),
    with fill-in and memory diagnostics and a SuperLU/MUMPS backend switch.
  - `ShiftInvertEigs` -- the k eigenpairs nearest a complex shift (sparse
    shift-invert Arnoldi on top of `SparseLU`, reusing its symbolic analysis
    across shifts).
  - `default_backend` / `set_default_backend` / `get_default_backend` --
    process-wide override that `SparseLU(backend="auto")` resolves against,
    for forcing an internal-`SparseLU` computation onto one factorization
    engine (a backend-equivalence check).

See `docs/physics/nd-tensor-hamiltonian.md`.
"""

from __future__ import annotations

from .eigs import ShiftInvertEigs
from .inner import c_product
from .kron import kron_sum
from .sparse_lu import (
    SparseLU,
    default_backend,
    get_default_backend,
    set_default_backend,
)

__all__ = [
    "kron_sum",
    "c_product",
    "SparseLU",
    "ShiftInvertEigs",
    "default_backend",
    "set_default_backend",
    "get_default_backend",
]
