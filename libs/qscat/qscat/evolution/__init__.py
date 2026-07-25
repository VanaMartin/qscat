"""Time propagation (Crank-Nicolson, ...).

Public API:
  - `make_cn_stepper` -- Crank-Nicolson propagator for the time-dependent
    Schrodinger equation, `d/dt psi = -i H psi`, for a general (possibly
    non-Hermitian) complex Hamiltonian matrix `H`. See
    `docs/physics/n2-td-cross-section.md` for the N2 resonance application.
  - `make_sparse_cn_stepper` -- the sparse sibling for large sparse `H`,
    factoring once with `SparseLU`.
"""

from __future__ import annotations

from .crank_nicolson import make_cn_stepper, make_sparse_cn_stepper

__all__ = ["make_cn_stepper", "make_sparse_cn_stepper"]
