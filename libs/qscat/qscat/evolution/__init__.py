"""Time propagation (Crank-Nicolson, ...).

Public API:
  - `make_cn_stepper` -- Crank-Nicolson propagator for the time-dependent
    Schrodinger equation, `d/dt psi = -i H psi`, for a general (possibly
    non-Hermitian) complex Hamiltonian matrix `H`. See
    `docs/physics/n2-td-cross-section.md` for the N2 resonance application.
  - `make_sparse_cn_stepper` -- the sparse sibling for large sparse `H`,
    factoring once with `SparseLU`.
  - `make_pade_stepper` -- order-N diagonal-Pade generalization of the sparse
    CN stepper (order 1 == Crank-Nicolson); `O(dt^(2N+1))` per step, for the
    higher-order accuracy the TD cross section needs to converge to the TI
    oracle. `pade_roots` exposes the denominator roots.
"""

from __future__ import annotations

from .crank_nicolson import make_cn_stepper, make_sparse_cn_stepper
from .pade import make_pade_stepper, pade_roots

__all__ = ["make_cn_stepper", "make_pade_stepper", "make_sparse_cn_stepper", "pade_roots"]
