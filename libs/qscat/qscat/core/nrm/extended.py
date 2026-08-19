"""The extended-space (time-local) form of the nonlocal resonance model.

Gertitschke & Domcke, Phys. Rev. A 47, 1031 (1993), Eq. (2.1) is a memory
integral. It does not have to be evaluated as one. The kernel this package
builds (`nonlocal_potential.nonlocal_operator`, PRA 77 Eq. 60-61) is a sum of
resolvents,

    F(E) = sum_n diag(V_dn) (E - H_n)^-1 diag(V_dn),
    H_n  = T_R + V_0(R) + E_n(R)

so introducing one auxiliary nuclear packet `phi_n` per projected electronic
state turns Eq. (2.1) into a time-LOCAL coupled system,

    i d/dt Psi_d = (T_N + V_d) Psi_d + sum_n V_dn phi_n
    i d/dt phi_n = H_n phi_n + V_dn Psi_d

i.e. propagation under one arrow-shaped block Hamiltonian. Eliminating the
`phi_n` from the time-independent version returns PRA 77 Eq. (52) with exactly
the `F(E)` above -- this is a RESUMMATION of Eq. (2.1), not an approximation
of it, and `test_nrm_extended.py` gates that identity.

SYMMETRY: under ECS the matrix is complex SYMMETRIC, not Hermitian. Every
transpose here is `.T`; a `.conj().T` would silently produce a different
operator.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, cast

import numpy as np
import scipy.sparse as sp

from qscat.dvr import FemDvrEcsGrid, kinetic_sparse

from .ingredients import NrmIngredients
from .nonlocal_potential import (
    check_nodes_coincide,
    check_tail_coupling,
    continue_to_tail,
)

if TYPE_CHECKING:
    from qscat.model import ResonanceModel

__all__ = ["extended_hamiltonian"]


def extended_hamiltonian(
    ing: NrmIngredients,
    nuclear_grid: FemDvrEcsGrid,
    model: ResonanceModel,
    *,
    n_states: int | None = None,
) -> sp.csr_matrix:
    """The arrow block Hamiltonian of the module docstring.

    Block 0 (rows/columns `0:N_R`) is the discrete state; block `k+1` is
    projected electronic state `k`, in `ing`'s own adiabatic-tracking order.

    Raises
    ------
    ValueError
        If `n_states` is outside `[0, ing.E_n.shape[1]]`, if `ing.R` does not
        coincide with the grid's real DVR nodes, or if a projected state
        carries appreciable coupling into the ECS tail (the same Eq. (67)
        guard `nonlocal_operator` applies, for the same reason).
    """
    n_avail = ing.E_n.shape[1]
    n_use = n_avail if n_states is None else int(n_states)
    if n_use > n_avail or n_use < 0:
        raise ValueError(f"n_states={n_states} outside the available range [0, {n_avail}]")
    check_nodes_coincide(ing.R, nuclear_grid)

    t_nuc = kinetic_sparse(nuclear_grid, model.mu)
    v0 = np.asarray(model.v0(nuclear_grid.points), dtype=np.complex128)
    v_d = continue_to_tail(ing.v_d_discrete, ing.R, nuclear_grid)
    tail = nuclear_grid.points.imag != 0.0

    h_d = (t_nuc + sp.diags(v_d)).tocsr()
    if n_use == 0:
        return h_d

    arms: list[sp.csr_matrix] = []
    couplings: list[sp.csr_matrix] = []
    for n in range(n_use):
        v_dn = continue_to_tail(ing.V_dn[:, n], ing.R, nuclear_grid)
        check_tail_coupling(v_dn, tail, n)
        e_n = continue_to_tail(ing.E_n[:, n], ing.R, nuclear_grid)
        arms.append((t_nuc + sp.diags(v0 + e_n)).tocsr())
        couplings.append(sp.diags(v_dn).tocsr())

    # (N_R, n_use * N_R). scipy-stubs' bmat/block_diag overloads don't resolve
    # cleanly for a mixed csr_matrix/csr_array block list (same interaction
    # noted in linalg/sparse_lu.py); `cast` sidesteps it. Every block is
    # complex128: `coup`/`diag_arms` via the explicit `dtype=` above, `h_d`
    # because `continue_to_tail` and `kinetic_sparse` already are.
    coup: sp.csr_matrix = sp.hstack(couplings, format="csr", dtype=np.complex128)
    off_diag: sp.csr_matrix = coup.T.tocsr()
    diag_arms = cast("sp.csr_matrix", sp.block_diag(arms, format="csr", dtype=np.complex128))
    blocks: list[list[sp.csr_matrix | None]] = [[h_d, coup], [off_diag, diag_arms]]
    out = cast("sp.csr_matrix", sp.bmat(blocks, format="csr"))
    return out
