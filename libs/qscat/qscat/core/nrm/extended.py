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

`initial_packet`/`LaunchBasis` build the OTHER half of the propagation this
Hamiltonian drives: PRA 47 Eq. (2.5)'s launch state `V_dk_i(R) chi_v(R)`,
factorized by an SVD along the energy axis so a whole energy sweep is
propagated from a handful of columns instead of one propagation per energy.

WHY THE SVD IS LEGITIMATE. PRA 47 Eq. (2.17) removes the launch state's
energy dependence analytically, but only under Eq. (2.16)'s separability,
which these models do not satisfy exactly -- they satisfy it *numerically*.
Measured 2026-08-19 (`AsymptoticDiscreteState`, i.e. choice B, `n_states=3`,
9 energies per window), the ratios `sigma_j / sigma_1` of the launch
matrix's singular values are

    F2 DA  (0.010-0.050 Ha): 1, 5.69e-3, 2.36e-4, 5.30e-7
    N2 VE  (0.060-0.160 Ha): 1, 9.78e-4, 1.23e-6, 1.0e-8

i.e. `r=2-3` columns reconstruct the whole sweep to ~1e-4-1e-7, not `r=1`
(Eq. 2.17's own claim) nor `r=n_E`. `PhysicalDiscreteState` (choice A) is
NOT this low-rank: at `rank_tol=1e-6` it needs rank 7 (F2 DA) / rank 5 (N2
VE) against choice B's 3 / 1 at the same tolerance -- its `R`-dependent
`phi_d` breaks Eq. (2.16)'s separability more than choice B's `R`-independent
one does, so these ratios are NOT a universal property of the model, only of
the discrete-state choice actually measured. Measured 2026-08-19 on the F2 and
N2 fixtures of `libs/qscat/tests/test_nrm_extended.py`, over 9 energies per
window; `docs/physics/nonlocal-resonance-model.md` carries the full table.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, cast

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp

from qscat.dvr import FemDvrEcsGrid, kinetic_sparse

from .coupling import v_dk_plus
from .discrete_state import DiscreteState
from .ingredients import NrmIngredients
from .nonlocal_potential import (
    check_nodes_coincide,
    check_tail_coupling,
    continue_to_tail,
)

if TYPE_CHECKING:
    from qscat.model import ResonanceModel

__all__ = ["LaunchBasis", "extended_hamiltonian", "initial_packet"]


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


@dataclass(frozen=True)
class LaunchBasis:
    """A low-rank factorization of the energy-dependent launch state.

    PRA 47 Eq. (2.5)'s `Psi_d(R,0) = V_dk_i(R) chi_v(R)` at a batch of
    incident energies, SVD'd along the energy axis and truncated at
    `rank_tol`. `vectors @ coeffs` reconstructs the per-energy launch state
    (embedded in the extended-space block Hamiltonian's `d`-block, with the
    arm blocks zero) to the accuracy the truncation allows -- exactly, at
    `rank_tol=0.0`.

    `truncation_error` is `sigma_{r+1}/sigma_1` -- a bound relative to the
    LARGEST column's norm, not to each column's own. It does not bound the
    per-energy reconstruction error at a WEAK column (one whose norm is far
    below `sigma_1`): measured on F2's DA window at `rank_tol=1e-6`,
    `truncation_error=5.26e-7` while the worst per-column relative error is
    `1.15e-6`, about 2.2x larger, because that window's column norms differ
    by only ~2.5x. On a window where the launch state's magnitude varies by
    orders across energy (dissociative attachment near threshold does
    exactly that), the gap between `truncation_error` and the worst
    per-column error can be far larger -- callers who need a per-energy
    guarantee should compute it directly (`np.linalg.norm(vectors @ coeffs
    - M, axis=0) / np.linalg.norm(M, axis=0)`), not read `truncation_error`
    as one.
    """

    vectors: npt.NDArray[np.complex128]  # (N_ext, r), arm blocks zero
    coeffs: npt.NDArray[np.complex128]  # (r, n_E)
    energies: npt.NDArray[np.float64]  # (n_E,) incident electron kinetic energies
    e_total: npt.NDArray[np.float64]  # (n_E,) total energies (transform frequencies)
    truncation_error: float  # sigma_{r+1}/sigma_1 -- sigma_1-relative, NOT per-energy (see above)

    @property
    def rank(self) -> int:
        return int(self.vectors.shape[1])


def initial_packet(
    nuclear_grid: FemDvrEcsGrid,
    elec_grid: FemDvrEcsGrid,
    model: ResonanceModel,
    phi_d: DiscreteState,
    ing: NrmIngredients,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    energies: npt.ArrayLike,
    *,
    n_states: int | None = None,
    rank_tol: float = 1e-6,
) -> LaunchBasis:
    """PRA 47 Eq. (2.5) `Psi_d(R,0) = V_dk_i(R) chi_v(R)`, factorized.

    `V_dk_i` is `qscat.core.nrm.coupling.v_dk_plus` evaluated at each real
    incident electron kinetic energy in `energies` and continued onto the
    full (real + ECS-tail) nuclear grid, so this is the exact right-hand
    side `nrm_da_cross_section`/`nrm_ve_cross_section` build per energy --
    `propagation.propagate_nrm` only saves work if this launch state matches
    theirs exactly, which `test_launch_state_drives_the_time_independent_solve`
    checks by feeding the reconstruction straight into `solve_nuclear`.

    Eq. (2.17) removes the energy dependence analytically, but only under
    Eq. (2.16)'s separability, which these models do not satisfy exactly --
    they satisfy it *numerically* (module docstring's measured singular
    values). The launch matrix `M[R, j] = V_dk(R; E_j) chi_{v_init}(R)` --
    one column per incident energy `E_j` -- is SVD'd and truncated at
    `rank_tol`; the kept left singular vectors are what `propagate_nrm` steps
    (embedded in the `d`-block, arms zero, since the launch state lives
    entirely in the discrete-state block), and `coeffs` reconstructs every
    energy from them by linearity of the resolvent that eventually solves
    for `Psi_d`. Rank 1 IS Eq. (2.17); higher ranks are the controlled
    generalization the paper's approximation does not make available.
    `rank_tol=0.0` keeps every mode (`min(N_R, n_E)` of them) -- the
    truncation is then the only approximation this factorization
    introduces, to round-off.

    The SVD is of the RAW launch matrix, not a normalized one: truncating by
    relative singular value (`sv > rank_tol * sv[0]`) already makes
    `rank_tol` scale-free, so normalizing first would only add rounding.

    Parameters
    ----------
    nuclear_grid, elec_grid : FemDvrEcsGrid
        The nuclear and electronic radial grids.
    model : ResonanceModel
        Supplies `surface`, `v0`, and `ell` (via `v_dk_plus`).
    phi_d : DiscreteState
        The discrete-state choice under test.
    ing : NrmIngredients
        Energy-independent ingredients; only `ing.R` and `ing.E_n`'s state
        count are used (the latter bounds `n_states`).
    eps, chi : ndarray
        Neutral vibrational energies and states (`qscat.core.vibrational`).
        `eps` sits beside `chi` exactly as it does on the TI entry points,
        and is used only for `e_total = energies + eps[v_init]`.
    v_init : int
        Initial vibrational level.
    energies : array_like
        Incident electron kinetic energies (hartree), all strictly positive.
    n_states : int, optional
        The number of projected electronic arm blocks `vectors` is sized
        for (the extended-space layout `extended_hamiltonian` builds).
        `None` (default) uses every state `ing` carries.
    rank_tol : float, optional
        Truncate singular values below `rank_tol * sv[0]`. `0.0` keeps
        every mode (exact reconstruction to round-off); default `1e-6`.

    Returns
    -------
    LaunchBasis

    Raises
    ------
    ValueError
        If any energy is non-positive, if `n_states` is outside
        `[0, ing.E_n.shape[1]]`, or if `ing.R` does not coincide with the
        grid's real DVR nodes (the same guard `extended_hamiltonian` and
        `nonlocal_operator` apply, for the same reason -- a mismatched `ing`
        would otherwise give a silently wrong launch state).
    """
    e = np.atleast_1d(np.asarray(energies, dtype=np.float64))
    if np.any(e <= 0.0):
        raise ValueError("incident energies must be positive")
    n_avail = ing.E_n.shape[1]
    n_use = n_avail if n_states is None else int(n_states)
    if n_use > n_avail or n_use < 0:
        raise ValueError(f"n_states={n_states} outside the available range [0, {n_avail}]")
    check_nodes_coincide(ing.R, nuclear_grid)
    n_r = nuclear_grid.n

    m = np.empty((n_r, e.size), dtype=np.complex128)
    for j, e_kin in enumerate(e):
        v_dk_i = continue_to_tail(
            v_dk_plus(elec_grid, model, phi_d, ing.R, float(e_kin)), ing.R, nuclear_grid
        )
        m[:, j] = v_dk_i * chi[v_init]

    u, sv, vh = np.linalg.svd(m, full_matrices=False)
    keep = sv.size if rank_tol <= 0.0 else max(1, int(np.sum(sv > rank_tol * sv[0])))
    truncation_error = float(sv[keep] / sv[0]) if keep < sv.size else 0.0

    vectors = np.zeros(((1 + n_use) * n_r, keep), dtype=np.complex128)
    vectors[:n_r, :] = u[:, :keep]
    coeffs = sv[:keep, None] * vh[:keep, :]
    return LaunchBasis(
        vectors=vectors,
        coeffs=coeffs,
        energies=e,
        e_total=e + float(eps[v_init]),
        truncation_error=truncation_error,
    )
