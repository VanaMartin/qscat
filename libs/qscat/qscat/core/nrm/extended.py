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

__all__ = [
    "LaunchBasis",
    "extended_hamiltonian",
    "initial_packet",
    "lcp_initial_packet",
    "lcp_limit_hamiltonian",
]


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


def lcp_limit_hamiltonian(
    nuclear_grid: FemDvrEcsGrid,
    model: ResonanceModel,
    v_res: npt.NDArray[np.complex128],
) -> sp.csr_matrix:
    """PRA 47 Eq. (2.15): the `d`-block alone, with a LOCAL complex potential.

        i d/dt Psi_d = [T_N + V_d(R) + Delta_L(R) - (i/2) Gamma_L(R)] Psi_d

    Eq. (2.11) is the statement this builds on: the memory kernel of Eq. (2.1)
    collapses to `i[Delta_L(R) - (i/2)Gamma_L(R)] delta(R-R') delta(t)`, i.e.
    the local complex potential IS the Markovian approximation to the nonlocal
    model. Localizing the kernel removes the arms entirely, so this is
    `extended_hamiltonian` with `n_states=0` and a different diagonal -- and
    `N_R` square instead of `(1 + n_states) * N_R`, which is why the local
    route costs seconds where the nonlocal one costs tens of minutes.

    WHICH `V_d` GOES IN `v_res`, AND WHY IT IS NOT THIS PACKAGE'S. Eq. (2.14),
    `E_res(R) - V_d(R) + V_0(R) - Delta(E_res(R),R) = 0`, rearranges to
    `V_d + Delta_L = E_res + V_0`, so Eq. (2.15)'s bracket is the RESONANCE
    POSITION plus the neutral curve -- exactly `qscat.core.lcp`'s `Vd` (whose
    `model.surface` already carries `v0`), NOT this package's
    `NrmIngredients.v_d_discrete` (PRA 77 Eq. 20, which is `V_d` WITHOUT the
    level shift `Delta_L`). The two differ by 0.0053 Ha on F2 and 0.0229 Ha on
    NO -- but WHERE that is measured is worth a factor of 70, so read the
    numbers with their definition attached. `docs/physics/nonlocal-resonance-
    model.md` sec. 4's 0.0053 Ha (F2) is at the NRM doorway `argmax |chi_0|`
    (R = 2.745); at the LCP doorway `argmax sqrt(Gamma/2pi)|chi_0|`
    (R = 2.486, 0.26 bohr further in) it is 0.0423 Ha instead, and over
    the region where `|chi_0|` exceeds 5% of its maximum the difference sweeps
    0.00095 -> 0.067 Ha. The DECK is worth 1.3x by comparison (0.04231 on a
    55-point electronic grid against 0.03283 on the 132-point production one,
    same nuclear deck; at the NRM doorway 0.00535 against 0.00534, i.e.
    deck-stable to three figures).

    The difference is not cosmetic: substituting `v_d_discrete` for
    `qscat.core.lcp`'s `Vd` here does NOT reproduce `lcp_da_cross_section`
    (measured ratios in `docs/physics/nrm-time-dependent.md` sec. 6). It decays
    over the same `R`-range as `Gamma` -- but NOT only where `Gamma` is nonzero:
    `Gamma_L(R)` is `Gamma` at ONE energy while `Delta_L(R) = P Int (dE'/2pi)
    Gamma(E',R)/(E_res - E')` (Eq. 2.12a/2.13b) integrates over ALL `E'`, so
    `Delta_L` is free to be nonzero where `Gamma_L` vanishes, and measurably is
    (9.5e-4 Ha at R = 3.01, where `Gamma = 0`).

    Callers therefore build `v_res = Vd - (i/2) Gamma` from
    `qscat.core.lcp.local_complex_potential`. Passing `v_res` in rather than
    computing it keeps `qscat.core.nrm` from depending on `qscat.core.lcp`, and
    keeps that measurement reproducible from outside.

    Parameters
    ----------
    nuclear_grid : FemDvrEcsGrid
        The nuclear radial grid.
    model : ResonanceModel
        Consumed only for `model.mu`.
    v_res : ndarray
        `V_d + Delta_L - (i/2) Gamma_L` on the FULL nuclear grid (real nodes
        and ECS tail), `nuclear_grid.n` entries.

    Returns
    -------
    sp.csr_matrix
        `T_N + diag(v_res)`, `(N_R, N_R)`, complex symmetric.

    Raises
    ------
    ValueError
        If `v_res` does not have one entry per nuclear grid point.
    """
    v = np.asarray(v_res, dtype=np.complex128)
    if v.shape != (nuclear_grid.n,):
        raise ValueError(
            f"v_res has shape {v.shape}, expected ({nuclear_grid.n},) -- one entry per nuclear node"
        )
    out: sp.csr_matrix = (kinetic_sparse(nuclear_grid, model.mu) + sp.diags(v)).tocsr()
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
        """Number of separable-expansion vectors (columns of `vectors`)."""
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


def lcp_initial_packet(
    nuclear_grid: FemDvrEcsGrid,
    gamma: npt.NDArray[np.float64],
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    energies: npt.ArrayLike,
) -> LaunchBasis:
    """The LOCAL launch state `sqrt(Gamma_L(R)/2pi) chi_v(R)`, as a rank-1 basis.

    The Markovian counterpart of `initial_packet`, and the doorway
    `qscat.core.lcp.lcp_da_cross_section` already solves with. It is
    ENERGY-INDEPENDENT -- `Gamma_L(R) = Gamma(E_res(R), R)` (PRA 47 Eq. 2.12a)
    is evaluated at the R-dependent resonance position, not at the incident
    energy -- so the launch matrix is EXACTLY rank 1 and no SVD is taken: the
    single column is the normalized doorway and `coeffs` is its norm repeated
    across energies. Rank 1 is what Eq. (2.17) claims for the nonlocal launch
    state under Eq. (2.16)'s separability; here it is not an approximation but
    the shape of the local model.

    WHY THIS DOORWAY AND NOT `initial_packet`'s. Eq. (2.11) localizes the
    KERNEL; read strictly, Eq. (2.5)'s launch state `V_dk_i(R) chi_v(R)` is
    untouched by it, and `V_dk_i` carries the INCIDENT energy `E_i` where
    `Gamma_L` carries `E_res(R)`. The two coincide only under Eq. (2.16)'s
    separability with an `E`-independent `Gamma(E)`, which these models do not
    have. `qscat.core.lcp` -- the local model this repository ships, measures,
    and calls "the LCP" -- uses `sqrt(Gamma_L/2pi) chi_v`, so that is what the
    Markovian route must use to be comparable to it at all. Keeping Eq. (2.5)'s
    launch against the local kernel is a THIRD model, neither the shipped LCP
    nor the nonlocal one; it is measured and reported in
    `docs/physics/nrm-time-dependent.md` sec. 6 rather than offered as an
    option here.

    Parameters
    ----------
    nuclear_grid : FemDvrEcsGrid
        The nuclear radial grid.
    gamma : ndarray
        `Gamma_L(R) >= 0` on the full nuclear grid (`qscat.core.lcp.
        local_complex_potential`'s second return), `nuclear_grid.n` entries.
    eps, chi : ndarray
        Neutral vibrational energies and states (`qscat.core.vibrational`),
        exactly as on `initial_packet`.
    v_init : int
        Initial vibrational level.
    energies : array_like
        Incident electron kinetic energies (hartree), all strictly positive.
        They enter only through `e_total`; the launch column itself is the
        same for every one of them.

    Returns
    -------
    LaunchBasis
        `rank == 1`, `truncation_error == 0.0` (exact, not truncated).

    Raises
    ------
    ValueError
        If any energy is non-positive, if `gamma` is not one entry per
        nuclear node, if any entry of `gamma` is negative, or if the doorway
        vanishes identically (a `Gamma` that is zero everywhere describes no
        resonance at all).
    """
    e = np.atleast_1d(np.asarray(energies, dtype=np.float64))
    if np.any(e <= 0.0):
        raise ValueError("incident energies must be positive")
    g = np.asarray(gamma, dtype=np.float64)
    if g.shape != (nuclear_grid.n,):
        raise ValueError(
            f"gamma has shape {g.shape}, expected ({nuclear_grid.n},) -- one entry per nuclear node"
        )
    # A width is non-negative, and `local_complex_potential` clamps it. Left
    # unchecked, one negative entry makes `sqrt` NaN, the NaN survives the
    # propagation and the transform, and `sigma_DA` comes back NaN with
    # nothing saying where it started.
    if np.any(g < 0.0):
        raise ValueError(
            f"gamma has {int(np.sum(g < 0.0))} negative entries (min {float(g.min()):.3g}); "
            "a width cannot be negative and sqrt(Gamma/2pi) would be NaN"
        )

    doorway = np.sqrt(g / (2.0 * np.pi)).astype(np.complex128) * chi[v_init]
    norm = float(np.linalg.norm(doorway))
    if norm == 0.0:
        raise ValueError("the LCP doorway sqrt(Gamma/2pi)*chi is identically zero")

    return LaunchBasis(
        vectors=(doorway / norm)[:, None],
        coeffs=np.full((1, e.size), norm, dtype=np.complex128),
        energies=e,
        e_total=e + float(eps[v_init]),
        truncation_error=0.0,
    )
