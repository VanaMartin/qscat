"""The NRM's per-R ingredients: `E_n(R)`, `V_dn(R)`, `V_d(R)`.

PRA 77 Sec. IV, Eq. (56)-(59). At each nuclear node the electronic Hilbert
space is split by the Feshbach projectors and the P-space is diagonalized:

    Q = |phi_d><phi_d|,  P = 1 - Q                        Eq. (57)-(58)
    (P H_el P) phi_n = E_n(R) phi_n                       Eq. (56)
    V_dn(R) = <phi_d| H_el |phi_n>                        Eq. (59)
    V_d(R)  = V_0(R) + <phi_d| H_el |phi_d>               Eq. (20)

`P H_el P` is complex SYMMETRIC, not Hermitian, so eigenvectors are normalized
under the c-product (the paper: "we have to use for the wave functions the
scalar product defined without complex conjugation", p. 012710-6).

This is the expensive, energy-INDEPENDENT half of the calculation: one dense
electronic eigenproblem per nuclear node, computed once and reused across the
whole energy sweep.

**Sign, not just order, must be tracked.** `np.linalg.eig` returns an
arbitrary phase per node, and c-normalizing (`c_product(v, v) = 1`) fixes it
only up to `c^2 = 1`, i.e. an overall `+-1` (the bilinear c-product, unlike a
Hermitian norm, has no other unit-modulus solution). `_track` matches state
identity across `R` by nearest eigenvalue, which does nothing about that
sign: measured on F2 (`electronic_grid(r_max=16, order=8, n_complex=6)`,
`R = linspace(6.0, 1.8, 25)`), `c_product(phi_n(R_j), phi_n(R_{j+1}))` is
`-0.995` to `-1.000` at several `(n, R)` pairs with `|V_dn| ~ 0.1`, not tail
noise. Eq. (60) is bilinear in `V_dn(R_i)` and `V_dn(R_j)` for the SAME `n`,
so an untracked sign flips the sign of that state's contribution to `F` --
wrong, and not caught by any shape/finiteness check. `_sign_align` fixes this
by orienting each tracked eigenvector against its predecessor.

**The kernel is genuinely discontinuous in `R` for the physical discrete
state (choice A).** Independent of the sign bug above: at F2's bound/
resonance crossing (`R_c ~ 2.59`), `V_dn` jumps from `~1e-13` (bound branch,
identically zero by the bound eigenrelation, see below) to `~0.26` (scattering
branch) across a `dR = 0.02` step -- a real, physical step, not a numerical
artifact. Choice A switches which STATE `phi_d` even is (a bound eigenvector
vs. a truncated scattering function) at the crossing, and the paper's
`exp[-i delta(R)]` phase-fixing (p. 012710-6) -- meant "to obtain a discrete
state that varies smoothly at the crossing point" -- cannot smooth over that,
since one side is an eigenvector and the other is not. Choice B
(`AsymptoticDiscreteState`, R-independent by construction) has no such step
and varies smoothly across the same span. This is the per-`R` ingredient
layer's instance of the paper's own choice-A breakdown diagnosis
(p. 012710-8); Task 9 measures whether it drives the DA discrepancy.

**On the bound branch, `V_dn == 0` is exact, not a numerical artifact.** A
bound `phi_d` is (to numerical precision) an eigenvector of the FULL `H_el`,
so `H_el phi_d = E_d phi_d` and, since `P` annihilates `phi_d`
(`d @ P == 0`, `P`'s defining property), `V_dn = <phi_d|H_el|phi_n> =
E_d <phi_d|phi_n> = E_d (d @ P @ ...) = 0` identically -- for EVERY bound-
branch `R`, including right at the leakiest point of the crossing (F2's
`R ~ 2.6`, where 1.08% of `phi_d`'s c-norm sits on the ECS contour; see
`discrete_state.PhysicalDiscreteState`'s docstring). This is not incidental:
the paper states it directly (p. 012710-7, "Because the bound state is the
eigenfunction of `H_el`, the discrete-state-continuum coupling `V_dk(R)`
goes to zero" -- a nonzero coupling at large `R` "has no physical meaning"),
it is Eq. (67)'s decoupling holding EXACTLY rather than asymptotically, and
it matches `qscat.core.lcp`'s own `Gamma` support condition (no width where
there is no open channel).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

from qscat.dvr import FemDvrEcsGrid, eigen
from qscat.exceptions import ConvergenceError
from qscat.linalg import c_product

from .discrete_state import DiscreteState, electronic_hamiltonian

if TYPE_CHECKING:
    from qscat.model import ResonanceModel

__all__ = ["NrmIngredients", "nrm_ingredients"]

# A P H_el P eigenvector with |c_product(v, v)| below this is self-orthogonal
# to numerical precision and cannot be c-normalized.
_MIN_NORM2 = 1e-12
# The null mode's |E| must be this much smaller than the next-smallest, and
# its |<phi_d|v>| this much larger than any other state's, or the
# identification is ambiguous and we refuse to guess.
_NULL_ENERGY_RATIO = 1e-3
_NULL_OVERLAP_MIN = 0.5


@dataclass(frozen=True)
class NrmIngredients:
    """The energy-independent inputs to the nonlocal potential.

    Attributes
    ----------
    R : ndarray
        The nuclear nodes, descending, as supplied to `nrm_ingredients`.
    v_d_discrete : ndarray
        `V_d(R)` of Eq. (20), complex, INCLUDING `V_0(R)`. This is the paper's
        discrete-state potential, NOT `qscat.core.lcp`'s `Vd`.
    E_n : ndarray
        `(n_R, n_states)` projected electronic energies of Eq. (56),
        EXCLUDING `V_0(R)` (which Eq. 61 adds separately). Adiabatically
        tracked across `R`; the null mode is dropped.
    V_dn : ndarray
        `(n_R, n_states)` discrete-continuum couplings of Eq. (59), aligned
        with `E_n`.
    """

    R: npt.NDArray[np.float64]
    v_d_discrete: npt.NDArray[np.complex128]
    E_n: npt.NDArray[np.complex128]
    V_dn: npt.NDArray[np.complex128]


def _c_normalize_columns(
    vecs: npt.NDArray[np.complex128],
) -> npt.NDArray[np.complex128]:
    # c_product per column -- the bilinear (no-conjugation) norm the paper's
    # convention requires, not numpy's Hermitian v^dagger v.
    norms2 = np.array(
        [c_product(vecs[:, i], vecs[:, i]) for i in range(vecs.shape[1])],
        dtype=np.complex128,
    )
    if np.any(np.abs(norms2) < _MIN_NORM2):
        bad = int(np.argmin(np.abs(norms2)))
        raise ConvergenceError(
            f"P H_el P eigenvector {bad} is c-product self-orthogonal "
            f"(norm^2={norms2[bad]!r}); it cannot be c-normalized, so the "
            "P-space basis is unusable at this R"
        )
    out: npt.NDArray[np.complex128] = vecs / np.sqrt(norms2)
    return out


def _drop_null_mode(
    evals: npt.NDArray[np.complex128],
    vecs: npt.NDArray[np.complex128],
    d: npt.NDArray[np.complex128],
) -> tuple[npt.NDArray[np.complex128], npt.NDArray[np.complex128]]:
    """Remove the spurious zero mode `P` introduces (its eigenvector is phi_d).

    Identified by BOTH criteria -- smallest `|E|` and largest `|<phi_d|v>|` --
    and refuses rather than guesses when they disagree. Index-based selection
    would be wrong: the eigensolver's ordering is not stable across `R`.
    """
    by_energy = int(np.argmin(np.abs(evals)))
    overlaps = np.abs(d @ vecs)
    by_overlap = int(np.argmax(overlaps))
    if by_energy != by_overlap:
        raise ConvergenceError(
            f"ambiguous P H_el P null mode: smallest |E| at index {by_energy} "
            f"(E={evals[by_energy]!r}) but largest <phi_d|v> at index "
            f"{by_overlap} (overlap={overlaps[by_overlap]:.3g})"
        )
    if overlaps[by_overlap] < _NULL_OVERLAP_MIN:
        raise ConvergenceError(
            f"the candidate null mode overlaps phi_d by only "
            f"{overlaps[by_overlap]:.3g}; expected ~1"
        )
    keep = np.ones(evals.size, dtype=bool)
    keep[by_energy] = False
    rest = np.abs(evals[keep])
    if rest.size and np.abs(evals[by_energy]) > _NULL_ENERGY_RATIO * rest.min():
        raise ConvergenceError(
            f"the candidate null mode's |E|={np.abs(evals[by_energy]):.3g} is "
            f"not clearly separated from the next-smallest {rest.min():.3g}"
        )
    return evals[keep], vecs[:, keep]


def _track(
    prev: npt.NDArray[np.complex128], evals: npt.NDArray[np.complex128]
) -> npt.NDArray[np.intp]:
    """Greedy nearest-energy matching of `evals` onto `prev`'s ordering.

    Each previous state claims its nearest unclaimed new eigenvalue, in order.
    Without this the eigensolver's ordering permutes states between adjacent
    `R` and `E_n(R)`/`V_dn(R)` acquire discontinuities that corrupt `F`.
    """
    order = np.empty(prev.size, dtype=np.intp)
    taken = np.zeros(evals.size, dtype=bool)
    for i in range(prev.size):
        dist = np.abs(evals - prev[i])
        dist[taken] = np.inf
        j = int(np.argmin(dist))
        order[i] = j
        taken[j] = True
    return order


def _sign_align(
    prev_vecs: npt.NDArray[np.complex128], vecs: npt.NDArray[np.complex128]
) -> npt.NDArray[np.complex128]:
    """Orient each (already order-tracked) eigenvector against its predecessor.

    c-normalization (`c_product(v, v) = 1`) fixes an eigenvector only up to
    `c` with `c^2 = 1`, i.e. `c = +-1` -- the bilinear c-product, unlike a
    Hermitian inner product, admits no other unit-modulus phase. `_track`
    matches state IDENTITY across `R` by nearest eigenvalue and says nothing
    about that residual sign, so `np.linalg.eig` is free to flip it from one
    `R` to the next. Left unaligned, a flipped sign flips `V_dn`'s sign for
    that state at that `R`, which corrupts Eq. (60): `F(E,R_i,R_j)` is
    bilinear in `V_dn(R_i)` and `V_dn(R_j)` for the SAME `n`. Columns must
    already be in `prev_vecs`'s state order (post-`_track`) before calling
    this -- it aligns sign only, not identity.
    """
    overlaps = np.array(
        [c_product(prev_vecs[:, i], vecs[:, i]) for i in range(vecs.shape[1])],
        dtype=np.complex128,
    )
    # The residual freedom is exactly +-1 (real), so a genuine sign flip
    # shows up as overlap ~= -1 and a preserved sign as overlap ~= +1; taking
    # the sign of the real part is robust to the small complex noise
    # `np.linalg.eig` introduces around that +-1.
    signs = np.sign(overlaps.real)
    signs[signs == 0.0] = 1.0
    out: npt.NDArray[np.complex128] = vecs * signs
    return out


def nrm_ingredients(
    elec_grid: FemDvrEcsGrid,
    model: ResonanceModel,
    phi_d: DiscreteState,
    R_values: npt.NDArray[np.float64],
) -> NrmIngredients:
    """Build `E_n(R)`, `V_dn(R)` and `V_d(R)` on the given nuclear nodes.

    Parameters
    ----------
    elec_grid : FemDvrEcsGrid
        The electronic radial grid.
    model : ResonanceModel
        Supplies `surface`, `v0` and `ell`.
    phi_d : DiscreteState
        The discrete state (choice A or B).
    R_values : ndarray
        Nuclear nodes, strictly DESCENDING. The tracking walk starts at the
        largest `R`, where `phi_d -> phi_b` and the P-space states are cleanly
        ordered, and continues inward.

    Returns
    -------
    NrmIngredients

    Raises
    ------
    ValueError
        If `R_values` is not strictly descending.
    ConvergenceError
        If the `P H_el P` null mode cannot be identified unambiguously, or an
        eigenvector is c-product self-orthogonal.
    """
    R = np.asarray(R_values, dtype=np.float64)
    if R.size > 1 and np.any(np.diff(R) >= 0.0):
        raise ValueError("R_values must be strictly descending")

    n_states = elec_grid.n - 1
    v_d = np.empty(R.size, dtype=np.complex128)
    e_n = np.empty((R.size, n_states), dtype=np.complex128)
    v_dn = np.empty((R.size, n_states), dtype=np.complex128)
    ident = np.eye(elec_grid.n, dtype=np.complex128)
    prev_evals: npt.NDArray[np.complex128] | None = None
    prev_vecs: npt.NDArray[np.complex128] | None = None

    for j in range(R.size):
        d = phi_d.phi_d(float(R[j]))
        h_el = electronic_hamiltonian(elec_grid, model, float(R[j]))
        p = ident - np.outer(d, d)  # Eq. (57)-(58)
        php = p @ h_el @ p
        evals, vecs = eigen(php)  # complex symmetric; Hermitian-normalized
        vecs = _c_normalize_columns(vecs)
        evals, vecs = _drop_null_mode(evals, vecs, d)
        if prev_evals is not None and prev_vecs is not None:
            order = _track(prev_evals, evals)
            evals, vecs = evals[order], vecs[:, order]
            vecs = _sign_align(prev_vecs, vecs)  # identity tracked; sign not
        prev_evals, prev_vecs = evals, vecs

        v_d[j] = complex(model.v0(float(R[j]))) + d @ (h_el @ d)  # Eq. (20)
        e_n[j] = evals  # Eq. (56)
        v_dn[j] = d @ (h_el @ vecs)  # Eq. (59)

    return NrmIngredients(R=R, v_d_discrete=v_d, E_n=e_n, V_dn=v_dn)
