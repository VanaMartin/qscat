"""The nonlocal, complex, energy-dependent potential `F(E,R,R')`.

PRA 77 Eq. (53), evaluated by the ECS+DVR scheme of Eq. (60)-(61):

    F(E,R_i,R_j) = sum_n sqrt(W_i) V_dn(R_i) M(n)^-1_ij V_dn(R_j) sqrt(W_j)
    M(n)_ij      = [E - T_R - V_0(R) - E_n(R)]_ij

This is what distinguishes the nonlocal model from the LCP: `F` couples
different `R`, and depends on the total energy. Its local limit,
`F -> -(i/2) Gamma(E,R) delta(R-R')`, collapses Eq. (52) to the LCP nuclear
equation `qscat.core.lcp` solves -- the bridge `test_nrm_dissociation.py`
gates on.

WEIGHTS. Eq. (60)'s `sqrt(W)` factors are the DVR function-value to
coefficient conversion. This package works entirely in coefficient space (as
`lcp_da_cross_section` does), where `M(n)^-1` already IS the coefficient-space
Green's function and multiplication by `V_dn(R)` is `diag(V_dn)`. The weights
are therefore already absorbed and must NOT be applied again.

No singularity arises in the sum: `E_n(R)` is complex while `E` is real
(p. 012710-6). This is exactly what the ECS discretization buys over the
analytic treatment of the singular `k`-integral in Eq. (53).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

from qscat.dvr import FemDvrEcsGrid, kinetic

from .ingredients import NrmIngredients

if TYPE_CHECKING:
    from qscat.model import ResonanceModel

__all__ = [
    "TAIL_COUPLING_MAX",
    "check_nodes_coincide",
    "check_tail_coupling",
    "continue_to_tail",
    "nonlocal_operator",
]

# The ECS tail must carry no appreciable coupling. Measured on F2/NO with
# AsymptoticDiscreteState built at R_inf = R0 (the production choice):
# max|V_dn(R0)| is ~1e-12-2e-12 and FLAT in R0 (F2: 1.8e-12 at R0=6.0, 2.4e-12
# at R0=10.7; NO: 6.8e-13 at R0=6.0, 1.0e-12 at R0=10.7) -- because a discrete
# state built at R_inf = R0 is then an exact eigenvector of H_el(R0), so
# V_dn(R0) is ZERO IDENTICALLY regardless of how large or small the box is.
# Shrinking the box does NOT walk this guard toward tripping. What it actually
# catches is a TRUNCATED or otherwise malformed ingredient set -- one whose
# outermost node isn't genuinely the phi_d's R_inf (e.g. phi_d pinned at
# R_inf=10.7 but the ingredient range only sampled out to R=6.0): measured
# 7.2e-7 there, within a factor of ~1.4 of this threshold and ~1e6x the
# genuine R0 value above. See task-6-report.md for the full measurement.
TAIL_COUPLING_MAX = 1e-6


def continue_to_tail(
    values: npt.NDArray[np.complex128],
    ing_R: npt.NDArray[np.float64],
    nuclear_grid: FemDvrEcsGrid,
) -> npt.NDArray[np.complex128]:
    """Map an ingredient defined on `ing_R` onto the full nuclear grid.

    The ingredients are electronic quantities evaluated at real `R`; the ECS
    tail is filled with the outermost-real value (the design note's "ECS
    tail" rule). For a discrete state satisfying Eq. (67) that value is ~0.
    """
    pts = nuclear_grid.points
    real = pts.imag == 0.0
    out = np.empty(nuclear_grid.n, dtype=np.complex128)
    idx = np.argmin(np.abs(ing_R[None, :] - pts[real].real[:, None]), axis=1)
    out[real] = values[idx]
    out[~real] = values[int(np.argmax(ing_R))]  # outermost real R
    return out


def check_nodes_coincide(ing_R: npt.NDArray[np.float64], nuclear_grid: FemDvrEcsGrid) -> None:
    """Require `ing_R` to be exactly the nuclear grid's real DVR nodes.

    `continue_to_tail` maps by NEAREST real R, which degrades silently rather
    than erroring on a mismatched ingredient set: e.g. an every-4th-node
    subsample still produces *some* output -- a piecewise-constant
    continuation with node-to-node jumps up to the full coupling magnitude --
    with no shape or finiteness signal that anything is wrong. `V_dn` is
    partly self-limiting (a large enough mismatch eventually trips the tail
    guard), but `E_n` is not: it is O(0.005-70) across the real region and
    would be silently flattened onto whichever handful of ingredient nodes
    happen to be nearest. Checked once per call, on `ing.R` against the full
    node set (independent of `n_states`, since both `E_n` and `V_dn` share it).
    """
    real_R = np.sort(nuclear_grid.points[nuclear_grid.points.imag == 0.0].real)
    ing_sorted = np.sort(ing_R)
    if ing_sorted.size != real_R.size or not np.allclose(ing_sorted, real_R, rtol=0.0, atol=1e-9):
        raise ValueError(
            f"ingredient nodes (n={ing_sorted.size}) do not coincide with "
            f"nuclear_grid's real nodes (n={real_R.size}); nonlocal_operator "
            "requires ing.R to be exactly the grid's real DVR nodes -- a "
            "subsampled or otherwise mismatched ingredient set would "
            "otherwise be silently piecewise-constant-continued rather than "
            "rejected (see nrm_ingredients, which must be called with "
            "R_values = this same nuclear_grid's real points)"
        )


def check_tail_coupling(
    v_dn: npt.NDArray[np.complex128], tail: npt.NDArray[np.bool_], n: int
) -> None:
    """Require state `n`'s continued coupling `v_dn` to vanish on the ECS tail.

    Eq. (67): a discrete state that has genuinely decoupled by `R0` leaves
    `V_dn` ~0 there (see `continue_to_tail`'s "outermost-real value" rule and
    the `TAIL_COUPLING_MAX` measurements above it). A nonzero tail value means
    either `phi_d` does not satisfy Eq. (67) on this grid, or `ing.R`'s
    outermost node is not actually where `phi_d` has decoupled (e.g. a
    truncated ingredient range) -- shared by `nonlocal_operator` and
    `extended_hamiltonian` so the same threshold and diagnosis apply to both.
    """
    tail_coupling = np.max(np.abs(v_dn[tail])) if np.any(tail) else 0.0
    if tail_coupling > TAIL_COUPLING_MAX:
        raise ValueError(
            f"state {n} carries coupling {tail_coupling:.3g} into the "
            "ECS tail; either the discrete state does not satisfy Eq. "
            "(67) at ing.R's outermost node, or that node is not "
            "actually where phi_d has decoupled (e.g. a truncated "
            "ingredient range) -- F would pick up an unphysical tail "
            "contribution either way"
        )


def nonlocal_operator(
    ing: NrmIngredients,
    nuclear_grid: FemDvrEcsGrid,
    model: ResonanceModel,
    e_total: float,
    *,
    n_states: int | None = None,
) -> npt.NDArray[np.complex128]:
    """`F(E)` in nuclear DVR coefficient space -- Eq. (60)-(61).

    Parameters
    ----------
    ing : NrmIngredients
        The energy-independent ingredients from `nrm_ingredients`.
    nuclear_grid : FemDvrEcsGrid
        The nuclear radial grid (exterior-complex-scaled).
    model : ResonanceModel
        Supplies `mu` and `v0`.
    e_total : float
        The TOTAL energy `E` (electron kinetic + initial vibrational),
        hartree. Eq. (61)'s `E`.
    n_states : int or None, optional
        Truncate the sum over `n` to the lowest `n_states` states, in the
        order `ing.E_n`/`ing.V_dn` already carry -- `nrm_ingredients`'s
        adiabatic tracking, seeded by `qscat.dvr.eigen`'s ordering at the
        LARGEST `R` in `R_values`. That seed order is ascending `Re(E_n)`
        only (not `|E_n|` or any width-aware criterion), so "lowest
        `n_states`" means lowest real part at the outermost node, not
        necessarily lowest at every `R`. `None` (default) sums over all of
        them. The convergence knob Task 8 sweeps.

    Returns
    -------
    ndarray
        The `(n, n)` complex-symmetric nonlocal potential matrix.

    Raises
    ------
    ValueError
        If `n_states` exceeds the number of available states, or is negative.
        If `ing.R` does not coincide with `nuclear_grid`'s real DVR nodes.
        If a projected state carries non-negligible coupling into the ECS
        tail (the discrete state does not satisfy Eq. (67) on this grid).
    """
    n_avail = ing.E_n.shape[1]
    n_use = n_avail if n_states is None else int(n_states)
    if n_use > n_avail or n_use < 0:
        raise ValueError(f"n_states={n_states} outside the available range [0, {n_avail}]")
    check_nodes_coincide(ing.R, nuclear_grid)

    t_nuc = kinetic(nuclear_grid, model.mu)
    v0 = np.asarray(model.v0(nuclear_grid.points), dtype=np.complex128)
    ident = np.eye(nuclear_grid.n, dtype=np.complex128)
    tail = nuclear_grid.points.imag != 0.0

    out = np.zeros((nuclear_grid.n, nuclear_grid.n), dtype=np.complex128)
    for n in range(n_use):
        v_dn = continue_to_tail(ing.V_dn[:, n], ing.R, nuclear_grid)
        check_tail_coupling(v_dn, tail, n)
        e_n = continue_to_tail(ing.E_n[:, n], ing.R, nuclear_grid)
        m = e_total * ident - t_nuc - np.diag(v0 + e_n)  # Eq. (61)
        g = np.linalg.inv(m)
        out += v_dn[:, None] * g * v_dn[None, :]  # Eq. (60), coefficient space
    return out
