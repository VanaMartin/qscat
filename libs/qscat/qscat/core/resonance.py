"""Exact (non-Born-Oppenheimer) resonance states of the 2-D model.

A resonance of the full electronic x nuclear problem is a pole of the 2-D
S-matrix: an eigenvalue `E_r - i*Gamma/2` of the complex-scaled 2-D Hamiltonian
with outgoing boundary conditions in BOTH coordinates. No Born-Oppenheimer
separation, no discrete-state choice, no local approximation -- so its
eigenvector `Psi(r, R)` does not in general factorize into `phi(r) chi(R)`.

That is what distinguishes these states from `qscat.core.lcp.resonance_levels`,
which builds the Born-Oppenheimer approximation to the same quantities by
diagonalizing the nuclear problem inside the local complex potential curve.
Comparing the two measures the non-adiabatic error directly, in the pole
positions and widths, rather than through a cross section that mixes many poles
with a background.

Identification is by ECS angle stability, generalized to two angles. Under
complex scaling a resonance eigenvalue is (nearly) independent of the rotation
angle while the discretized continuum rotates with it -- but a 2-D problem has
two continua, electronic (`r -> inf`) and nuclear (`R -> inf`), each rotating
with its own angle. A state is accepted only if it survives BOTH tests: it must
be stable when `theta_r` moves and when `theta_R` moves. Running them separately
rather than moving both angles at once costs one extra spectrum and buys the
diagnosis of WHICH continuum a rejected state belonged to.

See `docs/physics/exact-2d-resonances.md`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

from qscat.dvr import TensorGrid
from qscat.ecs import match_angle_stable
from qscat.linalg import ShiftInvertEigs, c_product

if TYPE_CHECKING:
    from qscat.model import ResonanceModel

__all__ = ["ExactResonanceStates", "exact_resonance_states"]

# Two eigenvalues closer than this (relative to |E|) are the same state seen
# from two shifts, not two states: shift-invert windows around neighbouring
# seeds overlap, so the pooled spectrum contains duplicates by construction.
_DEDUPE_RTOL = 1e-9


@dataclass(frozen=True)
class ExactResonanceStates:
    """Exact 2-D resonance states and the evidence that they are resonances.

    Attributes
    ----------
    energies : ndarray of complex128, shape (m,)
        `E_r - i*Gamma/2`, ascending in `Re E`, taken from the base grid.
    widths : ndarray of float64, shape (m,)
        `Gamma = -2 Im E`.
    states : ndarray of complex128, shape (n, m)
        `states[:, i]` is the 2-D eigenvector of `energies[i]` on the base
        grid, flattened in `TensorGrid` order and c-product normalized.
    residual_electronic : ndarray of float64, shape (m,)
        `|E_base - E_theta_r|`: how far the eigenvalue moved when the
        ELECTRONIC ECS angle changed. Small means the state does not live in
        the electronic continuum.
    residual_nuclear : ndarray of float64, shape (m,)
        `|E_base - E_theta_R|`: the same for the NUCLEAR angle.

    Notes
    -----
    Both residuals are reported rather than a single midpoint energy. With two
    independent angle partners there is no one midpoint, and the pair carries
    strictly more information: a state can be solidly bound in one coordinate
    and marginal in the other, which one number would hide.

    A small residual is necessary evidence, not sufficient proof. It says the
    eigenvalue did not move when the contour did; a grid too coarse to resolve
    the state can produce that for the wrong reason. Convergence in the grid is
    a separate question, answered by refining it.
    """

    energies: npt.NDArray[np.complex128]
    widths: npt.NDArray[np.float64]
    states: npt.NDArray[np.complex128]
    residual_electronic: npt.NDArray[np.float64]
    residual_nuclear: npt.NDArray[np.float64]


def _pooled_spectrum(
    model: ResonanceModel,
    tgrid: TensorGrid,
    shifts: npt.ArrayLike,
    k: int,
) -> tuple[npt.NDArray[np.complex128], npt.NDArray[np.complex128]]:
    """Eigenpairs near every shift, on one grid, pooled and de-duplicated.

    The Hamiltonian is factored once per shift through a single
    `ShiftInvertEigs`, so the symbolic analysis is reused across the sweep --
    the reason a multi-seed search is not a multiple of a single-seed one.
    """
    H = model.hamiltonian(tgrid)
    solver = ShiftInvertEigs(H, k=k)
    vals: list[npt.NDArray[np.complex128]] = []
    vecs: list[npt.NDArray[np.complex128]] = []
    for sigma in np.atleast_1d(np.asarray(shifts, dtype=np.complex128)):
        v, w = solver.near(complex(sigma))
        vals.append(v)
        vecs.append(w)
    all_vals = np.concatenate(vals)
    all_vecs = np.concatenate(vecs, axis=1)

    keep: list[int] = []
    for i, e in enumerate(all_vals):
        if not any(abs(e - all_vals[j]) <= _DEDUPE_RTOL * max(abs(e), 1.0) for j in keep):
            keep.append(i)
    idx = np.asarray(keep, dtype=np.intp)
    return all_vals[idx], all_vecs[:, idx]


def exact_resonance_states(
    model: ResonanceModel,
    grid_base: TensorGrid,
    grid_electronic: TensorGrid,
    grid_nuclear: TensorGrid,
    *,
    shifts: npt.ArrayLike,
    window: tuple[float, float, float, float],
    k: int = 8,
    rel_tol: float = 1e-4,
    atol: float = 1e-8,
) -> ExactResonanceStates:
    """Exact 2-D resonance states near `shifts`, by two-angle ECS stability.

    Parameters
    ----------
    model : ResonanceModel
        Supplies the sparse 2-D Hamiltonian via `model.hamiltonian(tgrid)`.
    grid_base : TensorGrid
        The tensor-product grid the returned energies and states belong to.
    grid_electronic : TensorGrid
        The same grid with the ELECTRONIC ECS angle changed, and only that.
    grid_nuclear : TensorGrid
        The same grid with the NUCLEAR ECS angle changed, and only that.
    shifts : array_like of complex
        Seed shifts to search near. Physically motivated seeds are the
        Born-Oppenheimer levels from `qscat.core.lcp.resonance_levels`; they are
        passed in rather than computed here so that the exact solver never
        depends on the approximation it exists to measure.
    window : tuple of float
        `(re_lo, re_hi, im_lo, im_hi)`, the complex-plane box the angle-stability
        match is restricted to.
    k : int, optional
        Eigenpairs per shift per grid (default 8).
    rel_tol, atol : float, optional
        Angle-stability acceptance: a state is kept when its partner satisfies
        `|E_a - E_b| < max(rel_tol*|E_a|, atol)`, for BOTH partners.

    Returns
    -------
    ExactResonanceStates
        Possibly empty. An empty result is a normal outcome -- no angle-stable
        state near these shifts -- and not an error.

    Raises
    ------
    ValueError
        If `window` catches no eigenvalue at all in one of the three spectra,
        which means the window or the shifts are misplaced rather than that
        nothing is there.

    Notes
    -----
    Cost is one sparse factorization per (grid, shift): three grids times the
    number of seeds. On a production 2-D deck that factorization dominates
    everything else, so keep the seed list short and prefer the MUMPS backend.
    """
    vals_a, vecs_a = _pooled_spectrum(model, grid_base, shifts, k)
    vals_b, _ = _pooled_spectrum(model, grid_electronic, shifts, k)
    vals_c, _ = _pooled_spectrum(model, grid_nuclear, shifts, k)

    _, res_ab, idx_ab = match_angle_stable(vals_a, vals_b, window, rel_tol=rel_tol, atol=atol)
    _, res_ac, idx_ac = match_angle_stable(vals_a, vals_c, window, rel_tol=rel_tol, atol=atol)

    # Stable under BOTH angle changes. Intersecting the index sets (rather than
    # the energies) keeps each accepted state tied to its own eigenvector and to
    # both of its residuals.
    res_by_index_ab = dict(zip(idx_ab.tolist(), res_ab.tolist(), strict=True))
    res_by_index_ac = dict(zip(idx_ac.tolist(), res_ac.tolist(), strict=True))
    accepted = sorted(set(res_by_index_ab) & set(res_by_index_ac))

    if not accepted:
        empty_c = np.empty(0, dtype=np.complex128)
        empty_f = np.empty(0, dtype=np.float64)
        return ExactResonanceStates(
            energies=empty_c,
            widths=empty_f,
            states=np.empty((vecs_a.shape[0], 0), dtype=np.complex128),
            residual_electronic=empty_f,
            residual_nuclear=empty_f,
        )

    idx = np.asarray(accepted, dtype=np.intp)
    energies = vals_a[idx]
    order = np.argsort(energies.real)
    idx = idx[order]
    energies = energies[order]

    states = vecs_a[:, idx]
    for i in range(states.shape[1]):
        states[:, i] /= np.sqrt(c_product(states[:, i], states[:, i]))

    return ExactResonanceStates(
        energies=np.asarray(energies, dtype=np.complex128),
        widths=np.asarray(-2.0 * energies.imag, dtype=np.float64),
        states=np.asarray(states, dtype=np.complex128),
        residual_electronic=np.asarray([res_by_index_ab[int(i)] for i in idx], dtype=np.float64),
        residual_nuclear=np.asarray([res_by_index_ac[int(i)] for i in idx], dtype=np.float64),
    )
