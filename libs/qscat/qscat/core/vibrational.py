"""Model-independent neutral-molecule vibrational states solver.

Promoted from `projects/n2_ti_cross_section/vibrational.py`: diagonalize
`T_nuc(mu) + diag(v0(R))` -> (eps_v, chi_v). ONE change from the original:
`v0` is a caller-supplied callable (the neutral molecule's potential-energy
curve) rather than a hardcoded import -- this is what keeps `qscat.core`
model-independent (it must never import `qscat.model` or any `projects.*`).

For a diatomic like N2 whose neutral potential well is entirely inside the
real region of the grid (bound states, R ~ 1.5-3 bohr), the eigenvalues are
real and angle-independent on the ECS grid -- the complex tail only affects
the discretized continuum/dissociative states, not these low-lying bound
ones. We select the `n` eigenpairs with the smallest |Im(E)| (the
bound-state signature) among the lowest-Re(E) eigenvalues, ordered ascending
in Re(E).

Normalization convention: `qscat.dvr.eigen()` returns eigenvectors with
numpy's Hermitian normalization (`v^dagger v = 1`); c-product
renormalization (`v^T v = 1`, the convention appropriate for ECS-basis
observables) is a *separate* step that `eigen()`'s own docstring flags
callers must apply themselves when it matters. For real bound vibrational
eigenvectors built here, `chi` is real, so `chi_v^dagger chi_v` and
`chi_v^T chi_v` are literally the same sum-of-squares -- they coincide to
machine precision in the Hermitian-normalized output `eigen()` already hands
back, so no extra c-product step is needed and `chi` is usable directly
here, matching the DVR/ECS convention used throughout this repo: the DVR
basis is pre-normalized by `1/sqrt(w)`, so the inner product is a plain dot
of coefficient vectors. This coincidence does NOT generalize: downstream code
handling genuinely complex vectors (e.g. resonance eigenvectors, or the
doorway/driven-equation solution) must apply the c-product convention
explicitly per `eigen()`'s own note, since Hermitian and c-product norms
differ once the vector is complex. Downstream, `chi_v(R_j)` is this raw
coefficient `chi[v, j]` -- no extra `sqrt(weight)` factor is needed because
the basis functions themselves already absorb it.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import NamedTuple

import numpy as np
import numpy.typing as npt

from qscat.dvr import FemDvrEcsGrid, eigen, kinetic
from qscat.exceptions import GridError

__all__ = ["VibrationalBasis", "vibrational_states"]


class VibrationalBasis(NamedTuple):
    """The vibrational eigenbasis returned by `vibrational_states`.

    A `NamedTuple`, so the historical `eps, chi = vibrational_states(...)`
    unpacking keeps working unchanged, while `.eps` / `.chi` give named access.
    """

    eps: npt.NDArray[np.float64]  # (n,) real vibrational energies (Ha), ascending
    chi: npt.NDArray[np.complex128]  # (n, grid.n) eigenvectors, one row per level

# Bound-state signature: true bound levels have |Im(E)| ~ 1e-15 on this ECS
# grid, while the discretized continuum/dissociative states jump to
# |Im(E)| ~ 1e-7 or larger. This tolerance sits comfortably between the two.
_IM_TOL_HA = 1e-6


def vibrational_states(
    grid: FemDvrEcsGrid,
    mu: float,
    n: int,
    v0: Callable[[npt.ArrayLike], npt.NDArray[np.complexfloating]],
) -> VibrationalBasis:
    """The `n` lowest bound eigenpairs of `T_nuc(mu) + diag(v0(R))`.

    `v0` is the neutral molecule's potential-energy curve, evaluated pointwise
    on `grid.points` -- passed in by the caller (e.g. `qscat.model.N2.v0`)
    rather than hardcoded, so this solver stays model-independent. The
    eigenvectors are Hermitian-normalized (``v^dagger v = 1``) by
    `qscat.dvr.eigen()`, which coincides with the c-product norm
    (``v^T v = 1``) for these real bound-state vectors.

    Parameters
    ----------
    grid : FemDvrEcsGrid
        The nuclear radial grid.
    mu : float
        Nuclear reduced mass (atomic units).
    n : int
        Number of lowest bound vibrational states to return.
    v0 : callable
        The neutral potential-energy curve ``v0(R)``, evaluated pointwise on
        ``grid.points``.

    Returns
    -------
    VibrationalBasis
        A ``(eps, chi)`` named tuple: ``eps`` the real part of the ``n``
        lowest-Re(E) eigenvalues (Hartree), ascending; ``chi`` the
        ``(n, grid.n)`` array of eigenvectors, one row per level.

    Raises
    ------
    GridError
        If any of the ``n`` lowest-Re(E) eigenvalues has
        ``|Im(E)| > _IM_TOL_HA`` -- i.e. ``n`` reached past the true bound
        spectrum into quasi-continuum/ECS states, which are not valid
        vibrational levels.
    """
    T = kinetic(grid, mu)
    H0 = T + np.diag(v0(grid.points))
    E, V = eigen(H0)  # already sorted ascending by Re(E)

    E_n = E[:n]
    bad = np.abs(E_n.imag) > _IM_TOL_HA
    if np.any(bad):
        bad_idx = np.flatnonzero(bad).tolist()
        raise GridError(
            f"vibrational_states(n={n}) requested more states than there "
            f"are bound levels: eigenvalue index/indices {bad_idx} (within "
            f"the n={n} lowest-Re(E) selection) have |Im(E)| > "
            f"{_IM_TOL_HA} Ha, i.e. they are quasi-continuum/ECS states, "
            "not true bound vibrational levels. Reduce n to stay within "
            "the true bound spectrum."
        )

    eps = E_n.real
    chi = V[:, :n].T
    return VibrationalBasis(eps, chi)
