"""Neutral N2 vibrational states on the nuclear FEM-DVR-ECS grid
(sub-project #3, Task 1).

`.superpowers/sdd/ti-cross-section-extraction.md` section 2: diagonalize
`T_nuc(mu) + diag(V0(R))` -> (eps_v, chi_v). The neutral N2 potential well is
entirely inside the real region of the grid (bound states, R ~ 1.5-3 bohr),
so its eigenvalues are real and angle-independent on the ECS grid -- the
complex tail only affects the discretized continuum/dissociative states, not
these low-lying bound ones. We select the `n` eigenpairs with the smallest
|Im(E)| (the bound-state signature) among the lowest-Re(E) eigenvalues,
ordered ascending in Re(E).

Normalization convention: `qscat.dvr.eigen()` returns eigenvectors with
numpy's Hermitian normalization (`v^dagger v = 1`); c-product
renormalization (`v^T v = 1`, the convention appropriate for ECS-basis
observables) is a *separate* step that `eigen()`'s own docstring flags
callers must apply themselves when it matters. For the real bound
vibrational eigenvectors built here, `chi` is real, so `chi_v^dagger chi_v`
and `chi_v^T chi_v` are literally the same sum-of-squares -- they coincide
to machine precision in the Hermitian-normalized output `eigen()` already
hands back, so no extra c-product step is needed and `chi` is usable
directly here, matching the DVR/ECS convention used throughout this repo
(see `.superpowers/sdd/ti-cross-section-extraction.md` section 2's
S-matrix formula "DVR basis pre-normalized by 1/sqrt(w) => inner product is
plain dot of coeff vectors"). This coincidence does NOT generalize:
downstream code handling genuinely complex vectors (e.g. resonance
eigenvectors, or the driven-equation/doorway solution in Task 3) must apply
the c-product convention explicitly per `eigen()`'s own note, since
Hermitian and c-product norms differ once the vector is complex. Downstream
(the doorway/driven-equation formulas), `chi_v(R_j)` is this raw coefficient
`chi[v, j]` -- no extra `sqrt(weight)` factor is needed because the basis
functions themselves already absorb it.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import numpy.typing as npt
from qscat.dvr import FemDvrEcsGrid, eigen, kinetic

# Reuse the already-validated N2 neutral Morse potential from sub-project #2
# (cross-import via sys.path insert, like `projects/n2_resonance/test_potential.py`
# does for `validation/n2/model.py`).
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "n2_resonance"))
from potential import v0  # noqa: E402

__all__ = ["vibrational_states"]

# Bound-state signature: true bound levels have |Im(E)| ~ 1e-15 on this ECS
# grid, while the discretized continuum/dissociative states (starting
# ~index 111 for the grid built by `nuclear_grid.n2_nuclear_grid()`) jump to
# |Im(E)| ~ 1e-7 or larger. This tolerance sits comfortably between the two.
_IM_TOL_HA = 1e-6


def vibrational_states(
    grid: FemDvrEcsGrid, mu: float, n: int
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.complex128]]:
    """The `n` lowest bound eigenpairs of `T_nuc(mu) + diag(V0(R))`.

    Returns `(eps, chi)`: `eps` is the real part of the `n` lowest-Re(E)
    eigenvalues (Hartree), ascending; `chi` is the (n, grid.n) array of the
    corresponding eigenvectors, Hermitian-normalized (`v^dagger v = 1`) by
    `qscat.dvr.eigen()` -- which coincides with the c-product norm
    (`v^T v = 1`) for these real bound-state vectors (see module docstring
    for the normalization convention).

    Raises `ValueError` if any of the `n` lowest-Re(E) eigenvalues has
    `|Im(E)| > _IM_TOL_HA`: that signals `n` reached past the true bound
    spectrum into quasi-continuum/ECS states, which are not valid
    "vibrational levels" and must not be silently returned as such.
    """
    T = kinetic(grid, mu)
    H0 = T + np.diag(v0(grid.points))
    E, V = eigen(H0)  # already sorted ascending by Re(E)

    E_n = E[:n]
    bad = np.abs(E_n.imag) > _IM_TOL_HA
    if np.any(bad):
        bad_idx = np.flatnonzero(bad).tolist()
        raise ValueError(
            f"vibrational_states(n={n}) requested more states than there "
            f"are bound levels: eigenvalue index/indices {bad_idx} (within "
            f"the n={n} lowest-Re(E) selection) have |Im(E)| > "
            f"{_IM_TOL_HA} Ha, i.e. they are quasi-continuum/ECS states, "
            "not true bound vibrational levels. Reduce n to stay within "
            "the true bound spectrum."
        )

    eps = E_n.real
    chi = V[:, :n].T
    return eps, chi
