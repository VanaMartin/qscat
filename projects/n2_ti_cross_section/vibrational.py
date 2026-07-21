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

Normalization convention (c-product-free, matches the DVR/ECS convention
used throughout this repo -- see `operators.eigen` and
`.superpowers/sdd/ti-cross-section-extraction.md` section 2's S-matrix
formula "DVR basis pre-normalized by 1/sqrt(w) => inner product is plain dot
of coeff vectors"): `chi` is returned as `numpy.linalg.eig`'s raw
eigenvectors, i.e. `chi_v^T chi_v = 1` (complex-symmetric/c-product norm,
not Hermitian `chi_v^dagger chi_v = 1`) built directly on the
weight-pre-normalized FEM-DVR basis. Downstream (the doorway/driven-equation
formulas), `chi_v(R_j)` is this raw coefficient `chi[v, j]` -- no extra
`sqrt(weight)` factor is needed because the basis functions themselves
already absorb it.
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


def vibrational_states(
    grid: FemDvrEcsGrid, mu: float, n: int
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.complex128]]:
    """The `n` lowest bound eigenpairs of `T_nuc(mu) + diag(V0(R))`.

    Returns `(eps, chi)`: `eps` is the real part of the `n` lowest-Re(E)
    eigenvalues (Hartree), ascending; `chi` is the (n, grid.n) array of the
    corresponding raw (`v^T v = 1`) eigenvectors (see module docstring for
    the normalization convention).
    """
    T = kinetic(grid, mu)
    H0 = T + np.diag(v0(grid.points))
    E, V = eigen(H0)  # already sorted ascending by Re(E)

    eps = E[:n].real
    chi = V[:, :n].T
    return eps, chi
