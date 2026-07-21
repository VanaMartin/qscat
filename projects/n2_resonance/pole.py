"""Two-angle ECS resonance pole finder for the N2 electronic problem
(sub-project #2, Task 2 -- the crux).

Builds the fixed-R electronic Hamiltonian H_el(R) = T + diag(v_eff_el(r, R))
on a `qscat.dvr.FemDvrEcsGrid` and locates the resonance pole by two-angle
matching: the eigenvalue of the ECS-rotated spectrum that is (nearly)
independent of the ECS rotation angle is the physical resonance pole; the
discretized continuum eigenvalues rotate with the angle and do not match
between two grids built at different angles.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from potential import v_eff_el
from qscat.dvr import FemDvrEcsGrid, eigen, hamiltonian

__all__ = ["electronic_hamiltonian", "find_pole"]


def electronic_hamiltonian(R: float, grid: FemDvrEcsGrid) -> npt.NDArray[np.complex128]:
    """Assemble H_el(R) = kinetic(grid, 1.0) + diag(v_eff_el(grid.points, R))."""
    return hamiltonian(grid, lambda z: v_eff_el(z, R), 1.0)


def _filter_window(
    E: npt.NDArray[np.complex128], window: tuple[float, float, float, float]
) -> npt.NDArray[np.complex128]:
    re_lo, re_hi, im_lo, im_hi = window
    mask = (
        (E.real >= re_lo)
        & (E.real <= re_hi)
        & (E.imag >= im_lo)
        & (E.imag <= im_hi)
    )
    return E[mask]


def find_pole(
    R: float,
    grid_a: FemDvrEcsGrid,
    grid_b: FemDvrEcsGrid,
    window: tuple[float, float, float, float],
) -> tuple[complex, float]:
    """Locate the resonance pole common to two different-angle ECS grids.

    Diagonalizes H_el(R) on `grid_a` and `grid_b`, restricts each spectrum to
    the search `window = (re_lo, re_hi, im_lo, im_hi)`, and returns the
    (ea, eb) pair -- one eigenvalue from each grid -- with the smallest
    |ea - eb|. That pair is the angle-stable pole; discretized-continuum
    eigenvalues rotate with the ECS angle and do not match this closely.

    Returns `(E_pole, residual)` where `E_pole = 0.5*(ea+eb)` and
    `residual = |ea-eb|`.
    """
    Ea, _ = eigen(electronic_hamiltonian(R, grid_a))
    Eb, _ = eigen(electronic_hamiltonian(R, grid_b))

    fa = _filter_window(Ea, window)
    fb = _filter_window(Eb, window)

    if fa.size == 0 or fb.size == 0:
        raise ValueError(
            f"find_pole: window {window} contains no eigenvalues on "
            f"{'grid_a' if fa.size == 0 else 'grid_b'} "
            f"(found {fa.size} in A, {fb.size} in B) -- window too tight "
            "or grid too coarse."
        )

    # For each candidate in A, distance to nearest candidate in B.
    diffs = np.abs(fa[:, None] - fb[None, :])
    i, j = np.unravel_index(np.argmin(diffs), diffs.shape)
    ea, eb = fa[i], fb[j]
    residual = float(np.abs(ea - eb))
    E_pole = complex(0.5 * (ea + eb))
    return E_pole, residual
