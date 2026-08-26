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
from qscat.dvr import FemDvrEcsGrid, eigen, hamiltonian
from qscat.ecs import find_resonance_pole

from projects.n2_resonance.potential import v0, v_eff_el

__all__ = ["electronic_hamiltonian", "find_pole", "resonance_curve"]


def electronic_hamiltonian(R: float, grid: FemDvrEcsGrid) -> npt.NDArray[np.complex128]:
    """Assemble H_el(R) = kinetic(grid, 1.0) + diag(v_eff_el(grid.points, R))."""
    return hamiltonian(grid, lambda z: v_eff_el(z, R), 1.0)


def find_pole(
    R: float,
    grid_a: FemDvrEcsGrid,
    grid_b: FemDvrEcsGrid,
    window: tuple[float, float, float, float],
) -> tuple[complex, float]:
    """Locate the resonance pole common to two different-angle ECS grids.

    Diagonalizes H_el(R) on `grid_a` and `grid_b` and delegates the
    eigenvalue-matching itself to the general `qscat.ecs.find_resonance_pole`
    (promoted from this function -- see Task 3). Returns `(E_pole,
    residual)` where `E_pole = 0.5*(ea+eb)` and `residual = |ea-eb|` for the
    closest-matching pair restricted to `window = (re_lo, re_hi, im_lo,
    im_hi)`.
    """
    Ea, _ = eigen(electronic_hamiltonian(R, grid_a))
    Eb, _ = eigen(electronic_hamiltonian(R, grid_b))
    return find_resonance_pole(Ea, Eb, window)


def resonance_curve(
    R_grid: npt.ArrayLike,
    grid_a: FemDvrEcsGrid,
    grid_b: FemDvrEcsGrid,
    *,
    window0: tuple[float, float, float, float] = (0.04, 0.16, -0.05, 0.0),
    seed_R: float = 2.01943,
    re_half_width: float = 0.05,
    im_half_width: float = 0.05,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Trace the resonance pole E_pole(R) across `R_grid` by continuation.

    The pole drifts a *lot* over a wide R range: at short R it is a genuine
    complex shape resonance (as at `R0`), but as R stretches past ~2.2-2.3
    Bohr it crosses into a real, angle-independent bound state (the familiar
    dissociative-attachment picture -- the anion curve dips below the
    neutral one). A single fixed search window can't span both regimes, so
    this walks `R_grid` by continuation instead: starting from the index
    nearest `seed_R` (the R0-equilibrium pole region the fixed `window0` --
    Re in [0.04, 0.16] Ha, see `test_pole.WINDOW` -- is seeded for), it walks
    outward in both directions (increasing then decreasing R), each step
    calling `find_pole(R, grid_a, grid_b, window)` with `window` a
    +/-`re_half_width`/+/-`im_half_width` box recentered on the *previous*
    step's matched pole. This tracks the pole smoothly through the
    resonance-to-bound-state crossing and avoids mode-hopping onto an
    unrelated (continuum or other-pole) eigenvalue; if hops are observed, fix
    via `re_half_width`/`im_half_width` (or `seed_R`/`window0`), not by
    loosening the physics window.

    Returns `(E_res, Gamma, V_d)` arrays aligned with `R_grid`:
      - `E_res = Re(E_pole)` (Hartree).
      - `Gamma = max(0, -2*Im(E_pole))` (Hartree; clipped at 0 so the real,
        angle-independent bound-state branch reports zero width rather than
        numerical noise).
      - `V_d = v0(R) + E_res` (Hartree) -- the resonance/anion curve measured
        from the neutral Morse potential `v0(R)`.
    """
    Rs = np.asarray(R_grid, dtype=np.float64)
    n = Rs.size
    E_res = np.empty(n, dtype=np.float64)
    Gamma = np.empty(n, dtype=np.float64)
    V_d = np.empty(n, dtype=np.float64)

    def _walk(indices: range) -> None:
        window = window0
        for idx in indices:
            R = Rs[idx]
            E_pole, _residual = find_pole(float(R), grid_a, grid_b, window)
            E_res[idx] = E_pole.real
            Gamma[idx] = max(0.0, -2.0 * E_pole.imag)
            # v0 is qscat.model.N2.v0, which computes in complex128 (ECS-
            # safe); float() on a 0-d complex raises TypeError even though R
            # is real here and v0(R) is real to round-off.
            V_d[idx] = float(np.real(v0(R))) + E_pole.real
            window = (
                E_pole.real - re_half_width,
                E_pole.real + re_half_width,
                E_pole.imag - im_half_width,
                E_pole.imag + im_half_width,
            )

    seed_idx = int(np.argmin(np.abs(Rs - seed_R)))
    _walk(range(seed_idx, n))  # seed_idx, seed_idx+1, ... n-1 (increasing R)
    _walk(range(seed_idx - 1, -1, -1))  # seed_idx-1, ... 0 (decreasing R)
    return E_res, Gamma, V_d
