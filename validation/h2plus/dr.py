"""H2+ dissociative-recombination (DR) driver.

`compute_dr` wraps `qscat.core.dr_cross_section` (the model-independent DR
T-matrix engine, sub-project #A Task 4) with `qscat.model.H2P` (the H2+
ionic model, sub-project D). `main()` is a Docker/MUMPS smoke run on the
full 1300-bohr deck (`validation.h2plus.config.full_grid`) -- NOT part of
the test suite (that grid is too heavy for a laptop; see
`validation/h2plus/test_dr.py` for the laptop-feasible small-proxy gate).
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from qscat.core import dr_cross_section
from qscat.core.vibrational import vibrational_states
from qscat.dvr import TensorGrid
from qscat.model import H2P

from validation.h2plus.config import N_CHANNELS, energy_grid, full_grid

__all__ = ["compute_dr"]


def compute_dr(
    grid: TensorGrid,
    *,
    energies: npt.ArrayLike,
    n_channels: int = N_CHANNELS,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """sigma_DR(E) for H2+ on `grid`.

    Diagonalizes the neutral H2's vibrational states (`v0 == 0`, the H2+
    target's own ground vibrational level) up to `max(3, n_channels + 1)`
    states -- one more than `n_channels` so `dr_cross_section`'s Rydberg
    exit-channel loop always has a properly-converged incident vibrational
    state to draw on -- then calls `qscat.core.dr_cross_section`.

    Returns `(energies, sigma)`: `energies` as passed in (as an array),
    `sigma` shape `(len(energies), n_channels)`.
    """
    eps, chi = vibrational_states(grid.grids[1], H2P.mu, max(3, n_channels + 1), H2P.v0)
    sigma = dr_cross_section(grid, H2P, eps, chi, 0, energies, n_channels=n_channels)
    return np.asarray(energies), np.asarray(sigma)


def main() -> None:
    """Docker/MUMPS smoke run: the full 1300-bohr H2+ deck, a couple of
    energies. NOT run in the test suite -- this grid is ~1.15M unknowns.
    """
    grid = full_grid()
    energies = energy_grid()[:3]
    print(f"H2+ DR smoke run: grid.size={grid.size}, energies={energies}")
    e_out, sigma = compute_dr(grid, energies=energies, n_channels=N_CHANNELS)
    out_path = "h2plus_dr_smoke.npz"
    np.savez(out_path, energies=e_out, sigma=sigma)
    print(f"sigma_DR shape={sigma.shape}")
    print(sigma)
    print(f"saved to {out_path}")


if __name__ == "__main__":
    main()
