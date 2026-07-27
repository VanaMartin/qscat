"""Exact-2D TI VE cross-section oracle for the diatomic molecules (NO, F2, ...).

The primary deliverable of the model port: for each molecule, the exact-2D
time-independent σ_{0→v'}(E) computed through the promoted library
(`qscat.core` + `qscat.model`) -- the ORACLE the LCP/TD approximations are
tested against (there is no independent golden data for NO/F2, unlike N2's
Houfek `CSVE.V00.J00`; the exact solver IS the reference, per the research
program). This module builds the grid from `config.MoleculeConfig`, computes
the curve via `qscat.core.driven.ve_cross_section`, and (`main`) saves the
figure + `.npz`.

Run: `uv run python -m validation.diatomic.curves` (both molecules; a few
minutes each on the working grid).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import numpy.typing as npt
from qscat.core.driven import ve_cross_section
from qscat.core.grids import electronic_grid, nuclear_grid
from qscat.core.plot import plot_cross_sections
from qscat.core.vibrational import vibrational_states
from qscat.dvr import TensorGrid

from validation.diatomic.config import CONFIGS, MoleculeConfig

__all__ = ["build_grid", "compute_ti_curve", "FIGURE_DIR", "main"]

FIGURE_DIR = Path(__file__).resolve().parents[2] / "docs" / "physics" / "figures"
VPRIMES = [0, 1, 2]


def build_grid(cfg: MoleculeConfig) -> TensorGrid:
    """The (electronic × nuclear) FEM-DVR-ECS `TensorGrid` for `cfg`."""
    return TensorGrid(
        [
            electronic_grid(r_max=cfg.e_r_max, order=cfg.e_order, n_complex=cfg.e_n_complex),
            nuclear_grid(r_max=cfg.n_r_max, quadrature=cfg.n_quadrature, n_complex=cfg.n_n_complex),
        ]
    )


def compute_ti_curve(
    cfg: MoleculeConfig, E_grid: npt.NDArray[np.float64] | None = None
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """`(E_grid, sigma[E, v'], eps)` — the exact-2D TI σ_{0→v'}(E) for `cfg`.

    `sigma` is bohr², channels `VPRIMES = [0, 1, 2]` (elastic + first two
    excitations); `eps` are the neutral vibrational energies. `E_grid` defaults
    to `cfg.energy_grid()`.
    """
    tgrid = build_grid(cfg)
    eps, chi = vibrational_states(tgrid.grids[1], cfg.model.mu, cfg.n_vib, cfg.model.v0)
    E = cfg.energy_grid() if E_grid is None else np.asarray(E_grid, dtype=np.float64)
    sigma = ve_cross_section(tgrid, cfg.model, eps, chi, 0, VPRIMES, E)
    return E, np.asarray(sigma, dtype=np.float64), eps


def main() -> None:
    """Compute + save the committed σ(E) figure and `.npz` for every molecule."""
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    for name, cfg in CONFIGS.items():
        E, sigma, eps = compute_ti_curve(cfg)
        stem = f"{name.lower()}-2d-ti-cross-section"
        np.savez(FIGURE_DIR / f"{stem}.npz", E=E, sigma=sigma, eps=eps, channels=np.array(VPRIMES))
        plot_cross_sections(
            E,
            sigma,
            channels=VPRIMES,
            title=f"{name} exact-2D VE cross section sigma_{{0->v'}}(E) (qscat.core)",
            path=FIGURE_DIR / f"{stem}.png",
        )
        print(
            f"{name}: wrote {stem}.png/.npz; vib spacing eps[1]-eps[0] = {eps[1] - eps[0]:.5f} Ha"
        )


if __name__ == "__main__":
    main()
