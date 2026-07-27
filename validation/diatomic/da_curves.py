"""Exact-2D TI dissociative-attachment sigma_DA(E) for F2 and NO (the oracle).

Computed on eMoScat's per-molecule NUCLEAR grid (`cfg.da_grid()`) -- the fine
discretisation the fast K_R~58 dissociation wave needs. No independent DA data
exists (only N2 VE has Houfek's); the exact-2D TI solver IS the reference. N2's
DA channel is closed in the window (+0.5 Ha), so only F2 (exothermic) and NO
(~0.17 Ha) are shown. Run via `uv run python -m validation.diatomic.da_curves`.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from qscat.core.dissociation import da_cross_section
from qscat.core.plot import plot_cross_sections
from qscat.core.vibrational import vibrational_states

from validation.diatomic.config import CONFIGS, MoleculeConfig
from validation.diatomic.curves import FIGURE_DIR


def compute_da_curve(
    cfg: MoleculeConfig, E_grid: npt.NDArray[np.float64]
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """`(E_grid, sigma_DA[E, 0])` on the eMoScat per-molecule nuclear grid."""
    tg = cfg.da_grid()
    eps, chi = vibrational_states(tg.grids[1], cfg.model.mu, cfg.n_vib, cfg.model.v0)
    E = np.asarray(E_grid, dtype=np.float64)
    sigma = da_cross_section(tg, cfg.model, eps, chi, 0, E)
    return E, np.asarray(sigma, dtype=np.float64)


def main() -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    specs = {
        "F2": (np.linspace(0.01, 0.05, 13), "f2-2d-ti-da-cross-section"),
        "NO": (np.linspace(0.15, 0.30, 12), "no-2d-ti-da-cross-section"),
    }
    for name, (E, stem) in specs.items():
        _, sigma = compute_da_curve(CONFIGS[name], E)
        np.savez(FIGURE_DIR / f"{stem}.npz", E=E, sigma=sigma)
        plot_cross_sections(
            E,
            sigma,
            channels=None,  # single DA channel; title carries the meaning
            title=f"{name} exact-2D TI dissociative attachment sigma_DA(E) (qscat.core)",
            path=FIGURE_DIR / f"{stem}.png",
        )
        print(f"{name}: wrote {stem}.png/.npz")


if __name__ == "__main__":
    main()
