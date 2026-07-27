"""LCP dissociative-attachment sigma_DA(E) vs the exact-2D oracle, for F2/NO.

Sub-project B's scientific deliverable: `qscat.core.lcp.local_complex_potential`
reduces the exact 2-D (electronic r x nuclear R) resonance problem to a 1-D
local complex potential `(V_d(R), Gamma(R))`; `lcp_da_cross_section` solves
the resulting 1-D driven equation for sigma_DA(E). This module computes that
curve on the SAME fine per-molecule nuclear grid as the exact-2D oracle
(`validation.diatomic.da_curves.compute_da_curve`) and overlays the two, so
the approximation is judged against the exact solver it is meant to
approximate -- not against itself. Run via
`uv run python -m validation.diatomic.lcp_da_curves`.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from qscat.core.lcp import lcp_da_cross_section, local_complex_potential
from qscat.core.plot import plot_cross_sections
from qscat.core.vibrational import vibrational_states

from validation.diatomic.config import CONFIGS, MoleculeConfig
from validation.diatomic.curves import FIGURE_DIR


def compute_lcp_da_curve(
    cfg: MoleculeConfig, E_grid: npt.NDArray[np.float64]
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """`(E_grid, sigma_DA_LCP[E])` on the eMoScat per-molecule nuclear grid."""
    g_R = cfg.lcp_nuclear_grid()
    Vd, Gamma = local_complex_potential(cfg.model, g_R, *cfg.lcp_elec_grids())
    eps, chi = vibrational_states(g_R, cfg.model.mu, cfg.n_vib, cfg.model.v0)
    E = np.asarray(E_grid, dtype=np.float64)
    sigma = lcp_da_cross_section(g_R, cfg.model.mu, Vd, Gamma, eps, chi, 0, E)
    return E, np.asarray(sigma, dtype=np.float64)


def main() -> None:
    from validation.diatomic.da_curves import compute_da_curve

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    specs = {
        "F2": (np.linspace(0.01, 0.05, 13), "f2-2d-da-lcp-vs-exact"),
        "NO": (np.linspace(0.15, 0.30, 12), "no-2d-da-lcp-vs-exact"),
    }
    for name, (E, stem) in specs.items():
        cfg = CONFIGS[name]
        _, sigma_lcp = compute_lcp_da_curve(cfg, E)
        E_exact, sigma_exact = compute_da_curve(cfg, E)
        plot_cross_sections(
            E,
            sigma_lcp[:, None],
            channels=None,  # single DA channel; title carries the meaning
            reference=(E_exact, sigma_exact),
            title=f"{name} sigma_DA(E): LCP approximation vs exact-2D oracle (qscat.core)",
            path=FIGURE_DIR / f"{stem}.png",
        )
        print(f"{name}: wrote {stem}.png")


if __name__ == "__main__":
    main()
