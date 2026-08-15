"""Reproducible H2+ dissociative-recombination sigma_DR(E) curves + figures.

This is the generator behind the committed DR figure(s) under
`docs/physics/figures/`. It runs the FULL 1300-bohr H2+ deck
(`validation.h2plus.config.full_grid`, ~1.15M unknowns) under MUMPS, so -- like
`validation.h2plus.dr.main` -- it is **Docker/MUMPS-only, NOT part of the test
suite** (that grid is far too heavy for a laptop; the laptop-feasible gate is
`validation/h2plus/test_dr.py` on the small proxy grid). The exact-2D TI solver
IS the oracle here: no independent golden DR data ships (eMoScat's
`output/H2+/sigma.txt` is absent from the snapshot).

Run inside the `qmodeling-base`/`test` container (system MUMPS present) with the
repo mounted, e.g.::

    uv run python -m validation.h2plus.dr_curves

Two curves are produced:

* **full-range** (`full_range_curve`): the coarse `config.energy_grid()`
  (0.001..0.050 Ha, step 0.001), first `N_CHANNELS` Rydberg exit channels ->
  `h2plus-dr-cross-section.png` (linear, one line per channel).
* **short-range** (`short_range_curve`): 200 log-spaced energies across the
  DR1 resonance in [0.005, 0.007] Ha -> `h2plus-dr-cross-section-shortrange.png`
  (log-log DR1 + DR2), the committed accuracy figure.

Each curve also writes a sidecar `.npz` (and the short-range a `.csv`) next to
its figure so the numbers are recoverable without re-solving.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import numpy.typing as npt
from qscat.linalg import default_backend

from validation.h2plus.config import N_CHANNELS, energy_grid, full_grid
from validation.h2plus.dr import compute_dr

__all__ = ["FIGURE_DIR", "full_range_curve", "short_range_curve", "main"]

FIGURE_DIR = Path(__file__).resolve().parents[2] / "docs" / "physics" / "figures"


def _sweep(
    energies: npt.NDArray[np.float64],
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """`(energies, sigma)` on the full deck under MUMPS.

    Forces the ``"mumps"`` backend so ``"auto"`` cannot silently fall back to
    SuperLU (far larger RSS/time on the 1.15M complex-symmetric factor).
    """
    grid = full_grid()
    with default_backend("mumps"):
        return compute_dr(grid, energies=energies, n_channels=N_CHANNELS)


def full_range_curve() -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Coarse `config.energy_grid()` sweep -> linear multi-channel figure."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    energies = energy_grid()
    energies, sigma = _sweep(energies)

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    np.savez(FIGURE_DIR / "h2plus-dr-cross-section.npz", energies=energies, sigma=sigma)

    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    for n in range(sigma.shape[1]):
        ax.plot(energies, sigma[:, n], marker="o", ms=3.0, lw=1.3, label=f"Rydberg channel n={n}")
    ax.set_xlabel(r"collision energy $E$ (Ha)")
    ax.set_ylabel(r"$\sigma_{\mathrm{DR}}$ (bohr$^2$)")
    ax.set_title(
        rf"H$_2^+$ TI dissociative recombination — first {sigma.shape[1]} Rydberg channels"
    )
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "h2plus-dr-cross-section.png", dpi=130)
    plt.close(fig)
    print("wrote h2plus-dr-cross-section.png/.npz")
    return energies, sigma


def short_range_curve() -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """200 log-spaced energies across the DR1 resonance -> log-log figure.

    Reproduces the committed `h2plus-dr-cross-section-shortrange.png`: DR1 (n=0)
    and DR2 (n=1) on log-log axes (positive values only -- log axes cannot show
    closed-channel zeros).
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    energies = np.logspace(np.log10(0.005), np.log10(0.007), 200)
    energies, sigma = _sweep(energies)

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    stem = "h2plus-dr-cross-section-shortrange"
    header = "E_Ha," + ",".join(f"sigma_DR{c + 1}_n{c}" for c in range(sigma.shape[1]))
    np.savetxt(
        FIGURE_DIR / f"{stem}.csv",
        np.column_stack([energies, sigma]),
        delimiter=",",
        header=header,
        comments="",
    )
    np.savez(FIGURE_DIR / f"{stem}.npz", energies=energies, sigma=sigma)

    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    series = [(0, "DR1 (n=0)", "C0"), (1, "DR2 (n=1)", "C1")]
    for col, label, c in series:
        if col >= sigma.shape[1]:
            continue
        y = sigma[:, col]
        m = y > 0  # log scale cannot show closed-channel zeros
        ax.plot(energies[m], y[m], marker="o", ms=2.6, lw=1.0, color=c, label=label)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"collision energy $E$ (Ha, log)")
    ax.set_ylabel(r"$\sigma_{\mathrm{DR}}$ (bohr$^2$, log)")
    ax.set_title(r"H$_2^+$ DR cross section — short range 0.005–0.007 Ha (log–log, 200 pts)")
    ax.grid(alpha=0.3, which="both")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / f"{stem}.png", dpi=140)
    plt.close(fig)

    i = int(np.argmax(sigma[:, 0]))
    print(f"wrote {stem}.png/.csv/.npz")
    print(f"DR1 peak at E={energies[i]:.6e} Ha, sigma={sigma[i, 0]:.4e} bohr^2")
    return energies, sigma


def main() -> None:
    full_range_curve()
    short_range_curve()


if __name__ == "__main__":
    main()
