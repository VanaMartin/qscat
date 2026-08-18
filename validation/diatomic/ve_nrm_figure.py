"""The committed N2 vibrational-excitation figure: exact vs LCP vs the NRM.

    uv run python -m validation.diatomic.ve_nrm_figure

Renders `docs/physics/figures/n2-ve-nrm-vs-exact.png` (and the `.npz` of the
curves behind it) from `validation.diatomic.ve_nrm.compare`, framed so it can
be laid directly next to PRA 77's own panels:

* **rows** are the transitions Fig. 4 plots for the N2-like model, `0->0` on
  top and `0->1` below;
* **columns** are the two discrete-state choices -- choice A (the R-dependent
  "physical" state) is what **Fig. 4** shows, choice B (the R-independent
  asymptotic state) is what **Fig. 8** shows;
* the **four curves per panel** are the four Fig. 4's caption defines (exact
  2-D solid, LCP short-dashed, nonlocal without background long-dashed,
  nonlocal with background crosses), in those line STYLES -- the colours are
  NOT the paper's (it draws the LCP blue, the no-background curve green and
  the crosses red; `_STYLES` below is red/blue/green), so pair the curves by
  dash pattern and marker, not by colour;
* the axes are **linear**, `E` in hartree over 0.06-0.16 and `sigma` in
  `a0^2` with Fig. 4's own ticks (0...50 for `0->0`, 0...14 for `0->1`)
  -- the scale on which the paper's "practically the same" is asserted
  (p. 012710-8, Fig. 4; p. 012710-10, Fig. 8 caption).

Fig. 8 has no N2 `0->1` panel, so the bottom-right panel is this repo's
result with no published counterpart to check it against -- read it against
the exact 2-D curve in the same panel, which is the oracle throughout.

The energy step is 0.001 Ha, 101 points across the window: N2's boomerang
oscillations are ~0.010 Ha apart here, so that is ten points per oscillation
-- at 0.002 the peaks are visibly clipped and angular. Measured cost is 10.5 s
of fixed setup (grids, the LCP pole walk, both ingredient sets) plus 8.0 s per
energy for all six curves, so ~14 min on a 12-core laptop (SuperLU,
`OMP_NUM_THREADS=8`); the dominant term is the exact 2-D driven solve, one
factorization per energy. Every energy is solved independently, so splitting
the window into chunks and concatenating gives the same curves as one call --
which is how the committed run was executed.

`validation/` may import `qscat` and `projects`; the reverse is forbidden.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import numpy.typing as npt
from qscat.core.plot import ComparisonPanel, plot_route_comparison

from .ve_nrm import VeComparison, compare

__all__ = ["ENERGIES", "FIGURE_PATH", "N_STATES", "VPRIMES", "panels", "render"]

# PRA 77 Fig. 4's own x range for the N2-like model.
ENERGIES = np.round(np.arange(0.060, 0.160 + 1e-9, 0.001), 6)

VPRIMES = [0, 1]

# The measured Eq. (60) state-sum truncation (see test_ve_nrm.py's ladders).
N_STATES = 100

_REPO_ROOT = Path(__file__).resolve().parents[2]
FIGURE_PATH = _REPO_ROOT / "docs" / "physics" / "figures" / "n2-ve-nrm-vs-exact.png"

# Fig. 4's four curves, in Fig. 4's own line STYLES (p. 012710-8, Fig. 4
# caption) -- but not its colours; see the module docstring.
_STYLES: dict[str, dict[str, object]] = {
    "exact 2-D": {"color": "black", "linestyle": "-", "linewidth": 1.6},
    "LCP": {"color": "tab:red", "linestyle": (0, (3, 2)), "linewidth": 1.3},
    "nonlocal, no background": {
        "color": "tab:blue",
        "linestyle": (0, (7, 3)),
        "linewidth": 1.3,
    },
    "nonlocal + background": {
        "color": "tab:green",
        "linestyle": "none",
        "marker": "x",
        "markersize": 4.0,
        "markeredgewidth": 1.0,
    },
}

# Fig. 4's y ticks per transition.
_YTICKS = {0: np.arange(0.0, 50.0 + 1e-9, 10.0), 1: np.arange(0.0, 14.0 + 1e-9, 2.0)}


def panels(c: VeComparison) -> list[list[ComparisonPanel]]:
    """`c` as the 2x2 panel grid described in the module docstring.

    Rows are `VPRIMES` (the transitions), columns are discrete-state choice A
    then B. Every panel carries the same four series labels, so `_STYLES`
    gives one route one appearance across the whole figure.
    """
    rows: list[list[ComparisonPanel]] = []
    for j, vp in enumerate(c.vprimes):
        row: list[ComparisonPanel] = []
        for choice, sigma_bg, sigma_nobg in (
            ("A (R-dependent, Fig. 4)", c.sigma_nrm_a, c.sigma_nrm_a_nobg),
            ("B (R-independent, Fig. 8)", c.sigma_nrm_b, c.sigma_nrm_b_nobg),
        ):
            row.append(
                ComparisonPanel(
                    series={
                        "exact 2-D": c.sigma_exact[:, j],
                        "LCP": c.sigma_lcp[:, j],
                        "nonlocal, no background": sigma_nobg[:, j],
                        "nonlocal + background": sigma_bg[:, j],
                    },
                    title=rf"N$_2$  $0 \to {vp}$   —   discrete state {choice}",
                    ylim=(0.0, float(_YTICKS[vp][-1])),
                    yticks=_YTICKS[vp],
                )
            )
        rows.append(row)
    return rows


def render(
    energies: npt.NDArray[np.float64] | None = None,
    *,
    path: Path | None = None,
    n_states: int | None = N_STATES,
) -> VeComparison:
    """Run the comparison and write the figure (plus its `.npz`).

    Returns the `VeComparison` so a caller can assert on the numbers rather
    than on the file having appeared.
    """
    e = ENERGIES if energies is None else np.asarray(energies, dtype=np.float64)
    out = FIGURE_PATH if path is None else Path(path)
    c = compare("N2", e, VPRIMES, n_states=n_states)

    out.parent.mkdir(parents=True, exist_ok=True)
    plot_route_comparison(
        c.energies,
        panels(c),
        styles=_STYLES,
        path=out,
        xlim=(float(e[0]), float(e[-1])),
        xlabel="electron energy (hartree)",
        ylabel=r"$\sigma$ ($a_0^2$)",
        suptitle=(
            "N$_2$ vibrational excitation: exact 2-D vs LCP vs the nonlocal "
            "resonance model\n(compare PRA 77 Fig. 4, left column; Fig. 8, right column)"
        ),
    )
    np.savez(
        out.with_suffix(".npz"),
        energies=c.energies,
        vprimes=np.asarray(c.vprimes),
        n_states=np.asarray(n_states if n_states is not None else -1),
        sigma_exact=c.sigma_exact,
        sigma_lcp=c.sigma_lcp,
        sigma_nrm_a=c.sigma_nrm_a,
        sigma_nrm_a_nobg=c.sigma_nrm_a_nobg,
        sigma_nrm_b=c.sigma_nrm_b,
        sigma_nrm_b_nobg=c.sigma_nrm_b_nobg,
    )
    return c


if __name__ == "__main__":  # pragma: no cover - a driver, not a test
    comparison = render()
    print(f"wrote {FIGURE_PATH} ({comparison.energies.size} energies)")
