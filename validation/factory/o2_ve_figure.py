"""The O2 VE overlay: the factory model's exact 2-D cross sections against the
paper's own nonlocal (NRM) and local (LCP) curves -- theory against theory.

`python -m validation.factory.o2_ve_figure --run runs/o2-ve [--out PNG]`
(default `docs/physics/figures/o2-2d-ti-ve-vs-alt-houfek.png`)

Reads `cross_section.csv` written by `qscat-run` for `apps/qscat-run/examples/
o2-ve.yaml` (columns `energy`, `ti:ve:v0->0` ...), applies the statistical
factor `g(2Pi_g) = 2/3` of Alt & Houfek, PRA 103, 032829 (2021), p. 032829-4
(the model is one electronic symmetry; the paper's total sums two spin-orbit
components of 1/3 each, or takes 2/3 without the splitting), and overlays
each 0 -> v' panel on the vector-extracted Fig. 5 curves
(`validation/factory/data/o2/fig5_ve_0{v}_{nrm,lcp}.csv`).

Read it as two questions, not one: do the exact 2-D peaks of the fitted
potential land where the paper's nonlocal comb has them (positions are the
fit, to the spectral check's +-7 meV), and does the paper's LCP fail the
same way this repository's own LCP does. Nothing here is experiment: the
paper's measured traces (its Figs. 7-9) are not extracted.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from qscat.units import HARTREE_TO_EV

__all__ = ["G_STAT", "load_run", "figure", "main"]

G_STAT = 2.0 / 3.0  # p. 032829-4
DATA = Path(__file__).parent / "data" / "o2"
FIGURE = Path("docs/physics/figures/o2-2d-ti-ve-vs-alt-houfek.png")
N_PANELS = 6


def load_run(run_dir: Path) -> tuple[np.ndarray, dict[int, np.ndarray]]:
    """`(E_eV, {v': sigma_a0^2 * g})` from `run_dir/cross_section.csv`."""
    with (run_dir / "cross_section.csv").open() as f:
        rows = list(csv.reader(f))
    head, body = rows[0], np.array(rows[1:], dtype=float)
    E = body[:, 0] * HARTREE_TO_EV
    out: dict[int, np.ndarray] = {}
    for j, key in enumerate(head[1:], start=1):
        if key.startswith("ti:ve:v0->"):
            out[int(key.split("->")[1])] = G_STAT * body[:, j]
    if not out:
        raise ValueError(f"no `ti:ve:v0->v'` columns in {run_dir / 'cross_section.csv'}")
    return E, out


def figure(run_dir: Path, out: Path = FIGURE) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    E, sig = load_run(run_dir)
    fig, axes = plt.subplots(3, 2, figsize=(11, 12))
    for v, ax in enumerate(axes.ravel()[:N_PANELS]):
        for name, colour, ls in (("nrm", "tab:blue", "-"), ("lcp", "tab:green", "--")):
            d = np.loadtxt(DATA / f"fig5_ve_0{v}_{name}.csv", delimiter=",")
            ax.plot(d[:, 0], d[:, 1], ls, color=colour, lw=0.9, label=f"A&H {name.upper()}")
        unsplit = DATA / f"fig7_ve_0{v}_nrm.csv"
        if unsplit.exists():  # Fig. 7: the NRM without spin-orbit splitting
            d = np.loadtxt(unsplit, delimiter=",")
            ax.plot(d[:, 0], d[:, 1], ":", color="k", lw=0.9, label="A&H NRM, no spin-orbit")
        if v in sig:
            ax.plot(E, sig[v], "-", color="tab:red", lw=0.9, label="exact 2-D, factory ($g=2/3$)")
        nrm = np.loadtxt(DATA / f"fig5_ve_0{v}_nrm.csv", delimiter=",")
        ax.set(xlim=(nrm[0, 0], nrm[-1, 0]), xlabel="E (eV)", ylabel="$\\sigma$ ($a_0^2$)")
        ax.set_title(f"VE 0 $\\to$ {v}")
        ax.legend(fontsize=7)
    fig.suptitle(
        "O$_2$ vibrational excitation: the factory model's exact 2-D cross section "
        "vs the paper's NRM and LCP (Fig. 5)",
        fontsize=10,
    )
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", type=Path, required=True, help="qscat-run output directory")
    ap.add_argument("--out", type=Path, default=FIGURE)
    a = ap.parse_args()
    print(f"[O2] wrote {figure(a.run, a.out)}")


if __name__ == "__main__":
    main()
