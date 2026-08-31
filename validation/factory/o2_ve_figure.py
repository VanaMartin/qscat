"""The O2 VE overlay: the factory model's spin-orbit-resolved exact 2-D cross
section against the paper's own nonlocal (NRM) and local (LCP) curves --
theory against theory.

`python -m validation.factory.o2_ve_figure --so12 runs/o2-so12-ve --so32 runs/o2-so32-ve`
(`--out`, default `docs/physics/figures/o2-2d-ti-ve-spin-orbit-vs-alt-houfek.png`)

Reads the two components' `cross_section.csv` written by `qscat-run` for
`apps/qscat-run/examples/o2-so{12,32}-ve.yaml` (columns `energy`,
`ti:ve:v0->0` ...), sums them with the statistical factor 1/3 each (Alt &
Houfek, PRA 103, 032829 (2021), p. 032829-4 -- the same composition as the
paper's own curves, whose every peak is a doublet), and overlays each
0 -> v' panel on the vector-extracted Fig. 5 curves
(`validation/factory/data/o2/fig5_ve_0{v}_{nrm,lcp}.csv`). It also reports
the separation of the 0 -> 1 doublet near 1.04 eV, Allan's v' = 9 pair.

Read it as two questions, not one: do the exact 2-D doublets of the fitted
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

__all__ = ["G_STAT", "figure", "load_run", "main"]

G_STAT = 2.0 / 3.0  # p. 032829-4
DATA = Path(__file__).parent / "data" / "o2"
FIGURE = Path("docs/physics/figures/o2-2d-ti-ve-spin-orbit-vs-alt-houfek.png")
N_PANELS = 6


def load_run(run_dir: Path) -> tuple[np.ndarray, dict[int, np.ndarray]]:
    """`(E_eV, {v': sigma_a0^2 * g})` from `run_dir/cross_section.csv`.

    The sweeps behind this figure are ~390 kB each and are published rather
    than committed, so a fresh clone has the pointer but not the numbers.
    Run `qscat-run fetch <run_dir>` first; the error below says so rather
    than reporting a missing file, because the file is not missing so much
    as not yet downloaded.
    """
    csv_path = run_dir / "cross_section.csv"
    if not csv_path.is_file() and (run_dir / "artifacts.json").is_file():
        raise FileNotFoundError(
            f"{csv_path} is published, not committed. Download it with:\n"
            f"    qscat-run fetch {run_dir}"
        )
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


def load_split_run(run_so12: Path, run_so32: Path) -> tuple[np.ndarray, dict[int, np.ndarray]]:
    """The spin-orbit-resolved cross section: `1/3 sigma(2Pi_1/2) + 1/3
    sigma(2Pi_3/2)` (p. 032829-4), each component on its own energy mesh
    (its levels sit +-Delta_SO/2 apart), summed on the union of the meshes
    by linear interpolation -- fine against 121 points per peak."""
    E1, s1 = load_run(run_so12)
    E2, s2 = load_run(run_so32)
    E = np.unique(np.concatenate([E1, E2]))
    out = {}
    for v in sorted(set(s1) & set(s2)):
        # load_run applied g = 2/3 to each; a component carries 1/3
        out[v] = 0.5 * (np.interp(E, E1, s1[v]) + np.interp(E, E2, s2[v]))
    return E, out


def doublet_separation(E: np.ndarray, sigma: np.ndarray, e0: float, half: float = 0.025) -> float:
    """Distance (eV) between the two tallest local maxima within `e0 +- half`."""
    m = (E > e0 - half) & (E < e0 + half)
    Ew, sw = E[m], sigma[m]
    peaks = [i for i in range(1, sw.size - 1) if sw[i] > sw[i - 1] and sw[i] > sw[i + 1]]
    top = sorted(peaks, key=lambda i: -sw[i])[:2]
    if len(top) < 2:
        return float("nan")
    return float(abs(Ew[top[0]] - Ew[top[1]]))


def figure(run_so12: Path, run_so32: Path, out: Path = FIGURE) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    E_so, sig_so = load_split_run(run_so12, run_so32)
    fig, axes = plt.subplots(3, 2, figsize=(11, 12))
    for v, ax in enumerate(axes.ravel()[:N_PANELS]):
        for name, colour, ls in (("nrm", "tab:blue", "-"), ("lcp", "tab:green", "--")):
            d = np.loadtxt(DATA / f"fig5_ve_0{v}_{name}.csv", delimiter=",")
            ax.plot(d[:, 0], d[:, 1], ls, color=colour, lw=0.9, label=f"A&H {name.upper()}")
        if v in sig_so:
            ax.plot(
                E_so,
                sig_so[v],
                "-",
                color="tab:red",
                lw=0.9,
                label="exact 2-D, factory, $\\frac{1}{3}\\Pi_{1/2}+\\frac{1}{3}\\Pi_{3/2}$",
            )
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
    ap.add_argument("--so12", type=Path, required=True, help="qscat-run output dir, O2_SO12")
    ap.add_argument("--so32", type=Path, required=True, help="qscat-run output dir, O2_SO32")
    ap.add_argument("--out", type=Path, default=FIGURE)
    a = ap.parse_args()
    print(f"[O2] wrote {figure(a.so12, a.so32, a.out)}")
    E, s = load_split_run(a.so12, a.so32)
    # Allan's v'=9 doublet in 0->1 sits at ~1.05 eV (p. 032829-5); the
    # paper's model gives 17.8 meV, the measurement 19.6 +- 1.0.
    sep = doublet_separation(E, s[1], 1.037)
    print(f"[O2] 0->1 doublet separation near 1.04 eV: {sep * 1e3:.1f} meV")


if __name__ == "__main__":
    main()
