"""The committed exact-2D / LCP dissociative-attachment figures for NO and F2.

    uv run python -m validation.diatomic.da_figure NO
    uv run python -m validation.diatomic.da_figure F2 --outdir docs/physics/figures

Writes, per molecule, into `docs/physics/figures/`:

    {mol}-2d-ti-da-cross-section.{png,npz}   exact 2-D TI sigma_DA(E)
    {mol}-2d-da-lcp-vs-exact.{png,npz}       the same, with the LCP overlaid
                                             and an LCP/exact ratio panel

WHY THIS IS A DRIVER AND NOT A qscat-run CONFIG. The per-molecule curve drivers
were deliberately retired into `apps/qscat-run` (see CLAUDE.md), and the plain
sigma(E) curves genuinely are configs -- `apps/qscat-run/examples/figures/
no-da-dense.yaml` computes exactly the numbers below. What is not reachable
from a config is the PRESENTATION these two figures need: a second panel
carrying the LCP/exact ratio (the point of the comparison is the ratio's energy
dependence, which is invisible when two curves five orders apart share one log
axis) and, for NO, an overlay of a published reference on a DA rather than a VE
observable. `qscat_run.artifacts._write_cross_section_png` draws one panel and
keys its reference series as `ref:ve:*`. This is the same deliberate exception
`validation/diatomic/ve_nrm_figure.py` already is, for the same reason.

The energy grid and both decks match the previously committed figures exactly,
so a diff of the output isolates the one thing that changed: `da_cross_section`
now reads the outgoing dissociation flux instead of the post-form volume
T-matrix (commit `d476516`). On F2 that moves sigma_DA by 0.05%; on NO it moves
it by four to seven orders, because NO's sigma_DA is ~1e-9 bohr^2 and the volume
integral delivered it as the residue of a ~1e6-fold cancellation. See
`validation/diatomic/test_no_da_thesis.py`.

Production decks: ~78k unknowns for NO at 151 energies. Run it in Docker with
MUMPS (`docker/run.sh`), not on a laptop.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Final

import numpy as np
import numpy.typing as npt
from qscat.core.dissociation import da_cross_section
from qscat.core.lcp import lcp_da_cross_section, local_complex_potential
from qscat.core.vibrational import vibrational_states
from qscat.dvr import TensorGrid

from validation.diatomic.config import CONFIGS, MoleculeConfig
from validation.diatomic.nrm import setup

# The energy grids the committed figures use, kept byte-identical so that
# regenerating isolates the extraction change. NO starts below its +0.1719 Ha
# threshold on purpose: the leading zeros show where the channel opens.
ENERGY_GRIDS: Final[dict[str, tuple[float, float, float]]] = {
    # molecule: (min, max, step)
    "NO": (0.150, 0.300, 0.001),
    "F2": (0.010, 0.050, 0.001),
}

# Published reference overlaid on the NO panels; absent for F2, for which no
# DA cross section has ever been published for this model.
REFERENCE: Final[dict[str, tuple[str, str]]] = {
    "NO": (
        "data/vana-2017-fig3.14-no-da.dat",
        "Vana 2017 Fig. 3.14 (read off)",
    ),
}

# Electronic real-region extents the `--ladder` DIAGNOSTIC recomputes the LCP
# over. These curves are NOT PLOTTED -- see LCP_RMAX_SPREAD below for why.
#
# The exact route is edge-insensitive (it reads a boundary flux), but the LCP's
# `V_d`/`Gamma` come from an ECS resonance-pole walk, and on NO that walk does
# NOT converge in `r_max`: measured at E=0.175 it gives 1.30e-4, 9.36e-5,
# 4.30e-2, 1.04e-2, 7.16e-3, 9.33e-4 across the ladder below, non-monotone, with
# the asymptotic width Gamma(R->R_inf) shrinking as the box grows instead of
# settling.
#
# F2 is measured on the same ladder rather than assumed converged, and it fails
# DIFFERENTLY: five of the six walks agree to ~1.8% (sigma_LCP spans
# [0.48945, 1.5794] at r_max=16 against [0.49841, 1.5514] at 80), and only the
# r_max=96 walk breaks, to [0.027386, 3.5123]. So F2's documented 0.263 -> 1.734
# sweep is a real property of its shipped deck; NO has no such property to quote.
LCP_RMAX_LADDER: Final[tuple[float, ...]] = (16.0, 32.0, 48.0, 64.0, 80.0, 96.0)

# The deck the rest of the repo runs on, and the ONLY LCP curve these figures
# draw.
SHIPPED_RMAX: Final[float] = 16.0

# Worst per-energy `hi/lo` over LCP_RMAX_LADDER, measured 2026-08-24 with
# `--ladder`. Recorded rather than recomputed because the figure does not draw
# the ladder and each rung costs a pole walk; rerun with `--ladder` to refresh.
#
# The two numbers mean different things, which is why neither figure plots the
# ladder and why an aggregated band would have been worse than useless: NO's
# walk is non-monotone across the whole range, so its LCP is UNDETERMINED; F2's
# agrees to ~1.8% over 16-80 bohr and breaks only at 96, so its documented
# 0.263 -> 1.734 sweep is a real property of the shipped deck and the 45x is one
# bad rung, not an uncertainty.
LCP_RMAX_SPREAD: Final[dict[str, float]] = {"NO": 3.98e4, "F2": 45.0}

_FIGURE_DIR = Path("docs/physics/figures")


def energies_for(molecule: str) -> npt.NDArray[np.float64]:
    lo, hi, step = ENERGY_GRIDS[molecule]
    return np.round(np.arange(lo, hi + 0.5 * step, step), 10)


def reference_for(molecule: str) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], str]:
    """`(E, sigma, label)` for `molecule`'s published anchor, or empty arrays."""
    spec = REFERENCE.get(molecule)
    if spec is None:
        empty = np.empty(0, dtype=np.float64)
        return empty, empty, ""
    rel, label = spec
    table = np.loadtxt(Path(__file__).parent / rel)
    return table[:, 0].astype(np.float64), table[:, 1].astype(np.float64), label


def compute(
    molecule: str, energies: npt.NDArray[np.float64], *, ladder: bool = False
) -> tuple[npt.NDArray[np.float64], dict[float, npt.NDArray[np.float64]]]:
    """`(sigma_exact, {r_max: sigma_lcp})` on `molecule`'s eMoScat nuclear deck.

    The exact route runs once, on the shipped deck: it reads a boundary flux, so
    it is invariant under the electronic box to 4 digits over 16 -> 96 bohr
    (`validation/diatomic/test_no_da_thesis.py` gates that).

    `ladder=False` (the default) runs the LCP once, on the shipped deck -- the
    only curve the figure draws. `ladder=True` additionally runs it at every
    `LCP_RMAX_LADDER` entry, which is a DIAGNOSTIC: it answers whether the pole
    walk determines `V_d`/`Gamma` on this molecule at all. It costs one extra
    pole walk per entry and nothing in the figure depends on it, so it is opt-in.

    Every route shares the same nuclear deck and the same vibrational basis, so
    the comparison is differential.
    """
    cfg: MoleculeConfig = CONFIGS[molecule]
    tgrid: TensorGrid = cfg.da_grid()
    nuc = tgrid.grids[1]
    eps, chi = vibrational_states(nuc, cfg.model.mu, cfg.n_vib, cfg.model.v0)

    sigma_exact = np.asarray(
        da_cross_section(tgrid, cfg.model, eps, chi, 0, energies), dtype=np.float64
    )[:, 0]

    # `setup` is the one place the LCP's paired two-ECS-angle electronic decks
    # are defined; `_ANGLE_B_DEG` rides along with whatever `r_max` it is given.
    sigma_lcp: dict[float, npt.NDArray[np.float64]] = {}
    rungs = LCP_RMAX_LADDER if ladder else (SHIPPED_RMAX,)
    for r_max in rungs:
        s = setup(molecule, e_r_max=r_max)
        vd, gamma = local_complex_potential(s.model, s.nuc, s.elec, s.elec_b)
        sigma_lcp[r_max] = np.asarray(
            lcp_da_cross_section(s.nuc, s.model.mu, vd, gamma, s.eps, s.chi, 0, energies),
            dtype=np.float64,
        )
    return sigma_exact, sigma_lcp


def lcp_envelope(
    sigma_lcp: dict[float, npt.NDArray[np.float64]],
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], float]:
    """`(lo, hi, worst_spread)` over the `r_max` ladder, per energy.

    `worst_spread` is the largest `hi/lo` at any energy where the channel is
    open -- the single number that says whether the LCP curve is determined.
    """
    stack = np.vstack([sigma_lcp[r] for r in sorted(sigma_lcp)])
    lo = stack.min(axis=0)
    hi = stack.max(axis=0)
    open_ = lo > 0.0
    worst = float((hi[open_] / lo[open_]).max()) if open_.any() else float("nan")
    return lo, hi, worst


def _title_suffix(energies: npt.NDArray[np.float64], step: float) -> str:
    return f"({energies.size} energies, step {step:.4f} Ha)"


def write_figures(
    molecule: str,
    energies: npt.NDArray[np.float64],
    sigma_exact: npt.NDArray[np.float64],
    sigma_lcp: dict[float, npt.NDArray[np.float64]],
    outdir: Path,
) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    step = ENERGY_GRIDS[molecule][2]
    suffix = _title_suffix(energies, step)
    ref_e, ref_s, ref_label = reference_for(molecule)
    # A zero sigma means a closed channel, not a small one: drop it rather than
    # letting a log axis clip it to the bottom of the frame.
    ex = np.where(sigma_exact > 0.0, sigma_exact, np.nan)
    _lo, _hi, worst = lcp_envelope(sigma_lcp)
    shipped = np.where(sigma_lcp[SHIPPED_RMAX] > 0.0, sigma_lcp[SHIPPED_RMAX], np.nan)
    outdir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    # --- exact only -------------------------------------------------------
    # With a reference, this becomes two panels. The left one is the full
    # window on a log axis (sigma_DA falls 19 orders across it, so nothing else
    # renders); the right one reproduces the SOURCE FIGURE'S OWN AXES -- its
    # energy range, on a linear 1e-9 bohr^2 scale -- because a comparison
    # against a published panel is made by eye, and six points squeezed into
    # 3% of a 19-decade log axis cannot be read. Same reasoning as
    # `ve_nrm_figure.py`'s layout against PRA 77 Fig. 4.
    stem = outdir / f"{molecule.lower()}-2d-ti-da-cross-section"
    if ref_e.size:
        fig, (ax, axz) = plt.subplots(1, 2, figsize=(13, 5.5))
    else:
        fig, ax = plt.subplots(figsize=(9, 6))
        axz = None
    ax.plot(energies, ex, "-", color="tab:blue", label=r"exact 2-D TI $\sigma_{DA}$")
    if ref_e.size:
        ax.plot(ref_e, ref_s, "o--", color="k", markersize=5, linewidth=1.0, label=ref_label)
    ax.set_xlabel("E (Hartree)")
    ax.set_ylabel(r"$\sigma_{DA}$ (bohr$^2$)")
    ax.set_yscale("log")
    ax.set_title(f"{molecule} dissociative attachment, exact 2-D  {suffix}")
    ax.legend(fontsize="small")
    ax.grid(True, which="both", alpha=0.2)

    if axz is not None:
        zlo, zhi = float(ref_e.min()) - 0.002, float(ref_e.max()) + 0.010
        sel = (energies >= zlo) & (energies <= zhi)
        axz.plot(energies[sel], ex[sel] * 1e9, "-", color="tab:blue", label="exact 2-D TI")
        axz.plot(ref_e, ref_s * 1e9, "o--", color="k", markersize=5, linewidth=1.0, label=ref_label)
        axz.set_xlim(zlo, zhi)
        axz.set_ylim(0.0, 2.0)
        axz.set_xlabel("E (Hartree)")
        axz.set_ylabel(r"$\sigma_{DA}$ ($10^{-9}$ bohr$^2$)")
        axz.set_title("on the source panel's own axes")
        axz.legend(fontsize="small")
        axz.grid(True, alpha=0.2)
    fig.tight_layout()
    fig.savefig(stem.with_suffix(".png"), dpi=120)
    plt.close(fig)
    np.savez(stem.with_suffix(".npz"), E=energies, sigma=sigma_exact)
    written += [stem.with_suffix(".png"), stem.with_suffix(".npz")]

    # --- LCP vs exact, with the ratio panel -------------------------------
    stem = outdir / f"{molecule.lower()}-2d-da-lcp-vs-exact"
    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=(9, 8), sharex=True, gridspec_kw={"height_ratios": [3, 1]}
    )
    ax.plot(energies, ex, "-", color="tab:blue", label="exact 2-D TI (oracle)")
    # ONLY the shipped deck's LCP is drawn. The `r_max` ladder is COMPUTED (see
    # `lcp_rmax_spread`, reported in the caption and in
    # `docs/physics/diatomic-ve-cross-sections.md`) but deliberately NOT PLOTTED:
    # on NO the pole walk does not converge, so those curves are not alternative
    # estimates of the same quantity, they are failed computations. Drawing them
    # beside a real curve invites a reader to average them or to read a spread as
    # an uncertainty band, when what the ladder actually establishes is a yes/no
    # -- whether the method determines `V_d`/`Gamma` on this molecule at all.
    # That belongs in prose, and it is there.
    ax.plot(
        energies,
        shipped,
        "--",
        color="tab:red",
        linewidth=1.6,
        label=(
            rf"LCP, shipped deck ($r_{{max}}$ = {int(SHIPPED_RMAX)} $a_0$)"
            + (
                rf" — undetermined: {LCP_RMAX_SPREAD[molecule]:.3g}$\times$ over "
                rf"$r_{{max}}$ = {int(min(LCP_RMAX_LADDER))}-{int(max(LCP_RMAX_LADDER))} $a_0$"
                if LCP_RMAX_SPREAD.get(molecule, 1.0) > 1.05
                else ""
            )
        ),
    )
    if ref_e.size:
        ax.plot(ref_e, ref_s, "o--", color="k", markersize=5, linewidth=1.0, label=ref_label)
    ax.set_ylabel(r"$\sigma_{DA}$ (bohr$^2$)")
    ax.set_yscale("log")
    ax.set_title(rf"{molecule} $\sigma_{{DA}}$: LCP vs exact-2D oracle  {suffix}")
    ax.legend(fontsize="small")
    ax.grid(True, which="both", alpha=0.2)

    good = np.isfinite(ex) & (ex > 0.0)
    axr.plot(energies, np.where(good, shipped / ex, np.nan), "-", color="tab:green", linewidth=1.6)
    axr.axhline(1.0, color="k", linestyle=":", linewidth=1.0)
    axr.set_yscale("log")
    axr.set_xlabel("E (Hartree)")
    axr.set_ylabel("LCP / exact")
    axr.grid(True, which="both", alpha=0.2)
    fig.tight_layout()
    fig.savefig(stem.with_suffix(".png"), dpi=120)
    plt.close(fig)
    # `sigma_lcp` keeps its shipped-deck meaning for any existing reader; the
    # per-r_max columns are what the envelope on the figure is drawn from.
    columns: dict[str, npt.NDArray[np.float64]] = {
        "E": energies,
        "sigma_exact": sigma_exact,
        "sigma_lcp": sigma_lcp[SHIPPED_RMAX],
    }
    columns.update({f"sigma_lcp_rmax{int(r)}": v for r, v in sigma_lcp.items()})
    # numpy types savez's second positional as `allow_pickle`, so a **kwargs
    # splat of arrays trips the stub; same ignore `qscat_run.artifacts` uses.
    np.savez(stem.with_suffix(".npz"), **columns)  # type: ignore[arg-type]
    written += [stem.with_suffix(".png"), stem.with_suffix(".npz")]
    return written


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("molecule", choices=sorted(ENERGY_GRIDS), help="which molecule to plot")
    p.add_argument("--outdir", type=Path, default=_FIGURE_DIR)
    p.add_argument(
        "--ladder",
        action="store_true",
        help=(
            "also run the LCP at every LCP_RMAX_LADDER entry and report the spread. "
            "A diagnostic -- it says whether the pole walk determines V_d/Gamma at "
            "all -- and it is NOT drawn: on NO those walks do not converge, so they "
            "are failed computations rather than alternative estimates."
        ),
    )
    args = p.parse_args()

    energies = energies_for(args.molecule)
    sigma_exact, sigma_lcp = compute(args.molecule, energies, ladder=args.ladder)
    open_ = sigma_exact > 0.0
    _lo, _hi, worst = lcp_envelope(sigma_lcp)
    shipped = sigma_lcp[SHIPPED_RMAX]
    print(
        f"{args.molecule}: {energies.size} energies, {int(open_.sum())} above the DA threshold; "
        f"sigma_exact in [{sigma_exact[open_].min():.4e}, {sigma_exact[open_].max():.4e}]"
    )
    print(
        f"  LCP at the shipped deck (r_max={SHIPPED_RMAX:g}): sigma/exact in "
        f"[{(shipped[open_] / sigma_exact[open_]).min():.4g}, "
        f"{(shipped[open_] / sigma_exact[open_]).max():.4g}]"
    )
    print(
        f"  LCP r_max spread over {list(LCP_RMAX_LADDER)}: worst hi/lo = {worst:.4g} "
        f"({'NOT converged -- reported, not plotted' if worst > 1.05 else 'converged'})"
    )
    for r in LCP_RMAX_LADDER:
        v = sigma_lcp[r][open_]
        print(f"    r_max={r:5g}: sigma_LCP in [{v.min():.4e}, {v.max():.4e}]")
    for path in write_figures(args.molecule, energies, sigma_exact, sigma_lcp, args.outdir):
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
