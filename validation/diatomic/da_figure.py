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
    molecule: str, energies: npt.NDArray[np.float64]
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """`(sigma_exact, sigma_lcp)` on `molecule`'s eMoScat production deck.

    Both routes run on the SAME nuclear deck and the same electronic grid the
    committed figures used, so the comparison is differential.
    """
    cfg: MoleculeConfig = CONFIGS[molecule]
    tgrid: TensorGrid = cfg.da_grid()
    elec, nuc = tgrid.grids
    eps, chi = vibrational_states(nuc, cfg.model.mu, cfg.n_vib, cfg.model.v0)

    sigma_exact = np.asarray(
        da_cross_section(tgrid, cfg.model, eps, chi, 0, energies), dtype=np.float64
    )[:, 0]

    # The LCP needs the second ECS-angle electronic grid for its pole walk;
    # `validation.diatomic.nrm.setup` is the one place that pairing is defined.
    s = setup(molecule)
    vd, gamma = local_complex_potential(s.model, s.nuc, s.elec, s.elec_b)
    sigma_lcp = np.asarray(
        lcp_da_cross_section(s.nuc, s.model.mu, vd, gamma, s.eps, s.chi, 0, energies),
        dtype=np.float64,
    )
    assert elec.n == s.elec.n, "the LCP and exact routes must share the electronic deck"
    return sigma_exact, sigma_lcp


def _title_suffix(energies: npt.NDArray[np.float64], step: float) -> str:
    return f"({energies.size} energies, step {step:.4f} Ha)"


def write_figures(
    molecule: str,
    energies: npt.NDArray[np.float64],
    sigma_exact: npt.NDArray[np.float64],
    sigma_lcp: npt.NDArray[np.float64],
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
    lcp = np.where(sigma_lcp > 0.0, sigma_lcp, np.nan)
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
        lo, hi = float(ref_e.min()) - 0.002, float(ref_e.max()) + 0.010
        sel = (energies >= lo) & (energies <= hi)
        axz.plot(energies[sel], ex[sel] * 1e9, "-", color="tab:blue", label="exact 2-D TI")
        axz.plot(ref_e, ref_s * 1e9, "o--", color="k", markersize=5, linewidth=1.0, label=ref_label)
        axz.set_xlim(lo, hi)
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
    ax.plot(energies, lcp, "--", color="tab:red", label="LCP approximation")
    if ref_e.size:
        ax.plot(ref_e, ref_s, "o--", color="k", markersize=5, linewidth=1.0, label=ref_label)
    ax.set_ylabel(r"$\sigma_{DA}$ (bohr$^2$)")
    ax.set_yscale("log")
    ax.set_title(rf"{molecule} $\sigma_{{DA}}$: LCP vs exact-2D oracle  {suffix}")
    ax.legend(fontsize="small")
    ax.grid(True, which="both", alpha=0.2)

    ratio = np.where(np.isfinite(ex) & (ex > 0.0), lcp / ex, np.nan)
    axr.plot(energies, ratio, "-", color="tab:green")
    axr.axhline(1.0, color="k", linestyle=":", linewidth=1.0)
    axr.set_yscale("log")
    axr.set_xlabel("E (Hartree)")
    axr.set_ylabel("LCP / exact")
    axr.grid(True, which="both", alpha=0.2)
    fig.tight_layout()
    fig.savefig(stem.with_suffix(".png"), dpi=120)
    plt.close(fig)
    np.savez(stem.with_suffix(".npz"), E=energies, sigma_exact=sigma_exact, sigma_lcp=sigma_lcp)
    written += [stem.with_suffix(".png"), stem.with_suffix(".npz")]
    return written


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("molecule", choices=sorted(ENERGY_GRIDS), help="which molecule to plot")
    p.add_argument("--outdir", type=Path, default=_FIGURE_DIR)
    args = p.parse_args()

    energies = energies_for(args.molecule)
    sigma_exact, sigma_lcp = compute(args.molecule, energies)
    open_ = sigma_exact > 0.0
    print(
        f"{args.molecule}: {energies.size} energies, {int(open_.sum())} above the DA threshold; "
        f"sigma_exact in [{sigma_exact[open_].min():.4e}, {sigma_exact[open_].max():.4e}], "
        f"LCP/exact in [{(sigma_lcp[open_] / sigma_exact[open_]).min():.4g}, "
        f"{(sigma_lcp[open_] / sigma_exact[open_]).max():.4g}]"
    )
    for path in write_figures(args.molecule, energies, sigma_exact, sigma_lcp, args.outdir):
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
