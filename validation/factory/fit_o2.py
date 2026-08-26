"""Fit the O2 target: `python -m validation.factory.fit_o2 [--out DIR] [--n-nodes N]`.

Phase 1 of the O2 work: the curves are the vector extraction of the published
figure (precision ~0.02 eV), so this measures what the two-dimensional form can
do on a real target, not the final O2 model. Writes the FitReport as JSON and an
overlay figure of the target against the fitted model's own curves.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np

from projects.potential_factory.extract import walk_t1
from projects.potential_factory.fit import fit, model_gamma_tilde
from projects.potential_factory.report import FitReport, Tolerances
from projects.potential_factory.tracker import ElectronicPair
from validation.factory.targets.o2 import o2_seed, o2_target
from validation.factory.targets.o2_data import EV, PRECISION_HA, load_o2

GRID = {"r_max": 32.0, "order": 11, "n_complex": 8}
# Image-match tolerances: the extraction's own precision floor (~0.02 eV) on the
# curve tiers, twice that on E_res; Gamma to 20 %. Not the sensitivity budget.
IMAGE_TOL = Tolerances(
    v0_rms=PRECISION_HA,
    omega_e_rel=0.02,
    e_res_rms=2 * PRECISION_HA,
    gamma_rel=0.20,
    gamma_floor=2e-3,
    coupling_log_rms=0.3,
)


def run(
    out: Path,
    *,
    n_nodes: int = 40,
    polish_nfev: int = 400,
    grid: dict | None = None,
    R_max: float = 6.0,
    lam_coeffs: int = 9,
    alpha_coeffs: int = 3,
) -> FitReport:
    g = grid or GRID
    pair = ElectronicPair(
        angles=(35.0, 44.0),
        r_max=float(g["r_max"]),
        order=int(g["order"]),
        n_complex=int(g["n_complex"]),
    )
    target = o2_target(R_range=(1.85, R_max))
    t0 = time.perf_counter()
    model, report = fit(
        target,
        o2_seed(),
        pair=pair,
        tol=IMAGE_TOL,
        n_beta=5,
        n_nodes=n_nodes,
        lam_coeffs=lam_coeffs,
        alpha_coeffs=alpha_coeffs,
        continue_on_miss=True,  # image match: every tier reports, even after a miss
        polish_nfev=polish_nfev,
    )
    out.mkdir(parents=True, exist_ok=True)
    report.to_json(out / "o2-fit-report.json")
    for t in report.tiers:
        print(f"[O2] {t.name}: {t.status}  rms={t.rms:.3e} max={t.max:.3e}  {t.detail}", flush=True)
    print(
        f"[O2] crossing_R = {report.crossing_R}, DA sign = {report.da_threshold_sign}, "
        f"{time.perf_counter() - t0:.0f} s",
        flush=True,
    )
    _figure(model, target, pair, out)
    return report


def _figure(model, target, pair, out: Path) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    c = load_o2()
    lo, hi = target.resonance.R_range  # the model is only constrained inside the fitted range
    R_desc = np.linspace(hi, lo, 60)
    s0 = float(target.resonance.v_ion(hi) - model.v0(hi).real)
    R_ok, shift, gamma = walk_t1(model, pair, R_desc, seed_energy=complex(s0, 0.0))
    v0_fit = model.v0(R_ok).real
    fig, ax = plt.subplots(2, 2, figsize=(11, 8))
    ax[0, 0].plot(c.R, c.v0 * EV, "k.", ms=3, label="target $V_0$ (Fig. 2)")
    ax[0, 0].plot(c.R, model.v0(c.R).real * EV, "-", label="factory $V_0$")
    ax[0, 0].plot(c.R, c.v_ion * EV, "b.", ms=3, label="target $V_{ion}$")
    ax[0, 0].plot(R_ok, (v0_fit + shift) * EV, "-", color="tab:orange", label="factory $V_{ion}$")
    ax[0, 0].axvline(c.R_c, color="grey", ls=":")
    ax[0, 0].set(xlabel="R (bohr)", ylabel="E (eV)", title="O$_2$ curves: target vs factory")
    ax[0, 0].legend(fontsize=7)
    ax[0, 1].plot(c.R, (c.v_ion - c.v0) * EV, "b.", ms=3, label="target $E_{res}$")
    ax[0, 1].plot(R_ok, shift * EV, "-", color="tab:orange", label="factory $E_{res}$")
    ax[0, 1].axhline(0, color="grey", lw=0.6)
    ax[0, 1].set(xlabel="R (bohr)", ylabel="$E_{res}$ (eV)", xlim=(1.8, 3.0), ylim=(-2, 3))
    ax[0, 1].legend(fontsize=7)
    ax[1, 0].plot(
        c.R, c.gamma * EV, "k.", ms=3, label="target $\\Gamma$ (Fig. 2, $\\times 2$ halved)"
    )
    ax[1, 0].plot(R_ok, gamma * EV, "-", color="tab:orange", label="factory $\\Gamma$")
    ax[1, 0].set(xlabel="R (bohr)", ylabel="$\\Gamma$ (eV)", xlim=(1.8, 3.0))
    ax[1, 0].legend(fontsize=7)
    eps = np.geomspace(*target.coupling.eps_window, 40)
    try:  # needs a bound anion at R_inf; a fit that misses the asymptote has none
        gt = {
            R: model_gamma_tilde(model, pair, eps, np.array([R]), target.resonance.R_inf)[:, 0]
            for R in (1.9, 2.1)
        }
    except ValueError as err:
        gt = {}
        ax[1, 1].set_title(f"no factory width: {str(err)[:60]}...", fontsize=7)
    for R in (1.9, 2.1):
        (ln,) = ax[1, 1].loglog(
            eps * EV, target.coupling.gamma_tilde(eps, R) * EV, "-", label=f"Table II, R={R}"
        )
        if R in gt:
            ax[1, 1].loglog(
                eps * EV, gt[R] * EV, "--", color=ln.get_color(), label=f"factory, R={R}"
            )
    ax[1, 1].set(xlabel="$\\epsilon$ (eV)", ylabel="$\\tilde\\Gamma(\\epsilon, R)$ (eV)")
    ax[1, 1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out / "o2-factory-fit.png", dpi=130)
    plt.close(fig)
    print(f"[O2] wrote {out / 'o2-factory-fit.png'}", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=Path("runs/factory-o2"))
    ap.add_argument("--n-nodes", type=int, default=40)
    ap.add_argument("--polish-nfev", type=int, default=400)
    ap.add_argument("--r-max", type=float, default=6.0, help="upper end of the fitted R range")
    ap.add_argument("--lam-coeffs", type=int, default=9)
    ap.add_argument("--alpha-coeffs", type=int, default=3)
    a = ap.parse_args()
    run(
        a.out,
        n_nodes=a.n_nodes,
        polish_nfev=a.polish_nfev,
        R_max=a.r_max,
        lam_coeffs=a.lam_coeffs,
        alpha_coeffs=a.alpha_coeffs,
    )


if __name__ == "__main__":
    main()
