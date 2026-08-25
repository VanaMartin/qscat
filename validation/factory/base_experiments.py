"""Base experiments of the potential factory on the published N2/NO/F2 models.

Three stages, run as ``python -m validation.factory.base_experiments
--molecule N2 --stage all``:

``curves``
    The resonant curves ``E_res(R)``, ``Gamma(R)`` and ``V_ion(R) = V_0 +
    E_res`` of the published model, extracted with the factory's gated pole
    walk on a ladder of electronic grids so that convergence is SHOWN, not
    assumed. Writes one CSV per grid and a three-panel figure (the curves, the
    width, the grid-to-grid differences).
``fit``
    The round trip at production quality: the factory refits the model from
    its own converged curves (T0, T1, T3) and writes the ``FitReport``.
``xs``
    Exact 2-D cross sections on the molecule's ``emoscat`` production deck
    for the PUBLISHED model and the REFITTED model: vibrational excitation
    ``sigma_{0->v'}`` for ``v' = 0, 1, 2`` and, where the channel exists,
    dissociative attachment. The figure overlays both with a ratio panel.
    This is the observable-level answer to "what does curve-level agreement
    buy?" -- nothing here is compared with experiment.

Heavy stages (``xs``) are meant for the MUMPS-equipped container; ``curves``
and ``fit`` run on a laptop.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from qscat.core.dissociation import da_cross_section
from qscat.core.driven import ve_cross_section
from qscat.core.vibrational import vibrational_states
from qscat.model import F2, N2, NO, DiatomicResonanceModel
from qscat_run.presets import PRESETS

from projects.potential_factory.ansatz import (
    FlexibleDiatomicModel,
    SmoothR,
    from_diatomic,
    with_params,
)
from projects.potential_factory.extract import extract_target, walk_t1
from projects.potential_factory.fit import fit
from projects.potential_factory.report import FitReport
from projects.potential_factory.tracker import ElectronicPair

EV = 27.211386

# (model, (R_hi, R_lo) of the curve scan in bohr, the published crossing R_c)
MOLECULES: dict[str, tuple[DiatomicResonanceModel, tuple[float, float], float]] = {
    "N2": (N2, (3.2, 1.5), N2.R_c),
    "NO": (NO, (3.4, 1.6), NO.R_c),
    "F2": (F2, (4.2, 1.9), F2.R_c),
}

# Electronic-grid ladder: real-region extent, DVR order, complex tail elements.
GRID_LADDER: list[dict[str, float | int]] = [
    {"r_max": 16.0, "order": 7, "n_complex": 6},
    {"r_max": 24.0, "order": 9, "n_complex": 8},
    {"r_max": 32.0, "order": 11, "n_complex": 8},
    {"r_max": 40.0, "order": 13, "n_complex": 10},
]

VPRIMES = [0, 1, 2]


def _grid_label(g: dict[str, float | int]) -> str:
    return f"r{int(g['r_max'])}_o{g['order']}_c{g['n_complex']}"


def _energies(spec) -> np.ndarray:  # EnergySpec(min, max, step)
    return np.arange(spec.min, spec.max + 0.5 * spec.step, spec.step)


# --------------------------------------------------------------------------- curves


def run_curves(mol: str, out: Path, *, n_nodes: int = 45, ladder: list[dict] | None = None) -> dict:
    model, (R_hi, R_lo), R_c = MOLECULES[mol]
    ladder = ladder or GRID_LADDER
    R_desc = np.linspace(R_hi, R_lo, n_nodes)
    results = []
    for g in ladder:
        t0 = time.perf_counter()
        pair = ElectronicPair(
            angles=(35.0, 44.0),
            r_max=float(g["r_max"]),
            order=int(g["order"]),
            n_complex=int(g["n_complex"]),
        )
        R_ok, shift, gamma = walk_t1(model, pair, R_desc)
        v0 = model.v0(R_ok).real
        rows = np.column_stack([R_ok, v0, v0 + shift, shift, gamma])
        label = _grid_label(g)
        path = out / f"{mol}-curves-{label}.csv"
        np.savetxt(
            path, rows, delimiter=",", header="R_bohr,V0_Ha,V_ion_Ha,E_res_Ha,Gamma_Ha", comments=""
        )
        results.append(
            {
                "grid": g,
                "label": label,
                "R": R_ok,
                "shift": shift,
                "gamma": gamma,
                "v0": v0,
                "seconds": time.perf_counter() - t0,
                "n_points": pair.grid_a.n,
            }
        )
        print(
            f"[{mol}] curves {label}: {R_ok.size}/{n_nodes} nodes, "
            f"{pair.grid_a.n} DVR pts, {results[-1]['seconds']:.1f} s",
            flush=True,
        )

    # grid-to-grid differences on the common nodes
    conv = []
    for a, b in zip(results[:-1], results[1:], strict=False):
        common = np.intersect1d(a["R"], b["R"])
        ia = np.searchsorted(a["R"][::-1], common)
        ib = np.searchsorted(b["R"][::-1], common)
        ds = np.abs(a["shift"][::-1][ia] - b["shift"][::-1][ib])
        dg = np.abs(a["gamma"][::-1][ia] - b["gamma"][::-1][ib])
        conv.append(
            {
                "from": a["label"],
                "to": b["label"],
                "R": common,
                "dE": ds,
                "dG": dg,
                "max_dE": float(ds.max()),
                "max_dG": float(dg.max()),
            }
        )
        print(
            f"[{mol}] {a['label']} -> {b['label']}: "
            f"max|dE_res| = {ds.max():.2e} Ha, max|dGamma| = {dg.max():.2e} Ha",
            flush=True,
        )

    summary = {
        "molecule": mol,
        "R_c_published": R_c,
        "grids": [
            {
                "label": r["label"],
                **r["grid"],
                "n_points": r["n_points"],
                "nodes": int(r["R"].size),
                "seconds": r["seconds"],
            }
            for r in results
        ],
        "convergence": [
            {
                "from": c["from"],
                "to": c["to"],
                "max_dE_res_Ha": c["max_dE"],
                "max_dGamma_Ha": c["max_dG"],
            }
            for c in conv
        ],
    }
    (out / f"{mol}-curves-summary.json").write_text(json.dumps(summary, indent=2))
    _plot_curves(mol, results, conv, R_c, out / f"{mol}-factory-curves.png")
    return summary


def _plot_curves(mol: str, results: list[dict], conv: list[dict], R_c: float, path: Path) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping figure", flush=True)
        return
    fig, ax = plt.subplots(3, 1, figsize=(7.5, 10), sharex=False)
    fin = results[-1]
    ax[0].plot(fin["R"], fin["v0"] * EV, "k-", label="$V_0(R)$")
    for r in results:
        ax[0].plot(
            r["R"],
            (r["v0"] + r["shift"]) * EV,
            "-",
            lw=1.0,
            alpha=0.8,
            label=f"$V_{{ion}}$ {r['label']}",
        )
    ax[0].axvline(R_c, color="grey", ls=":", label=f"published $R_c$ = {R_c}")
    ax[0].set(
        xlabel="R (bohr)",
        ylabel="energy (eV)",
        title=f"{mol}-like model: neutral and resonance curves",
    )
    ax[0].legend(fontsize=7)
    for r in results:
        ax[1].plot(r["R"], r["gamma"] * EV, "-", lw=1.0, alpha=0.8, label=r["label"])
    ax[1].axvline(R_c, color="grey", ls=":")
    ax[1].set(xlabel="R (bohr)", ylabel="$\\Gamma(R)$ (eV)", title="autodetachment width")
    ax[1].legend(fontsize=7)
    for c in conv:
        ax[2].semilogy(
            c["R"],
            np.maximum(c["dE"], 1e-16),
            "o-",
            ms=3,
            lw=0.8,
            label=f"|ΔE_res| {c['from']}→{c['to']}",
        )
        ax[2].semilogy(
            c["R"],
            np.maximum(c["dG"], 1e-16),
            "s--",
            ms=3,
            lw=0.8,
            label=f"|ΔΓ| {c['from']}→{c['to']}",
        )
    ax[2].set(
        xlabel="R (bohr)",
        ylabel="grid-to-grid difference (Ha)",
        title="convergence with the electronic grid",
    )
    ax[2].legend(fontsize=6, ncol=2)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)
    print(f"[{mol}] wrote {path}", flush=True)


# --------------------------------------------------------------------------- fit


def _seed(model: DiatomicResonanceModel) -> FlexibleDiatomicModel:
    """The round-trip test's perturbed seed: every fitted constant moved off its value."""
    return with_params(
        from_diatomic(model),
        {
            "D_e": model.D0 * 1.2,
            "R_e": model.R0 * 1.05,
            "beta0": model.alpha0 * 0.9,
            "lam.f_inf": model.lambda_inf * 1.1,
            "alpha.f_inf": model.alpha_c * 1.3,
        },
    )


def run_fit(
    mol: str, out: Path, *, grid: dict | None = None, n_nodes: int = 24, n_eps: int = 8
) -> FitReport:
    model, (R_hi, R_lo), _ = MOLECULES[mol]
    g = grid or GRID_LADDER[-1]
    pair = ElectronicPair(
        angles=(35.0, 44.0),
        r_max=float(g["r_max"]),
        order=int(g["order"]),
        n_complex=int(g["n_complex"]),
    )
    t0 = time.perf_counter()
    target = extract_target(
        model,
        pair=pair,
        R_desc=np.linspace(R_hi, R_lo, n_nodes),
        n_eps=n_eps,
        name=f"{mol} (published model, {_grid_label(g)})",
    )
    fitted, report = fit(target, _seed(model), pair=pair, n_nodes=n_nodes)
    report.to_json(out / f"{mol}-fit-report.json")
    for t in report.tiers:
        print(
            f"[{mol}] {t.name}: {t.status}  rms={t.rms:.3e} max={t.max:.3e}  {t.detail}", flush=True
        )
    print(
        f"[{mol}] fit: crossing_R = {report.crossing_R}, "
        f"DA sign = {report.da_threshold_sign}, {time.perf_counter() - t0:.0f} s",
        flush=True,
    )
    return report


def refitted_model(mol: str, out: Path) -> FlexibleDiatomicModel:
    """Rebuild the refitted model from the committed FitReport parameters."""
    model = MOLECULES[mol][0]
    report = FitReport.from_json(out / f"{mol}-fit-report.json")
    base = from_diatomic(model)
    if "shell.f_inf" in report.parameters:
        base = base.with_shell(
            SmoothR(f_inf=0.0, f_0=0.0, f_1=1.0, R_f=model.R0, R_e=model.R0),
            report.parameters["alpha_b"],
            report.parameters["r_b"],
        )
    return with_params(base, report.parameters)


# --------------------------------------------------------------------------- cross sections

# DA energy window (E_min, E_max, n) in Hartree where the preset VE window does
# not reach the DA threshold: NO's channel opens at +0.1719 Ha, above its
# 0.004-0.120 Ha VE sweep. F2's DA is exothermic and shares the VE window.
DA_WINDOW: dict[str, tuple[float, float, int]] = {"NO": (0.175, 0.30, 26)}


def _ratio(num: np.ndarray, den: np.ndarray) -> np.ndarray:
    """`num/den` where `den > 0`, NaN elsewhere (closed channels give 0/0)."""
    out = np.full_like(den, np.nan, dtype=np.float64)
    np.divide(num, den, out=out, where=den > 0)
    return out


def _max_dev(num: np.ndarray, den: np.ndarray, floor: float = 0.0) -> float | None:
    mask = den > floor
    if not mask.any():
        return None
    return float(np.max(np.abs(num[mask] / den[mask] - 1.0)))


def run_xs(mol: str, out: Path, *, backend: str = "auto", n_energies: int | None = None) -> dict:
    model = MOLECULES[mol][0]
    refit = refitted_model(mol, out)
    preset = PRESETS[f"{mol}:emoscat"]
    if backend != "auto":
        from qscat.linalg import set_default_backend

        set_default_backend(backend)  # type: ignore[arg-type]
    tgrid = preset.ti_grid()
    E = _energies(preset.default_energies)
    has_da = "da" in preset.valid_observables and mol != "N2"
    E_da = np.linspace(*DA_WINDOW[mol]) if mol in DA_WINDOW else E
    if n_energies is not None:
        if n_energies < E.size:
            E = E[np.linspace(0, E.size - 1, n_energies).astype(int)]
        if n_energies < E_da.size:
            E_da = E_da[np.linspace(0, E_da.size - 1, n_energies).astype(int)]
    print(
        f"[{mol}] xs on {tgrid.grids[0].n} x {tgrid.grids[1].n} = {tgrid.size} unknowns, "
        f"{E.size} VE energies, {E_da.size if has_da else 0} DA energies, backend={backend}",
        flush=True,
    )
    results: dict[str, dict] = {}
    for name, m in (("published", model), ("refit", refit)):
        t0 = time.perf_counter()
        eps, chi = vibrational_states(tgrid.grids[1], m.mu, preset.n_vib, m.v0)
        sig_ve = np.asarray(ve_cross_section(tgrid, m, eps, chi, 0, VPRIMES, E), dtype=np.float64)
        sig_da = None
        if has_da:
            sig_da = np.asarray(
                da_cross_section(tgrid, m, eps, chi, 0, E_da, n_channels=1), dtype=np.float64
            )[:, 0]
        results[name] = {
            "eps": eps,
            "sig_ve": sig_ve,
            "sig_da": sig_da,
            "seconds": time.perf_counter() - t0,
        }
        print(
            f"[{mol}] {name}: VE+DA sweep {results[name]['seconds']:.0f} s; "
            f"eps_0 = {eps[0]:.6f} Ha",
            flush=True,
        )
    header = "E_Ha," + ",".join(
        f"sigma_ve_0to{v}_{n}" for n in ("published", "refit") for v in VPRIMES
    )
    cols = [E] + [
        results[n]["sig_ve"][:, i] for n in ("published", "refit") for i in range(len(VPRIMES))
    ]
    np.savetxt(
        out / f"{mol}-xs-ve.csv", np.column_stack(cols), delimiter=",", header=header, comments=""
    )
    if has_da:
        np.savetxt(
            out / f"{mol}-xs-da.csv",
            np.column_stack([E_da, results["published"]["sig_da"], results["refit"]["sig_da"]]),
            delimiter=",",
            header="E_Ha,sigma_da_published,sigma_da_refit",
            comments="",
        )
    summary = {
        "molecule": mol,
        "unknowns": int(tgrid.size),
        "ve_energies": int(E.size),
        "da_energies": int(E_da.size) if has_da else 0,
        "seconds": {n: results[n]["seconds"] for n in results},
        "eps_v": {n: results[n]["eps"].tolist() for n in results},
        "ve_ratio_max_dev": {
            f"0->{v}": _max_dev(
                results["refit"]["sig_ve"][:, i], results["published"]["sig_ve"][:, i]
            )
            for i, v in enumerate(VPRIMES)
        },
    }
    if has_da:
        p = results["published"]["sig_da"]
        summary["da_ratio_max_dev"] = _max_dev(results["refit"]["sig_da"], p, 1e-3 * float(p.max()))
        summary["da_sigma_max_published"] = float(p.max())
    (out / f"{mol}-xs-summary.json").write_text(json.dumps(summary, indent=2))
    _plot_xs(mol, E, E_da, results, has_da, out)
    return summary


def _plot_xs(
    mol: str, E: np.ndarray, E_da: np.ndarray, results: dict, has_da: bool, out: Path
) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping figure", flush=True)
        return
    p, r = results["published"], results["refit"]
    fig, ax = plt.subplots(
        2, 1, figsize=(7.5, 7.5), sharex=True, gridspec_kw={"height_ratios": [3, 1]}
    )
    for i, v in enumerate(VPRIMES):
        (line,) = ax[0].semilogy(E * EV, p["sig_ve"][:, i], "-", label=f"published $0\\to{v}$")
        ax[0].semilogy(
            E * EV, r["sig_ve"][:, i], "--", color=line.get_color(), label=f"refit $0\\to{v}$"
        )
        ax[1].plot(
            E * EV, _ratio(r["sig_ve"][:, i], p["sig_ve"][:, i]), "-", color=line.get_color()
        )
    ax[0].set(
        ylabel="$\\sigma_{0\\to v'}$ (bohr$^2$)",
        title=f"{mol}-like model: exact 2-D VE, published vs factory refit",
    )
    ax[0].legend(fontsize=7, ncol=2)
    ax[1].axhline(1.0, color="grey", lw=0.8)
    ax[1].set(xlabel="E (eV)", ylabel="refit / published")
    fig.tight_layout()
    fig.savefig(out / f"{mol}-factory-ve.png", dpi=130)
    plt.close(fig)
    if has_da:
        fig, ax = plt.subplots(
            2, 1, figsize=(7.5, 6.5), sharex=True, gridspec_kw={"height_ratios": [3, 1]}
        )
        ax[0].semilogy(E_da * EV, np.maximum(p["sig_da"], 1e-300), "-", label="published")
        ax[0].semilogy(E_da * EV, np.maximum(r["sig_da"], 1e-300), "--", label="refit")
        ax[0].set(
            ylabel="$\\sigma_{DA}$ (bohr$^2$)",
            title=f"{mol}-like model: exact 2-D DA, published vs factory refit",
        )
        ax[0].legend(fontsize=8)
        ax[1].plot(E_da * EV, _ratio(r["sig_da"], p["sig_da"]), "-")
        ax[1].axhline(1.0, color="grey", lw=0.8)
        ax[1].set(xlabel="E (eV)", ylabel="refit / published")
        fig.tight_layout()
        fig.savefig(out / f"{mol}-factory-da.png", dpi=130)
        plt.close(fig)
    print(f"[{mol}] wrote cross-section figures", flush=True)


# --------------------------------------------------------------------------- CLI


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--molecule", choices=list(MOLECULES), required=True)
    ap.add_argument("--stage", choices=["curves", "fit", "xs", "all"], default="all")
    ap.add_argument("--out", type=Path, default=Path("runs/factory-base"))
    ap.add_argument(
        "--backend", default="auto", help="SparseLU backend for the xs stage: auto|mumps|scipy"
    )
    ap.add_argument(
        "--n-energies", type=int, default=None, help="subsample the preset energy sweep (xs stage)"
    )
    ap.add_argument(
        "--quick",
        action="store_true",
        help="smoke settings: 2-grid ladder, 12 curve nodes, 6 energies",
    )
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    stages = ["curves", "fit", "xs"] if a.stage == "all" else [a.stage]
    ladder = GRID_LADDER[:2] if a.quick else None
    if "curves" in stages:
        run_curves(a.molecule, a.out, n_nodes=12 if a.quick else 45, ladder=ladder)
    if "fit" in stages:
        run_fit(
            a.molecule,
            a.out,
            grid=(ladder or GRID_LADDER)[-1],
            n_nodes=10 if a.quick else 24,
            n_eps=4 if a.quick else 8,
        )
    if "xs" in stages:
        run_xs(a.molecule, a.out, backend=a.backend, n_energies=6 if a.quick else a.n_energies)


if __name__ == "__main__":
    main()
