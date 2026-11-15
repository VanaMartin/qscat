"""Saved-only local CO plots, with development/validation roles kept separate."""

from pathlib import Path

import numpy as np
from qscat.units import hartree_to_ev

from projects.ukrmol_co.local_representation import DOMAIN_BOHR, LocalRepresentation


def plot_saved_curves(
    curve: LocalRepresentation, evidence: dict, output_dir: str | Path
) -> list[Path]:
    """Render curves and residuals into a new directory outside the checkout.

    No fitting or inference from the plotted samples. Call load_plot_evidence
    first when using the accepted external delivery instead of synthetic inputs.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output = Path(output_dir).resolve()
    if output.is_relative_to(Path(__file__).resolve().parents[2]):
        raise ValueError("Generated experiment plots belong outside the checkout")
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite plot directory: {output}")
    radius = np.linspace(*DOMAIN_BOHR, 601)
    neutral = evidence["inputs/accepted-neutral-result.json"]["fits"]["5z"]["data"]
    pairs = evidence["ten-qualified-pairs.json"]["pairs"]
    figures = []
    paths = []
    try:
        fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
        figures.append(fig)
        for ax, values, title in [
            (axes[0, 0], curve.v0(radius), "5Z neutral: molecular zero at R=2.1323 bohr"),
            (axes[0, 1], curve.er(radius), "Native incident resonance position Er"),
            (axes[1, 0], curve.full_gamma(radius), "Native FULL width Gamma (not Gamma/2)"),
            (axes[1, 1], curve.vres_hybrid(radius), "Mixed-method V0_5Z + Er_native"),
        ]:
            ax.plot(radius, hartree_to_ev(values), label="Frozen representation")
            ax.set(title=title, xlabel="R (bohr)", ylabel="Energy (eV)", xlim=DOMAIN_BOHR)
            ax.grid(alpha=0.2)
        for role, marker, color in [
            ("historical_knot", "o", "black"),
            ("historical_midpoint", "x", "gray"),
        ]:
            rows = [r for r in neutral if r["historical_role"] == role]
            axes[0, 0].scatter(
                [r["R_bohr"] for r in rows],
                hartree_to_ev(np.array([r["actual_hartree"] for r in rows])),
                marker=marker,
                color=color,
                label=role.replace("_", " "),
                s=25,
            )
        for role, marker, color, label in [
            ("C2_development", "o", "black", "8 development nodes"),
            ("C2_prospective_validation", "D", "tab:orange", "2 prospective validation points"),
        ]:
            rows = [p for p in pairs if p["role"] == role]
            r = np.array([p["R_bohr"] for p in rows])
            er = np.array([p["Er_ev"] for p in rows])
            gamma = np.array([p["full_Gamma_ev"] for p in rows])
            for ax, values in [
                (axes[0, 1], er),
                (axes[1, 0], gamma),
                (axes[1, 1], er + hartree_to_ev(curve.v0(r))),
            ]:
                ax.scatter(r, values, marker=marker, color=color, label=label, s=30, zorder=3)
        for ax in axes.flat:
            ax.legend(fontsize=8)
        fig.suptitle(
            "CO local finite model, 1.9–2.5 bohr; no extrapolation\n"
            "Hybrid is NOT a correlated anion total; no electronic/global accuracy claim"
        )

        residual, axs = plt.subplots(1, 3, figsize=(13, 4.5), constrained_layout=True)
        figures.append(residual)
        axs[0].scatter([r["R_bohr"] for r in neutral], [r["error_mev"] for r in neutral], s=20)
        axs[0].set(title="5Z historical fit residuals\nNot prospective physical holdouts")
        validation_r = []
        width_budgets = []
        for label in ["lower", "upper"]:
            score = evidence[f"inputs/{label}-frozen-score.json"]
            r = score["geometry_bohr"]
            validation_r.append(r)
            width_budgets.append(max(0.01 * hartree_to_ev(curve.full_gamma(r)), 0.001) * 1000)
            for sector, marker, color in [("B1", "D", "tab:orange"), ("B2", "+", "black")]:
                row = score["sectors"][sector]
                for ax, actual, predicted in [
                    (axs[1], row["position_ev"], hartree_to_ev(curve.er(r))),
                    (axs[2], row["full_width_ev"], hartree_to_ev(curve.full_gamma(r))),
                ]:
                    ax.scatter(
                        r,
                        (actual - predicted) * 1000,
                        marker=marker,
                        color=color,
                        label=sector if label == "lower" else None,
                        s=45,
                    )
        axs[1].axhline(5, ls="--", color="tab:red", label="Unchanged pointwise gates")
        axs[1].axhline(-5, ls="--", color="tab:red")
        for sign in [-1, 1]:
            axs[2].scatter(
                validation_r,
                sign * np.array(width_budgets),
                marker="_",
                color="tab:red",
                s=150,
                label="Unchanged pointwise gates" if sign == 1 else None,
            )
        axs[1].set(title="Prospective Er: actual − frozen prediction")
        axs[2].set(title="Prospective FULL Gamma: actual − prediction")
        for ax in axs:
            ax.axhline(0, color="gray", lw=0.7)
            ax.set(xlabel="R (bohr)", ylabel="Residual (meV)", xlim=DOMAIN_BOHR)
            ax.grid(alpha=0.2)
        axs[1].legend(fontsize=8)
        axs[2].legend(fontsize=8)
        residual.suptitle("Saved residuals only; no refit and no uniform confidence bound")
        output.mkdir(parents=True, exist_ok=False)
        for figure, name in [(fig, "local-curves"), (residual, "saved-residuals")]:
            for extension in ["png", "svg"]:
                path = output / f"{name}.{extension}"
                figure.savefig(path, dpi=170)
                paths.append(path)
    finally:
        for figure in figures:
            plt.close(figure)
    return paths
