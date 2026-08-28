"""Two panels from the screen report: where the pole goes as the anisotropy
is turned on, and how far the width moves against the gate line."""

from __future__ import annotations

import json

import numpy as np

from validation.coupled.observable import GAMMA_TOL
from validation.coupled.screen import RESULTS

FIGURE = "docs/physics/figures/no-coupled-pole-trajectory.png"
R_MARK = 2.4  # bohr, inside NO's resonant region


def main() -> str:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    report = json.loads((RESULTS / "screen.json").read_text())
    R = np.asarray(report["R"], dtype=np.float64)
    j = int(np.argmin(np.abs(R - R_MARK)))

    fig, (ax_traj, ax_gam) = plt.subplots(1, 2, figsize=(11, 4.5))

    # Not every curve's walk reaches every value in the declared `S_VALUES`
    # ladder -- the walk stops where Gamma/eps exceeds 1, and stops earlier for
    # some n_channels than others (measured: every walk here stops by s = 0.5 or
    # 0.6). Iterate the s-values each curve actually has, not the full ladder,
    # or a shorter walk KeyErrors on the s the longer ones reached.
    s_max = 0.0
    for n_ch, curves in sorted(report["s_curves"].items(), key=lambda kv: int(kv[0])):
        s_here = sorted(curves.keys(), key=float)
        s_max = max(s_max, float(s_here[-1]))
        re = [curves[s]["v_d"][j] for s in s_here]
        im = [-0.5 * curves[s]["gamma"][j] for s in s_here]
        ax_traj.plot(re, im, "-o", ms=3, lw=1.0, label=f"$N_l$ = {n_ch}")
    ax_traj.plot(
        report["s_curves"]["1"]["0.0"]["v_d"][j],
        -0.5 * report["s_curves"]["1"]["0.0"]["gamma"][j],
        "k*",
        ms=12,
        label="$s = 0$ (shipped model)",
    )
    ax_traj.set(
        xlabel="Re $E$ (Ha)",
        ylabel="$-\\Gamma/2$ (Ha)",
        # `s_max` is the LONGEST walk's endpoint (whichever N_l reached it),
        # not a value every curve here reaches -- the shorter walks stop
        # earlier because the pole has already left the resonant window for
        # that channel count; see the per-curve legend.
        title=f"pole as $s$: 0 $\\to$ {s_max:g} (longest walk) at "
        f"$\\kappa$ = 0.3, $R$ = {R[j]:.2f} bohr",
    )
    ax_traj.legend(fontsize=8)

    # `kappa_curves[n_ch]` is `{kappa: {s: payload}}` -- one more level than the
    # `s_curves[n_ch]` walked above -- so selecting kappa = 0.5 still leaves an
    # {s: payload} dict, not a payload. The full and fixed-l walks need not stop
    # at the same s, so the comparison point is the largest s BOTH reached, the
    # same matched point `observable._summarize` uses for the gate.
    n_max = str(max(report["n_channels"]))
    full_walk = report["kappa_curves"][n_max]["0.5"]
    fixed_walk = report["kappa_curves"]["1"]["0.5"]
    s_common = sorted(set(full_walk) & set(fixed_walk), key=float)[-1]
    g_full = np.asarray(full_walk[s_common]["gamma"])
    g_fixed = np.asarray(fixed_walk[s_common]["gamma"])
    ax_gam.plot(R, g_full, "-", lw=1.2, label=f"full, $N_l$ = {n_max}")
    ax_gam.plot(R, g_fixed, "--", lw=1.2, label="fixed-$l$")
    ax_gam.set(
        xlabel="$R$ (bohr)",
        ylabel="$\\Gamma$ (Ha)",
        title=f"$(s, \\kappa) = ({s_common}, 0.5)$",
    )
    ax_gam.legend(fontsize=8, loc="upper right")

    ax_rel = ax_gam.twinx()
    rel = np.abs(g_full - g_fixed) / np.maximum(g_fixed, 1e-12)
    ax_rel.plot(R, rel, ":", color="tab:red", lw=1.0)
    ax_rel.axhline(GAMMA_TOL, color="tab:red", lw=0.8, alpha=0.5)
    ax_rel.set_ylabel("relative shift (dotted); gate line", color="tab:red")

    fig.tight_layout()
    fig.savefig(FIGURE, dpi=130)
    plt.close(fig)
    print(f"[coupled] wrote {FIGURE}")
    return FIGURE


if __name__ == "__main__":
    main()
