"""Does the renormalisation restore the bound anion? Counted, not assumed.

Splitting the well without rescaling `lam` leaves the deeper centre only
`(1+kappa)/2` of it, and the anion unbinds beyond `R ~ 0.7/s` for any
`s > 0`: at `s = 0.3` NONE of the 41 grid points binds, against 34 in the
shipped model. If the per-`R` rescaling in `renormalise.py` does what it
claims, that count has to come back.

Measured, it does: **33 bound of the 39 points that solve**. The two that do
not, `R = 2.26` and `2.37`, straddle the crossing, where neither the
bound-state filter nor the pole walk classifies the change of character
cleanly; one of them is bound in the shipped model and one is not, so a
crossing treatment recovers 34 exactly. The widths agree alongside: at
`R = 2.15` the shipped model gives `Gamma = 0.011064`, the renormalised model
0.011244 (+1.6 %), the unrenormalised one 0.055500 -- five times too wide,
and worst approaching the crossing where `Gamma` is smallest.

Run: `python -m validation.coupled.bound_count`
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from validation.coupled.renormalise import (
    RESULTS,
    e_res,
    e_res_pole,
    lam_factor,
    lam_factor_pole,
)

__all__ = ["CROSSING", "S_RUN", "main"]

S_RUN = 0.3
# Below this the anion is a resonance and needs the two-angle pole walk; above
# it a bound state, reachable by the cheaper shift-invert filter. The value is
# the shipped model's own crossing, where E_res changes sign.
CROSSING = 2.25


def main(results: Path = RESULTS) -> dict[str, object]:
    """Solve `f(R)` across the whole grid and count what binds."""
    from qscat.model import NO

    d = json.loads((results / "screen.json").read_text())
    R_all = np.asarray(d["R"], dtype=float)
    v0 = np.real(np.asarray(NO.v0(R_all)))
    target = np.asarray(d["s_curves"]["1"]["0.0"]["v_d"], dtype=float) - v0
    g_ship = np.asarray(d["s_curves"]["1"]["0.0"]["gamma"], dtype=float)
    g_unren = np.asarray(d["s_curves"]["1"]["0.3"]["gamma"], dtype=float)

    rows: list[dict[str, float | None]] = []
    print("   R      target      f(R)      Gamma(renorm)   Gamma(shipped)  Gamma(unrenorm)")
    for i, R in enumerate(R_all):
        t = float(target[i])
        n_ch = max(4, int(np.ceil(7 * S_RUN * R / 2)))
        if R < CROSSING:
            f = lam_factor_pole(S_RUN, float(R), t, 4)
            got = e_res_pole(S_RUN, float(R), f, 4, complex(t)) if f else None
        else:
            f = lam_factor(S_RUN, float(R), t, n_ch)
            got = e_res(S_RUN, float(R), f, n_ch, complex(t)) if f else None
        if f is None or got is None:
            print(f"  {R:5.2f}  {t:+.6f}    no root")
            rows.append({"R": float(R), "f": None, "gamma": None})
            continue
        g = -2 * got.imag
        print(
            f"  {R:5.2f}  {t:+.6f}   {f:.5f}    {g:12.6f}   {g_ship[i]:12.6f}   {g_unren[i]:12.6f}"
        )
        rows.append({"R": float(R), "f": f, "gamma": g, "target": t})

    solved = [r for r in rows if r["gamma"] is not None]
    bound = sum(1 for r in solved if abs(float(r["gamma"] or 0.0)) < 1e-6)
    print()
    print(f"solved at {len(solved)}/41 points")
    print(f"BOUND (|Gamma| < 1e-6) after renormalisation : {bound}/41")
    print(f"BOUND in the shipped model                   : {int((g_ship < 1e-6).sum())}/41")
    print(f"BOUND unrenormalised at s=0.3                : {int((g_unren < 1e-6).sum())}/41")

    results.mkdir(parents=True, exist_ok=True)
    (results / "renormalise_s03_grid.json").write_text(json.dumps(rows, indent=1))
    return {"rows": rows, "bound": bound, "solved": len(solved)}


if __name__ == "__main__":
    main()
