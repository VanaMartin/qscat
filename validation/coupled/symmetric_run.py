"""The symmetric comparison: BOTH models pinned to the shipped anion curve.

`renormalised_campaign` solves one `f(R)` -- the one that pins the COUPLED
model -- and runs every channel count on it. That answers "given an
interaction tuned so the coupled model is right, how wrong is the truncation",
which conflates truncation error with a `lam` mismatch the truncated model
never asked for.

Here `N_l = 1` gets its OWN `f_1(R)`, solved so it too reproduces the shipped
`E_res(R)`. Both models then carry the correct resonance curve and the only
difference left is the partial-wave coupling itself -- the same discipline
that turned the 58 % width claim into 0.56 %.

One asymmetry cannot be removed and is itself a finding: `f_1` does NOT
converge to the analytic `2/(1+kappa)`. It grows without bound (1.60 at
R = 6, 2.66 at R = 9), because the asymptotic state is an electron localised
on ONE nucleus, which a single partial wave about the MOLECULAR centre cannot
represent at any well depth. So the pinning is honest only in the interaction
region -- which is where the vibrational states live (weight beyond R = 4 is
1e-29), hence where this observable is decided.

Run: `python -m validation.coupled.symmetric_run`
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from projects.no_coupled_channels.renormalised import ScaleTable
from validation.coupled.cross_section import V_INIT, VPRIMES
from validation.coupled.energies import sweep_energies
from validation.coupled.renormalised_campaign import (
    KAPPA_RUN,
    R_ECS_NUCLEAR,
    RESULTS,
    run_one,
)
from validation.diatomic.config import CONFIGS

__all__ = ["main"]


def main(results: Path = RESULTS) -> dict[str, object]:
    d = np.load(results / "f_table_s03_nl1.npz")
    table = ScaleTable.for_ecs_grid(
        np.asarray(d["R"], dtype=float),
        np.asarray(d["f"], dtype=float),
        r_ecs=R_ECS_NUCLEAR,
        # NOT 2/(1+kappa): N_l=1 never reaches it. Hold the last solved value,
        # which is the honest statement of where this model actually is.
        asymptote=float(np.asarray(d["f"], dtype=float)[-1]),
    )
    E = sweep_energies()
    tgrid = CONFIGS["NO"].da_grid()
    t0 = time.perf_counter()
    out = run_one(1, E, tgrid, table)
    wall = time.perf_counter() - t0
    print(f"[sym] N_l=1 on its OWN f_1(R): {E.size} energies in {wall / 60:.1f} min")

    report = {
        "s": 0.3,
        "kappa": KAPPA_RUN,
        "n_channels": 1,
        "v_init": V_INIT,
        "vprimes": VPRIMES,
        "own_table": True,
        "sigma": {
            "E": out.E.tolist(),
            "total": out.total.tolist(),
            "restricted": out.restricted.tolist(),
            "wall_clock_s": wall,
        },
    }
    (results / "cross_section_nl1_own.json").write_text(json.dumps(report, indent=1))
    print(f"[sym] wrote {results / 'cross_section_nl1_own.json'}")
    return report


if __name__ == "__main__":
    main()
