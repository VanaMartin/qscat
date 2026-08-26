"""The two spin-orbit components of the O2 model: 2Pi_{1/2} and 2Pi_{3/2}.

`python -m validation.factory.fit_o2_so [--out validation/factory/results]`

Alt & Houfek (Sec. III A, Fig. 1) place the two components symmetrically at
`V_ion -/+ Delta_SO(R)/2` around the 2Pi_g curve with the same width and run
the nuclear dynamics for each; the VE cross section is the sum with the
statistical factor 1/3 each (p. 032829-4). Here the same is done to the
potential: from `qscat.model.O2`, `lam(R)` (and `alpha(R)`) are re-polished
against each shifted anion curve (`refine_resonance`: polish only, no
re-tracking -- the T1 tier alone, the neutral curve is untouched) and the
asymptote tier re-verified with the component's own atom + ion limit. Two
reports, `o2-so12-fit-report.json` and `o2-so32-fit-report.json`, whose
constants become `qscat.model.O2_SO12` / `O2_SO32`.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
from qscat.model import O2
from qscat.model.flexible import params

from projects.potential_factory.fit import check_asymptote, refine_resonance
from projects.potential_factory.report import FitReport, ecs_bounded
from projects.potential_factory.tracker import ElectronicPair
from validation.factory.fit_o2 import GRID, IMAGE_TOL
from validation.factory.targets.o2 import o2_target

__all__ = ["COMPONENTS", "fit_component", "main"]

COMPONENTS = {"so12": -1, "so32": +1}


def fit_component(so: int, *, n_nodes: int = 40, polish_nfev: int = 400) -> FitReport:
    pair = ElectronicPair(
        angles=(35.0, 44.0),
        r_max=float(GRID["r_max"]),
        order=int(GRID["order"]),
        n_complex=int(GRID["n_complex"]),
    )
    target = o2_target(so=so)
    assert target.resonance is not None
    t0 = time.perf_counter()
    # polish only: the component's curve is the parent's moved by +-10 meV,
    # far under the fit's residual -- re-tracking from a seed re-fits the
    # smooth forms from scratch and landed 2Pi_1/2 in a wrong basin (E_res
    # rms 31 mHa, measured); `refine_resonance` moves O2's coefficients
    # from where they are.
    model, t1 = refine_resonance(
        target.resonance, O2, pair=pair, n_nodes=n_nodes, tol=IMAGE_TOL, polish_nfev=polish_nfev
    )
    asym = check_asymptote(target.resonance, model, pair, IMAGE_TOL)
    nuclear_deg = 35.0
    R_tail = 12.0 + np.linspace(0.1, 6.0, 8) * np.exp(1j * np.deg2rad(nuclear_deg))
    report = FitReport(
        target_name=target.name,
        parameters=params(model),
        tiers=[t1, asym],
        ecs_bounds_deg=ecs_bounded(model, pair, R_tail, nuclear_deg=nuclear_deg),
        crossing_R=None,
        da_threshold_sign=None,
        provenance={
            k: {"source": v.source, "locator": v.locator} for k, v in target.provenance.items()
        },
    )
    for t in report.tiers:
        head = f"[O2 so={so:+d}] {t.name}: {t.status} rms={t.rms:.3e} max={t.max:.3e}"
        print(f"{head} {t.detail[:220]}")
    print(f"[O2 so={so:+d}] {time.perf_counter() - t0:.0f} s")
    return report


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", type=Path, default=Path("validation/factory/results"))
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    for name, so in COMPONENTS.items():
        fit_component(so).to_json(a.out / f"o2-{name}-fit-report.json")


if __name__ == "__main__":
    main()
