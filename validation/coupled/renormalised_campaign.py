"""The cross section on the RENORMALISED model -- the experiment the earlier
campaigns were clearing the ground for.

Every earlier run used the bare `TwoCentreWell`, whose split `lam` leaves the
anion unbound beyond `R ~ 0.7/s` and misplaces `E_res` by hundreds of meV
right through the Franck-Condon region. Those runs measure a model that is not
NO. Here `lam(R)` is scaled by the tabulated `f(R)` that pins the shipped
model's `E_res(R)` at every `R`, so the anion curve, the crossing and the
dissociation limit are all reproduced by construction, `Gamma(R)` is within
1.6 % of the shipped width, and the ONLY thing differing between `N_l = 1`
and `N_l = 4` is the partial-wave coupling itself.

Run: `python -m validation.coupled.renormalised_campaign`
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import numpy.typing as npt
from qscat.core.vibrational import vibrational_states
from qscat.dvr import TensorGrid
from qscat.model import NO

from projects.no_coupled_channels.model import CoupledModel
from projects.no_coupled_channels.renormalised import RenormalisedTwoCentreWell, ScaleTable
from projects.no_coupled_channels.scattering import CoupledSigma, coupled_ve_cross_section
from validation.coupled.cross_section import N_VIB, V_INIT, VPRIMES
from validation.coupled.energies import sweep_energies
from validation.diatomic.config import CONFIGS

__all__ = ["RESULTS", "load_table", "main", "run_one"]

RESULTS = Path("validation/coupled/results")
S_RUN = 0.3
KAPPA_RUN = 0.5
N_CHANNEL_VALUES = (1, 3, 4)


def _channels_from_argv(default: tuple[int, ...]) -> tuple[int, ...]:
    """Allow `--channels 6` so a convergence rung can be added without a rerun."""
    import sys

    if "--channels" in sys.argv:
        return tuple(int(x) for x in sys.argv[sys.argv.index("--channels") + 1].split(","))
    return default


# The nuclear ECS contour of NO's DA deck turns complex here (measured: the
# last real point is 9.000, the first rotated one 9.0035+0.0035i). `f` must be
# constant from this radius outward or it is not an analytic function of R on
# the tail.
R_ECS_NUCLEAR = 9.0


def load_table(results: Path = RESULTS) -> ScaleTable:
    """The tabulated `f(R)` from `renormalise.py`, truncated at the ECS radius."""
    d = np.load(results / "f_table_s03.npz")
    return ScaleTable.for_ecs_grid(
        np.asarray(d["R"], dtype=float),
        np.asarray(d["f"], dtype=float),
        r_ecs=R_ECS_NUCLEAR,
        asymptote=2.0 / (1.0 + KAPPA_RUN),
    )


def run_one(
    n_channels: int,
    energies: npt.NDArray[np.float64],
    tgrid: TensorGrid,
    table: ScaleTable,
) -> CoupledSigma:
    """One model's sweep on the renormalised well."""
    well = RenormalisedTwoCentreWell(base=NO, table=table, s=S_RUN, kappa=KAPPA_RUN)
    model = CoupledModel(well=well, n_channels=n_channels)
    eps, chi = vibrational_states(tgrid.grids[1], NO.mu, N_VIB, NO.v0)
    return coupled_ve_cross_section(tgrid, model, eps, chi, V_INIT, VPRIMES, energies)


def main(
    results: Path = RESULTS, channels: tuple[int, ...] | None = None, filename: str | None = None
) -> dict[str, object]:
    channels = channels or N_CHANNEL_VALUES
    E = sweep_energies()
    tgrid = CONFIGS["NO"].da_grid()
    table = load_table(results)
    report: dict[str, object] = {
        "s": S_RUN,
        "kappa": KAPPA_RUN,
        "v_init": V_INIT,
        "vprimes": VPRIMES,
        "n_channels": list(channels),
        "renormalised": True,
        "sigma": {},
    }
    for n_ch in channels:
        t0 = time.perf_counter()
        out = run_one(n_ch, E, tgrid, table)
        elapsed = time.perf_counter() - t0
        report["sigma"][str(n_ch)] = {  # type: ignore[index]
            "E": out.E.tolist(),
            "total": out.total.tolist(),
            "restricted": out.restricted.tolist(),
            "wall_clock_s": elapsed,
        }
        print(f"[renorm] N_l={n_ch}: {E.size} energies in {elapsed / 60:.1f} min")

    results.mkdir(parents=True, exist_ok=True)
    name = filename or "cross_section_renormalised.json"
    (results / name).write_text(json.dumps(report, indent=1))
    print(f"[renorm] wrote {results / name}")
    return report


if __name__ == "__main__":
    ch = _channels_from_argv(N_CHANNEL_VALUES)
    tag = "" if ch == N_CHANNEL_VALUES else "_nl" + "-".join(str(c) for c in ch)
    main(channels=ch, filename=f"cross_section_renormalised{tag}.json")
