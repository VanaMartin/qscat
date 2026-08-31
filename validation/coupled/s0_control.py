"""The control the cross-section campaign never ran: does the coupled solver,
on the campaign's own deck and mesh, reproduce the SHIPPED NO model?

At `s = 0` the two-centre well collapses to the shipped isotropic Gaussian for
any `kappa`, so `n_channels = 1` here IS `qscat.model.NO` -- the same model
behind `docs/physics/figures/no-2d-ti-cross-section.png`. That published curve
oscillates over three decades; the campaign's curves at `s = 0.3` are smooth
single humps. This run decides whether that difference is physics (the
anisotropy destroying the bound anion) or an artefact of the deck and mesh the
campaign chose.

It is deliberately cheap: one channel, so ~0.4 s per energy with MUMPS.

Run: `python -m validation.coupled.s0_control`
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import numpy.typing as npt
from qscat.core.vibrational import vibrational_states
from qscat.model import NO

from projects.no_coupled_channels.anisotropy import TwoCentreWell
from projects.no_coupled_channels.model import CoupledModel
from projects.no_coupled_channels.scattering import coupled_ve_cross_section
from validation.coupled.cross_section import N_VIB, V_INIT, VPRIMES
from validation.coupled.energies import sweep_energies
from validation.diatomic.config import CONFIGS

__all__ = ["RESULTS", "main"]

RESULTS = Path("validation/coupled/results")


def main(
    results: Path = RESULTS,
    energies: npt.NDArray[np.float64] | None = None,
    grid: object | None = None,
) -> dict[str, object]:
    """Run `s = 0`, `n_channels = 1` on the campaign's deck and mesh."""
    E = sweep_energies() if energies is None else np.asarray(energies, dtype=np.float64)
    tgrid = CONFIGS["NO"].da_grid() if grid is None else grid
    eps, chi = vibrational_states(tgrid.grids[1], NO.mu, N_VIB, NO.v0)

    model = CoupledModel(well=TwoCentreWell(base=NO, s=0.0, kappa=0.5), n_channels=1)
    t0 = time.perf_counter()
    out = coupled_ve_cross_section(tgrid, model, eps, chi, V_INIT, VPRIMES, E)
    wall = time.perf_counter() - t0
    print(f"[s0] {E.size} energies in {wall / 60:.1f} min")

    report = {
        "s": 0.0,
        "kappa": 0.5,
        "v_init": V_INIT,
        "vprimes": VPRIMES,
        "sigma": {
            "E": out.E.tolist(),
            "total": out.total.tolist(),
            "restricted": out.restricted.tolist(),
            "wall_clock_s": wall,
        },
    }
    results.mkdir(parents=True, exist_ok=True)
    (results / "s0_control.json").write_text(json.dumps(report, indent=1))
    print(f"[s0] wrote {results / 's0_control.json'}")
    return report


if __name__ == "__main__":
    main()
