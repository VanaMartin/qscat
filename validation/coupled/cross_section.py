"""The coupled VE cross section: does the fixed-l reduction change the observable?

Runs the exact 2-D driven solve at (s, kappa) = (0.3, 0.5) -- the anisotropy
where the preceding phase's channel truncation is converged to 0.2 % and the
width difference is 58 % over all 41 comparable R -- for N_l = 1 (fixed-l), 3
and 4, on one deck and one energy mesh.

N_l = 2 is deliberately absent from the sweep: it was measured 30 % from
converged against N_l = 4 at (s, kappa) = (0.3, 0.5), so it would be neither
the oracle nor a useful convergence step. It
appears only in the parity identity gate, which lives in
`projects/no_coupled_channels/test_scattering.py`.

The production run needs MUMPS and does not belong on a laptop: measured,
SuperLU factorises the N_l = 2 deck in 208 s and gives no analysis reuse. Run
`probe_one_energy` first -- the cost estimate behind this campaign is an
extrapolation, and one measurement is cheaper than a wrong night.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp
from qscat.core.vibrational import vibrational_states
from qscat.dvr import TensorGrid
from qscat.linalg import SparseLU
from qscat.model import NO

from projects.no_coupled_channels.anisotropy import TwoCentreWell
from projects.no_coupled_channels.model import CoupledModel
from projects.no_coupled_channels.scattering import CoupledSigma, coupled_ve_cross_section
from validation.coupled.energies import sweep_energies
from validation.diatomic.config import CONFIGS

__all__ = [
    "KAPPA_RUN",
    "N_CHANNEL_VALUES",
    "N_VIB",
    "RESULTS",
    "S_RUN",
    "VPRIMES",
    "V_INIT",
    "main",
    "probe_one_energy",
    "run_one",
]

S_RUN = 0.3
KAPPA_RUN = 0.5
V_INIT = 0
VPRIMES = [0, 1, 2, 3, 4]
N_CHANNEL_VALUES = (1, 3, 4)
N_VIB = 8
RESULTS = Path("validation/coupled/results")


def _model(n_channels: int, s: float = S_RUN, kappa: float = KAPPA_RUN) -> CoupledModel:
    return CoupledModel(well=TwoCentreWell(base=NO, s=s, kappa=kappa), n_channels=n_channels)


def probe_one_energy(n_channels: int) -> dict[str, float]:
    """Measure one energy on the production deck: build, factor, solve.

    The campaign's cost estimate is an extrapolation from SuperLU scaling and
    a published MUMPS ratio. This is the measurement that replaces it, and it
    costs one energy rather than a night.
    """
    tgrid = CONFIGS["NO"].da_grid()
    model = _model(n_channels)

    t0 = time.perf_counter()
    H = model.hamiltonian(tgrid)
    build_s = time.perf_counter() - t0

    ident = sp.identity(H.shape[0], dtype=np.complex128, format="csr")
    a = sp.csc_matrix(0.05 * ident - H)
    t0 = time.perf_counter()
    lu = SparseLU(a)
    factor_s = time.perf_counter() - t0

    b = np.zeros(H.shape[0], dtype=np.complex128)
    b[0] = 1.0
    t0 = time.perf_counter()
    lu.solve(b)
    solve_s = time.perf_counter() - t0

    return {
        "n_channels": float(n_channels),
        "unknowns": float(H.shape[0]),
        "nnz": float(H.nnz),
        "build_s": build_s,
        "factor_s": factor_s,
        "solve_s": solve_s,
    }


def run_one(
    n_channels: int,
    energies: npt.NDArray[np.float64],
    grid: TensorGrid | None = None,
    s: float = S_RUN,
    kappa: float = KAPPA_RUN,
) -> CoupledSigma:
    """One model's sweep, on a deck the caller supplies.

    `main` builds the deck ONCE and threads it into every model, so that "one
    deck" is structural rather than three deterministic reconstructions that
    happen to agree. The `None` default exists only so a caller running a
    single model on its own need not build one.
    """
    tgrid = CONFIGS["NO"].da_grid() if grid is None else grid
    eps, chi = vibrational_states(tgrid.grids[1], NO.mu, N_VIB, NO.v0)
    return coupled_ve_cross_section(
        tgrid, _model(n_channels, s, kappa), eps, chi, V_INIT, VPRIMES, energies
    )


def main(
    results: Path = RESULTS,
    energies: npt.NDArray[np.float64] | None = None,
    grid: TensorGrid | None = None,
    s: float = S_RUN,
    kappa: float = KAPPA_RUN,
    filename: str = "cross_section.json",
) -> dict[str, object]:
    """Every model on ONE mesh and ONE deck; writes `filename`.

    `s` is exposed because the anisotropy is not a free knob: past s ~ 0.15 the
    two-centre well unbinds the anion (see `s0_control.py`), and the campaign
    has to be run inside the window where the bound state survives.
    """
    E = sweep_energies() if energies is None else np.asarray(energies, dtype=np.float64)
    # ONE deck, built here and threaded into every model. Rebuilding it inside
    # the loop would give three deterministic reconstructions that agree by
    # luck of purity rather than by construction -- and on the production path
    # nothing would ever check that they had.
    tgrid = CONFIGS["NO"].da_grid() if grid is None else grid
    report: dict[str, object] = {
        "s": s,
        "kappa": kappa,
        "v_init": V_INIT,
        "vprimes": VPRIMES,
        "n_channels": list(N_CHANNEL_VALUES),
        "sigma": {},
    }
    for n_ch in N_CHANNEL_VALUES:
        t0 = time.perf_counter()
        out = run_one(n_ch, E, grid=tgrid, s=s, kappa=kappa)
        elapsed = time.perf_counter() - t0
        report["sigma"][str(n_ch)] = {  # type: ignore[index]
            "E": out.E.tolist(),
            "total": out.total.tolist(),
            "restricted": out.restricted.tolist(),
            "wall_clock_s": elapsed,
        }
        print(f"[coupled] N_l={n_ch}: {E.size} energies in {elapsed / 60:.1f} min")

    results.mkdir(parents=True, exist_ok=True)
    (results / filename).write_text(json.dumps(report, indent=1))
    print(f"[coupled] wrote {results / filename}")
    return report


if __name__ == "__main__":
    import sys

    if "--probe" in sys.argv:
        for n in N_CHANNEL_VALUES:
            print(probe_one_energy(n))
    elif "--s" in sys.argv:
        s_val = float(sys.argv[sys.argv.index("--s") + 1])
        main(s=s_val, filename=f"cross_section_s{s_val:g}.json")
    else:
        main()
