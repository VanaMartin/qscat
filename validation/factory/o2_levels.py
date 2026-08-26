"""The spectral check: O2's anion vibrational levels in the fitted model
against the same levels in the extracted target curves.

`python -m validation.factory.o2_levels [--report PATH] [--out DIR]`

The vibrational-excitation cross section of O2 is a comb of narrow peaks at
the quasi-bound levels of the anion, so the metric that predicts the figure
is not an rms over the curve but the levels themselves: `E_v - E_0(neutral)`
is where each peak sits, `Gamma_v` how wide it is. Both are 1-D nuclear
eigenproblems in the complex curve `V_ion(R) - i Gamma(R)/2` (Born-Oppenheimer,
no nuclear ECS tail -- every level of interest lies far below the O + O^-
limit, so the only decay is electronic), solved once in the target's curves
and once in the fitted model's (its `E_res(R)`/`Gamma(R)` from a seeded
`walk_t1`). Writes `o2-anion-levels.csv` and prints the table.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import numpy.typing as npt
from qscat.dvr import ElementSpec, FemDvrEcsGrid, GridSpec, eigen, kinetic

from projects.potential_factory.extract import walk_t1
from projects.potential_factory.report import FitReport
from projects.potential_factory.tracker import ElectronicPair
from validation.factory.fit_o2 import GRID
from validation.factory.targets.o2 import o2_model_from_report, o2_target
from validation.factory.targets.o2_data import EV, load_o2

__all__ = ["LevelTable", "anion_levels", "main"]

REPORT = Path("validation/factory/results/o2-fit-report.json")
_N_LEVELS = 30  # reaches ~2.7 eV above the neutral v=0, the top of the VE window
_R_GRID = (1.3, 7.5, 0.1, 12)  # real FEM-DVR: [x_min, x_max], element length, quadrature


@dataclass(frozen=True)
class LevelTable:
    v: npt.NDArray[np.int64]
    E_target: npt.NDArray[np.float64]  # E_v - E_0(neutral), Ha, from the extracted curves
    E_fit: npt.NDArray[np.float64]  # the same in the fitted model
    gamma_target: npt.NDArray[np.float64]  # Gamma_v, Ha
    gamma_fit: npt.NDArray[np.float64]


def _nuclear_grid() -> FemDvrEcsGrid:
    x_min, x_max, h, quad = _R_GRID
    n = round((x_max - x_min) / h)
    return FemDvrEcsGrid(GridSpec(quadrature=quad, elements=[ElementSpec(h)] * n, x_min=x_min))


def _levels(T: npt.NDArray[np.complex128], V: npt.NDArray[np.complex128], n: int) -> np.ndarray:
    E, _ = eigen(T + np.diag(V))
    return np.asarray(E[:n], dtype=np.complex128)


def anion_levels(report: FitReport, *, n_levels: int = _N_LEVELS, so: int = 0) -> LevelTable:
    model = o2_model_from_report(report)
    g = GRID
    pair = ElectronicPair(
        angles=(35.0, 44.0),
        r_max=float(g["r_max"]),
        order=int(g["order"]),
        n_complex=int(g["n_complex"]),
    )
    target = o2_target(so=so)
    assert target.resonance is not None
    c = load_o2()
    grid = _nuclear_grid()
    R = grid.points.real
    T = kinetic(grid, model.mu)

    # neutral v=0 in each curve (the peak positions are measured from it)
    E0_t = _levels(T, np.interp(R, c.R, c.v0).astype(np.complex128), 1)[0].real
    E0_f = _levels(T, model.v0(R), 1)[0].real

    # the target's complex curve, on its own table (clamped outside it: every
    # level below is far inside the table's turning points)
    lo, hi = target.resonance.R_range
    Rt = np.linspace(lo + 1e-6, hi - 1e-6, 600)
    v_t = np.array([float(target.resonance.v_ion(r)) for r in Rt])
    g_t = np.array([float(target.resonance.gamma(r)) for r in Rt])
    ok = np.isfinite(v_t) & np.isfinite(g_t)
    Vt = np.interp(R, Rt[ok], v_t[ok]) - 0.5j * np.interp(R, Rt[ok], g_t[ok])
    Et = _levels(T, Vt.astype(np.complex128), n_levels)

    # the model's curve from a seeded per-node walk, descending from the far end
    R_desc = np.linspace(hi, lo, 120)
    s0 = float(target.resonance.v_ion(hi) - model.v0(hi).real)
    R_ok, shift, gamma = walk_t1(model, pair, R_desc, seed_energy=complex(s0, 0.0))
    order = np.argsort(R_ok)
    R_ok, shift, gamma = R_ok[order], shift[order], gamma[order]
    v_f = model.v0(R_ok).real + shift
    Vf = np.interp(R, R_ok, v_f) - 0.5j * np.interp(R, R_ok, gamma)
    Ef = _levels(T, Vf.astype(np.complex128), n_levels)

    return LevelTable(
        v=np.arange(n_levels),
        E_target=Et.real - E0_t,
        E_fit=Ef.real - E0_f,
        gamma_target=-2.0 * Et.imag,
        gamma_fit=-2.0 * Ef.imag,
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--report", type=Path, default=REPORT)
    ap.add_argument("--out", type=Path, default=Path("validation/factory/results"))
    ap.add_argument(
        "--name", default="o2", help="output stem: o2, o2-so12, o2-so32 (the SO components)"
    )
    ap.add_argument("--so", type=int, default=0, help="target component, -1/0/+1, as the report")
    a = ap.parse_args()
    t = anion_levels(FitReport.from_json(a.report), so=a.so)
    a.out.mkdir(parents=True, exist_ok=True)
    np.savetxt(
        a.out / f"{a.name}-anion-levels.csv",
        np.column_stack(
            [t.v, t.E_target * EV, t.E_fit * EV, t.gamma_target * EV, t.gamma_fit * EV]
        ),
        delimiter=",",
        header="v,E_target_eV,E_fit_eV,Gamma_target_eV,Gamma_fit_eV",
        comments="",
        fmt=["%d", "%.5f", "%.5f", "%.6f", "%.6f"],
    )
    print(" v   E_v-E_0 target(eV)  fit(eV)   diff(meV)   Gamma_v target(meV)  fit(meV)")
    for i in range(t.v.size):
        print(
            f"{t.v[i]:2d}   {t.E_target[i] * EV:10.4f}   {t.E_fit[i] * EV:9.4f}   "
            f"{(t.E_fit[i] - t.E_target[i]) * EV * 1e3:8.1f}     "
            f"{t.gamma_target[i] * EV * 1e3:10.3f}   {t.gamma_fit[i] * EV * 1e3:10.3f}"
        )
    d = (t.E_fit - t.E_target) * EV * 1e3
    print(f"[O2] peak positions: max |diff| = {np.max(np.abs(d)):.1f} meV over v=0..{t.v[-1]}")


if __name__ == "__main__":
    main()
