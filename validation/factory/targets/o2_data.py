"""The O2 curves of Alt & Houfek (2021) as extracted from Fig. 2 (see extract_fig2.py).

`load_o2()` stitches the four extracted curves into what the factory's T1 tier
wants: `V_0(R)`, the resonance curve `V_ion(R) = V_0 + E_res` on one grid (the
dashed real part below the crossing, the bound anion curve above it), and the
local width `Gamma(R)` (the figure plots `2 Gamma`; zero beyond its last point).
Energies are converted to Hartree here; nothing else is altered.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import numpy.typing as npt

EV = 27.211386  # eV per Hartree
_DATA = Path(__file__).resolve().parents[1] / "data" / "o2"
PRECISION_HA = 0.02 / EV  # the extractor's stated vertical precision


@dataclass(frozen=True)
class O2Curves:
    R: npt.NDArray[np.float64]  # bohr, ascending
    v0: npt.NDArray[np.float64]  # Ha, zero at the O(3P)+O(3P) limit
    v_ion: npt.NDArray[np.float64]  # Ha, same zero
    gamma: npt.NDArray[np.float64]  # Ha, >= 0
    R_c: float  # bohr, where the dashed resonance curve meets the bound anion curve
    precision_Ha: float
    source: str


def _read(name: str) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    arr = np.loadtxt(_DATA / f"{name}.csv", delimiter=",", comments="#")
    order = np.argsort(arr[:, 0])
    return arr[order, 0], arr[order, 1] / EV


def load_o2() -> O2Curves:
    R0, v0 = _read("v0")
    Rb, vb = _read("v_ion_bound")
    Rd, vd = _read("e_res_dashed")
    Rg, g2 = _read("gamma_x2")
    R_c = float(Rd[-1])  # the dashed curve ends where the bound curve begins
    # one ascending R grid for V_ion: dashed part below R_c, bound part above
    Ri = np.concatenate([Rd[Rd < R_c], Rb[Rb >= R_c]])
    vi = np.concatenate([vd[Rd < R_c], vb[Rb >= R_c]])
    Ri, keep = np.unique(Ri, return_index=True)
    vi = vi[keep]
    # everything on the V0 grid, restricted to where V_ion is known
    lo, hi = max(R0.min(), Ri.min()), min(R0.max(), Ri.max())
    R = R0[(R0 >= lo) & (R0 <= hi)]
    v0_R = np.interp(R, R0, v0)
    vion_R = np.interp(R, Ri, vi)
    gamma_R = np.clip(np.interp(R, Rg, g2 / 2.0, left=g2[0] / 2.0, right=0.0), 0.0, None)
    return O2Curves(
        R=R,
        v0=v0_R,
        v_ion=vion_R,
        gamma=gamma_R,
        R_c=R_c,
        precision_Ha=PRECISION_HA,
        source="Alt & Houfek, PRA 103, 032829 (2021) Fig. 2, vector-extracted",
    )
