"""Load and integrity-check the Houfek time-independent N₂ VE cross-section data."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

DATA_PATH = Path(__file__).parent / "data" / "CSVE.V00.J00"


@dataclass(frozen=True)
class CrossSectionData:
    energy: np.ndarray  # (N,) Hartree
    sigma: np.ndarray   # (N, 31) bohr²; column j = v=0->j (j=0 elastic)

    @property
    def n_energies(self) -> int:
        return int(self.energy.shape[0])

    @property
    def n_channels(self) -> int:
        return int(self.sigma.shape[1])


def load(path: Path = DATA_PATH) -> CrossSectionData:
    raw = np.loadtxt(path)
    return CrossSectionData(energy=raw[:, 0].copy(), sigma=raw[:, 1:].copy())


def integrity_checks(path: Path = DATA_PATH) -> list[tuple[str, bool, str]]:
    d = load(path)
    out: list[tuple[str, bool, str]] = []
    out.append(("C1 shape 400x32", d.n_energies == 400 and d.n_channels == 31,
                f"{d.n_energies} energies, {d.n_channels} channels"))
    inc = bool(np.all(np.diff(d.energy) > 0))
    out.append(("C2 energy strictly increasing (Ha)", inc,
                f"[{d.energy[0]:.4g}, {d.energy[-1]:.4g}]"))
    out.append(("C3 cross sections non-negative (bohr^2)", bool(np.all(d.sigma >= 0.0)),
                f"min={d.sigma.min():.3e}"))

    def first_open(j: int) -> float:
        nz = np.nonzero(d.sigma[:, j] > 0)[0]
        return float(d.energy[nz[0]]) if nz.size else float("inf")

    opens = [first_open(j) for j in range(1, d.n_channels)]
    finite = [o for o in opens if np.isfinite(o)]
    out.append(("C4 channel thresholds ordered", finite == sorted(finite),
                f"{sum(np.isfinite(opens))} channels open in range"))
    return out
