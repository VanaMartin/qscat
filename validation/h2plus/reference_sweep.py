"""The published H2+ dissociative-recombination cross-section sweep.

`data/dr_cross_sections.txt` is the author's own time-independent sigma_DR(E)
sweep for this model (M. Vana, doctoral thesis, Charles University 2017,
Fig. 4.7): 5001 energies over [0.0001, 0.050] Ha at 1e-5 Ha spacing, four
columns -- E, then the three Rydberg exit channels DR_0, DR_1, DR_2. Stored
verbatim, including two duplicated rows the file ships with (identical
duplicates, harmless; `load` leaves them alone and sorts by energy).

**The stored values carry a 2*pi this repository's `dr_cross_section` does
not**, so `load(rescale=True)` (the default) divides them by `2*pi` to put both
on this repository's convention. That factor was measured, not assumed: across
twelve compared energies the ratio reference/computed clusters on 2*pi and the
DR_0 geometric mean lands on 1.001. See docs/physics/h2plus-dr.md. (The
corresponding VE sweep differs by `(2*pi)^2` and is not stored here.)

**Reduced mass.** The sweep was computed with eMoScat's `mu = 918.25`, not this
repository's corrected `918.076` -- substituting 918.25 drops the mean level
discrepancy from 2.4e-6 to 1.1e-7 Ha, which is how the value was identified.
Anything meant to line up with these curves at resonance precision must be
computed at `REFERENCE_MU`; see `mu_matched_model`.

**This sweep is a poor pointwise validation target and a good figure
background.** Its resonances have a median FWHM of 2e-5 Ha and 54 of its 80
prominent peaks are narrower than two of its own samples, so sampling it at a
point mostly measures resonance-position agreement. The repository gates on the
level table (`reference_levels.py`) instead.
"""

from __future__ import annotations

import dataclasses
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import numpy.typing as npt
from qscat.model import H2P, ResonanceModel

__all__ = ["REFERENCE_MU", "DrSweep", "load", "mu_matched_model", "DATA_PATH"]

DATA_PATH = Path(__file__).parent / "data" / "dr_cross_sections.txt"

# The reduced mass the published sweep and level table were computed with.
# NOT this repository's shipped value -- see the module docstring.
REFERENCE_MU = 918.25


@dataclass(frozen=True)
class DrSweep:
    energy: npt.NDArray[np.float64]  # (N,) electron energy, Ha
    sigma: npt.NDArray[np.float64]  # (N, 3) bohr^2; column n = Rydberg channel n
    rescaled: bool  # whether the 2*pi convention factor has been divided out


def load(path: Path = DATA_PATH, *, rescale: bool = True) -> DrSweep:
    """The sweep, sorted by energy, on this repository's sigma convention."""
    raw = np.loadtxt(path)
    raw = raw[np.argsort(raw[:, 0], kind="stable")]
    sigma = raw[:, 1:4].copy()
    if rescale:
        sigma /= 2.0 * math.pi
    return DrSweep(energy=raw[:, 0].copy(), sigma=sigma, rescaled=rescale)


def mu_matched_model(model: ResonanceModel = H2P) -> ResonanceModel:
    """`model` with the reduced mass the published data was computed with.

    For comparing computed level or pole positions against the stored sweep at
    resonance precision. The repository's own results keep the corrected
    `918.076`; this is strictly a comparability shim, and which mass produced a
    given number belongs in that number's caption.
    """
    return dataclasses.replace(model, mu=REFERENCE_MU)  # type: ignore[type-var]
