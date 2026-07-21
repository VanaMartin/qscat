"""C5 wiring: sigma at the 6 `reference.ANCHOR_COORDS` via the validated
time-independent (TI) VE cross-section solver (`projects/n2_ti_cross_section`),
GENERALLY classified as GATED (real PASS/FAIL) or DOCUMENTED-LIMITED (an
informational NOTE naming the mechanism, never a FAIL).

Cross-import note: unlike `resonance.py` (which reimplements the electronic
grid in-line so Group B stays independent of `projects/`), this module
imports the TI solver itself, package-absolute
(`projects.n2_ti_cross_section`) -- the object under test for C5 *is* that
project's resolvent/driven-equation solver, so there is nothing to keep
independent of it here. This mirrors the already-established reverse
cross-import: `projects/n2_ti_cross_section/test_cross_section.py` already
imports `validation.n2.reference` to get the anchor coordinates; this module
closes the loop from the other side.

Classification (per `.superpowers/sdd/task-4-brief.md`), decided GENERALLY
from the anchor's `(energy, channel)`, never by hardcoding which of the 6
coordinates is which:

- **GATED** (channel is a VE channel, `channel >= 1`, AND clear of its own
  threshold: `E_tot - eps[channel] > reference.ANCHOR_MARGIN_HA`, where
  `E_tot = E + eps[0]` since every anchor is a v_init=0 transition): a real
  PASS/FAIL gate at `1/reference.ANCHOR_FACTOR <= ratio <=
  reference.ANCHOR_FACTOR`.
- **DOCUMENTED-LIMITED** (elastic, `channel == 0`, OR failing the margin
  test): the LCP model is known -- on physical grounds established in
  Task 3, not empirically patched -- to diverge from Houfek's 2D
  calculation here:
    - elastic: the doorway/driven-equation formula is built purely from the
      resonance's `V_d(R)`/`Gamma(R)` and structurally omits the
      non-resonant *background* (direct/potential) scattering that
      dominates the elastic channel away from the resonance peak.
    - near-threshold: `Gamma(R)` has no explicit electron-energy dependence
      (a *local* CP), so every channel opens with the wrong (non-Wigner)
      threshold power law -- `sigma` diverges as `~1/E` toward each
      threshold purely from the model's structure.
  These anchors are reported (ratio always printed) but never fail the
  harness.

`vres_on_grid` costs ~7s (a two-angle electronic pole search continued
across ~300 nuclear grid points); `_build_system` is `lru_cache`d so it runs
at most once per process even though both `experiment.run_checks()` and its
own tests call into this module.
"""

from __future__ import annotations

import functools
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import numpy.typing as npt
from qscat.dvr import FemDvrEcsGrid

from projects.n2_ti_cross_section.cross_section import ve_cross_section
from projects.n2_ti_cross_section.nuclear_grid import n2_nuclear_grid
from projects.n2_ti_cross_section.vibrational import vibrational_states
from projects.n2_ti_cross_section.vres import vres_on_grid
from validation.n2 import reference

__all__ = ["AnchorResult", "compute_anchor_results"]

_CONFIG = json.loads((Path(__file__).parent / "config.json").read_text())
MU = _CONFIG["reduced_mass"]  # N2 nuclear reduced mass (a.u.); matches
# `projects/n2_ti_cross_section/test_cross_section.py`.
N_VIB = 6  # v=0..5 -- enough to cover vprimes up to 3 used by the anchors.


@dataclass(frozen=True)
class AnchorResult:
    """One C5 anchor's TI-solver-vs-Houfek comparison + gating classification."""

    energy_ha: float  # Houfek's tabulated energy nearest the anchor coordinate
    channel: int  # v' (0 = elastic)
    sigma_computed: float  # bohr^2, from the TI solver
    sigma_houfek: float  # bohr^2, from CSVE.V00.J00
    ratio: float  # sigma_computed / sigma_houfek
    gated: bool  # True: real PASS/FAIL gate; False: DOCUMENTED-LIMITED (NOTE)
    mechanism: str  # empty if gated; else names the LCP limitation


System = tuple[
    FemDvrEcsGrid,
    npt.NDArray[np.float64],
    npt.NDArray[np.complex128],
    npt.NDArray[np.complex128],
    npt.NDArray[np.float64],
]


@functools.lru_cache(maxsize=1)
def _build_system() -> System:
    """Build the shared grid / vibrational states / V_d(R),Gamma(R) ONCE.

    `vres_on_grid` walks a two-angle-matched electronic pole-finder
    continuation across ~300 nuclear grid points (~7s) -- cached so repeated
    calls within one process (e.g. `experiment.run_checks()` plus this
    module's own tests) pay that cost exactly once.
    """
    grid = n2_nuclear_grid()
    eps, chi = vibrational_states(grid, MU, N_VIB)
    Vd, Gamma = vres_on_grid(grid)
    return grid, eps, chi, Vd, Gamma


def classify(e_ha: float, channel: int, eps: npt.NDArray[np.float64]) -> tuple[bool, str]:
    """GATED iff `channel` is a VE channel clear of its own threshold.

    General rule (not keyed off specific coordinates): elastic (`channel ==
    0`) is always DOCUMENTED-LIMITED; a VE channel is GATED only if
    `E_tot - eps[channel] > reference.ANCHOR_MARGIN_HA`, with `E_tot = e_ha +
    eps[0]` (every anchor is a v_init=0 transition, matching
    `ve_cross_section`'s own `E_tot = E + eps[v_init]`).
    """
    if channel == 0:
        return False, "LCP-limited: elastic omits non-resonant background scattering"
    e_tot = e_ha + eps[0]
    excess = e_tot - eps[channel]
    if excess <= reference.ANCHOR_MARGIN_HA:
        return False, (
            "LCP-limited: local width (no electron-energy dependence) => ~1/E "
            f"near-threshold divergence (only {excess:.4f} Ha above the "
            f"v'={channel} threshold, margin={reference.ANCHOR_MARGIN_HA} Ha)"
        )
    return True, ""


def compute_anchor_results() -> list[AnchorResult]:
    """Compute sigma at all 6 C5 anchors and classify each GATED/DOCUMENTED-LIMITED.

    Uses `reference.anchors()` to resolve each `(energy_ha, channel)` anchor
    coordinate to Houfek's nearest tabulated energy row and sigma value (the
    golden-data lookup, never hardcoded), then evaluates the TI solver at
    that same row energy so the comparison is apples-to-apples.
    """
    grid, eps, chi, Vd, Gamma = _build_system()

    results: list[AnchorResult] = []
    for e_row, channel, sigma_houfek in reference.anchors():
        sigma_computed = float(
            ve_cross_section(grid, MU, Vd, Gamma, eps, chi, 0, [channel], e_row)[0]
        )
        ratio = sigma_computed / sigma_houfek if sigma_houfek != 0 else float("inf")
        gated, mechanism = classify(e_row, channel, eps)
        results.append(
            AnchorResult(
                energy_ha=e_row,
                channel=channel,
                sigma_computed=sigma_computed,
                sigma_houfek=sigma_houfek,
                ratio=ratio,
                gated=gated,
                mechanism=mechanism,
            )
        )
    return results
