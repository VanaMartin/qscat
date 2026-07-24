"""Task 5: sigma at the 6 `reference.ANCHOR_COORDS` via the exact 2-D
driven-Lippmann-Schwinger solver (`projects/n2_2d_cross_section`), compared
THREE ways against both Houfek's golden data and the 1-D LCP solver's own
anchor table (`validation/n2/cross_section.py`, sub-project #3/#4).

Framing (do not lose this while reading numbers below): the model is a GIVEN
testbed, never tuned to match anything. Houfek's `CSVE.V00.J00` is the GATE
(an independent implementation of the same model + method) -- agreement
there certifies THIS solver, not the model's realism. Once gated, the exact
2-D result is the ORACLE and the 1-D LCP is the thing UNDER TEST, so
`ratio_lcp_vs_exact` (`sigma_lcp / sigma_exact`) is the primary scientific
output; a large, well-characterized LCP discrepancy is a successful result
of this benchmark, not an error to minimize.

Structure mirrors `cross_section.py` (C5) exactly: a `lru_cache`d
`_build_system()` for the expensive one-time setup, GATED/DOCUMENTED-LIMITED
classification decided from the anchor's `(energy, channel)` -- never
hardcoded coordinates -- and a `compute_*_results()` entry point that
resolves anchors via `reference.anchors()`.

Cross-import note: `validation/` importing `projects/` is allowed (the
reverse is forbidden); the object under test for this benchmark *is*
`projects.n2_2d_cross_section`'s exact solver, so there is nothing to keep
independent of it here -- same rationale `cross_section.py` gives for its
own import of `projects.n2_ti_cross_section`.

Classification: rather than re-derive GATED-vs-DOCUMENTED-LIMITED from
scratch (which would risk silently drifting from C5's rule on a borderline
anchor), this module reuses `AnchorResult.gated` / `.mechanism` from
`cross_section.compute_anchor_results()` verbatim for the corresponding
`(energy, channel)` anchor. This *is* "the same gating rule as C5" -- not a
reimplementation of it that could disagree, but literally the same computed
values -- so the two tables partition the 6 anchors identically by
construction. `mechanism` names the known LCP-model limitation (per
`Exact2dResult`'s own docstring), which is unchanged by which solver
recomputed `sigma_exact`: it describes why the LCP itself is untrustworthy
there, not a property of the exact solver.

Grouping: `reference.anchors()` yields 6 `(energy_row, channel, sigma_houfek)`
triples that collapse to 3 distinct `energy_row` values (four channels share
E=0.2 Ha); anchors are grouped by energy so `ve_cross_section_2d` is called
once per distinct energy with the full channel list, reusing one sparse LU
factorization across every channel at that energy -- 6 anchors, 3 solves.
"""

from __future__ import annotations

import functools
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from qscat.dvr import TensorGrid

from projects.n2_2d_cross_section.convergence import working_tgrid
from projects.n2_2d_cross_section.cross_section_2d import ve_cross_section_2d
from projects.n2_2d_cross_section.hamiltonian2d import MU
from projects.n2_ti_cross_section.vibrational import vibrational_states
from validation.n2 import reference
from validation.n2.cross_section import N_VIB, compute_anchor_results

__all__ = ["Exact2dResult", "build_system", "compute_exact2d_results"]


@dataclass(frozen=True)
class Exact2dResult:
    """One anchor, compared three ways.

    The exact 2-D solver is the ORACLE here; `ratio_lcp_vs_exact` is the
    scientific deliverable, and `ratio_exact_vs_houfek` is the gate that
    certifies the oracle.
    """

    energy_ha: float
    channel: int  # v' (0 = elastic)
    sigma_exact: float  # bohr^2, this sub-project
    sigma_lcp: float  # bohr^2, sub-project #3
    sigma_houfek: float  # bohr^2, CSVE.V00.J00
    ratio_exact_vs_houfek: float  # V4 -- the GATE
    ratio_lcp_vs_exact: float  # V5 -- the DELIVERABLE
    ratio_lcp_vs_houfek: float  # context, already known
    gated: bool  # same classification rule as C5
    mechanism: str  # empty if gated; else the known LCP limitation


System = tuple[TensorGrid, npt.NDArray[np.float64], npt.NDArray[np.complex128]]


@functools.lru_cache(maxsize=1)
def _build_system() -> System:
    """Build the working 2-D tensor grid and neutral-N2 vibrational states ONCE.

    `working_tgrid()` is Task 4's converged grid (N~=26857, ~3s per sparse
    LU solve); vibrational states use the SAME `N_VIB` as the LCP solver
    (`cross_section.N_VIB`, imported rather than restated) so `eps`/`chi`
    cover every channel index the anchors use.
    """
    tgrid = working_tgrid()
    eps, chi = vibrational_states(tgrid.grids[1], MU, N_VIB)
    return tgrid, eps, chi


build_system = _build_system


@functools.lru_cache(maxsize=1)
def compute_exact2d_results() -> list[Exact2dResult]:
    """Compute sigma at all 6 anchors with the exact 2-D solver and compare
    against both Houfek (`reference.anchors()`) and the LCP solver
    (`cross_section.compute_anchor_results()`, never recomputed here).

    Anchors are grouped by their (already-resolved) Houfek row energy so
    each distinct energy costs exactly one sparse LU factorization, reused
    across every channel sharing that energy via one
    `ve_cross_section_2d(..., vprimes=[...], E)` call.
    """
    tgrid, eps, chi = _build_system()
    anchor_rows = reference.anchors()  # [(e_row, channel, sigma_houfek), ...]
    lcp_by_key = {(r.energy_ha, r.channel): r for r in compute_anchor_results()}

    # Group anchor channels by their shared Houfek row energy, preserving
    # first-seen order so distinct energies are solved in a stable sequence.
    channels_by_energy: dict[float, list[int]] = {}
    for e_row, channel, _sigma_houfek in anchor_rows:
        channels_by_energy.setdefault(e_row, []).append(channel)

    sigma_exact_by_key: dict[tuple[float, int], float] = {}
    for e_row, channels in channels_by_energy.items():
        sigma = ve_cross_section_2d(tgrid, eps, chi, 0, channels, e_row)
        for channel, s in zip(channels, sigma, strict=True):
            sigma_exact_by_key[(e_row, channel)] = float(s)

    results: list[Exact2dResult] = []
    for e_row, channel, sigma_houfek in anchor_rows:
        sigma_exact = sigma_exact_by_key[(e_row, channel)]
        lcp = lcp_by_key[(e_row, channel)]
        ratio_exact_vs_houfek = (
            sigma_exact / sigma_houfek if sigma_houfek != 0 else float("inf")
        )
        ratio_lcp_vs_exact = (
            lcp.sigma_computed / sigma_exact if sigma_exact != 0 else float("inf")
        )
        results.append(
            Exact2dResult(
                energy_ha=e_row,
                channel=channel,
                sigma_exact=sigma_exact,
                sigma_lcp=lcp.sigma_computed,
                sigma_houfek=sigma_houfek,
                ratio_exact_vs_houfek=ratio_exact_vs_houfek,
                ratio_lcp_vs_exact=ratio_lcp_vs_exact,
                ratio_lcp_vs_houfek=lcp.ratio,
                gated=lcp.gated,
                mechanism=lcp.mechanism,
            )
        )
    return results
