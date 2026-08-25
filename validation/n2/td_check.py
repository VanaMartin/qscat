"""Group D wiring: sigma_TD at the GATED C5 anchors via the time-dependent
(Crank-Nicolson propagation + energy transform) solver
(`projects/n2_td_cross_section`), cross-checked against both the TI solver
(the exact differential oracle TD converges to, see
`docs/physics/n2-td-cross-section.md`) and the Houfek golden data (the same
factor-3 cross-model bound `cross_section.py`/Group C5 already uses).

Only the 4 GATED `validation/n2/reference.ANCHOR_COORDS` anchors are
checked here -- the 2 DOCUMENTED-LIMITED anchors (elastic, near-threshold)
are known LCP-model limitations that apply equally to the TD solver (it
shares the same `V_d(R)`/`Gamma(R)`/doorway machinery as the TI solver), so
re-checking them under Group D would just re-litigate Group C5's already-
documented exclusions.

The correlation functions `c_v'(t)` needed for the energy transform are
E-independent, so ONE Crank-Nicolson trajectory (`v_init=0`) computed with
`vprimes` = every distinct GATED channel and `E` = every distinct GATED
energy yields sigma_TD at all 4 gated anchors -- not 4 separate
propagations. `cross_section.build_system` (a public alias for the same
`functools.lru_cache`d builder) is the one `cross_section.py` uses for
Group C5, so the ~7s `vres_on_grid` cost is paid at most once per process
even though Groups C5 and D both need it.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from projects.n2_td_cross_section.td_cross_section import td_ve_cross_section
from validation.n2 import cross_section, reference

__all__ = ["DT", "N_STEPS", "TDResult", "compute_td_results"]

# Propagation config: T = n_steps*dt = 1500 a.u., dt = 0.025 a.u. -- the same
# (dt, n_steps) tuned and validated in
# `projects/n2_td_cross_section/test_td_cross_section.py` (V1/V2), where it
# gives sigma_TD within ~1-3.5% of the TI oracle at these energies and
# depletes the resonance (||psi(T)||/||psi(0)|| ~ 1e-2 << 0.1) so the
# finite-time energy transform is not truncated. Costs ~9s (one propagation,
# amortized across all 4 gated anchors).
DT = 0.025
N_STEPS = 60000


@dataclass(frozen=True)
class TDResult:
    """One GATED anchor's sigma_TD vs. sigma_TI and sigma_Houfek."""

    energy_ha: float
    channel: int
    sigma_td: float
    sigma_ti: float
    sigma_houfek: float
    ratio_td_ti: float  # sigma_td / sigma_ti -- TD-vs-TI agreement (the real gate)
    ratio_td_houfek: float  # sigma_td / sigma_houfek -- same factor-3 cross-model bound as C5
    ok: bool  # True iff both ratios pass their bounds


def compute_td_results() -> list[TDResult]:
    """sigma_TD at the GATED C5 anchors, from a single CN propagation."""
    grid, eps, chi, Vd, Gamma = cross_section.build_system()
    anchors = cross_section.compute_anchor_results()  # reuses the cached system
    gated = [r for r in anchors if r.gated]
    if not gated:
        return []

    e_list = sorted({r.energy_ha for r in gated})
    vp_list = sorted({r.channel for r in gated})
    e_index = {e: i for i, e in enumerate(e_list)}
    vp_index = {vp: i for i, vp in enumerate(vp_list)}

    sigma_td = td_ve_cross_section(
        grid,
        cross_section.MU,
        Vd,
        Gamma,
        eps,
        chi,
        0,
        vp_list,
        np.array(e_list),
        dt=DT,
        n_steps=N_STEPS,
    )

    results: list[TDResult] = []
    for r in gated:
        td = float(sigma_td[e_index[r.energy_ha], vp_index[r.channel]])
        ratio_td_ti = td / r.sigma_computed if r.sigma_computed != 0 else float("inf")
        ratio_td_houfek = td / r.sigma_houfek if r.sigma_houfek != 0 else float("inf")
        ok = (abs(ratio_td_ti - 1.0) <= 0.10) and (
            1.0 / reference.ANCHOR_FACTOR <= ratio_td_houfek <= reference.ANCHOR_FACTOR
        )
        results.append(
            TDResult(
                energy_ha=r.energy_ha,
                channel=r.channel,
                sigma_td=td,
                sigma_ti=r.sigma_computed,
                sigma_houfek=r.sigma_houfek,
                ratio_td_ti=ratio_td_ti,
                ratio_td_houfek=ratio_td_houfek,
                ok=ok,
            )
        )
    return results
