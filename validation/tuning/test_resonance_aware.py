"""Task 5 gate: the resonance-aware DA nuclear mesh (`propose_grid(...,
channel="dissociation")`) closes Genuine finding #3 (`test_emoscat_decks.py`
/ `docs/physics/discretisation-tuning.md`) -- the 1-D-probe-passing but
2-D-UNCONVERGED F2 nuclear grid. This module is the SIZE + CONVERGENCE gate
for that fix, on both F2 (the case finding #3 diagnosed) and H2+ (the other
resonant/dissociative model in the repo).

**Controller-verified numbers (2026-07-28) -- used here verbatim, not
re-derived (the resonance scan + 2-D solves behind them are expensive; see
the per-test docstrings for the cost budget each one paid):**

- F2: `propose_grid(F2, "nuclear", (0.01, 0.05), channel="dissociation")` ->
  n=1000 (order 14). The eMoScat deck
  (`validation.tuning.test_emoscat_decks.CONFIGS["F2"].da_grid().grids[1]`)
  = 974 points (order 14). Ratio 1000/974 = 1.027x -- DECK-PARITY, not the
  "10-20% smaller" figure this sub-project originally hoped for (see below).
  sigma_DA(E=0.03) = 1.6562 bohr^2 on this grid, CONVERGED (the eMoScat
  deck's own reference value is 1.66; finding #3's refine^2 value was
  1.658). The OLD v0-only grid (`channel="ve"`) was smaller (609 points)
  but gave sigma_DA ~ 0.31 -- 5x too low (finding #3's bug).
- H2+: resonant nuclear grid -> n=489 (order 8). The proxy deck
  (`validation.h2plus.config.proxy_grid().grids[1]`) = 510 points (order
  8). Ratio 489/510 = 0.959x -- ~4% SMALLER. Convergence is NOT
  laptop-verifiable (the full 2-D DR problem is ~1.15M unknowns --
  Docker/MUMPS-sized, consistent with `validation/h2plus/test_dr.py`'s
  existing small-proxy-only laptop gate); size + successful build is the
  laptop-verifiable part here.

**The honest finding (state plainly, not spun):** the "10-20% smaller than
the hand deck" expectation this sub-project set out with does NOT hold for
F2. eMoScat's F2 DA deck is a near-optimal expert hand-tuning (it already
hand-places extra-fine sub-0.1-bohr elements at the R~2.5-2.7 bohr
resonance crossing -- see finding #3). The OLD `propose_grid` (`channel=
"ve"`, v0-alone) was smaller than that deck (609 vs 974) ONLY because it
was under-converged, not because it found a genuinely cheaper discretisation
-- that was finding #3 itself. To reach 2-D convergence you need
approximately deck-sized resolution; the resonance-aware tuner reaches it
AUTOMATICALLY, at deck-parity (F2, 1.027x) or a few percent under (H2+,
0.959x). The deliverable this gate certifies is **convergence + automation
at deck-competitive size**, not a point-count win.
"""

from __future__ import annotations

import numpy as np
import pytest
from qscat.core.dissociation import da_cross_section
from qscat.core.vibrational import vibrational_states
from qscat.dvr import TensorGrid
from qscat.model import F2, H2P
from qscat.tuning import propose_grid, refine

from validation.diatomic.config import CONFIGS
from validation.h2plus.config import proxy_grid

# Deck-parity margin for F2's resonance-aware grid vs the eMoScat deck --
# 1000/974 = 1.027x measured; 5% headroom around that.
_F2_DECK_PARITY_MARGIN = 1.05

# Minimum DVR order each resonant grid must reach -- sanity that the
# exit-wave sizing (`order_for_wavenumber` against K_exit) actually fired,
# rather than falling back to some low default order.
_F2_MIN_ORDER = 14
_H2P_MIN_ORDER = 8

_F2_ENERGY_RANGE = (0.01, 0.05)
_H2P_ENERGY_RANGE = (0.0, 0.05)


@pytest.mark.slow
def test_f2_resonant_nuclear_grid_is_deck_parity_sized() -> None:
    """SIZE gate: the resonance-aware F2 nuclear grid is deck-PARITY sized
    (<= 1.05x the 974-point eMoScat deck), not smaller -- see the module
    docstring's honest finding. The assertion itself is cheap; building the
    resonant grid pays a ~60-90s adiabatic resonance-curve scan (two-angle
    ECS pole match over the R range) -- @slow.
    """
    deck = CONFIGS["F2"].da_grid().grids[1]

    g = propose_grid(F2, "nuclear", _F2_ENERGY_RANGE, channel="dissociation")

    assert np.all(np.isfinite(g.points))
    assert g.spec.quadrature >= _F2_MIN_ORDER
    assert g.n <= _F2_DECK_PARITY_MARGIN * deck.n, (g.n, deck.n)


@pytest.mark.slow
def test_h2p_resonant_nuclear_grid_is_no_larger_than_proxy_deck() -> None:
    """SIZE gate: the resonance-aware H2+ nuclear grid is no larger than the
    510-point proxy deck (in fact ~4% smaller, 489 vs 510 -- see the module
    docstring). Same ~60-90s resonance-scan cost as the F2 grid -- @slow.
    """
    deck = proxy_grid().grids[1]

    g = propose_grid(H2P, "nuclear", _H2P_ENERGY_RANGE, channel="dissociation")

    assert np.all(np.isfinite(g.points))
    assert g.spec.quadrature >= _H2P_MIN_ORDER
    assert g.n <= deck.n, (g.n, deck.n)


# rtol for sigma_DA(base resonance-aware grid) vs sigma_DA(once h-refined) --
# the same convergence bar `test_tuning_propose.py::
# test_resonant_nuclear_grid_converges_f2_da_cross_section` (qscat-level)
# uses; here re-run as the validation-level deck-comparison gate, alongside
# the eMoScat-deck-anchored context in the module docstring.
_2D_CONVERGENCE_RTOL = 0.15

# Below this, the resonant grid has not meaningfully lifted off the old
# v0-only grid's ~0.31 unconverged value toward the converged ~1.6-1.66.
_CONVERGED_FLOOR = 1.0


@pytest.mark.slow
def test_f2_resonant_nuclear_grid_2d_da_cross_section_converges() -> None:
    """CONVERGENCE gate (finding #3, closed): copies the
    `test_emoscat_decks.py::test_f2_2d_da_cross_section_spot_check` harness
    (electronic grid via `propose_grid(F2, "electronic", ...)`, vibrational
    eps/chi via `vibrational_states`, `da_cross_section` on
    `TensorGrid([g_elec, g_R])`), but builds `g_R` via `channel=
    "dissociation"` (the resonance-aware mesh) instead of the OLD `channel=
    "ve"` (v0-alone) path that spot-check used to diagnose the bug.

    Controller-verified: sigma_base = 1.6562 bohr^2 on this exact grid
    (E=0.03), agreeing with the eMoScat deck's own reference (1.66) and
    finding #3's refine^2 value (1.658). NOT re-run here -- a full pass
    costs ~8 minutes (two driven 2-D solves, up to ~1000-pt nuclear x the
    F2 electronic grid, on SuperLU) -- far over the fast-suite budget, so
    this lives behind `@slow` exactly as its sibling spot-check does.
    """
    e_probe = np.array([0.03])

    g_elec = propose_grid(F2, "electronic", _F2_ENERGY_RANGE)
    g_nuc = propose_grid(F2, "nuclear", _F2_ENERGY_RANGE, channel="dissociation")
    g_nuc_refined = refine(g_nuc)

    sigmas = []
    for g_n in (g_nuc, g_nuc_refined):
        eps, chi = vibrational_states(g_n, F2.mu, 4, F2.v0)
        tg = TensorGrid([g_elec, g_n])
        sigma = float(da_cross_section(tg, F2, eps, chi, 0, e_probe)[0, 0])
        assert np.isfinite(sigma) and sigma >= 0.0
        sigmas.append(sigma)
    sigma_base, sigma_refined = sigmas

    rel = abs(sigma_base - sigma_refined) / abs(sigma_refined)
    assert rel < _2D_CONVERGENCE_RTOL, (sigma_base, sigma_refined, rel)
    assert sigma_base > _CONVERGED_FLOOR, (sigma_base, sigma_refined)


# H2+ 2-D convergence is DEFERRED, not gated here: the full DR problem is
# ~1.15M unknowns at production size (electronic real region to 1300 bohr --
# see `validation.h2plus.config.full_grid`), Docker/MUMPS-sized. Even the
# laptop-feasible `proxy_grid` is too heavy for a routine gate (192k
# unknowns; `validation/h2plus/test_dr.py` runs its own @slow tests on a
# SMALLER ad hoc grid instead, not `proxy_grid`, for exactly this reason).
# There is no laptop-feasible way to 2-D-converge-check the resonance-aware
# H2+ nuclear grid against a real DR cross section without that same
# Docker/MUMPS budget; the SIZE gate above (
# `test_h2p_resonant_nuclear_grid_is_no_larger_than_proxy_deck`) is the
# laptop-verifiable half of this molecule's story. A Docker-scale 2-D H2+
# DR convergence check is a follow-on, consistent with the existing H2+
# 2-D handling (`validation/h2plus/test_dr.py`, `docs/physics/h2plus-dr.md`).
