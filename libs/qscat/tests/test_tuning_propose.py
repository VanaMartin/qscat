from __future__ import annotations

import math

import numpy as np
import pytest
from qscat.core.dissociation import anion_electronic_states
from qscat.core.grids import electronic_grid, nuclear_grid
from qscat.model import F2, N2
from qscat.tuning import interaction_region, probe_channel_representation, propose_grid


def test_propose_grid_nuclear_n2_is_a_valid_grid_in_a_sane_range():
    g = propose_grid(N2, "nuclear", (0.04, 0.18))
    assert np.all(np.isfinite(g.points))
    assert np.all(np.isfinite(g.weights))
    assert 10.0 <= g.R0 <= 22.0

    deck_n = nuclear_grid().n  # the committed N2 nuclear deck's point count
    assert g.n < 2 * deck_n


def test_propose_grid_rejects_empty_energy_range():
    with pytest.raises(ValueError):
        propose_grid(N2, "nuclear", (0.18, 0.04))


def test_propose_grid_electronic_n2_is_a_valid_grid():
    g = propose_grid(N2, "electronic", (0.04, 0.18))
    assert np.all(np.isfinite(g.points))
    assert g.R0 > 0.0
    assert g.n > 0


def test_propose_grid_nuclear_f2_resolves_the_da_channel_wave():
    # F2's actual dissociative-attachment (DA) wavenumber at the top of this
    # energy range: K = sqrt(2*mu*(E_max - eps_e)), eps_e the anion bound
    # electronic state at the dissociation limit (docs/physics/
    # diatomic-ve-cross-sections.md quotes the rough "K_R ~ 58" at E~0.03;
    # this is the precise value at E_max=0.05, ~78 -- Task 8's
    # validation.tuning.calibrate calibrated the tuner's phase constant
    # against exactly this wave). The coarse shared N2-style grid famously
    # failed to resolve it (sigma_DA off by ~36 orders); the a-priori
    # proposed grid must not repeat that failure.
    e_max = 0.05
    elec = electronic_grid(r_max=16.0, order=8, n_complex=6)
    eps_e, _ = anion_electronic_states(elec, F2, R_inf=10.7, n_states=1)
    K = math.sqrt(2.0 * F2.mu * (e_max - float(eps_e[0])))

    g = propose_grid(F2, "nuclear", (0.01, e_max))
    result = probe_channel_representation(g, k=K, l=0, mass=F2.mu, rtol=1e-3)
    assert result.converged, result.detail


def test_ve_channel_default_unchanged():
    # channel="ve" (the default) must be BYTE-IDENTICAL to omitting channel
    # entirely -- the whole tuner depends on this path staying v0-only.
    g_default = propose_grid(F2, "nuclear", (0.01, 0.05))
    g_ve = propose_grid(F2, "nuclear", (0.01, 0.05), channel="ve")
    assert g_default.n == g_ve.n
    assert np.allclose(g_default.points, g_ve.points)


def test_propose_grid_rejects_unknown_channel():
    with pytest.raises(ValueError):
        propose_grid(F2, "nuclear", (0.01, 0.05), channel="bogus")


def test_propose_grid_rejects_dissociation_electronic():
    with pytest.raises(ValueError):
        propose_grid(F2, "electronic", (0.01, 0.05), channel="dissociation")


def test_resonant_nuclear_mesh_refines_the_crossing():
    # Small electronic grids injected via elec_grids + a coarse
    # resonance_n_dense keep this real-eigensolve test fast (well under 1s
    # per grid; no @pytest.mark.slow needed).
    #
    # Energy range: NOT the (0.01, 0.05) used elsewhere in this file. At
    # E_max=0.05, F2's V_d(R) asymptote (~-0.13 Ha -- the ground anion
    # state's real, R-independent electron-affinity binding relative to
    # v0's ~0 dissociation limit) is already classically ALLOWED
    # (E_max - V_d_asymptote ~ 0.18 > 0) with a wavenumber (~80) that swamps
    # v0's own outer wavenumber (~42) *everywhere* past R_hi, not just at
    # the crossing -- combined_profile's worst-case merge then over-
    # resolves the whole outer region, and the crossing's SHARE of the
    # (larger) total point budget actually falls relative to the v0-only
    # grid's. At E_max=0.2 the two asymptotic wavenumbers are close enough
    # that this global effect no longer dominates and the crossing-local
    # refinement (from combined_profile's higher in-region k, plus the
    # Gamma-peak turning point) shows through as intended.
    small_elec = (
        electronic_grid(r_max=12.0, order=6, n_complex=4, angle_deg=35.0),
        electronic_grid(r_max=12.0, order=6, n_complex=4, angle_deg=44.0),
    )
    g_ve = propose_grid(F2, "nuclear", (0.01, 0.2))  # v0-only
    g_res = propose_grid(
        F2,
        "nuclear",
        (0.01, 0.2),
        channel="dissociation",
        elec_grids=small_elec,
        resonance_n_dense=10,
    )
    R_lo, R_hi = interaction_region(F2)

    def frac_in(g: object) -> float:
        rp = g.real_points[g.real_points < g.R0]  # type: ignore[attr-defined]
        return float(((rp >= R_lo) & (rp <= R_hi)).sum() / max(rp.size, 1))

    # the resonance-aware grid packs MORE real points into [R_lo, R_hi]
    # than the v0-only one.
    assert frac_in(g_res) > frac_in(g_ve)
