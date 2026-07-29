from __future__ import annotations

import math

import numpy as np
import pytest
from qscat.core.dissociation import anion_electronic_states, da_cross_section
from qscat.core.grids import electronic_grid, nuclear_grid
from qscat.core.vibrational import vibrational_states
from qscat.dvr import TensorGrid
from qscat.model import F2, N2
from qscat.tuning import (
    order_for_wavenumber,
    probe_channel_representation,
    propose_grid,
    refine,
)
from qscat.tuning.resonance import resonance_curve


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


def _small_resonance_elec_grids() -> tuple[object, object]:
    # Small electronic grids injected via elec_grids + a coarse
    # resonance_n_dense keep the crossing/order test fast (well under 1s
    # per grid; no @pytest.mark.slow needed) -- see F2's resonance-aware
    # nuclear-mesh design anchors (docs/physics/discretisation-tuning.md).
    return (
        electronic_grid(r_max=12.0, order=6, n_complex=4, angle_deg=35.0),
        electronic_grid(r_max=12.0, order=6, n_complex=4, angle_deg=44.0),
    )


def test_resonant_nuclear_mesh_order_and_crossing_refinement():
    # The crossing/order test (task-3-brief Step 3): the resonant nuclear
    # mesh's quadrature order matches order_for_wavenumber(K_exit, min_len)
    # (independently recomputed here from the SAME small elec_grids/
    # resonance_n_dense the grid itself used), and the finest real element
    # sits inside the super-refined crossing window around R* ~ 2.6 -- NOT
    # at the inner wall (the exact failure the prior argmax(Gamma) design
    # had).
    small_elec = _small_resonance_elec_grids()
    e_max = 0.05
    # resonance_n_dense must be large enough for resonance_pole_walk to
    # track the pole all the way to the crossing without breaking (and
    # freezing early) -- empirically ~18-20 for F2 with these small grids;
    # still well under a second, no @pytest.mark.slow needed.
    resonance_n_dense = 20

    g_res = propose_grid(
        F2,
        "nuclear",
        (0.01, e_max),
        channel="dissociation",
        elec_grids=small_elec,
        resonance_n_dense=resonance_n_dense,
    )

    # Independently reproduce K_exit/order the same way propose_grid does
    # internally, from the SAME (small, fast) resonance curve.
    R, Vd, _Gamma = resonance_curve(F2, *small_elec, R_max=12.0, n_dense=resonance_n_dense)
    Vd_real = np.real(Vd)
    vd_asym = float(Vd_real[-1])
    K_exit = math.sqrt(2.0 * F2.mu * max(e_max - vd_asym, e_max))
    min_len = 0.15  # qscat.tuning.propose._NUCLEAR_MIN_LEN
    expected_order = order_for_wavenumber(K_exit, min_len)

    assert g_res.spec.quadrature == expected_order
    assert expected_order >= 10

    # The outermost Re(V_d) - v0 sign change, independently located.
    v0_at_R = np.real(F2.v0(R))
    diff = Vd_real - v0_at_R
    sign = np.sign(diff)
    crossings = np.nonzero(np.diff(sign))[0]
    assert crossings.size > 0
    i = int(crossings[-1])
    x0, x1 = R[i], R[i + 1]
    f0, f1 = diff[i], diff[i + 1]
    r_star = float(x0 + (-f0 / (f1 - f0)) * (x1 - x0))
    assert abs(r_star - 2.6) < 0.2  # the design anchor: R* ~ 2.598

    real_elements = [el for el in g_res.spec.elements if el.angle_deg == 0.0]
    real_lengths = [el.length for el in real_elements]
    finest = min(real_lengths)
    assert finest <= 0.05

    boundaries = np.concatenate([[0.0], np.cumsum(real_lengths)])
    idx = real_lengths.index(finest)
    mid = 0.5 * (boundaries[idx] + boundaries[idx + 1])

    # The clamp bounds the crossing window's half-width to _CROSSING_DELTA_MAX
    # (0.18 bohr); allow one extra base-element length of slack, since
    # `refine_elements_in_window` refines any element OVERLAPPING the window in
    # full (so an element straddling the window's edge is refined even though
    # part of it sits just outside).
    assert r_star - 0.18 - 0.2 <= mid <= r_star + 0.18 + 0.2
    # NOT at the inner wall (the argmax(Gamma) failure mode this design
    # supersedes puts the "refined" region at R~0.06 instead).
    assert mid > 1.5


# rtol for the resonant nuclear grid's sigma_DA vs its own once-refined
# solve -- the load-bearing convergence claim (task-3-brief Step 4): the
# resonant a-priori grid should be converged on the FIRST pass, unlike the
# plain v0-only grid (~5x off, see validation/tuning/
# test_emoscat_decks.py::test_f2_2d_da_cross_section_spot_check).
_RESONANT_2D_CONVERGENCE_RTOL = 0.15

# Below this, the resonant grid has not meaningfully lifted off the v0-only
# grid's ~0.31 bohr^2 (the unconverged value) toward the converged ~1.6.
_RESONANT_2D_LIFT_FLOOR = 1.0


@pytest.mark.slow
def test_resonant_nuclear_grid_converges_f2_da_cross_section():
    """The LOAD-BEARING 2-D convergence gate (task-3-brief Step 4): the
    resonance-aware nuclear mesh (crossing super-refine + exit-wave order)
    must give a sigma_DA CONVERGED on the first (unrefined) pass, in
    contrast to the plain v0-only grid's ~5x-off ~0.31 bohr^2
    (`validation/tuning/test_emoscat_decks.py::
    test_f2_2d_da_cross_section_spot_check`).

    Harness copied from that same spot-check: electronic grid via
    `propose_grid(F2, "electronic", ...)`, anion eps/chi via
    `vibrational_states`, `TensorGrid([g_elec, g_R])`,
    `da_cross_section(tg, F2, eps, chi, 0, e_probe)`. Heavy (~2.5 min per
    2-D solve on SuperLU) -- @slow.
    """
    e_max = 0.05
    energy_range = (0.01, e_max)
    e_probe = np.array([0.03])

    g_elec = propose_grid(F2, "electronic", energy_range)
    g_R = propose_grid(F2, "nuclear", energy_range, channel="dissociation")
    g_R_refined = refine(g_R)

    sigmas = []
    for g_n in (g_R, g_R_refined):
        eps, chi = vibrational_states(g_n, F2.mu, 4, F2.v0)
        tg = TensorGrid([g_elec, g_n])
        sigma = float(da_cross_section(tg, F2, eps, chi, 0, e_probe)[0, 0])
        assert np.isfinite(sigma) and sigma >= 0.0
        sigmas.append(sigma)
    sigma_base, sigma_ref = sigmas

    rel = abs(sigma_base - sigma_ref) / abs(sigma_ref)
    assert rel < _RESONANT_2D_CONVERGENCE_RTOL, (sigma_base, sigma_ref, rel)
    assert sigma_base > _RESONANT_2D_LIFT_FLOOR, (sigma_base, sigma_ref)
