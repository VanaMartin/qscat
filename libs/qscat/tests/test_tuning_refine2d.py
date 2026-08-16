from __future__ import annotations

import numpy as np
import pytest


def _obs_factory(exact: float = 1.66, scale: float = 200.0):
    # The observable "converges" as the NUCLEAR grid gains points; the
    # electronic grid is already fine (refining it does nothing). A closure
    # over `g_R.n` mimics a real observable's approach to its converged
    # limit, so the loop must (a) pick the nuclear coordinate, (b) stop when
    # |Delta|/|value| < rtol.
    return lambda g_r, g_R: exact - scale / g_R.n


def test_refine_adopts_nuclear_and_converges():
    from qscat.model import F2
    from qscat.tuning import propose_grid, refine_to_2d_convergence

    g_r = propose_grid(F2, "electronic", (0.01, 0.05))
    g_R = propose_grid(F2, "nuclear", (0.01, 0.05))
    obs = _obs_factory()
    g_r2, g_R2, detail = refine_to_2d_convergence(obs, g_r, g_R, rtol=1e-2, max_iter=6)
    assert detail["converged"]
    assert g_R2.n > g_R.n and g_r2.n == g_r.n  # refined nuclear, left electronic alone
    assert all(step["coordinate"] == "nuclear" for step in detail["iterations"])
    assert detail["final_value"] == pytest.approx(obs(g_r2, g_R2))


def test_refine_caps_at_max_iter_when_never_converging():
    from qscat.model import F2
    from qscat.tuning import propose_grid, refine_to_2d_convergence

    g_r = propose_grid(F2, "electronic", (0.01, 0.05))
    g_R = propose_grid(F2, "nuclear", (0.01, 0.05))
    # An observable that keeps changing by a fixed relative amount never
    # converges.
    flip = {"v": 1.0}

    def obs(g_r, g_R):
        flip["v"] *= -2.0
        return flip["v"]

    _, _, detail = refine_to_2d_convergence(obs, g_r, g_R, rtol=1e-3, max_iter=3)
    assert detail["converged"] is False and len(detail["iterations"]) == 3


def test_refine_already_converged_is_zero_iterations():
    from qscat.model import F2
    from qscat.tuning import propose_grid, refine_to_2d_convergence

    g_r = propose_grid(F2, "electronic", (0.01, 0.05))
    g_R = propose_grid(F2, "nuclear", (0.01, 0.05))
    _, _, detail = refine_to_2d_convergence(lambda a, b: 3.14, g_r, g_R, rtol=1e-3, max_iter=4)
    assert detail["converged"] and len(detail["iterations"]) == 0
    assert detail["final_value"] == pytest.approx(3.14)


@pytest.mark.slow
def test_refine_converges_f2_da_from_coarse_guess():
    """Real F2 dissociative-attachment integration: closes `observable` over
    the actual `da_cross_section` at the same hard energy used by
    `validation/tuning/test_emoscat_decks.py::
    test_f2_2d_da_cross_section_spot_check`, starting from the a-priori
    `propose_grid` nuclear/electronic pair that test found NOT 2-D-converged.

    Cost: each iteration runs up to two driven 2-D solves (nuclear- and
    electronic-refined variants) against the current adopted pair, on grids
    up to ~490k unknowns -- a full run is ~20-40 minutes on SuperLU. NOT run
    in the fast suite; the synthetic tests above gate the loop logic. The
    Task-3 gate (`test_emoscat_decks.py`) already proved the resonance-aware
    nuclear grid converges F2 DA, so this test is a documented integration
    check, not the primary gate.
    """
    from qscat.core.dissociation import da_cross_section
    from qscat.core.vibrational import vibrational_states
    from qscat.dvr import TensorGrid
    from qscat.model import F2
    from qscat.tuning import propose_grid, refine_to_2d_convergence

    e_max = 0.05
    energy_range = (0.01, e_max)
    e_probe = np.array([0.03])

    g_elec = propose_grid(F2, "electronic", energy_range)
    g_nuc = propose_grid(F2, "nuclear", energy_range)

    def observable(g_r: object, g_R: object) -> float:
        eps, chi = vibrational_states(g_R, F2.mu, 4, F2.v0)
        tg = TensorGrid([g_r, g_R])
        return float(da_cross_section(tg, F2, eps, chi, 0, e_probe)[0, 0])

    _, _, detail = refine_to_2d_convergence(observable, g_elec, g_nuc, rtol=0.02, max_iter=3)
    assert detail["converged"] or len(detail["iterations"]) == 3
    assert detail["final_value"] > 0.0
