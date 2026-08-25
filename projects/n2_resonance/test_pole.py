"""Failing-first tests for the N2 two-angle ECS resonance pole finder
(Task 2, sub-project #2 -- the crux). Builds on the validated
`potential.v_eff_el` (Task 0/1) and `grid_n2.n2_electronic_grid` (Task 1).
"""

from projects.n2_resonance import pole
from projects.n2_resonance.grid_n2 import n2_electronic_grid

HARTREE_TO_EV = 27.211386245988
R0 = 2.01943
WINDOW = (0.04, 0.16, -0.05, 0.0)  # Re in [0.04,0.16] Ha, Im in [-0.05,0]


def test_V1_resonance_at_equilibrium():
    ga, gb = n2_electronic_grid(35.0), n2_electronic_grid(44.0)
    E, _resid = pole.find_pole(R0, ga, gb, WINDOW)
    Eres_eV = E.real * HARTREE_TO_EV
    Gamma_eV = max(0.0, -2 * E.imag) * HARTREE_TO_EV
    assert 2.3 <= Eres_eV <= 2.5, Eres_eV
    assert 0.35 <= Gamma_eV <= 0.55, Gamma_eV


def test_V2_pole_is_stable():
    ga, gb = n2_electronic_grid(35.0), n2_electronic_grid(44.0)
    E, resid = pole.find_pole(R0, ga, gb, WINDOW)
    assert resid < 1e-3, resid  # angle-stable (residual << Gamma)
    # resolution stability: coarser grid gives ~same pole (few %)
    ga2, gb2 = (
        n2_electronic_grid(35.0, n_real=6, n_complex=6),
        n2_electronic_grid(44.0, n_real=6, n_complex=6),
    )
    E2, _ = pole.find_pole(R0, ga2, gb2, WINDOW)
    assert abs(E - E2) / abs(E) < 0.05, (E, E2)
