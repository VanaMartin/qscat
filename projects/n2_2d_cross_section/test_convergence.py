"""The chosen working grid must actually be converged (V3).

Tolerances below are set just above the MEASURED spreads recorded in
`.superpowers/sdd/task-4-convergence-table.md`, not the ~1% acceptance
criterion from the spec: sigma on this problem is converged to ~1e-6
relative (four to six orders of magnitude tighter than 1%), so a 1%
tolerance here would be a tolerance that can never fail and would not
actually test anything.
"""

from __future__ import annotations

import numpy as np
import pytest
from qscat.dvr import TensorGrid

from projects.n2_2d_cross_section.convergence import WORKING_GRID, working_tgrid
from projects.n2_2d_cross_section.cross_section_2d import ve_cross_section_2d
from projects.n2_2d_cross_section.electronic_grid import n2_electronic_grid
from projects.n2_2d_cross_section.hamiltonian2d import MU
from projects.n2_ti_cross_section.nuclear_grid import n2_nuclear_grid
from projects.n2_ti_cross_section.vibrational import vibrational_states


def _sigma(tgrid: TensorGrid, E: float = 0.2, vp: int = 1) -> float:
    eps, chi = vibrational_states(tgrid.grids[1], MU, 4)
    return float(ve_cross_section_2d(tgrid, eps, chi, 0, [vp], E)[0])


@pytest.mark.slow
def test_working_grid_is_theta_independent() -> None:
    """THE decisive ECS check: bound-state/scattering results must not move
    when the complex-scaling contour rotates.

    Measured directly at WORKING_GRID's own settings (r_max=16, order=7,
    n_complex=5, nuc_quadrature=10, nuc_n_complex=5), N=26857, relative to
    the theta=35 deg value (sigma=1.256450927036e-01): theta=30 deg deviates
    by 5.521e-07, theta=40 deg by 1.911e-06 -- spread over {30, 35, 40} is
    1.911e-06. theta=25 deg deviates by 6.752e-05 (~35x worse) at this same
    n_complex=5 grid, confirming `.superpowers/sdd/task-4-convergence-table
    .md`'s note that a shallow 25 deg contour combined with few complex
    tail elements under-resolves the rotated continuum -- 25 deg is
    deliberately excluded from this sweep for that documented reason, not
    asserted "converged" here. 1e-4 sits ~50x above the measured 1.9e-6
    spread over {30, 35, 40}: tight relative to the spec's ~1% bar, with
    headroom against run-to-run solver noise.
    """
    base = dict(WORKING_GRID)
    sigmas = []
    for theta in (30.0, 35.0, 40.0):
        params = {**base, "angle_deg": theta}
        tg = TensorGrid(
            [
                n2_electronic_grid(
                    r_max=params["r_max"], angle_deg=theta,
                    order=params["order"], n_complex=params["n_complex"],
                ),
                n2_nuclear_grid(
                    quadrature=params["nuc_quadrature"], r_max=params["nuc_r_max"],
                    n_complex=params["nuc_n_complex"], angle_deg=theta,
                ),
            ]
        )
        sigmas.append(_sigma(tg))
    spread = (max(sigmas) - min(sigmas)) / np.mean(sigmas)
    assert spread < 1e-4, f"theta-dependence {spread:.2%}: grid is NOT converged"


@pytest.mark.slow
def test_working_grid_is_stable_under_refinement() -> None:
    """Refining past the working grid must not move sigma appreciably.

    Measured: WORKING_GRID's sigma (N=26857) deviates from the richer
    BASELINE grid's sigma (N=71476, ~2.7x larger) by 2.368e-06 relative
    (`.superpowers/sdd/task-4-convergence-table.md`). This refinement test
    goes further still (r_max*1.5, order+1, n_complex+3, nuc_quadrature+2,
    nuc_r_max+10, nuc_n_complex+3), so 1e-4 sits comfortably above the
    measured ~2e-6-scale deviation while remaining four orders of magnitude
    inside the spec's ~1% bar.
    """
    coarse = _sigma(working_tgrid())
    fine_tg = TensorGrid(
        [
            n2_electronic_grid(
                r_max=WORKING_GRID["r_max"] * 1.5,
                angle_deg=WORKING_GRID["angle_deg"],
                order=WORKING_GRID["order"] + 1,
                n_complex=WORKING_GRID["n_complex"] + 3,
            ),
            n2_nuclear_grid(
                quadrature=WORKING_GRID["nuc_quadrature"] + 2,
                r_max=WORKING_GRID["nuc_r_max"] + 10.0,
                n_complex=WORKING_GRID["nuc_n_complex"] + 3,
            ),
        ]
    )
    assert abs(_sigma(fine_tg) - coarse) / coarse < 1e-4
