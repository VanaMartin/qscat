"""Small-grid smoke gate for the H2+ TD-DR experiment (`validation.h2plus.td_dr`).

`@pytest.mark.slow`: even a tiny grid propagation with the mpmath Coulomb
outgoing wave is minutes-scale, and the FULL deck (`config.full_grid`, ~1.15M
unknowns) is Docker/MUMPS-only. This gate only asserts the charged (H2P) TD-DR
extraction RUNS end-to-end and returns finite, correctly-shaped sigma_DR for all
three nuclear extractors from ONE shared propagation -- NOT a converged cross
section (that needs the full deck + a long propagation; see `td_dr.main`).

`n_channels=2`: this tiny electronic grid (r_max=20) resolves only ~2 bound
Rydberg exit states -- the count the exit series can support scales with the
electronic grid extent (the full 1300-bohr deck supports eMoScat's 3, hence the
huge grid). Requesting more raises a clear `anion_electronic_states` ValueError.
"""

from __future__ import annotations

import numpy as np
import pytest
from qscat.core.grids import electronic_grid, nuclear_grid
from qscat.dvr import TensorGrid

from validation.h2plus.td_dr import compute_td_dr


def _small() -> TensorGrid:
    return TensorGrid(
        [
            electronic_grid(r_max=20.0, order=6, n_complex=3),
            nuclear_grid(r_max=14.0, quadrature=6, n_complex=3),
        ]
    )


@pytest.mark.slow
def test_td_dr_runs_and_is_finite_all_three_methods() -> None:
    tg = _small()
    energies = np.array([0.01, 0.02])
    e_out, sigmas = compute_td_dr(
        tg, energies=energies, n_steps=5, r0_incident=12.0, n_channels=2
    )

    assert np.array_equal(e_out, energies)
    assert set(sigmas) == {"flow", "delta", "tw"}
    for method, sigma in sigmas.items():
        assert sigma.shape == (2, 2), (method, sigma.shape)
        assert np.all(np.isfinite(sigma)), method
        assert np.all(sigma >= 0.0), method
