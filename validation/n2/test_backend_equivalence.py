"""V2 (Task 4): the factorization backend must not change the physics.

The exact-2D anchor tests (`validation/n2/test_anchors.py`, Group E) and the
#7 TD tests already run end-to-end through MUMPS in Docker, because
`SparseLU(backend="auto")` resolves to MUMPS there -- so "physics unchanged
through MUMPS" is largely validated by the whole 213-passed suite already.

This module adds ONE explicit, targeted differential check: recompute a #6
exact VE cross section TWICE on the SAME modest grid -- once with every
internal `SparseLU` forced through SuperLU, once through MUMPS -- and assert
the two cross sections agree to a tight rtol. `ve_cross_section_2d` builds its
`SparseLU` internally and exposes no `backend=` kwarg, so the forcing is done
with `qscat.linalg.default_backend`, the process-wide override that
`SparseLU(backend="auto")` consults (an explicit backend at a call site would
still win; here every internal site is `"auto"`).

Skipped unless MUMPS is importable, so it SKIPS on the Mac and RUNS in the
Docker `test` image. The grid is deliberately small (electronic r_max=16,
nuclear quadrature=10 -- the same shape as
`projects/n2_2d_cross_section/test_cross_section_2d.py`'s module grid) so this
is a few sparse solves, not a 250s converged-grid run.
"""

from __future__ import annotations

import numpy as np
import pytest
from qscat.dvr import TensorGrid
from qscat.linalg import default_backend
from qscat.linalg._mumps_backend import mumps_available

from projects.n2_2d_cross_section.cross_section_2d import ve_cross_section_2d
from projects.n2_2d_cross_section.electronic_grid import n2_electronic_grid
from projects.n2_2d_cross_section.hamiltonian2d import MU
from projects.n2_ti_cross_section.nuclear_grid import n2_nuclear_grid
from projects.n2_ti_cross_section.vibrational import vibrational_states

pytestmark = pytest.mark.skipif(
    not mumps_available(), reason="system MUMPS / qscat[mumps] not installed"
)


def _modest_system() -> tuple[TensorGrid, np.ndarray, np.ndarray]:
    tgrid = TensorGrid(
        [
            n2_electronic_grid(r_max=16.0, order=7, n_complex=5),
            n2_nuclear_grid(quadrature=10, r_max=22.0, n_complex=5),
        ]
    )
    eps, chi = vibrational_states(tgrid.grids[1], MU, 4)
    return tgrid, eps, chi


def test_scipy_and_mumps_give_the_same_ve_cross_section() -> None:
    """A #6 exact VE cross section is invariant under the factorization backend.

    Forces the internal `SparseLU` through SuperLU, then through MUMPS's
    complex-symmetric (SYM=2) path, via `default_backend`; both are exact
    solves of the SAME driven Lippmann-Schwinger system, so they must agree
    to solver round-off, NOT merely to a loose cross-model tolerance. A wrong
    MUMPS wiring (dropped lower triangle, wrong symmetry flag, transposed
    solve) would move sigma well outside this rtol.
    """
    tgrid, eps, chi = _modest_system()
    v_init, vprimes, energy = 0, [0, 1, 2], 0.2

    with default_backend("scipy"):
        sigma_scipy = ve_cross_section_2d(tgrid, eps, chi, v_init, vprimes, energy)
    with default_backend("mumps"):
        sigma_mumps = ve_cross_section_2d(tgrid, eps, chi, v_init, vprimes, energy)

    # Both are real, non-negative cross sections; the elastic channel is the
    # largest, so a relative check is meaningful on every open channel.
    np.testing.assert_allclose(sigma_mumps, sigma_scipy, rtol=1e-9, atol=0.0)
