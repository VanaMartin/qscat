"""TEMPORARY parity gate for the exp-M1 graduation (deleted with the old copy).

Pins `projects/n2_ti_cross_section/cross_section.py::ve_cross_section` (dense
np.linalg.solve) against `qscat.core.lcp.lcp_ve_cross_section` (sparse
SparseLU + refactor) on the real N2 system, before the projects copy is
deleted. Dense LAPACK vs SuperLU on the same ~300x300 complex-symmetric
matrix: the target is 1e-12 relative; if the measured cross-solver floor is
above that, gate at 10x the measured maximum and record it below (the
test_anchors.py tolerance-derivation pattern; cf. the ci-test-portability
rule against pinning sparse-solve outputs at 1e-12 cross-arch).

MEASURED (2026-08-26): max |sigma_new/sigma_old - 1| = 3.574e-13 over the
grid below (dense np.linalg.solve vs SuperLU, well under the 1e-12 target,
so RTOL is kept at 1e-12 rather than loosened).
"""

from __future__ import annotations

import numpy as np
from qscat.core.lcp import lcp_ve_cross_section
from qscat.model import N2

from projects.n2_ti_cross_section.cross_section import ve_cross_section
from validation.n2.cross_section import build_system

RTOL = 1e-12  # raise to 10x the measured floor if dense-vs-sparse exceeds it


def test_graduated_solver_reproduces_the_projects_copy():
    grid, eps, chi, Vd, Gamma = build_system()
    E = np.array([0.02, 0.05, 0.1, 0.15, 0.2])
    vprimes = [0, 1, 2, 3]
    old = ve_cross_section(grid, N2.mu, Vd, Gamma, eps, chi, 0, vprimes, E)
    new = lcp_ve_cross_section(grid, N2.mu, Vd, Gamma, eps, chi, 0, vprimes, E)
    assert old.shape == new.shape == (5, 4)
    dev = np.abs(new - old) / np.maximum(np.abs(old), 1e-300)
    print(f"max relative deviation old-vs-new: {dev.max():.3e}")
    np.testing.assert_allclose(new, old, rtol=RTOL)
