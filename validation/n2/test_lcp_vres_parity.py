"""moved from libs/qscat/tests/test_lcp.py so the library suite runs from the
sdist; the projects pole walk is the independent oracle for
`local_complex_potential`, and validation is the layer allowed to import
both.
"""

from __future__ import annotations

import numpy as np
import pytest
from qscat.core.grids import electronic_grid
from qscat.core.lcp import local_complex_potential
from qscat.model import N2

from projects.n2_ti_cross_section.nuclear_grid import n2_nuclear_grid
from projects.n2_ti_cross_section.vres import vres_on_grid

pytestmark = pytest.mark.slow


def _elec_grids():
    return (
        electronic_grid(r_max=16.0, order=7, n_complex=6, angle_deg=35.0),
        electronic_grid(r_max=16.0, order=7, n_complex=6, angle_deg=44.0),
    )


def test_matches_n2_vres_oracle():
    g_R = n2_nuclear_grid()
    ga, gb = _elec_grids()
    Vd, Gamma = local_complex_potential(N2, g_R, ga, gb)
    Vd_ref, Gamma_ref = vres_on_grid(g_R)
    real = g_R.points.imag == 0.0
    # compare on the resonance region where both are well-defined (R in [1.5,3.5])
    R = g_R.points.real
    band = real & (R > 1.5) & (R < 3.5)
    assert np.allclose(Vd[band].real, Vd_ref[band].real, atol=5e-3)
    assert np.allclose(Gamma[band], Gamma_ref[band], atol=5e-3)
