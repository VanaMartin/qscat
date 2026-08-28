"""Gate 3 of the spec: the quadrature route for the Legendre components
must agree with the closed form, over the range the campaign uses and up to
lambda = 10 (the largest l + l' reachable at N_l = 5)."""

from __future__ import annotations

import numpy as np
import pytest
from qscat.model import NO

from projects.no_coupled_channels.anisotropy import TwoCentreWell, v_lambda_closed_form

R_ELEC = np.linspace(0.2, 16.0, 61)[:, None]  # electronic r
R_NUC = np.array([1.8, 2.3, 2.9, 4.0, 6.0])[None, :]  # nuclear R


@pytest.mark.parametrize("lam", list(range(11)))
@pytest.mark.parametrize(("s", "kappa"), [(0.5, 0.0), (1.0, 0.3), (1.0, 0.5)])
def test_v_lambda_matches_the_closed_form(lam: int, s: float, kappa: float) -> None:
    well = TwoCentreWell(base=NO, s=s, kappa=kappa)
    quad = well.v_lambda(lam, R_ELEC, R_NUC)
    exact = v_lambda_closed_form(well, lam, R_ELEC, R_NUC)
    # Scale on the MONOPOLE, not on this component. A component that vanishes
    # by symmetry (every odd lambda at kappa = 0) has the closed form returning
    # exactly zero and the quadrature returning its round-off, so scaling on the
    # component itself divides 5e-16 by 0 and fails a correct implementation.
    # The physical question is whether the component is small compared with the
    # leading one -- which is what this measures.
    scale = float(np.max(np.abs(v_lambda_closed_form(well, 0, R_ELEC, R_NUC))))
    assert float(np.max(np.abs(quad - exact))) / scale < 1e-10


def test_odd_components_vanish_in_the_symmetric_well() -> None:
    well = TwoCentreWell(base=NO, s=1.0, kappa=0.0)
    for lam in (1, 3, 5):
        assert float(np.max(np.abs(well.v_lambda(lam, R_ELEC, R_NUC)))) < 1e-13
