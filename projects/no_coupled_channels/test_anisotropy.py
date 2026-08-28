"""The two-centre well's three load-bearing properties: it collapses to the
shipped isotropic model at s = 0, it produces no odd Legendre components at
kappa = 0, and its quadrature is converged."""

from __future__ import annotations

import numpy as np
import pytest
from qscat.model import NO

from projects.no_coupled_channels.anisotropy import TwoCentreWell

R_ELEC = np.linspace(0.3, 12.0, 41)[:, None]  # electronic r
R_NUC = np.array([1.8, 2.2, 2.6, 3.4])[None, :]  # nuclear R


def test_s0_collapses_to_the_shipped_isotropic_interaction() -> None:
    well = TwoCentreWell(base=NO, s=0.0, kappa=0.3)
    ref = NO.v_int(R_ELEC, R_NUC)
    np.testing.assert_allclose(well.v_block(1, 1, R_ELEC, R_NUC), ref, rtol=0, atol=1e-14)


@pytest.mark.parametrize(("l", "lp"), [(1, 2), (2, 3), (1, 4)])
def test_s0_gives_no_inter_channel_coupling(l: int, lp: int) -> None:
    well = TwoCentreWell(base=NO, s=0.0, kappa=0.3)
    off = well.v_block(l, lp, R_ELEC, R_NUC)
    assert np.max(np.abs(off)) < 1e-14


def test_kappa0_kills_the_odd_lambda_coupling() -> None:
    """A SYMMETRIC two-centre well is the homonuclear case: only even Legendre
    components survive, so l = 1 cannot reach l = 2 however large s is."""
    well = TwoCentreWell(base=NO, s=1.0, kappa=0.0)
    assert np.max(np.abs(well.v_block(1, 2, R_ELEC, R_NUC))) < 1e-13
    # ... but it CAN reach l = 3, across the even lambda = 2.
    assert np.max(np.abs(well.v_block(1, 3, R_ELEC, R_NUC))) > 1e-4


def test_kappa_opens_the_delta_l_equals_one_channel() -> None:
    well = TwoCentreWell(base=NO, s=1.0, kappa=0.3)
    assert np.max(np.abs(well.v_block(1, 2, R_ELEC, R_NUC))) > 1e-4


def test_quadrature_is_converged_at_64_nodes() -> None:
    coarse = TwoCentreWell(base=NO, s=1.0, kappa=0.3, n_nodes=64)
    fine = TwoCentreWell(base=NO, s=1.0, kappa=0.3, n_nodes=128)
    for l, lp in ((1, 1), (1, 2), (2, 3), (4, 4)):
        a = coarse.v_block(l, lp, R_ELEC, R_NUC)
        b = fine.v_block(l, lp, R_ELEC, R_NUC)
        scale = max(float(np.max(np.abs(b))), 1e-30)
        assert float(np.max(np.abs(a - b))) / scale < 1e-12


def test_complex_r_is_not_silently_made_real() -> None:
    """The ECS tail carries complex r; a well that coerces to float would
    destroy the analytic continuation without raising."""
    well = TwoCentreWell(base=NO, s=1.0, kappa=0.3)
    r = np.array([[3.0 + 1.5j]])
    out = well.v_block(1, 1, r, np.array([[2.2]]))
    assert np.iscomplexobj(out)
    assert abs(out.imag).max() > 1e-12
