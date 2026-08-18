"""Tests for the NRM vibrational-excitation background overlap (Eq. 38)."""

from __future__ import annotations

import numpy as np
import pytest
from qscat.core.grids import electronic_grid
from qscat.core.nrm.discrete_state import AsymptoticDiscreteState
from qscat.core.nrm.scattering import incident_coefficients
from qscat.core.nrm.vibrational_excitation import j_dk
from qscat.linalg import c_product
from qscat.model import F2


@pytest.fixture(scope="module")
def elec():
    return electronic_grid(r_max=16.0, order=8, n_complex=6)


def test_j_dk_is_the_c_product_of_phi_d_and_the_incident_wave(elec):
    """Eq. (38) in coefficient space is exactly a c-product -- no weights."""
    ds = AsymptoticDiscreteState(elec, F2, R_inf=10.7)
    R = np.array([3.0, 2.5, 2.0])
    got = j_dk(elec, ds, R, energy=0.05, ell=F2.ell)
    k = float(np.sqrt(2.0 * 0.05))
    inc = incident_coefficients(elec, k, F2.ell)
    for i, r in enumerate(R):
        assert abs(got[i] - c_product(ds.phi_d(float(r)), inc)) < 1e-14


def test_j_dk_is_r_independent_for_an_r_independent_discrete_state(elec):
    """Choice B's phi_d does not vary with R, so neither can its overlap."""
    ds = AsymptoticDiscreteState(elec, F2, R_inf=10.7)
    got = j_dk(elec, ds, np.array([6.0, 3.0, 2.0]), energy=0.05, ell=F2.ell)
    assert np.allclose(got, got[0], rtol=1e-14, atol=0.0)


def test_j_dk_scales_with_the_incident_normalization(elec):
    """A doubled incident wave doubles the overlap -- it is LINEAR in J_k,
    which a squared-quantity test could not detect.
    """
    ds = AsymptoticDiscreteState(elec, F2, R_inf=10.7)
    R = np.array([2.5])
    a = j_dk(elec, ds, R, energy=0.05, ell=F2.ell)[0]
    b = j_dk(elec, ds, R, energy=0.20, ell=F2.ell)[0]
    assert a != b and np.isfinite(a) and np.isfinite(b)


def test_j_dk_rejects_non_positive_energy(elec):
    ds = AsymptoticDiscreteState(elec, F2, R_inf=10.7)
    with pytest.raises(ValueError, match="positive"):
        j_dk(elec, ds, np.array([2.5]), energy=0.0, ell=F2.ell)
