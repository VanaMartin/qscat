"""Tests for the NRM vibrational-excitation background overlap (Eq. 38), the
resonant T-matrix (Eq. 31), and the background T-matrix (Eq. 37)."""

from __future__ import annotations

import numpy as np
import pytest
from qscat.core.grids import electronic_grid, segmented_grid
from qscat.core.nrm.coupling import v_dk_plus
from qscat.core.nrm.discrete_state import AsymptoticDiscreteState
from qscat.core.nrm.nonlocal_potential import continue_to_tail
from qscat.core.nrm.scattering import incident_coefficients
from qscat.core.nrm.vibrational_excitation import (
    _t_background_term2,
    j_dk,
    t_background,
    t_resonant,
)
from qscat.core.vibrational import vibrational_states
from qscat.linalg import c_product
from qscat.model import F2


class _ZeroInteraction:
    """Test double: same molecule as `model`, but `V_int == 0` everywhere.

    Forwards `mu`/`ell`/`charge`/`v0` unchanged so the underlying electronic
    Hamiltonian and vibrational structure are untouched; only the
    interaction that couples the electron to the nuclear coordinate is
    zeroed. `surface` is derived by peeling the real model's own `v_int` back
    off its own `surface`, rather than re-deriving the centrifugal term, so
    it stays correct for any `model.ell`.
    """

    def __init__(self, model) -> None:  # noqa: ANN001 -- test double, any ResonanceModel
        self._model = model

    @property
    def mu(self) -> float:
        return self._model.mu

    @property
    def ell(self) -> int:
        return self._model.ell

    @property
    def charge(self) -> int:
        return self._model.charge

    def v0(self, R):  # noqa: ANN001, ANN201 -- test double
        return self._model.v0(R)

    def v_int(self, r, R):  # noqa: ANN001, ANN201 -- test double
        r = np.asarray(r)
        return np.zeros_like(r, dtype=np.complex128)

    def surface(self, r, R):  # noqa: ANN001, ANN201 -- test double
        return self._model.surface(r, R) - self._model.v_int(r, R)

    def hamiltonian(self, tgrid):  # noqa: ANN001, ANN201 -- test double
        raise NotImplementedError

    def interaction_diag(self, tgrid):  # noqa: ANN001, ANN201 -- test double
        raise NotImplementedError


@pytest.fixture(scope="module")
def elec():
    return electronic_grid(r_max=16.0, order=8, n_complex=6)


@pytest.fixture(scope="module")
def nuc():
    return segmented_grid(
        ((9, 1.8), (1, 2.0), (5, 2.5), (4, 2.7), (20, 10.7)),
        ((1, 11.0), (1, 12.5), (1, 14.0), (3, 30.0)),
        angle_deg=45.0,
        quadrature=14,
    )


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


def test_t_resonant_is_the_weighted_c_product():
    """Eq. (31) with the paper's weights already absorbed by the coefficients."""
    rng = np.random.default_rng(0)
    n = 12
    chi_f = (rng.normal(size=n) + 1j * rng.normal(size=n)).astype(np.complex128)
    v = (rng.normal(size=n) + 1j * rng.normal(size=n)).astype(np.complex128)
    psi = (rng.normal(size=n) + 1j * rng.normal(size=n)).astype(np.complex128)
    assert abs(t_resonant(chi_f, v, psi) - c_product(chi_f, v * psi)) < 1e-14


def test_t_resonant_is_bilinear_not_sesquilinear():
    """A conjugating implementation would give a different answer for complex
    inputs -- this is the check that pins PRA 77's convention over Domcke's.
    """
    rng = np.random.default_rng(1)
    n = 8
    chi_f = (rng.normal(size=n) + 1j * rng.normal(size=n)).astype(np.complex128)
    v = np.ones(n, dtype=np.complex128)
    psi = (rng.normal(size=n) + 1j * rng.normal(size=n)).astype(np.complex128)
    bilinear = t_resonant(chi_f, v, psi)
    sesquilinear = np.vdot(chi_f, v * psi)
    assert abs(bilinear - sesquilinear) > 1e-6 * abs(bilinear)


def test_t_resonant_is_linear_in_psi():
    """T is linear in Psi_d -- doubling it doubles T. sigma ~ |T|^2 hides
    phase, so linearity is asserted on T itself, before squaring.
    """
    rng = np.random.default_rng(2)
    n = 10
    chi_f = rng.normal(size=n).astype(np.complex128)
    v = rng.normal(size=n).astype(np.complex128)
    psi = (rng.normal(size=n) + 1j * rng.normal(size=n)).astype(np.complex128)
    assert abs(t_resonant(chi_f, v, 2.0 * psi) - 2.0 * t_resonant(chi_f, v, psi)) < 1e-13


def test_t_resonant_rejects_mismatched_shapes():
    with pytest.raises(ValueError, match="same length"):
        t_resonant(
            np.ones(4, dtype=np.complex128),
            np.ones(5, dtype=np.complex128),
            np.ones(4, dtype=np.complex128),
        )


def test_t_background_vanishes_with_no_interaction(elec, nuc):
    """V_int = 0 kills term1; term2 is separately forced to zero by the
    literal zero `v_dk_f` passed in (not because `V_int = 0` makes `V_dk`
    vanish -- `v_dk_f` here is an input, not derived from `model`). Both
    zeros have independent sources, so the whole background T-matrix must
    be exactly zero -- the sharpest structural check available without an
    oracle.
    """
    ds = AsymptoticDiscreteState(elec, F2, R_inf=float(nuc.R0))
    R = np.sort(nuc.points[nuc.points.imag == 0.0].real)[::-1]
    eps, chi = vibrational_states(nuc, F2.mu, 2, F2.v0)
    zero_v = np.zeros(nuc.n, dtype=np.complex128)
    t = t_background(elec, nuc, _ZeroInteraction(F2), ds, R, chi[0], chi[0], zero_v, 0.05, 0.05)
    assert abs(t) < 1e-12


def test_t_background_term1_alone_is_nonzero_for_the_real_interaction(elec, nuc):
    """With `v_dk_f = 0` term 2 is forced to zero by construction (its
    formula multiplies by `v_dk_f` elementwise), regardless of the model --
    so with the REAL F2 interaction, `t_background` here returns term 1
    alone. A dropped or zeroed term 1 (e.g. `return 0.0 * term1 - term2`)
    would make this exactly zero; the real V_int is not identically zero on
    this grid, so a correct term 1 must not be either.
    """
    ds = AsymptoticDiscreteState(elec, F2, R_inf=float(nuc.R0))
    R = np.sort(nuc.points[nuc.points.imag == 0.0].real)[::-1]
    eps, chi = vibrational_states(nuc, F2.mu, 2, F2.v0)
    zero_v = np.zeros(nuc.n, dtype=np.complex128)
    t = t_background(elec, nuc, F2, ds, R, chi[0], chi[0], zero_v, 0.05, 0.05)
    assert abs(t) > 1e-3


def test_t_background_distinguishes_initial_from_final_energy(elec, nuc):
    """`e_kin_i` feeds the initial channel (`inc_i`, `j_dk`); `e_kin_f`
    feeds term 1's P-space continuum (`phi^+_{k_f}`). Swapping the two
    energies (holding chi_i, chi_f and v_dk_f fixed) must change the
    answer by more than floating-point noise -- otherwise the two roles
    could have been silently interchanged inside `t_background` without
    any test noticing.
    """
    ds = AsymptoticDiscreteState(elec, F2, R_inf=float(nuc.R0))
    R = np.sort(nuc.points[nuc.points.imag == 0.0].real)[::-1]
    eps, chi = vibrational_states(nuc, F2.mu, 2, F2.v0)
    v = continue_to_tail(v_dk_plus(elec, F2, ds, R, 0.04), R, nuc)
    a = t_background(elec, nuc, F2, ds, R, chi[0], chi[1], v, 0.05, 0.04)
    b = t_background(elec, nuc, F2, ds, R, chi[0], chi[1], v, 0.04, 0.05)
    assert abs(a - b) > 0.03 * abs(a)


def test_t_background_is_linear_in_the_final_channel(elec, nuc):
    """T^bg is linear in chi_f -- doubling it doubles T, asserted before
    squaring since sigma ~ |T|^2 is blind to it. `v_dk_f` is evaluated at
    `e_kin_f` to match `t_background`'s own contract.
    """
    ds = AsymptoticDiscreteState(elec, F2, R_inf=float(nuc.R0))
    R = np.sort(nuc.points[nuc.points.imag == 0.0].real)[::-1]
    eps, chi = vibrational_states(nuc, F2.mu, 2, F2.v0)
    v = v_dk_plus(elec, F2, ds, R, 0.04)
    v_full = continue_to_tail(v, R, nuc)
    a = t_background(elec, nuc, F2, ds, R, chi[0], chi[1], v_full, 0.05, 0.04)
    b = t_background(elec, nuc, F2, ds, R, chi[0], 2.0 * chi[1], v_full, 0.05, 0.04)
    assert abs(b - 2.0 * a) < 1e-10 * abs(a)


def test_t_background_term2_matches_an_independent_assembly(elec, nuc):
    """The second term of Eq. (37), built term by term from j_dk and v_dk."""
    ds = AsymptoticDiscreteState(elec, F2, R_inf=float(nuc.R0))
    R = np.sort(nuc.points[nuc.points.imag == 0.0].real)[::-1]
    eps, chi = vibrational_states(nuc, F2.mu, 2, F2.v0)
    j = continue_to_tail(j_dk(elec, ds, R, 0.05, F2.ell), R, nuc)
    v = continue_to_tail(v_dk_plus(elec, F2, ds, R, 0.04), R, nuc)
    expected_term2 = c_product(chi[1], v * j * chi[0])
    got = _t_background_term2(chi[0], chi[1], v, j)
    assert abs(got - expected_term2) < 1e-14


def test_t_background_rejects_non_positive_energies(elec, nuc):
    ds = AsymptoticDiscreteState(elec, F2, R_inf=float(nuc.R0))
    R = np.sort(nuc.points[nuc.points.imag == 0.0].real)[::-1]
    eps, chi = vibrational_states(nuc, F2.mu, 2, F2.v0)
    with pytest.raises(ValueError, match="positive"):
        t_background(
            elec,
            nuc,
            F2,
            ds,
            R,
            chi[0],
            chi[0],
            np.zeros(nuc.n, dtype=np.complex128),
            0.05,
            0.0,
        )
