from __future__ import annotations

import numpy as np
import pytest
from qscat.core.dissociation import anion_electronic_states
from qscat.core.grids import electronic_grid, nuclear_grid, segmented_grid
from qscat.core.lcp import lcp_da_cross_section, local_complex_potential
from qscat.core.vibrational import vibrational_states
from qscat.model import F2


def _elec_grids():
    return (
        electronic_grid(r_max=16.0, order=7, n_complex=6, angle_deg=35.0),
        electronic_grid(r_max=16.0, order=7, n_complex=6, angle_deg=44.0),
    )


# The F2 pole walk below is what every coarse-grid test here is built on, and
# it measured ~3.5 s each time it was repeated. It is the same walk on the
# same three grids in all of them, and no test mutates its output, so it runs
# once per module instead of once per test -- the repeated-recomputation
# anti-pattern docs/adr/0005 names. Read-only by contract: copy `Vd`/`Gamma`
# before modifying them. The two `@slow` tests below use their own, different
# grids and deliberately do not share these.
@pytest.fixture(scope="module")
def elec_grids():
    return _elec_grids()


@pytest.fixture(scope="module")
def coarse_nuc():
    return nuclear_grid(r_max=22.0, n_complex=6, quadrature=10)


@pytest.fixture(scope="module")
def coarse_lcp(coarse_nuc, elec_grids):
    ga, gb = elec_grids
    return local_complex_potential(F2, coarse_nuc, ga, gb)


def test_vd_gamma_shapes_and_gamma_nonneg(coarse_nuc, coarse_lcp):
    g_R = coarse_nuc
    Vd, Gamma = coarse_lcp
    assert Vd.shape == (g_R.n,) and Gamma.shape == (g_R.n,)
    assert Vd.dtype == np.complex128 and Gamma.dtype == np.float64
    assert np.all(Gamma >= 0.0)


def test_gamma_closes_and_vd_matches_anion_at_large_R(coarse_nuc, elec_grids, coarse_lcp):
    # As R -> R_inf the pole closes to the bound anion: Gamma -> ~0 and
    # V_d(R_inf) == eps_e (the exact-DA threshold from anion_electronic_states).
    g_R = coarse_nuc
    ga, _gb = elec_grids
    Vd, Gamma = coarse_lcp
    R = g_R.points
    real = R.imag == 0.0
    i_outer = np.flatnonzero(real)[np.argmax(R[real].real)]  # largest real R
    eps_e, _ = anion_electronic_states(ga, F2, g_R.R0, 1)
    assert Gamma[i_outer] < 1e-3  # closed at the edge
    assert abs(Vd[i_outer].real - eps_e[0]) < 5e-3  # == anion asymptote


def test_gamma_positive_in_resonance_region(coarse_nuc, coarse_lcp):
    # At smaller R (inside the crossing) the anion is a real resonance: Gamma>0.
    _Vd, Gamma = coarse_lcp
    R = coarse_nuc.points.real
    band = (R > 1.5) & (R < 2.5)
    assert Gamma[band].max() > 1e-4  # genuine width somewhere


# eMoScat F2 nuclear deck (verbatim from reference/eMoScat/input/F2/grids.txt, 2nd decl)
_F2_NUC_REAL = [(9, 1.8), (1, 2.0), (5, 2.5), (4, 2.596908), (4, 2.7), (40, 10.7)]
_F2_NUC_CPLX = [
    (1, 10.8),
    (1, 11.0),
    (1, 11.5),
    (1, 12.5),
    (1, 14.0),
    (1, 18.0),
    (4, 30.0),
    (2, 101.0),
]


def _f2_fine_grid():
    return segmented_grid(_F2_NUC_REAL, _F2_NUC_CPLX, angle_deg=35.0, quadrature=14)


def _lcp_inputs(g_R, n_vib=3):
    ga, gb = _elec_grids()
    Vd, Gamma = local_complex_potential(F2, g_R, ga, gb)
    eps, chi = vibrational_states(g_R, F2.mu, n_vib, F2.v0)
    return Vd, Gamma, eps, chi


@pytest.fixture(scope="module")
def coarse_lcp_inputs(coarse_nuc, coarse_lcp):
    """`_lcp_inputs` on the coarse grid, reusing the shared pole walk."""
    Vd, Gamma = coarse_lcp
    eps, chi = vibrational_states(coarse_nuc, F2.mu, 3, F2.v0)
    return Vd, Gamma, eps, chi


def test_lcp_da_shape_and_nonneg(coarse_nuc, coarse_lcp_inputs):
    g_R = coarse_nuc  # coarse ok: shape only
    Vd, Gamma, eps, chi = coarse_lcp_inputs
    s = lcp_da_cross_section(g_R, F2.mu, Vd, Gamma, eps, chi, 0, np.array([0.02, 0.03, 0.04]))
    assert s.shape == (3,) and np.all(np.isfinite(s)) and np.all(s >= 0.0)
    assert lcp_da_cross_section(g_R, F2.mu, Vd, Gamma, eps, chi, 0, 0.03).shape == ()  # scalar


def test_lcp_da_closed_channel_is_zero(coarse_nuc, coarse_lcp_inputs):
    # A below-threshold collision energy is closed -> sigma == 0 exactly.
    g_R = coarse_nuc
    Vd, Gamma, eps, chi = coarse_lcp_inputs
    eps_e = float(
        Vd[np.flatnonzero(g_R.points.imag == 0.0)][
            np.argmax(g_R.points.real[g_R.points.imag == 0.0])
        ].real
    )
    E_closed = (eps_e - eps[0]) - 0.05  # well below the DA threshold
    if E_closed > 0:  # F2 is exothermic -> threshold<0, so pick any tiny E:
        E_closed = None
    if E_closed is not None:
        s = lcp_da_cross_section(g_R, F2.mu, Vd, Gamma, eps, chi, 0, np.array([E_closed]))
        assert s[0] == 0.0
    # E<=0 is always closed:
    assert lcp_da_cross_section(g_R, F2.mu, Vd, Gamma, eps, chi, 0, np.array([-0.01]))[0] == 0.0


@pytest.mark.slow
def test_lcp_da_f2_magnitude_matches_exact_order():
    # THE gate: on the fine eMoScat grid with the value extraction, LCP sigma_DA(F2)
    # must land at the exact-2D oracle's ORDER (exact-2D ~1.66 bohr^2 at E=0.03);
    # the LCP is an approximation, so a ~2x band, not exact agreement.
    g_R = _f2_fine_grid()
    Vd, Gamma, eps, chi = _lcp_inputs(g_R)
    s = lcp_da_cross_section(g_R, F2.mu, Vd, Gamma, eps, chi, 0, np.array([0.02, 0.03, 0.04]))
    assert np.all(np.isfinite(s)) and np.all(s >= 0.0)
    assert 0.5 < s[1] < 5.0  # sigma_DA(0.03) ~ 1.47 (exact ~1.66); within ~2x


def test_lcp_da_return_wavefunction_parity_and_shape(coarse_nuc, coarse_lcp_inputs):
    # #2: return_wavefunction exposes the 1-D nuclear resolvent psi_sc(R) without
    # changing sigma; None for a closed channel.
    g_R = coarse_nuc
    Vd, Gamma, eps, chi = coarse_lcp_inputs
    E = np.array([0.02, 0.03])
    s_plain = lcp_da_cross_section(g_R, F2.mu, Vd, Gamma, eps, chi, 0, E)
    s2, psis = lcp_da_cross_section(g_R, F2.mu, Vd, Gamma, eps, chi, 0, E, return_wavefunction=True)
    assert np.array_equal(s_plain, s2)  # exact
    assert isinstance(psis, list) and len(psis) == 2
    for psi in psis:  # F2 exothermic -> both open
        assert psi is not None and psi.shape == (g_R.n,) and psi.dtype == np.complex128
    # scalar E -> a single array; E<=0 -> None (closed)
    _s1, psi1 = lcp_da_cross_section(
        g_R, F2.mu, Vd, Gamma, eps, chi, 0, 0.03, return_wavefunction=True
    )
    assert psi1 is not None and psi1.shape == (g_R.n,)
    _, psi0 = lcp_da_cross_section(
        g_R, F2.mu, Vd, Gamma, eps, chi, 0, -0.01, return_wavefunction=True
    )
    assert psi0 is None


def test_lcp_ve_matches_dense_reference(coarse_nuc, coarse_lcp_inputs):
    # Differential oracle: the same driven equation assembled densely and
    # solved with np.linalg.solve -- the projects toy model's formulation.
    from qscat.core.lcp import lcp_ve_cross_section
    from qscat.dvr import kinetic

    g_R = coarse_nuc
    Vd, Gamma, eps, chi = coarse_lcp_inputs
    E, vprimes = 0.05, [0, 1, 2]
    sigma = lcp_ve_cross_section(g_R, F2.mu, Vd, Gamma, eps, chi, 0, vprimes, E)

    doorway = np.sqrt(Gamma / (2.0 * np.pi))[None, :] * chi
    H = kinetic(g_R, F2.mu) + np.diag(Vd - 0.5j * Gamma)
    M = (E + eps[0]) * np.eye(g_R.n, dtype=np.complex128) - H
    xi = np.linalg.solve(M, doorway[0])
    e_tot = E + eps[0]
    expected = np.zeros(len(vprimes))
    for k, vp in enumerate(vprimes):
        if e_tot - eps[vp] > 0.0:
            S = np.dot(doorway[vp], xi)  # c-product: no conjugate
            expected[k] = 4.0 * np.pi**3 * np.abs(S) ** 2 / (2.0 * E)
    np.testing.assert_allclose(sigma, expected, rtol=1e-9)


def test_lcp_ve_shapes_closed_channels_and_sweep_reuse(coarse_nuc, coarse_lcp_inputs):
    from qscat.core.lcp import lcp_ve_cross_section

    g_R = coarse_nuc
    Vd, Gamma, eps, chi = coarse_lcp_inputs
    # scalar E -> (len(vprimes),); array E -> (len(E), len(vprimes))
    s1 = lcp_ve_cross_section(g_R, F2.mu, Vd, Gamma, eps, chi, 0, [0, 1], 0.05)
    assert s1.shape == (2,) and np.all(np.isfinite(s1)) and np.all(s1 >= 0.0)
    E = np.array([0.02, 0.05])
    s2 = lcp_ve_cross_section(g_R, F2.mu, Vd, Gamma, eps, chi, 0, [0, 1], E)
    assert s2.shape == (2, 2)
    # the refactor sweep must equal fresh per-energy solves
    for i, e in enumerate(E):
        fresh = lcp_ve_cross_section(g_R, F2.mu, Vd, Gamma, eps, chi, 0, [0, 1], float(e))
        np.testing.assert_allclose(s2[i], fresh, rtol=1e-12)
    # E <= 0 is closed -> exactly zero; a closed v' channel -> exactly zero
    assert np.all(lcp_ve_cross_section(g_R, F2.mu, Vd, Gamma, eps, chi, 0, [0, 1], -0.01) == 0.0)
    s3 = lcp_ve_cross_section(g_R, F2.mu, Vd, Gamma, eps, chi, 0, [2], 1e-4)
    assert s3[0] == 0.0  # eps[2]-eps[0] >> 1e-4: channel closed


def test_lcp_ve_return_wavefunction_parity(coarse_nuc, coarse_lcp_inputs):
    from qscat.core.lcp import lcp_ve_cross_section

    g_R = coarse_nuc
    Vd, Gamma, eps, chi = coarse_lcp_inputs
    E = np.array([0.02, 0.05])
    s_plain = lcp_ve_cross_section(g_R, F2.mu, Vd, Gamma, eps, chi, 0, [0, 1], E)
    s2, xis = lcp_ve_cross_section(
        g_R, F2.mu, Vd, Gamma, eps, chi, 0, [0, 1], E, return_wavefunction=True
    )
    assert np.array_equal(s_plain, s2)
    assert isinstance(xis, list) and len(xis) == 2
    for xi in xis:
        assert xi is not None and xi.shape == (g_R.n,) and xi.dtype == np.complex128
    _s, xi0 = lcp_ve_cross_section(
        g_R, F2.mu, Vd, Gamma, eps, chi, 0, [0], -0.01, return_wavefunction=True
    )
    assert xi0 is None


def test_resonance_eigenstate_at_peak_width(coarse_nuc, elec_grids, coarse_lcp):
    # #1: the resonance eigenstate at the width peak -- a genuine resonance
    # (Gamma>0), a c-product-normalized electronic eigenfunction, and Re(E_pole)
    # consistent with local_complex_potential's V_d at that R.
    from qscat.core.lcp import resonance_eigenstate_at_peak_width
    from qscat.linalg import c_product

    g_R = coarse_nuc
    ga, gb = elec_grids
    R_star, E_pole, phi = resonance_eigenstate_at_peak_width(F2, g_R, ga, gb)

    assert 1.0 < R_star < 3.0  # F2 resonance is inside the crossing
    assert phi.shape == (ga.n,) and phi.dtype == np.complex128
    assert -2.0 * E_pole.imag > 1e-4  # genuine width (Gamma>0)
    # Re(E_pole) reproduces local_complex_potential's V_d at R_star
    Vd, _ = coarse_lcp
    j = int(
        np.flatnonzero(g_R.points.imag == 0.0)[
            np.argmin(np.abs(g_R.points.real[g_R.points.imag == 0.0] - R_star))
        ]
    )
    assert abs(E_pole.real - float(Vd[j].real)) < 5e-3
    # c-product normalized over the electronic real region
    p = phi.copy()
    p[~(ga.real_points <= ga.R0)] = 0.0
    assert abs(complex(c_product(p, p)) - 1.0) < 1e-6
