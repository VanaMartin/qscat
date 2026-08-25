from __future__ import annotations

import numpy as np
import pytest
from qscat.core.dissociation import anion_electronic_states, da_cross_section, v_dr_diag
from qscat.core.grids import electronic_grid, nuclear_grid
from qscat.core.vibrational import vibrational_states
from qscat.dvr import TensorGrid
from qscat.model import F2, N2, NO

_eg = electronic_grid
_ng = nuclear_grid


def _eps0(model):
    g_R = nuclear_grid(r_max=22.0, n_complex=8, quadrature=12)
    eps, _ = vibrational_states(g_R, model.mu, 3, model.v0)
    return eps[0], g_R.R0


@pytest.mark.parametrize("model", [N2, NO, F2], ids=["N2", "NO", "F2"])
def test_one_bound_anion_state_real(model):
    g_r = electronic_grid(r_max=16.0, order=7, n_complex=6)
    _, R0 = _eps0(model)
    eps_e, phi = anion_electronic_states(g_r, model, R0, n_states=1)
    assert eps_e.shape == (1,) and phi.shape == (1, g_r.n)
    # c-product self-normalized over the real region ~ 1
    real = g_r.real_points <= g_r.R0
    p = phi[0].copy()
    p[~real] = 0.0
    assert abs(complex(p @ p) - 1.0) < 1e-6


def test_thresholds_have_correct_signs():
    # threshold(E_coll) = eps_e - eps[0]; F2 exothermic (<0), N2 closed (>0.3),
    # NO opens above its resonance (~0.17). No independent data -> sign/band gate.
    g_r = electronic_grid(r_max=16.0, order=7, n_complex=6)
    thr = {}
    for name, model in (("N2", N2), ("NO", NO), ("F2", F2)):
        eps0, R0 = _eps0(model)
        eps_e, _ = anion_electronic_states(g_r, model, R0, 1)
        thr[name] = float(eps_e[0]) - eps0
    assert thr["F2"] < 0.0  # exothermic: DA open at all E>0
    assert thr["N2"] > 0.3  # closed in the measurement window
    assert 0.10 < thr["NO"] < 0.25  # opens above the resonance


def test_raises_when_too_many_states_requested():
    g_r = electronic_grid(r_max=16.0, order=7, n_complex=6)
    _, R0 = _eps0(F2)
    with pytest.raises(ValueError):
        anion_electronic_states(g_r, F2, R0, n_states=50)


def _tgrid():
    return TensorGrid(
        [
            electronic_grid(r_max=14.0, order=6, n_complex=4),
            nuclear_grid(r_max=20.0, n_complex=4, quadrature=8),
        ]
    )


def test_v_dr_shape_and_dtype():
    tg = _tgrid()
    vdr = v_dr_diag(tg, F2)
    assert vdr.shape == (tg.size,) and vdr.dtype == np.complex128


def test_v_dr_equals_definition_pointwise():
    tg = _tgrid()
    model = F2
    R_inf = tg.grids[1].R0
    pts_r, pts_R = tg.points()  # (n_r,1), (1,n_R)
    expect = (
        model.interaction_diag(tg)
        + np.broadcast_to(model.v0(pts_R), tg.shape).ravel()
        - np.broadcast_to(model.v_int(pts_r, R_inf), tg.shape).ravel()
    )
    assert np.allclose(v_dr_diag(tg, model), expect, rtol=0, atol=1e-14)


def test_v_dr_tends_to_v0_at_large_R():
    # Where R is near R_inf, V_int(r,R) ~ V_int(r,R_inf), so V_DR ~ v0(R).
    tg = _tgrid()
    model = F2
    vdr = v_dr_diag(tg, model).reshape(tg.shape)  # (n_r, n_R)
    pts_R = tg.points()[1].ravel()
    j = int(np.argmin(np.abs(pts_R - tg.grids[1].R0)))  # column nearest R_inf
    v0_col = np.broadcast_to(model.v0(tg.points()[1]), tg.shape)[:, j]
    assert np.allclose(vdr[:, j], v0_col, rtol=0, atol=1e-10)


def _working():
    tg = TensorGrid(
        [_eg(r_max=16.0, order=8, n_complex=6), _ng(r_max=22.0, n_complex=8, quadrature=12)]
    )
    return tg


@pytest.mark.slow
def test_da_shape_scalar_and_array():
    tg = _working()
    eps, chi = vibrational_states(tg.grids[1], F2.mu, 3, F2.v0)
    s1 = da_cross_section(tg, F2, eps, chi, 0, 0.05)
    assert s1.shape == (1,)
    sN = da_cross_section(tg, F2, eps, chi, 0, np.array([0.05, 0.10]))
    assert sN.shape == (2, 1)
    assert np.all(sN >= 0.0) and np.all(np.isfinite(sN))


@pytest.mark.slow
def test_da_return_wavefunction_parity_and_shape():
    # return_wavefunction must not change sigma (byte-identical), and must hand
    # back the driven Psi+ (full tensor field) per energy -- the R6 hook that
    # lets qscat_run snapshot/animate the DA scattering state.
    tg = _working()
    eps, chi = vibrational_states(tg.grids[1], F2.mu, 3, F2.v0)
    E = np.array([0.05, 0.10])
    s_plain = da_cross_section(tg, F2, eps, chi, 0, E)
    s2, psis = da_cross_section(tg, F2, eps, chi, 0, E, return_wavefunction=True)
    assert np.array_equal(s_plain, s2)  # exact, not approx
    assert isinstance(psis, list) and len(psis) == 2
    for psi in psis:
        assert psi is not None and psi.shape == (tg.size,) and psi.dtype == np.complex128
    # scalar E -> a single array (None only below threshold)
    s1, psi1 = da_cross_section(tg, F2, eps, chi, 0, 0.05, return_wavefunction=True)
    assert psi1 is not None and psi1.shape == (tg.size,)
    # E <= 0 -> closed, Psi+ is None
    _, psi0 = da_cross_section(tg, F2, eps, chi, 0, -0.01, return_wavefunction=True)
    assert psi0 is None


@pytest.mark.slow
def test_n2_channel_closed_is_zero():
    # N2's DA threshold is +0.5 Ha -> sigma_DA == 0 across the whole VE window.
    tg = _working()
    eps, chi = vibrational_states(tg.grids[1], N2.mu, 3, N2.v0)
    E = np.array([0.04, 0.10, 0.18])
    s = da_cross_section(tg, N2, eps, chi, 0, E)
    assert np.all(s == 0.0)


@pytest.mark.slow
def test_f2_exothermic_da_is_positive():
    # F2 DA is open at all E>0; expect a nonzero, finite sigma in its resonance
    # window. No golden number (no independent DA data) -- positivity + soft
    # unitarity only.
    tg = _working()
    eps, chi = vibrational_states(tg.grids[1], F2.mu, 3, F2.v0)
    E = np.array([0.02, 0.03, 0.04])
    s = da_cross_section(tg, F2, eps, chi, 0, E)[:, 0]
    assert np.all(np.isfinite(s)) and np.all(s >= 0.0)
    assert s.max() > 0.0
    # soft unitarity: sigma_DA <= a few * pi/(2E) (partial-wave cap, generous
    # band for the under-resolved fast outgoing wave; see the convergence note)
    cap = np.pi / (2.0 * E)
    assert np.all(s < 50.0 * cap)


def _h2p_proxy():
    # small ionic proxy: electronic to ~60 bohr (holds a couple Rydberg states +
    # the incident), nuclear to ~14. Big enough for well-posedness, laptop-fast.
    from qscat.core.grids import electronic_grid, nuclear_grid
    from qscat.dvr import TensorGrid

    return TensorGrid(
        [
            electronic_grid(r_max=60.0, order=8, n_complex=6),
            nuclear_grid(r_max=22.0, n_complex=6, quadrature=10),
        ]
    )


@pytest.mark.slow
def test_dr_wellposed_and_threshold_ordered():
    from qscat.core.dissociation import dr_cross_section
    from qscat.core.vibrational import vibrational_states
    from qscat.model import H2P

    tg = _h2p_proxy()
    eps, chi = vibrational_states(tg.grids[1], H2P.mu, 3, H2P.v0)
    E = np.array([0.01, 0.03])
    s = dr_cross_section(tg, H2P, eps, chi, 0, E, n_channels=2)
    assert s.shape == (2, 2)
    assert np.all(np.isfinite(s)) and np.all(s >= 0.0)

    # R6: return_wavefunction hands back the driven Psi+ without changing sigma.
    s2, psis = dr_cross_section(tg, H2P, eps, chi, 0, E, n_channels=2, return_wavefunction=True)
    assert np.array_equal(s, s2)
    assert isinstance(psis, list) and len(psis) == 2
    for psi in psis:
        assert psi is not None and psi.shape == (tg.size,) and psi.dtype == np.complex128


def _dr_t_matrix_conjugated_channel0(tg, E: float) -> complex:
    """Mirrors `dr_cross_section`'s internals for ONE channel (n=0) but with a
    CONJUGATED T-matrix dot (eMoScat's `zdotc` convention) instead of the shipped
    c-product -- the reference the convention check compares against."""
    import scipy.sparse as sp
    from qscat.core.channels import channel_vector
    from qscat.linalg import SparseLU
    from qscat.model import H2P
    from qscat.special import riccati_bessel_en_mass

    model = H2P
    mu = model.mu
    g_R = tg.grids[1]
    eps, chi = vibrational_states(tg.grids[1], model.mu, 3, model.v0)
    eps_ryd, phi_ryd = anion_electronic_states(
        g_r=tg.grids[0], model=model, R_inf=g_R.R0, n_states=1
    )
    v_dr = v_dr_diag(tg, model)
    mask = tg.real_mask()
    sqrt_w_R = tg.sqrt_weights()[1].ravel()

    ident = sp.identity(tg.size, format="csc", dtype=np.complex128)
    e_tot = E + eps[0]
    lu = SparseLU((e_tot * ident - model.hamiltonian(tg)).tocsc())
    k = float(np.sqrt(2.0 * E))
    psi_i = channel_vector(tg, k, chi[0], model.ell, charge=model.charge)
    psi_plus = psi_i + lu.solve(model.interaction_diag(tg) * psi_i)
    v_psi = v_dr * psi_plus

    e_dr = e_tot - eps_ryd[0]
    assert e_dr > 0.0, "test picked a closed n=0 Rydberg channel"
    k_r = float(np.sqrt(2.0 * mu * e_dr))
    y_coeff = riccati_bessel_en_mass(g_R.real_points, k_r, 0, mu) * sqrt_w_R
    phi_f = tg.outer([phi_ryd[0], y_coeff])
    phi_f[~mask] = 0.0
    return complex(np.sum(np.conj(phi_f[mask]) * v_psi[mask]))


@pytest.mark.slow
def test_dr_cproduct_matches_conjugated_dot_on_proxy():
    """The CONVENTION check (promoted from the retired validation/h2plus TD-DR
    driver, removed in the qscat-run consolidation): `dr_cross_section`'s c-product T-matrix (no conjugate,
    the ECS-correct choice) agrees with eMoScat's conjugated-dot (`zdotc`)
    convention to <1e-2 relative on the proxy -- the rotated-nuclear-tail
    contribution is negligible there, so the convention question is settled."""
    from qscat.core.dissociation import dr_cross_section
    from qscat.model import H2P

    tg = _h2p_proxy()
    E = 0.03
    eps, chi = vibrational_states(tg.grids[1], H2P.mu, 3, H2P.v0)
    sigma_c0 = float(dr_cross_section(tg, H2P, eps, chi, 0, E, n_channels=1)[0])
    t_conj = _dr_t_matrix_conjugated_channel0(tg, E)
    sigma_conj0 = 4.0 * np.pi**3 * abs(t_conj) ** 2 / (2.0 * E)

    assert sigma_c0 > 0.0 and np.isfinite(sigma_c0)
    assert sigma_conj0 > 0.0 and np.isfinite(sigma_conj0)
    assert abs(sigma_c0 - sigma_conj0) / sigma_c0 < 1e-2


@pytest.mark.slow
def test_dr_amplitude_reproduces_the_returned_sigma() -> None:
    """The amplitude and sigma must not drift apart: sigma is 4pi^3|t|^2/2E."""
    from qscat.core.dissociation import dr_cross_section
    from qscat.model import H2P

    tg = _h2p_proxy()
    eps, chi = vibrational_states(tg.grids[1], H2P.mu, 3, H2P.v0)
    energies = np.array([0.012, 0.014])
    sigma, amp = dr_cross_section(
        tg, H2P, eps, chi, 0, energies, n_channels=2, return_amplitude=True
    )
    assert amp.shape == sigma.shape
    assert amp.dtype == np.complex128
    recomputed = 4.0 * np.pi**3 * np.abs(amp) ** 2 / (2.0 * energies[:, None])
    open_channels = sigma > 0.0
    assert np.allclose(recomputed[open_channels], sigma[open_channels], rtol=1e-12)
    # A closed channel contributes exactly zero to both.
    assert np.all(amp[~open_channels] == 0.0)


@pytest.mark.slow
def test_dr_amplitude_composes_with_the_wavefunction_return() -> None:
    from qscat.core.dissociation import dr_cross_section
    from qscat.model import H2P

    tg = _h2p_proxy()
    eps, chi = vibrational_states(tg.grids[1], H2P.mu, 3, H2P.v0)
    sigma, psi, amp = dr_cross_section(
        tg,
        H2P,
        eps,
        chi,
        0,
        0.012,
        n_channels=2,
        return_wavefunction=True,
        return_amplitude=True,
    )
    assert psi is not None and psi.shape == (tg.size,)
    assert amp.shape == sigma.shape


@pytest.mark.slow
def test_dr_amplitude_matches_conjugated_oracle_value_and_phase() -> None:
    """`return_amplitude`'s `t` must equal the file's INDEPENDENT oracle
    (`_dr_t_matrix_conjugated_channel0`, eMoScat's conjugated-dot `zdotc`
    convention) in VALUE, not just modulus -- this is a stronger, independent
    check than `test_dr_amplitude_reproduces_the_returned_sigma`, which only
    confirms `amp` and `sigma` are filled from the same `t` inside
    `dr_cross_section` and so cannot see a sign/phase bug in `t` itself.

    The two conventions differ, in general, by whether `phi_f` (the exit
    channel) is conjugated in the T-matrix dot. They coincide here because
    `phi_f` is REAL-valued in the surviving (non-ECS-tail) region masked by
    `real_mask()` -- it is built from a genuine bound Rydberg electronic
    state and a real Riccati-Bessel radial factor, so conjugating it is a
    no-op there. `test_dr_cproduct_matches_conjugated_dot_on_proxy` already
    established this numerically for sigma (to a loose 1e-2 bound, "the
    convention question is settled"); this test pins the SAME agreement for
    the complex amplitude itself at a tight tolerance (~1e-13 relative in
    practice, limited only by the sparse solve/normalization roundoff, not
    by any convention gap) -- tight enough that a spurious conjugation of
    `amp` (which would flip its imaginary part, moving it ~1e3x further from
    the oracle, per `test_dr_cproduct_matches_conjugated_dot_on_proxy`'s own
    `_dr_t_matrix_conjugated_channel0` reference) or a stray overall sign
    would fail it outright.
    """
    from qscat.core.dissociation import dr_cross_section
    from qscat.model import H2P

    tg = _h2p_proxy()
    E = 0.03
    eps, chi = vibrational_states(tg.grids[1], H2P.mu, 3, H2P.v0)
    sigma, amp = dr_cross_section(
        tg, H2P, eps, chi, 0, E, n_channels=1, return_amplitude=True
    )
    assert sigma[0] > 0.0, "test picked a closed channel"

    t = complex(amp[0])
    t_oracle = _dr_t_matrix_conjugated_channel0(tg, E)

    # VALUE and PHASE, not just modulus: a bare np.isclose on abs() would
    # pass under a conjugation or sign error, which is exactly what this
    # test exists to catch.
    assert np.isclose(t, t_oracle, rtol=1e-8, atol=0.0)
    # The conjugated relation must NOT also hold -- if it did, this test
    # could not distinguish the correct convention from its conjugate, i.e.
    # it would not be capable of failing on a spurious conjugation bug.
    assert not np.isclose(t, t_oracle.conjugate(), rtol=1e-3, atol=0.0)
