"""Small-proxy validation gate for the H2+ DR driver (`validation.h2plus.dr`).

Two `@pytest.mark.slow` tests on a laptop-feasible grid (NOT
`config.proxy_grid` -- 192k unknowns is too heavy for a routine gate; this
uses the same small size as `qscat`'s own `test_dr_wellposed_and_threshold_
ordered` proxy):

1. `test_dr_wellposed_and_threshold` -- well-posedness (finite, >=0, shape)
   plus a genuine below-threshold closed-channel check.
2. `test_cproduct_matches_conjugated_dot_on_proxy` -- the CONVENTION check:
   `dr_cross_section`'s c-product T-matrix (no conjugate) vs a locally
   reimplemented conjugated-dot T-matrix (eMoScat's `zdotc` convention)
   agree to <1e-2 relative on the proxy -- see `qscat.core.dissociation`'s
   module docstring and `docs/physics/n2-2d-cross-section.md` for the
   c-product-vs-conjugate-dot convention background.
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp
from qscat.core.channels import channel_vector
from qscat.core.dissociation import anion_electronic_states, v_dr_diag
from qscat.core.grids import electronic_grid, nuclear_grid
from qscat.core.vibrational import vibrational_states
from qscat.dvr import TensorGrid
from qscat.linalg import SparseLU
from qscat.model import H2P
from qscat.special import riccati_bessel_en_mass

from validation.h2plus.dr import compute_dr


def _small() -> TensorGrid:
    return TensorGrid(
        [
            electronic_grid(r_max=60.0, order=8, n_complex=6),
            nuclear_grid(r_max=22.0, n_complex=6, quadrature=10),
        ]
    )


@pytest.mark.slow
def test_dr_wellposed_and_threshold():
    tg = _small()
    energies = np.array([0.01, 0.03])
    energies_out, sigma = compute_dr(tg, energies=energies, n_channels=2)

    assert np.array_equal(energies_out, energies)
    assert sigma.shape == (2, 2)
    assert np.all(np.isfinite(sigma))
    assert np.all(sigma >= 0.0)

    # Genuine threshold check: force a below-threshold Rydberg channel.
    # eps_ryd are the anion/Rydberg electronic-state energies (same solver
    # DA uses); a collision energy E such that E_tot = E + eps[0] < eps_ryd[n]
    # must give sigma[n] == 0 exactly (dr_cross_section's `e_dr <= 0` guard).
    # On this proxy grid the n=0,1 Rydberg channels are open across the whole
    # [0.01, 0.03] window (eps_ryd - eps[0] ~ -1.29, -0.028 Ha) but the n=2
    # channel opens only above ~0.043 Ha -- verified numerically at proxy
    # build time (measured: eps_ryd - eps[0] = [-1.2875, -0.0276, 0.0426]) --
    # so requesting 3 channels gives a genuinely closed n=2 channel at both
    # probe energies without any faked/forced construction.
    eps, chi = vibrational_states(tg.grids[1], H2P.mu, 3, H2P.v0)
    eps_ryd, _ = anion_electronic_states(
        g_r=tg.grids[0], model=H2P, R_inf=tg.grids[1].R0, n_states=3
    )
    assert np.all(energies < eps_ryd[2] - eps[0]), (
        "test assumption broken: n=2 Rydberg channel expected closed across "
        f"{energies!r} (threshold {eps_ryd[2] - eps[0]!r} Ha)"
    )
    _, sigma3 = compute_dr(tg, energies=energies, n_channels=3)
    assert np.all(sigma3[:, 2] == 0.0)


def _t_matrix_conjugated_channel0(tg: TensorGrid, E: float) -> complex:
    """Mirrors `dr_cross_section`'s internals for ONE channel (n=0), but with
    a CONJUGATED T-matrix dot (`np.sum(np.conj(Phi[mask]) * (V_DR*psi_plus)
    [mask])`) instead of the shipped c-product. Returns the raw T-matrix
    (not sigma) so the caller can square/compare on equal footing with the
    c-product T.
    """
    model = H2P
    mu = model.mu
    g_R = tg.grids[1]
    R_inf = g_R.R0

    eps, chi = vibrational_states(tg.grids[1], model.mu, 3, model.v0)
    eps_ryd, phi_ryd = anion_electronic_states(
        g_r=tg.grids[0], model=model, R_inf=R_inf, n_states=1
    )
    v_dr = v_dr_diag(tg, model)
    mask = tg.real_mask()
    sqrt_w_R = tg.sqrt_weights()[1].ravel()

    H = model.hamiltonian(tg)
    v_diag = model.interaction_diag(tg)
    ident = sp.identity(tg.size, format="csc", dtype=np.complex128)

    e_tot = E + eps[0]
    a = (e_tot * ident - H).tocsc()
    lu = SparseLU(a)

    k = float(np.sqrt(2.0 * E))
    psi_i = channel_vector(tg, k, chi[0], model.ell, charge=model.charge)
    psi_plus = psi_i + lu.solve(v_diag * psi_i)
    v_psi = v_dr * psi_plus

    e_dr = e_tot - eps_ryd[0]
    assert e_dr > 0.0, "test picked a closed n=0 Rydberg channel"
    k_r = float(np.sqrt(2.0 * mu * e_dr))
    y_coeff = riccati_bessel_en_mass(g_R.real_points, k_r, 0, mu) * sqrt_w_R
    phi_f = tg.outer([phi_ryd[0], y_coeff])
    phi_f[~mask] = 0.0

    driven = v_psi
    t_conj = complex(np.sum(np.conj(phi_f[mask]) * driven[mask]))
    return t_conj


@pytest.mark.slow
def test_cproduct_matches_conjugated_dot_on_proxy():
    tg = _small()
    E = 0.03

    # (a) shipped c-product route, one channel.
    eps, chi = vibrational_states(tg.grids[1], H2P.mu, 3, H2P.v0)
    from qscat.core.dissociation import dr_cross_section

    sigma_c = dr_cross_section(tg, H2P, eps, chi, 0, E, n_channels=1)
    sigma_c0 = float(sigma_c[0])

    # (b) locally reimplemented conjugated-dot route, same channel.
    t_conj = _t_matrix_conjugated_channel0(tg, E)
    sigma_conj0 = 4.0 * np.pi**3 * abs(t_conj) ** 2 / (2.0 * E)

    assert sigma_c0 > 0.0 and np.isfinite(sigma_c0)
    assert sigma_conj0 > 0.0 and np.isfinite(sigma_conj0)

    rel_diff = abs(sigma_c0 - sigma_conj0) / sigma_c0
    assert rel_diff < 1e-2, (
        f"c-product sigma={sigma_c0!r} vs conjugated-dot sigma={sigma_conj0!r}, "
        f"rel_diff={rel_diff!r}"
    )
