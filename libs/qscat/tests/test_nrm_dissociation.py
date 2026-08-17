"""Tests for the NRM nuclear solve and sigma_DA (PRA 77 Eq. 52, 54)."""

from __future__ import annotations

import numpy as np
import pytest
from qscat.core.grids import electronic_grid, segmented_grid
from qscat.core.lcp import lcp_da_cross_section, local_complex_potential
from qscat.core.nrm.coupling import v_dk_plus
from qscat.core.nrm.discrete_state import AsymptoticDiscreteState, electronic_hamiltonian
from qscat.core.nrm.dissociation import (
    da_sigma_from_psi,
    nrm_da_cross_section,
    solve_nuclear,
)
from qscat.core.nrm.ingredients import nrm_ingredients
from qscat.core.nrm.nonlocal_potential import continue_to_tail, nonlocal_operator
from qscat.core.nrm.scattering import scattering_state
from qscat.core.vibrational import vibrational_states
from qscat.dvr import kinetic
from qscat.model import F2


@pytest.fixture(scope="module")
def nuc():
    return segmented_grid(
        ((9, 1.8), (1, 2.0), (5, 2.5), (4, 2.7), (20, 10.7)),
        ((1, 11.0), (1, 12.5), (1, 14.0), (3, 30.0)),
        angle_deg=45.0,
        quadrature=14,
    )


@pytest.fixture(scope="module")
def elec():
    return electronic_grid(r_max=16.0, order=8, n_complex=6)


def test_solve_nuclear_reduces_to_a_plain_resolvent(nuc):
    """With F = 0 the solve is (E - T - V_d)^-1 rhs, nothing more."""
    v_d = F2.v0(nuc.points).astype(np.complex128)
    rhs = np.zeros(nuc.n, dtype=np.complex128)
    rhs[5] = 1.0
    f = np.zeros((nuc.n, nuc.n), dtype=np.complex128)
    psi = solve_nuclear(nuc, F2.mu, v_d, f, rhs, e_total=0.05)
    a = 0.05 * np.eye(nuc.n, dtype=np.complex128) - kinetic(nuc, F2.mu) - np.diag(v_d)
    assert np.allclose(a @ psi, rhs, rtol=1e-8, atol=1e-12)


def test_da_sigma_matches_the_lcp_flux_formula(nuc):
    """Eq. (54) and lcp_da_cross_section's 4 pi^3 |S|^2/2E are the same number."""
    rng = np.random.default_rng(0)
    psi = (rng.normal(size=nuc.n) + 1j * rng.normal(size=nuc.n)).astype(np.complex128)
    e_kin, eps_e, e_tot = 0.03, -0.08, 0.03 - 0.05
    sigma = da_sigma_from_psi(nuc, F2.mu, psi, e_tot, eps_e, e_kin)

    real_idx = np.flatnonzero(nuc.points.imag == 0.0)
    b = int(real_idx[np.argmax(nuc.points[real_idx].real)])
    val = psi[b] / np.sqrt(complex(nuc.weights[b]))
    k_r = float(np.sqrt(2.0 * F2.mu * (e_tot - eps_e)))
    s_da = np.sqrt(k_r / (2.0 * np.pi * F2.mu)) * val
    assert abs(sigma - 4.0 * np.pi**3 * abs(s_da) ** 2 / (2.0 * e_kin)) < 1e-12


def test_closed_channel_gives_zero(nuc):
    psi = np.ones(nuc.n, dtype=np.complex128)
    assert da_sigma_from_psi(nuc, F2.mu, psi, -0.2, -0.08, 0.03) == 0.0


@pytest.mark.slow
def test_local_limit_reproduces_the_lcp(nuc, elec):
    """GATE (validation check 2): the LCP local-limit bridge.

    F -> -(i/2) Gamma(E,R) delta(R-R') collapses Eq. (52) onto the LCP nuclear
    equation. Driving solve_nuclear with that diagonal F and the LCP's own real
    doorway must reproduce lcp_da_cross_section to solver tolerance -- a
    differential oracle against already-validated code, isolating the nuclear
    solve and the flux extraction from the coupling's phase.
    """
    elec_b = electronic_grid(r_max=16.0, order=8, n_complex=6, angle_deg=40.0)
    vd, gamma = local_complex_potential(F2, nuc, elec, elec_b)
    eps, chi = vibrational_states(nuc, F2.mu, 4, F2.v0)
    v_init, e_kin = 0, 0.03
    e_tot = e_kin + eps[v_init]

    reference = float(lcp_da_cross_section(nuc, F2.mu, vd, gamma, eps, chi, v_init, e_kin))

    f_local = np.diag(-0.5j * gamma).astype(np.complex128)
    doorway = np.sqrt(gamma / (2.0 * np.pi)).astype(np.complex128) * chi[v_init]
    psi = solve_nuclear(nuc, F2.mu, vd, f_local, doorway, e_tot)
    eps_e = float(vd[_boundary(nuc)].real)
    got = da_sigma_from_psi(nuc, F2.mu, psi, e_tot, eps_e, e_kin)

    assert reference > 0.0, "the LCP reference is zero -- pick an open energy"
    assert abs(got - reference) < 1e-8 * reference


@pytest.mark.slow
def test_nrm_da_cross_section_is_within_a_loose_band_of_the_lcp(nuc, elec):
    """sigma_NRM lands within a factor of 10 of sigma_LCP, same grid/energy.

    "Finite and positive" alone is vacuous here: it survives e_total being
    silently replaced by e_kin (a 4.6e-32x error, still positive and finite).
    That is NOT exercised by the local-limit bridge gate, which drives
    solve_nuclear/da_sigma_from_psi directly and never calls
    nrm_da_cross_section's own e_total assembly or its nonlocal_operator/
    ing/rhs wiring at all. A loose 0.1-10x band against the (independently
    validated) LCP needs no oracle to justify -- the two models are not
    expected to agree closely (F is nonlocal and energy-dependent; Gamma is
    F's local, energy-independent-at-fixed-E limit) -- and it is tight enough
    to catch the e_total-class defect above (orders of magnitude) outright.
    Measured ratio ~0.56 (F2, AsymptoticDiscreteState, E=0.03, n_states=40).

    This band does NOT cover a zeroed/dropped F: measured, zeroing F changes
    sigma by 5.15x at this same (grid, energy), which still lands inside
    0.1-10x. Tightening the band to exclude it would be arbitrary (the
    honest margin at E=0.03 is ~0.15-1.8x) and fragile (the sigma_NRM/
    sigma_LCP ratio is itself strongly energy-dependent: 0.79/0.56/0.31 at
    E=0.02/0.03/0.05). `test_f_is_load_bearing_in_nrm_da_cross_section`
    below covers F directly instead.
    """
    ds = AsymptoticDiscreteState(elec, F2, R_inf=elec.R0)
    eps, chi = vibrational_states(nuc, F2.mu, 4, F2.v0)
    v_init, e_kin = 0, 0.03
    sigma_nrm = float(nrm_da_cross_section(nuc, elec, F2, ds, eps, chi, v_init, e_kin, n_states=40))
    assert np.isfinite(sigma_nrm) and sigma_nrm > 0.0

    elec_b = electronic_grid(r_max=16.0, order=8, n_complex=6, angle_deg=40.0)
    vd, gamma = local_complex_potential(F2, nuc, elec, elec_b)
    sigma_lcp = float(lcp_da_cross_section(nuc, F2.mu, vd, gamma, eps, chi, v_init, e_kin))

    ratio = sigma_nrm / sigma_lcp
    assert 0.1 < ratio < 10.0, f"sigma_NRM/sigma_LCP = {ratio:.3g}, outside the loose band"


@pytest.mark.slow
def test_f_is_load_bearing_in_nrm_da_cross_section(nuc, elec):
    """Zeroing F changes sigma_DA by a large factor -- F is not numerically inert.

    Rebuilds nrm_da_cross_section's own pipeline (ingredients, V_d, the RHS,
    F itself) explicitly and compares solve_nuclear/da_sigma_from_psi driven
    with the real F against the same call with F replaced by zeros. Needed
    because the loose LCP-comparison band above is blind to F being dropped
    (5.15x still lands inside its 0.1-10x window).

    MUST run at E >= 0.03 (do not "simplify" this to a cheaper energy): F's
    effect on sigma is itself energy-dependent and weak near threshold --
    only ~12% at E=0.02 -- which would make a >2x assertion fragile. At
    E=0.03 the measured effect is 5.15x, comfortably clear of the 2x bar.
    """
    ds = AsymptoticDiscreteState(elec, F2, R_inf=elec.R0)
    eps, chi = vibrational_states(nuc, F2.mu, 4, F2.v0)
    v_init, e_kin = 0, 0.03
    e_total = e_kin + eps[v_init]

    real = nuc.points.imag == 0.0
    R_desc = np.sort(nuc.points[real].real)[::-1]
    ing = nrm_ingredients(elec, F2, ds, R_desc)
    v_d_full = continue_to_tail(ing.v_d_discrete, ing.R, nuc)
    eps_e = float(v_d_full[_boundary(nuc)].real)
    v_dk = v_dk_plus(elec, F2, ds, ing.R, e_kin)
    rhs = continue_to_tail(v_dk, ing.R, nuc) * chi[v_init]

    f_true = nonlocal_operator(ing, nuc, F2, e_total, n_states=40)
    f_zero = np.zeros_like(f_true)

    psi_true = solve_nuclear(nuc, F2.mu, v_d_full, f_true, rhs, e_total)
    psi_zero = solve_nuclear(nuc, F2.mu, v_d_full, f_zero, rhs, e_total)
    sigma_true = da_sigma_from_psi(nuc, F2.mu, psi_true, e_total, eps_e, e_kin)
    sigma_zero = da_sigma_from_psi(nuc, F2.mu, psi_zero, e_total, eps_e, e_kin)

    ratio = sigma_zero / sigma_true
    assert ratio > 2.0 or ratio < 0.5, (
        f"zeroing F changed sigma by only {ratio:.3g}x at E={e_kin} -- F is not load-bearing"
    )


def test_coupling_phase_matches_the_background_scattering_phase(elec):
    """arg V_dk+(R) = delta_bg(R) + pi wherever |V_dk+(R)| is appreciable.

    The one independent check available anywhere in this package on V_dk+'s
    PHASE (Eq. 21) -- Task 4's gate on Eq. (68) checked only |V_dk+|^2. For
    real energy and a real short-range h, the P-space scattering solution
    phi_k+ (`scattering.scattering_state`) is a REAL function times a single
    overall complex constant on the INTERIOR of the real region (away from
    the R0 boundary node) -- measured phase spread ~2e-6 here. phi_d is real
    (Eq. 69/c-normalization), so under the bilinear c-product,
    V_dk+ = <phi_d|H_el|phi_k+> = c * (real number), i.e.
    arg(V_dk+) = arg(c) = delta_bg up to a systematic +-pi from that real
    number's sign -- measured +pi at every R checked here.

    WHY the R0 boundary node is excluded from delta_bg, not just convenient
    to exclude: the interior phase deviation from realness rises smoothly
    (single digits x1e-8 -> low 1e-7 over the last ten real nodes approaching
    R0) and then jumps ~1e6x higher AT R0 itself. The mechanism is in the DVR
    weights: every interior node's weight is real, but R0's is NOT -- it sums
    the last real element and the first ECS element (measured
    0.0975+0.0307j here) because it is the bridge node between them. A
    genuine matching bug would smear across several nodes rather than sit on
    exactly one at exactly the real/complex junction, so this is a real
    boundary-node artifact being excluded on principle, not noise being
    discarded to make the test pass.

    LIMITS (do not read more into this than it is). delta_bg is read from
    the SAME scattering_state call v_dk_plus itself makes internally, so this
    does NOT catch an error common to both sides -- e.g. substituting the
    phase-stripped standing wave for phi_k+ everywhere cancels exactly
    (verified: spread and residual both stay at their nominal ~1e-6 level,
    test still passes). What it DOES catch is a wrong pairing or conjugation
    INSIDE v_dk_plus specifically -- demonstrated by conjugating phi_k+ only
    in the d @ (h_el @ phi_k) pairing (see task-7-report.md's fix-round
    addendum), which pushes the residual to ~0.07-0.11 rad, far past the
    1.6e-3 rad tolerance used here. It does NOT catch a conjugated PAIRING on
    phi_d (phi_d is real to ~1e-8, so conjugating it is a no-op) and it does
    not pin delta_bg's own value -- a Breit-Wigner tie-in was tried and fails
    on this model (Gamma ~0.5 Ha is far too broad for the near-resonance
    approximation it relies on). The discriminating margin is set by
    2*delta_bg (~0.1 rad here, since a wrong-pairing residual lands near
    -2*delta_bg away from the correct +pi): this check weakens wherever the
    background phase itself is small.
    """
    ds = AsymptoticDiscreteState(elec, F2, R_inf=elec.R0)
    energy = 0.03
    ident = np.eye(elec.n, dtype=np.complex128)
    real = elec.points.imag == 0.0
    interior = real & (elec.points.real < elec.R0 - 1e-9)

    for R in (2.6, 3.0, 3.5):
        d = ds.phi_d(R)
        h_el = electronic_hamiltonian(elec, F2, R)
        p = ident - np.outer(d, d)
        php = p @ h_el @ p
        phi_k = scattering_state(php, elec, energy, F2.ell)

        vals = phi_k[interior]
        idx = int(np.argmax(np.abs(vals)))
        ref = vals[idx]
        mask = np.abs(vals) > 1e-2 * np.abs(ref)
        spread = float(np.max(np.abs((vals[mask] / ref).imag)))
        assert spread < 1e-3, f"phi_k+ is not real up to an overall phase at R={R}: {spread:.3g}"
        delta_bg = float(np.angle(ref))

        vdk = complex(v_dk_plus(elec, F2, ds, np.array([R]), energy)[0])
        assert abs(vdk) > 1e-3, f"|V_dk+|={abs(vdk):.3g} at R={R} is not appreciable"
        diff = (np.angle(vdk) - delta_bg - np.pi + np.pi) % (2.0 * np.pi) - np.pi
        assert abs(diff) < 1.6e-3, f"arg(V_dk+) - delta_bg - pi = {diff:.3g} rad at R={R}"


def _boundary(grid):
    real_idx = np.flatnonzero(grid.points.imag == 0.0)
    return int(real_idx[np.argmax(grid.points[real_idx].real)])
