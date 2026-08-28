"""Tests for the NRM's real-energy electronic scattering primitive."""

from __future__ import annotations

import numpy as np
import pytest
from qscat.core.grids import electronic_grid
from qscat.core.nrm.scattering import (
    electronic_free_hamiltonian,
    incident_coefficients,
    scattering_state,
    scattering_state_minus,
)
from qscat.dvr import kinetic
from qscat.linalg import c_product
from qscat.special import riccati_bessel_en, riccati_hankel_en


@pytest.fixture(scope="module")
def grid():
    return electronic_grid(r_max=16.0, order=8, n_complex=6)


def test_free_hamiltonian_is_kinetic_plus_centrifugal(grid):
    """H_free = T_r + l(l+1)/2r^2, nothing else."""
    h_free = electronic_free_hamiltonian(grid, ell=1)
    cent = 1.0 * 2.0 / (2.0 * grid.points**2)
    expected = kinetic(grid, 1.0) + np.diag(cent)
    assert np.allclose(h_free, expected, rtol=1e-12, atol=0.0)


def test_incident_coefficients_are_function_values_times_sqrt_w(grid):
    """The incident wave is a FUNCTION and picks up sqrt(w); the ECS tail is masked."""
    k = 0.5
    c = incident_coefficients(grid, k=k, ell=1)
    real = grid.points.imag == 0.0
    vals = riccati_bessel_en(grid.points[real].real, k, 1)
    expected = vals * np.sqrt(grid.weights[real])
    assert np.allclose(c[real], expected, rtol=1e-12, atol=0.0)
    assert np.all(c[~real] == 0.0)


def test_free_potential_gives_back_the_incident_wave(grid):
    """With h == H_free the scattered part vanishes: phi+ == the incident wave.

    NOTE this does NOT check the sign convention of the driven equation: with
    h == h_free the right-hand side `(h - h_free) @ inc` is identically zero
    regardless of the sign in front of either `(E*I - h)` or the source term,
    so `phi_sc = 0` (and this test passes) for every sign combination. What it
    does check is that `electronic_free_hamiltonian` and `incident_coefficients` are
    mutually consistent -- `H_free`'s regular solution really is
    `riccati_bessel_en` on this grid. See
    `test_scattering_is_unitary_for_a_real_potential` for the test that
    actually discriminates the sign.
    """
    h_free = electronic_free_hamiltonian(grid, ell=1)
    k = 0.5
    phi = scattering_state(h_free, grid, energy=0.5 * k**2, ell=1)
    inc = incident_coefficients(grid, k=k, ell=1)
    real = grid.points.imag == 0.0
    scale = np.max(np.abs(inc[real]))
    assert np.max(np.abs((phi - inc)[real])) < 1e-8 * scale


def test_attractive_well_scatters(grid):
    """A real attractive Gaussian must produce a non-trivial scattered wave.

    NOTE this does NOT discriminate the sign convention either: it thresholds
    `|phi - inc|`, and `|phi_sc|` is identical whichever sign phi_sc carries
    (`|x| == |-x|`). See `test_scattering_is_unitary_for_a_real_potential`.
    """
    h_free = electronic_free_hamiltonian(grid, ell=1)
    well = -5.0 * np.exp(-3.0 * grid.points**2)
    h = h_free + np.diag(well)
    k = 0.5
    phi = scattering_state(h, grid, energy=0.5 * k**2, ell=1)
    inc = incident_coefficients(grid, k=k, ell=1)
    real = grid.points.imag == 0.0
    scale = np.max(np.abs(inc[real]))
    assert np.max(np.abs((phi - inc)[real])) > 1e-3 * scale


def test_scattering_is_unitary_for_a_real_potential(grid):
    """|S| == 1 for a REAL potential is a sign-sensitive check the two tests
    above cannot make.

    A bug that flips the sign of either side of the driven equation
    (`(E*I - h)` or the RHS alone) leaves `|phi_sc|` unchanged everywhere --
    it just negates `phi_sc`, and `|x| == |-x|`, so no magnitude-threshold
    test on `phi` can catch it (confirmed empirically: both single-sided
    flips reproduce `test_attractive_well_scatters`'s scattered-wave
    magnitude exactly). Unitarity does not have this blind spot: with the
    elastic T-matrix `T = <J_k | (h - H_free) | phi+>` (c-product, the 1-D
    analogue of `qscat.core.driven`'s formula) and `S = 1 - 2 pi i T` (same
    convention as `qscat.core.driven`'s module docstring), a real potential
    must give `|S| = 1` exactly -- flipping the sign of `T` breaks that
    identity unless `T == 0`. Measured on this grid/well/k: `||S| - 1|` is
    ~6e-11 for the correct sign and ~1.4e-4 for either single-sided flip, a
    4-order-of-magnitude gap that `atol=1e-6` sits cleanly inside.
    """
    h_free = electronic_free_hamiltonian(grid, ell=1)
    well = -5.0 * np.exp(-3.0 * grid.points**2)
    h = h_free + np.diag(well)
    k = 0.5
    phi = scattering_state(h, grid, energy=0.5 * k**2, ell=1)
    inc = incident_coefficients(grid, k=k, ell=1)
    t = c_product(inc, (h - h_free) @ phi)
    s = 1.0 - 2.0j * np.pi * t
    assert abs(abs(s) - 1.0) < 1e-6


def test_minus_state_is_conjugate_on_the_real_region(grid):
    """Eq. (34): phi^- = (phi^+)^* where the ECS contour is real."""
    h_free = electronic_free_hamiltonian(grid, ell=1)
    well = -5.0 * np.exp(-3.0 * grid.points**2)
    h = h_free + np.diag(well)
    e = 0.125
    plus = scattering_state(h, grid, energy=e, ell=1)
    minus = scattering_state_minus(h, grid, energy=e, ell=1)
    real = grid.points.imag == 0.0
    assert np.allclose(minus[real], np.conjugate(plus[real]), rtol=1e-12, atol=0.0)


def test_minus_state_is_zero_on_the_ecs_tail(grid):
    """The identity holds only where the contour is real, so the tail is not
    claimed -- it is zeroed, and Eq. (37)'s integrand has no support there.
    """
    h_free = electronic_free_hamiltonian(grid, ell=1)
    h = h_free + np.diag(-5.0 * np.exp(-3.0 * grid.points**2))
    minus = scattering_state_minus(h, grid, energy=0.125, ell=1)
    assert np.all(minus[grid.points.imag != 0.0] == 0.0)


def _hankel_amplitudes(
    phi: np.ndarray, inc: np.ndarray, grid, k: float, ell: int, r_probe: tuple[float, float]
) -> np.ndarray:
    """Decompose a scattering state's scattered part onto outgoing/incoming
    Riccati-Hankel functions at two real probe points.

    `phi_sc(r) = phi(r) - J_k(r) = a * h1(r) + b * h2(r)`, `h1` the OUTGOING
    Riccati-Hankel function (`riccati_hankel_en`) and `h2 = conj(h1)` (valid
    for real `r`) its INCOMING counterpart. Solved as a 2x2 linear system
    from DVR coefficients converted to function values (`psi(r_j) =
    c_j / sqrt(w_j)`, the inverse of `incident_coefficients`' convention).

    This is the non-circular check Eq. (34) actually makes a claim about: a
    purely outgoing/incoming boundary condition beyond the potential's
    support, referenced against the analytic asymptotic forms directly
    rather than against the other state's own extraction.
    """
    idx = [int(np.argmin(np.abs(grid.points.real - r))) for r in r_probe]
    w = grid.weights
    sc_vals = np.array([(phi[i] - inc[i]) / np.sqrt(w[i]) for i in idx])
    r_arr = grid.points[idx].real
    h1 = riccati_hankel_en(r_arr, k, ell)
    h2 = np.conjugate(h1)
    m = np.array([[h1[0], h2[0]], [h1[1], h2[1]]])
    ab: np.ndarray = np.linalg.solve(m, sc_vals)
    return ab


def test_minus_state_is_purely_incoming_by_hankel_decomposition(grid):
    """GATE: boundary condition, not a re-derivation of the other state's S.

    Eq. (34)'s content is that beyond the potential's support, phi+'s
    scattered part is purely OUTGOING (`h1`) and phi-'s is purely INCOMING
    (`h2 = conj(h1)` for real r): decompose each scattered part onto both
    Hankel functions directly and assert the wrong one is (numerically)
    absent. This does not reference the other state's extraction at all, so
    it cannot be circular the way testing `s_minus == conj(s_plus)` via the
    SAME outgoing-referenced formula was (that reduced to `s_plus ==
    conj(s_plus)`, i.e. Im(T) == 0 -- see the removed prior version of this
    test / the task report).

    Measured on this grid/well/e at r_probe=(8.19, 13.61) (well outside the
    Gaussian well's support, inside R0=16): |b_plus/a_plus| ~ 2.95e-7 and
    |a_minus/b_minus| ~ 2.95e-7 (identical, as the algebra predicts:
    conjugating phi+'s scattered part swaps which Hankel function carries
    the small residual). `1e-5` sits ~34x above that residual. A mutant
    `scattering_state_minus` that omits the conjugate (returns phi+
    unchanged, tail-zeroed) gives |a/b| ~ 3.4e6 for the "minus" decomposition
    -- 11 orders of magnitude over the tolerance -- confirming the gate
    discriminates a real sign/conjugation bug.
    """
    h_free = electronic_free_hamiltonian(grid, ell=1)
    h = h_free + np.diag(-5.0 * np.exp(-3.0 * grid.points**2))
    e = 0.125
    ell = 1
    k = float(np.sqrt(2.0 * e))
    r_probe = (8.18605117, 13.61244973)

    plus = scattering_state(h, grid, e, ell)
    minus = scattering_state_minus(h, grid, e, ell)
    inc = incident_coefficients(grid, k, ell)

    a_plus, b_plus = _hankel_amplitudes(plus, inc, grid, k, ell, r_probe)
    a_minus, b_minus = _hankel_amplitudes(minus, inc, grid, k, ell, r_probe)

    tol = 1e-5
    assert abs(b_plus) < tol * abs(a_plus)
    assert abs(a_minus) < tol * abs(b_minus)


def test_rejects_non_positive_energy(grid):
    h_free = electronic_free_hamiltonian(grid, ell=1)
    with pytest.raises(ValueError, match="positive"):
        scattering_state(h_free, grid, energy=0.0, ell=1)
