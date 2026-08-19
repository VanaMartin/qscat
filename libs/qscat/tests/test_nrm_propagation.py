"""Propagation, the half-Fourier transform, and the packet diagnostics."""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp
from qscat.core.grids import nuclear_grid
from qscat.core.nrm.extended import LaunchBasis
from qscat.core.nrm.propagation import propagate_nrm


def _launch(psi0: np.ndarray, e_total: np.ndarray) -> LaunchBasis:
    """A raw-column `LaunchBasis` with `coeffs = I` -- exercises propagation
    and the transform, not the SVD factorization Task 2 already gates."""
    r = psi0.shape[1]
    return LaunchBasis(
        vectors=psi0.astype(np.complex128),
        coeffs=np.eye(r, dtype=np.complex128),
        energies=e_total,
        e_total=e_total,
        truncation_error=0.0,
    )


def test_transform_of_a_single_decaying_mode_is_the_resolvent():
    """A 1x1 'Hamiltonian' h has Psi(t) = e^{-iht}, whose transform is
    -i * i/(E-h) ... i.e. exactly (E-h)^-1. The smallest possible check that
    the -i prefactor and the quadrature weight are both right."""
    h = sp.csr_matrix(np.array([[0.20 - 0.05j]], dtype=np.complex128))
    psi0 = np.array([[1.0 + 0.0j]])
    e_total = np.array([0.10])
    res = propagate_nrm(h, _launch(psi0, e_total), nuclear_grid=None, dt=0.02, n_steps=40000)
    assert res.psi_d[0, 0] == pytest.approx(1.0 / (0.10 - (0.20 - 0.05j)), rel=1e-4)


def test_columns_do_not_talk_to_each_other():
    """Two energies propagated together == each propagated alone."""
    rng = np.random.default_rng(11)
    n = 12
    a = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    # symmetric, absorbing: -6j margin (checked against this seed's spectrum)
    # comfortably clears a+a.T's largest imaginary eigenvalue (~4.96), so
    # every mode of h actually decays under d/dt psi = -i h psi -- a smaller
    # shift (e.g. -3j) is NOT guaranteed to for an arbitrary complex-normal
    # a+a.T and, for this rng, does not (see test_survival_... below).
    h = sp.csr_matrix(a + a.T - 6j * np.eye(n))
    psi0 = (rng.normal(size=(n, 2)) + 1j * rng.normal(size=(n, 2))).astype(np.complex128)
    e = np.array([0.05, 0.09])
    both = propagate_nrm(h, _launch(psi0, e), nuclear_grid=None, dt=0.05, n_steps=400)
    for j in range(2):
        one = propagate_nrm(
            h, _launch(psi0[:, [j]], e[[j]]), nuclear_grid=None, dt=0.05, n_steps=400
        )
        assert np.allclose(both.psi_d[:, j], one.psi_d[:, 0], rtol=1e-10)


def test_survival_decays_and_unabsorbed_reports_it():
    rng = np.random.default_rng(3)
    n = 12
    a = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    # -3j is NOT enough margin for this seed: a+a.T's largest imaginary
    # eigenvalue is ~4.58, so a -3j shift leaves a genuinely GROWING mode
    # (verified against scipy.linalg.expm) and survival diverges rather than
    # decays. -6j clears it (largest shifted Im eigenvalue ~-1.42).
    h = sp.csr_matrix(a + a.T - 6j * np.eye(n))
    psi0 = np.ones((n, 1), dtype=np.complex128)
    res = propagate_nrm(h, _launch(psi0, np.array([0.05])), nuclear_grid=None, dt=0.05, n_steps=400)
    assert res.survival[0, 0] > res.survival[-1, 0]
    assert res.unabsorbed[0] == pytest.approx(res.survival[-1, 0])


def test_diagnostics_match_analytic_gaussian_at_t0():
    """A real-envelope Gaussian times a plane wave `exp(i*p0*(R-r0))` has
    EXACT `<R>_0 = r0` and `<P>_0 = p0`: the envelope's own contribution to
    `<-i d/dR>` is purely imaginary-suppressed by the c-product's Hermitian
    structure on the real region (it integrates to zero by symmetry), so
    only the phase gradient `p0` survives. This is the diagnostics-ON path
    (`nuclear_grid` given) that the other three tests never exercise --
    before the `sqrt(w)` fix in `_record`'s Eq. (4.6) term, this measured
    5.93 instead of 1.70 (a ~3.5x, grid-dependent error), not a rounding
    discrepancy.
    """
    grid = nuclear_grid()
    mask = grid.real_points <= grid.R0
    r = grid.real_points[mask]
    r0, sigma, p0 = 6.0, 0.7, 1.7
    psi_val = np.exp(-((r - r0) ** 2) / (2 * sigma**2)) * np.exp(1j * p0 * (r - r0))
    coeff = psi_val * np.sqrt(grid.weights[mask].real)

    psi0 = np.zeros((grid.n, 1), dtype=np.complex128)
    psi0[mask, 0] = coeff
    h_ext = sp.identity(grid.n, format="csr", dtype=np.complex128) * (-3j)

    res = propagate_nrm(
        h_ext, _launch(psi0, np.array([0.05])), nuclear_grid=grid, dt=0.01, n_steps=2
    )

    assert res.survival[0, 0] > 0.0
    assert res.centroid[0, 0] == pytest.approx(r0, abs=1e-9)
    assert res.momentum[0, 0] == pytest.approx(p0, abs=1e-6)


def test_low_rank_reconstruction_matches_dense_per_energy_propagation():
    """A rank-2 launch with random COMPLEX `coeffs`, gated against a dense
    per-energy reconstruction propagated on its own.

    `test_columns_do_not_talk_to_each_other` uses `coeffs = I`, which cannot
    distinguish `coeffs` from `coeffs.conj()` in the reconstruction
    `d = psi[:n_r, :] @ coeffs` (a c-product violation that would silently
    conjugate the SVD coefficients) -- with `coeffs = I` both give the same
    answer. Random complex, non-identity `coeffs` breaks that degeneracy.
    """
    rng = np.random.default_rng(7)
    n = 10
    a = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    h = sp.csr_matrix(a + a.T - 6j * np.eye(n))
    vectors = (rng.normal(size=(n, 2)) + 1j * rng.normal(size=(n, 2))).astype(np.complex128)
    coeffs = (rng.normal(size=(2, 3)) + 1j * rng.normal(size=(2, 3))).astype(np.complex128)
    e_total = np.array([0.05, 0.07, 0.09])
    launch = LaunchBasis(
        vectors=vectors, coeffs=coeffs, energies=e_total, e_total=e_total, truncation_error=0.0
    )
    batched = propagate_nrm(h, launch, nuclear_grid=None, dt=0.05, n_steps=300)
    for j in range(3):
        psi0_j = (vectors @ coeffs[:, [j]]).astype(np.complex128)
        one = propagate_nrm(
            h, _launch(psi0_j, e_total[[j]]), nuclear_grid=None, dt=0.05, n_steps=300
        )
        assert np.allclose(batched.psi_d[:, j], one.psi_d[:, 0], rtol=1e-8)
