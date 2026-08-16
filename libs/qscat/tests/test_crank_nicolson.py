"""Tests for `qscat.evolution.make_cn_stepper` (promoted from the N2
time-dependent cross-section sub-project's Task 1): exact-exp match for a
Hermitian H, unitarity (norm conservation) for a Hermitian H, and norm
decay for a non-Hermitian (absorbing) H.
"""

from __future__ import annotations

from typing import cast

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp
from qscat.evolution import make_cn_stepper, make_sparse_cn_stepper
from scipy.linalg import expm


def _random_sparse(n: int, density: float, rng: np.random.Generator) -> sp.csr_matrix:
    """`sp.random(..., format="csr")` typed as a concrete `csr_matrix`.

    scipy-stubs types `sp.random`'s return as the generic `spmatrix` mixin,
    which lacks the arithmetic/conversion overloads used by the tests below
    (the same incomplete-mixin quirk documented in
    `qscat.linalg.kron_sum`'s docstring). `format="csr"` already makes this a
    `csr_matrix` at runtime, so the cast changes no behavior.
    """
    return cast(
        "sp.csr_matrix",
        sp.random(n, n, density=density, format="csr", random_state=rng, dtype=float),
    )


def test_cn_matches_exact_exp_for_hermitian() -> None:
    # small Hermitian H: CN step ≈ exp(-i H dt) to O(dt^3) per step
    rng = np.random.default_rng(0)
    A = rng.standard_normal((5, 5)) + 1j * rng.standard_normal((5, 5))
    H = A + A.conj().T  # Hermitian
    dt = 1e-3
    step = make_cn_stepper(H, dt)
    psi0 = rng.standard_normal(5) + 1j * rng.standard_normal(5)
    exact = expm(-1j * H * dt) @ psi0
    # Local truncation error is O((||H|| dt)^3); for this random Hermitian H
    # (spectral norm ~5.5) the measured residual is ~2.2e-8, so the bound is
    # set at 5e-8 (comfortable margin) rather than the naive O(1)-norm 1e-8
    # estimate.
    assert np.linalg.norm(step(psi0) - exact) < 5e-8


def test_cn_unitary_for_hermitian() -> None:
    rng = np.random.default_rng(1)
    A = rng.standard_normal((6, 6)) + 1j * rng.standard_normal((6, 6))
    H = A + A.conj().T
    step = make_cn_stepper(H, 0.1)
    psi: npt.NDArray[np.complex128] = (rng.standard_normal(6) + 1j * rng.standard_normal(6)).astype(
        np.complex128
    )
    n0 = np.vdot(psi, psi).real
    for _ in range(50):
        psi = step(psi)
    assert abs(np.vdot(psi, psi).real - n0) < 1e-10  # CN preserves norm for Hermitian


def test_cn_decays_for_non_hermitian_decaying() -> None:
    # H with negative imaginary part (decaying) -> norm decreases
    H = np.diag([1.0 - 0.1j, 2.0 - 0.2j, 3.0 - 0.05j])
    step = make_cn_stepper(H, 0.05)
    psi = np.ones(3, dtype=complex)
    n0 = np.vdot(psi, psi).real
    for _ in range(100):
        psi = step(psi)
    assert np.vdot(psi, psi).real < n0  # decays


def test_sparse_cn_matches_dense_cn() -> None:
    """The sparse stepper must equal the dense one to round-off on the same H."""
    rng = np.random.default_rng(7)
    n = 40
    a = _random_sparse(n, 0.1, rng)
    b = _random_sparse(n, 0.1, rng)
    H_sp = (a + 1j * b).tocsr()
    H_dense = H_sp.toarray()
    dt = 0.05
    psi0 = (rng.standard_normal(n) + 1j * rng.standard_normal(n)).astype(np.complex128)
    dense_step = make_cn_stepper(H_dense, dt)
    sparse_step = make_sparse_cn_stepper(H_sp, dt)
    assert np.linalg.norm(sparse_step(psi0) - dense_step(psi0)) < 1e-11


def test_sparse_cn_matches_exact_exp_for_hermitian() -> None:
    rng = np.random.default_rng(8)
    n = 6
    a = _random_sparse(n, 0.3, rng)
    b = _random_sparse(n, 0.3, rng)
    m = (a + 1j * b).toarray()
    H = m + m.conj().T  # Hermitian
    dt = 1e-3
    step = make_sparse_cn_stepper(sp.csr_matrix(H), dt)
    psi0 = (rng.standard_normal(n) + 1j * rng.standard_normal(n)).astype(np.complex128)
    exact = expm(-1j * H * dt) @ psi0
    assert np.linalg.norm(step(psi0) - exact) < 5e-8


def test_sparse_cn_unitary_for_hermitian() -> None:
    rng = np.random.default_rng(9)
    n = 8
    a = _random_sparse(n, 0.4, rng)
    H_dense = a.toarray()
    H = sp.csr_matrix(H_dense + H_dense.T)  # real symmetric => Hermitian
    step = make_sparse_cn_stepper(H, 0.1)
    psi = (rng.standard_normal(n) + 1j * rng.standard_normal(n)).astype(np.complex128)
    n0 = np.vdot(psi, psi).real
    for _ in range(50):
        psi = step(psi)
    assert abs(np.vdot(psi, psi).real - n0) < 1e-10


def test_sparse_cn_decays_for_absorbing_H() -> None:
    """Negative-imaginary diagonal (ECS/optical) => norm strictly decreases."""
    H = sp.diags([1.0 - 0.1j, 2.0 - 0.2j, 3.0 - 0.05j], format="csr")
    step = make_sparse_cn_stepper(H, 0.05)
    psi = np.ones(3, dtype=np.complex128)
    n0 = np.vdot(psi, psi).real
    for _ in range(100):
        psi = step(psi)
    assert np.vdot(psi, psi).real < n0
