"""Tests for `qscat.linalg.kron_sum` (V1) and `qscat.linalg.c_product`.

`kron_sum` is checked against dense `np.kron` at D = 1, 2, 3, 4 with UNEQUAL
per-axis dimensions -- square cases hide transposed-index bugs.
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp
from qscat.linalg import c_product, kron_sum


def _dense_kron_sum(mats: list[np.ndarray]) -> np.ndarray:
    """Reference: sum_d I x ... x mats[d] x ... x I, built with dense np.kron."""
    sizes = [m.shape[0] for m in mats]
    total = int(np.prod(sizes))
    out = np.zeros((total, total), dtype=complex)
    for d, m in enumerate(mats):
        term = np.eye(1, dtype=complex)
        for e, n in enumerate(sizes):
            term = np.kron(term, m if e == d else np.eye(n, dtype=complex))
        out += term
    return out


def _random_mats(rng: np.random.Generator, sizes: list[int]) -> list[np.ndarray]:
    return [rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n)) for n in sizes]


@pytest.mark.parametrize("sizes", [[4], [3, 4], [2, 3, 4], [2, 3, 2, 3]])
def test_kron_sum_matches_dense_np_kron(sizes: list[int]) -> None:
    rng = np.random.default_rng(0)
    mats = _random_mats(rng, sizes)
    got = kron_sum([sp.csr_matrix(m) for m in mats]).toarray()
    want = _dense_kron_sum(mats)
    assert got.shape == (int(np.prod(sizes)),) * 2
    assert np.allclose(got, want, rtol=0, atol=1e-12)


def test_kron_sum_single_operator_is_identity_operation() -> None:
    rng = np.random.default_rng(1)
    (m,) = _random_mats(rng, [5])
    assert np.allclose(kron_sum([sp.csr_matrix(m)]).toarray(), m, rtol=0, atol=1e-14)


def test_kron_sum_acts_on_c_order_ravel() -> None:
    """LAST axis fastest: (A (x) I + I (x) B) vec(psi) == vec(A@psi + psi@B.T)."""
    rng = np.random.default_rng(2)
    A, B = _random_mats(rng, [3, 4])
    psi = rng.standard_normal((3, 4)) + 1j * rng.standard_normal((3, 4))
    got = kron_sum([sp.csr_matrix(A), sp.csr_matrix(B)]) @ psi.ravel()
    want = (A @ psi + psi @ B.T).ravel()
    assert np.allclose(got, want, rtol=0, atol=1e-12)


def test_kron_sum_rejects_empty() -> None:
    with pytest.raises(ValueError, match="at least one"):
        kron_sum([])


def test_kron_sum_rejects_non_square() -> None:
    with pytest.raises(ValueError, match="square"):
        kron_sum([sp.csr_matrix(np.zeros((2, 3)))])


def test_c_product_does_not_conjugate() -> None:
    """The whole point: c_product != vdot for complex vectors."""
    a = np.array([1j, 2.0])
    assert c_product(a, a) == pytest.approx(3.0 + 0j)  # (1j)^2 + 4 = 3
    assert np.vdot(a, a) == pytest.approx(5.0 + 0j)  # |1j|^2 + 4 = 5, NOT what we want


def test_c_product_is_symmetric() -> None:
    rng = np.random.default_rng(3)
    a = rng.standard_normal(6) + 1j * rng.standard_normal(6)
    b = rng.standard_normal(6) + 1j * rng.standard_normal(6)
    assert c_product(a, b) == pytest.approx(c_product(b, a))


def test_c_product_rejects_mismatched_shapes_even_with_equal_size() -> None:
    """`(n0, n1)` vs `(n1, n0)` with n0 != n1 has the SAME total element
    count, so a shape check performed only after `ravel()` would silently
    accept it and return a plausible-looking, transposed-axis-wrong number.
    The check must happen BEFORE flattening.
    """
    rng = np.random.default_rng(4)
    a = rng.standard_normal((2, 3)) + 1j * rng.standard_normal((2, 3))
    b = rng.standard_normal((3, 2)) + 1j * rng.standard_normal((3, 2))
    assert a.size == b.size  # same total size -- the trap a post-ravel check would miss
    with pytest.raises(ValueError, match="shape mismatch"):
        c_product(a, b)
