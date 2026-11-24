"""Separable analytic sums and independently projected dense-DVR oracle."""

import numpy as np
import pytest
import scipy.sparse as sp

from projects.potential_factory.nuclear_projection import projected_hamiltonian


def test_separable_spectrum_is_pairwise_sum():
    radial = np.array([[2, -1j], [-1j, 3]], complex)
    basis = np.eye(3)[:, :2]
    eps = np.array([0.2, 0.7])
    U = np.broadcast_to(np.array([0.1, -0.3])[:, None], (2, 3))
    H = projected_hamiltonian(radial, eps, basis, U)
    expected = (np.linalg.eigvals(radial + np.diag([0.1, -0.3]))[:, None] + eps).ravel()
    np.testing.assert_allclose(
        np.sort_complex(np.linalg.eigvals(H.toarray())),
        np.sort_complex(expected),
        rtol=1e-12,
        atol=1e-12,
    )


@pytest.mark.parametrize("count", [1, 2, 4])
def test_matches_dense_full_dvr_projection(count):
    rng = np.random.default_rng(739)
    nuclear = rng.normal(size=(4, 4))
    nuclear = nuclear + nuclear.T
    eps, vectors = np.linalg.eigh(nuclear)
    B = vectors[:, :count]
    radial = np.array([[2 - 0.1j, -0.3], [-0.3, 1 - 0.2j]])
    U = rng.normal(size=(2, 4)) + 1j * rng.normal(size=(2, 4))
    full = np.kron(radial, np.eye(4)) + np.kron(np.eye(2), nuclear) + np.diag(U.ravel())
    P = np.kron(np.eye(2), B)
    expected = P.T @ full @ P
    H = projected_hamiltonian(sp.csr_matrix(radial), eps[:count], B, U).toarray()
    np.testing.assert_allclose(H, expected, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(H, H.T, rtol=1e-12, atol=1e-12)
    if count == 4:
        np.testing.assert_allclose(
            np.sort_complex(np.linalg.eigvals(H)),
            np.sort_complex(np.linalg.eigvals(full)),
            rtol=1e-12,
            atol=1e-12,
        )


@pytest.mark.parametrize("which", ["weights", "complex_basis", "shape", "nan"])
def test_invalid_basis_or_samples_reject(which):
    B = np.eye(3)[:, :2]
    U = np.zeros((2, 3))
    if which == "weights":
        B = 2 * B
    elif which == "complex_basis":
        B = B.astype(complex)
    elif which == "shape":
        U = U.T
    else:
        U[0, 0] = np.nan
    with pytest.raises(ValueError):
        projected_hamiltonian(np.eye(2), [0.1, 0.2], B, U)
