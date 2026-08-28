"""Channel-outermost block assembly: the layout that makes every off-diagonal
block a plain diagonal matrix."""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp

from projects.no_coupled_channels.blocks import assemble_coupled

N = 5


def _diag_block(seed: int) -> sp.csr_matrix:
    rng = np.random.default_rng(seed)
    a = rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))
    return sp.csr_matrix(a + a.T)  # complex SYMMETRIC, as ECS Hamiltonians are


def test_no_coupling_gives_a_block_diagonal_matrix() -> None:
    blocks = [_diag_block(1), _diag_block(2)]
    H = assemble_coupled(blocks, [[None, None], [None, None]])
    assert H.shape == (2 * N, 2 * N)
    assert H[0:N, N : 2 * N].nnz == 0
    np.testing.assert_allclose(H[0:N, 0:N].toarray(), blocks[0].toarray())
    np.testing.assert_allclose(H[N : 2 * N, N : 2 * N].toarray(), blocks[1].toarray())


def test_coupling_lands_on_the_off_diagonal_block_as_a_diagonal() -> None:
    v = np.arange(1, N + 1).astype(np.complex128)
    H = assemble_coupled([_diag_block(1), _diag_block(2)], [[None, v], [v, None]])
    np.testing.assert_allclose(H[0:N, N : 2 * N].toarray(), np.diag(v))
    np.testing.assert_allclose(H[N : 2 * N, 0:N].toarray(), np.diag(v))


def test_the_result_stays_complex_symmetric() -> None:
    v = np.linspace(0.1, 0.5, N).astype(np.complex128) * (1 + 0.3j)
    H = assemble_coupled([_diag_block(1), _diag_block(2)], [[None, v], [v, None]]).toarray()
    np.testing.assert_allclose(H, H.T, atol=1e-14)


def test_a_diagonal_coupling_entry_is_refused() -> None:
    v = np.ones(N, dtype=np.complex128)
    with pytest.raises(ValueError, match="coupling\\[0\\]\\[0\\]"):
        assemble_coupled([_diag_block(1)], [[v]])


def test_mismatched_block_size_is_refused() -> None:
    with pytest.raises(ValueError, match="size"):
        assemble_coupled(
            [_diag_block(1), sp.csr_matrix(np.zeros((N + 1, N + 1), dtype=complex))],
            [[None, None], [None, None]],
        )
