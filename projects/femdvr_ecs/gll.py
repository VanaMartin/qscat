"""Gauss-Lobatto-Legendre nodes/weights and the collocation differentiation matrix.

Ported from the construction in eMoScat's FemDvrEcsGrid.cpp / DvrGrid.cpp, but built
with numpy.polynomial.legendre instead of the reference's hand-rolled QL eigensolver.
See .superpowers/sdd/femdvr-ecs-extraction.md sections 1-3.
"""

import numpy as np
from numpy.polynomial import legendre as L


def gll_nodes_weights(n: int) -> tuple[np.ndarray, np.ndarray]:
    """Gauss-Lobatto-Legendre nodes and weights on (-1, 1), n points incl. endpoints.

    Interior nodes are the roots of P'_{n-1}; endpoints are +/-1.
    Weights: w_i = 2 / (n(n-1) [P_{n-1}(x_i)]^2).
    """
    if n < 2:
        raise ValueError("n >= 2")
    # interior nodes = roots of P_{n-1}'
    coeff = np.zeros(n)          # P_{n-1}
    coeff[n - 1] = 1.0
    dcoeff = L.legder(coeff)     # P_{n-1}'
    # legroots returns complex128 (companion-matrix eigenvalues); roots of P_{n-1}'
    # are guaranteed real and simple in (-1, 1), so drop the numerically-zero
    # imaginary part to keep node/weight dtype real.
    interior = np.sort(L.legroots(dcoeff).real) if n > 2 else np.array([])
    x = np.concatenate(([-1.0], interior, [1.0]))
    Pn1 = L.legval(x, coeff)     # P_{n-1}(x_i)
    w = 2.0 / (n * (n - 1) * Pn1**2)
    return x, w


def diff_matrix(x: np.ndarray) -> np.ndarray:
    """Collocation differentiation matrix D with (D @ f)[j] = f'(x_j) for f sampled at x.

    D[j, i] = L_i'(x_j), via barycentric weights (robust for any node set).
    """
    n = x.size
    # barycentric weights
    bw = np.ones(n)
    for i in range(n):
        bw[i] = 1.0 / np.prod([x[i] - x[k] for k in range(n) if k != i])
    D = np.zeros((n, n))
    for j in range(n):
        for i in range(n):
            if i != j:
                D[j, i] = (bw[i] / bw[j]) / (x[j] - x[i])
        D[j, j] = -np.sum(D[j, :])   # negative sum trick
    return D
