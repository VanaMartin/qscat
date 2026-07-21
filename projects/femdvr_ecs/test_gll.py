import gll
import numpy as np


def test_gll_nodes_endpoints_and_count():
    x, w = gll.gll_nodes_weights(6)
    assert x.shape == (6,) and w.shape == (6,)
    assert np.isclose(x[0], -1.0) and np.isclose(x[-1], 1.0)
    assert np.all(np.diff(x) > 0)               # ascending
    assert np.isclose(w.sum(), 2.0)             # weights integrate 1 over (-1,1)


def test_gll_integrates_polynomials_exactly():
    # n-point GLL is exact for degree <= 2n-3
    x, w = gll.gll_nodes_weights(6)             # exact to degree 9
    for p in range(0, 10):
        approx = np.sum(w * x**p)
        exact = 0.0 if p % 2 else 2.0 / (p + 1)
        assert abs(approx - exact) < 1e-12, (p, approx, exact)


def test_diff_matrix_differentiates_exactly():
    x, _ = gll.gll_nodes_weights(7)
    D = gll.diff_matrix(x)                       # D[j,i] = L_i'(x_j); (D @ f)[j] = f'(x_j)
    for c in ([0, 1, 0, 0], [0, 0, 1, 0], [1, -2, 3, -1]):   # polynomials up to degree 3
        f = sum(ci * x**i for i, ci in enumerate(c))
        df = sum(i * ci * x**(i - 1) for i, ci in enumerate(c) if i >= 1)
        assert np.allclose(D @ f, df, atol=1e-10)
