import numpy as np
import qscat_kernels


def test_l2_norm_matches_numpy():
    v = [3.0, 4.0]
    assert qscat_kernels.l2_norm(v) == 5.0


def test_l2_norm_differential_vs_numpy():
    rng = np.random.default_rng(0)
    for _ in range(100):
        v = rng.standard_normal(rng.integers(1, 50)).tolist()
        assert abs(qscat_kernels.l2_norm(v) - float(np.linalg.norm(v))) < 1e-12
