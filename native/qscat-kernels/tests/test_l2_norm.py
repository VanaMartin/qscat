import numpy as np
import pytest

# CI builds the kernel only when native/** changes (docs/adr/0006); on jobs
# that skipped the build this module must skip, not fail collection.
qscat_kernels = pytest.importorskip("qscat_kernels")


def test_l2_norm_matches_numpy():
    v = [3.0, 4.0]
    assert qscat_kernels.l2_norm(v) == 5.0


def test_l2_norm_differential_vs_numpy():
    rng = np.random.default_rng(0)
    for _ in range(100):
        v = rng.standard_normal(rng.integers(1, 50)).tolist()
        assert abs(qscat_kernels.l2_norm(v) - float(np.linalg.norm(v))) < 1e-12
