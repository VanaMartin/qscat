"""The angular factor's only contract: orthonormality under the quadrature
this project actually integrates with."""

from __future__ import annotations

import numpy as np
import pytest

from projects.no_coupled_channels.angular import theta_lm

N_NODES = 64


@pytest.mark.parametrize("Lambda", [0, 1])
def test_theta_lm_is_orthonormal_under_gauss_legendre(Lambda: int) -> None:
    x, w = np.polynomial.legendre.leggauss(N_NODES)
    ells = list(range(Lambda, Lambda + 5))
    gram = np.array(
        [
            [float(np.sum(w * theta_lm(a, Lambda, x) * theta_lm(b, Lambda, x))) for b in ells]
            for a in ells
        ]
    )
    np.testing.assert_allclose(gram, np.eye(len(ells)), atol=1e-13)


def test_theta_lm_is_real_and_finite_at_the_poles() -> None:
    x = np.array([-1.0, 0.0, 1.0])
    out = theta_lm(2, 1, x)
    assert out.dtype == np.float64
    assert np.all(np.isfinite(out))
