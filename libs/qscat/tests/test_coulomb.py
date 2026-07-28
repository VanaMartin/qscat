from __future__ import annotations

import numpy as np
import pytest
from qscat.special import coulomb_f_en, coulomb_h1_en, riccati_bessel_en, riccati_hankel_en


def test_f_en_reduces_to_riccati_bessel_at_zero_charge():
    r = np.linspace(0.5, 40.0, 60)
    for l in (0, 1, 2):
        got = coulomb_f_en(r, 1.3, 0.0, 1.0, l)
        assert np.allclose(got.real, riccati_bessel_en(r, 1.3, l), rtol=1e-8, atol=1e-9)
        assert np.allclose(got.imag, 0.0, atol=1e-9)


def test_h1_en_reduces_to_riccati_hankel_at_zero_charge():
    # At eta=0: F_l(0,rho) = rho j_l(rho) but G_l(0,rho) = -rho y_l(rho) (note
    # the minus -- mpmath/Abramowitz-Stegun sign convention for the irregular
    # Coulomb function). So G+iF = rho(-y_l + i j_l) = i * rho(j_l + i y_l)
    # = i * h1_l(rho): coulomb_h1_en reduces to riccati_hankel_en up to an
    # overall phase of i, not identically. That i is consistent (mpmath
    # verified numerically), not a bug -- the regular part (coulomb_f_en,
    # the load-bearing one for the DR incident wave) matches exactly with no
    # such factor, per test_f_en_reduces_to_riccati_bessel_at_zero_charge.
    r = np.linspace(0.5, 30.0, 40)
    got = coulomb_h1_en(r, 1.0, 0.0, 1.0, 1)      # G + iF -> i * Riccati-Hankel h1
    assert np.allclose(got, 1j * riccati_hankel_en(r, 1.0, 1), rtol=1e-7, atol=1e-8)


def test_attractive_coulomb_known_value():
    # mpmath.coulombf(1, -0.5, 2.0) = 0.972687664241193 (eta = m z/k = -0.5)
    # coulomb_f_en(x, k, z, m, 1) with k*x=2, m z/k=-0.5: pick k=1, x=2, z=-0.5, m=1
    got = coulomb_f_en(np.array([2.0]), 1.0, -0.5, 1.0, 1)
    expect = np.sqrt(2.0 / np.pi) * 0.972687664241193
    assert abs(got[0] - expect) < 1e-9


def test_accepts_complex_ecs_argument():
    r = np.array([3.0 + 0.4j, 10.0 + 2.0j])       # ECS-rotated points
    got = coulomb_f_en(r, 1.0, -1.0, 1.0, 1)
    assert got.shape == (2,) and np.all(np.isfinite(got))


@pytest.mark.parametrize("bad", [(0.0, 1.0), (1.0, 0.0)])
def test_rejects_nonpositive_k_or_m(bad):
    with pytest.raises(ValueError):
        coulomb_f_en(np.array([1.0]), bad[0], -1.0, bad[1], 1)
