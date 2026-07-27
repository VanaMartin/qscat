from __future__ import annotations

import numpy as np
import pytest
from qscat.special import riccati_bessel_en, riccati_bessel_en_mass


def test_reduces_to_mass_one():
    r = np.linspace(0.1, 20.0, 200)
    for l in (0, 1, 2):
        got = riccati_bessel_en_mass(r, 1.3, l, 1.0)
        assert np.allclose(got, riccati_bessel_en(r, 1.3, l), rtol=0, atol=1e-14)


def test_l0_closed_form():
    # F_{E,0}(R) = sqrt(2 mu K / pi) R j_0(KR) = sqrt(2 mu / (pi K)) sin(KR)
    r = np.linspace(0.05, 15.0, 300)
    K, mu = 4.0, 918.25
    got = riccati_bessel_en_mass(r, K, 0, mu)
    expect = np.sqrt(2.0 * mu / (np.pi * K)) * np.sin(K * r)
    assert np.allclose(got, expect, rtol=1e-12, atol=1e-12)


def test_scales_as_sqrt_mu():
    r = np.linspace(0.1, 10.0, 50)
    a = riccati_bessel_en_mass(r, 2.0, 1, 1.0)
    b = riccati_bessel_en_mass(r, 2.0, 1, 4.0)
    assert np.allclose(b, 2.0 * a, rtol=1e-13, atol=1e-13)  # sqrt(4)=2


@pytest.mark.parametrize("bad", [0.0, -1.0])
def test_rejects_nonpositive_k(bad):
    with pytest.raises(ValueError):
        riccati_bessel_en_mass(np.array([1.0]), bad, 0, 918.25)
