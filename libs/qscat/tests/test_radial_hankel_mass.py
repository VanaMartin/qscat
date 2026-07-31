"""Tests for `riccati_hankel_en_mass` -- the mass-`mu` generalization of
`riccati_hankel_en` (analog of `riccati_bessel_en_mass` for the regular
function), used by the nuclear-axis `Flux` extractor's outgoing dissociation
wave (`qscat.core.td_extractors.Flux`, `axis="nuclear"`)."""

from __future__ import annotations

import numpy as np
import pytest
from qscat.special import riccati_bessel_en_mass, riccati_hankel_en, riccati_hankel_en_mass
from scipy.special import spherical_yn


def test_reduces_to_riccati_hankel_en_at_mass_one() -> None:
    r = np.linspace(0.1, 20.0, 200)
    for l in (0, 1, 2):
        got = riccati_hankel_en_mass(r, 1.3, l, 1.0)
        want = riccati_hankel_en(r, 1.3, l)
        np.testing.assert_allclose(got, want, rtol=1e-12, atol=1e-14)


def test_real_part_equals_riccati_bessel_en_mass() -> None:
    r = np.linspace(0.1, 20.0, 200)
    k, mu, l = 4.0, 918.25, 1
    got = riccati_hankel_en_mass(r, k, l, mu)
    want_real = riccati_bessel_en_mass(r, k, l, mu)
    np.testing.assert_allclose(got.real, want_real, rtol=1e-12, atol=1e-14)


def test_imag_part_is_riccati_neumann_en_mass() -> None:
    r = np.linspace(0.1, 20.0, 200)
    k, mu, l = 4.0, 918.25, 1
    got = riccati_hankel_en_mass(r, k, l, mu)
    want_imag = np.sqrt(2.0 * mu * k / np.pi) * r * spherical_yn(l, k * r)
    np.testing.assert_allclose(got.imag, want_imag, rtol=1e-12, atol=1e-14)


def test_scales_as_sqrt_mu() -> None:
    r = np.linspace(0.1, 10.0, 50)
    a = riccati_hankel_en_mass(r, 2.0, 1, 1.0)
    b = riccati_hankel_en_mass(r, 2.0, 1, 4.0)
    np.testing.assert_allclose(b, 2.0 * a, rtol=1e-13, atol=1e-13)  # sqrt(4)=2


@pytest.mark.parametrize("bad", [0.0, -1.0])
def test_rejects_nonpositive_k(bad: float) -> None:
    with pytest.raises(ValueError):
        riccati_hankel_en_mass(np.array([1.0]), bad, 0, 918.25)


@pytest.mark.parametrize("bad", [0.0, -1.0])
def test_rejects_nonpositive_mu(bad: float) -> None:
    with pytest.raises(ValueError):
        riccati_hankel_en_mass(np.array([1.0]), 1.0, 0, bad)
