"""Tests for the promoted energy-normalized free radial functions (`qscat.special`).

`riccati_bessel_en` is a pure move from `projects.n2_2d_cross_section.channels`
(pre-move values used here as the differential oracle); `riccati_hankel_en` is
the new outgoing sibling, `sqrt(2k/pi) r h_l^{(1)}(kr)` with
`h_l^{(1)} = j_l + i y_l`.
"""

from __future__ import annotations

import numpy as np
import pytest
from qscat.special import riccati_bessel_en, riccati_hankel_en
from scipy.special import spherical_jn, spherical_yn


def test_riccati_bessel_en_matches_pre_move_channels_values() -> None:
    from projects.n2_2d_cross_section.channels import (
        riccati_bessel_en as old_riccati_bessel_en,
    )

    r = np.linspace(0.1, 20.0, 25)
    k = 0.55
    l = 2

    got = riccati_bessel_en(r, k, l)
    want = old_riccati_bessel_en(r, k, l)

    np.testing.assert_allclose(got, want, atol=1e-14, rtol=0.0)


def test_riccati_hankel_en_real_part_equals_riccati_bessel_en() -> None:
    r = np.linspace(0.1, 20.0, 25)
    k = 0.7
    l = 1

    h = riccati_hankel_en(r, k, l)

    assert h.dtype == np.complex128
    np.testing.assert_allclose(h.real, riccati_bessel_en(r, k, l), atol=1e-14, rtol=0.0)


def test_riccati_hankel_en_imag_part_is_riccati_neumann_en() -> None:
    r = np.linspace(0.1, 20.0, 25)
    k = 0.7
    l = 1

    h = riccati_hankel_en(r, k, l)
    want_imag = np.sqrt(2.0 * k / np.pi) * r * spherical_yn(l, k * r)

    np.testing.assert_allclose(h.imag, want_imag, atol=1e-14, rtol=0.0)


def test_riccati_hankel_en_matches_direct_construction() -> None:
    r = np.linspace(0.1, 20.0, 25)
    k = 0.7
    l = 3

    got = riccati_hankel_en(r, k, l)
    want = np.sqrt(2.0 * k / np.pi) * r * (spherical_jn(l, k * r) + 1j * spherical_yn(l, k * r))

    np.testing.assert_allclose(got, want, atol=1e-14, rtol=0.0)


@pytest.mark.parametrize("k", [0.0, -1.0])
def test_riccati_bessel_en_rejects_nonpositive_k(k: float) -> None:
    with pytest.raises(ValueError):
        riccati_bessel_en(np.array([1.0, 2.0]), k, 0)


@pytest.mark.parametrize("k", [0.0, -1.0])
def test_riccati_hankel_en_rejects_nonpositive_k(k: float) -> None:
    with pytest.raises(ValueError):
        riccati_hankel_en(np.array([1.0, 2.0]), k, 0)
