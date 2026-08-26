"""`TailR`: the long-range-correct `lam(R)` form -- exact asymptote, power-law
approach, flat-parameter round trip, and an analytic gradient that matches
finite differences in `_smooth_reparam`'s layout."""

from __future__ import annotations

import numpy as np

from projects.potential_factory.ansatz import TailR, params, with_params
from projects.potential_factory.fit import _smooth_grad, _smooth_reparam
from validation.factory.targets.o2 import o2_seed


def test_tailr_reaches_f_inf_exactly_as_a_power_law():
    s = TailR(f_inf=5.3, coeffs=(1.0, 1.5, -0.4), R_e=2.28, p=3, q=4)
    far = np.array([50.0, 100.0, 200.0])
    d = (s(far) - 5.3).real
    # every term dies like 2 (R_e/R)^4 * P(1): ratio 16 per doubling
    assert np.all(np.abs(d) > 0)
    np.testing.assert_allclose(d[:-1] / d[1:], 16.0, rtol=2e-2)
    assert abs(complex(s(1e6)) - 5.3) < 1e-20


def test_tailr_params_round_trip_and_update():
    s = TailR(f_inf=5.3, coeffs=(1.0, 1.5), R_e=2.28)
    m = with_params(o2_seed(), {})
    m = m.__class__(**{**m.__dict__, "lam": s})
    p = params(m)
    assert p["lam.f_inf"] == 5.3 and p["lam.c1"] == 1.5 and "lam.f_1" not in p
    m2 = with_params(m, {"lam.f_inf": 5.0, "lam.c0": 2.0})
    assert isinstance(m2.lam, TailR) and m2.lam.f_inf == 5.0 and m2.lam.coeffs == (2.0, 1.5)


def test_tailr_gradient_matches_finite_differences():
    s = TailR(f_inf=5.3, coeffs=(1.0, 1.5, -0.4), R_e=2.28, p=3, q=4)
    x0, lo, hi, decode = _smooth_reparam(s, 1.8, 6.0)
    assert x0.size == 4 and all(np.isinf(b) for b in lo + hi)
    for R in (1.9, 2.5, 4.0, 14.0):
        g = _smooth_grad(s, R)
        h = 1e-6
        fd = []
        for i in range(x0.size):
            xp, xm = x0.copy(), x0.copy()
            xp[i] += h
            xm[i] -= h
            fd.append((float(decode(xp)(R).real) - float(decode(xm)(R).real)) / (2 * h))
        np.testing.assert_allclose(g, fd, rtol=1e-6, atol=1e-9)


def test_o2_seed_lam_is_a_tail_form_with_the_right_asymptote_sign():
    m = o2_seed()
    assert isinstance(m.lam, TailR) and m.lam.q == 4
    # deeper at the equilibrium than at infinity (the anion binds MORE in the
    # molecule than O + O^- at -EA), and monotone beyond the well
    R = np.array([2.3, 4.0, 8.0, 30.0])
    lam = m.lam(R).real
    assert lam[0] > lam[-1] and np.all(np.diff(lam[1:]) < 0)
