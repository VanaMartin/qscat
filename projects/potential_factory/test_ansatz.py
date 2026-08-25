from __future__ import annotations

import numpy as np
import pytest
from qscat.model import F2, N2, NO, ResonanceModel

from projects.potential_factory.ansatz import (
    FlexibleDiatomicModel,
    SmoothR,
    from_diatomic,
    pack,
    params,
    unpack,
    with_params,
    y_p,
)

R_REAL = np.linspace(1.2, 6.0, 25)
R_CPLX = 12.0 + np.linspace(0.1, 4.0, 9) * np.exp(1j * np.deg2rad(35.0))
r_REAL = np.linspace(0.05, 12.0, 40)
r_CPLX = 16.0 + np.linspace(0.1, 6.0, 7) * np.exp(1j * np.deg2rad(40.0))


@pytest.mark.parametrize("model", [N2, NO, F2], ids=["N2", "NO", "F2"])
def test_from_diatomic_reproduces_published_form_to_roundoff(model):
    flex = from_diatomic(model)
    for R in (R_REAL, R_CPLX):
        np.testing.assert_allclose(flex.v0(R), model.v0(R), rtol=0, atol=1e-14)
        np.testing.assert_allclose(flex.lam_R(R), model.lam(R), rtol=0, atol=1e-14)
        for r in (r_REAL, r_CPLX):
            rr, RR = np.meshgrid(r, R, indexing="ij")
            np.testing.assert_allclose(flex.v_int(rr, RR), model.v_int(rr, RR), rtol=0, atol=1e-14)
            np.testing.assert_allclose(
                flex.surface(rr, RR),
                model.surface(rr, RR),
                rtol=0,
                atol=1e-13,
            )


def test_flexible_model_is_a_resonance_model():
    assert isinstance(from_diatomic(N2), ResonanceModel)


def test_emo_with_one_beta_is_morse():
    flex = from_diatomic(N2)
    assert flex.betas == (N2.alpha0,)
    assert flex.D_e == N2.D0 and flex.R_e == N2.R0


def test_y_p_is_zero_at_R_e_and_tends_to_one():
    assert y_p(2.0, 2.0, 3) == 0.0
    assert abs(y_p(1e6, 2.0, 3) - 1.0) < 1e-12


def test_smooth_r_reduces_to_houfek_sigmoid():
    lam0 = (N2.lambda_c - N2.lambda_inf) * (1 + np.exp(N2.lambda_1 * (N2.R_c - N2.R_lambda)))
    s = SmoothR(f_inf=N2.lambda_inf, f_0=lam0, f_1=N2.lambda_1, R_f=N2.R_lambda)
    np.testing.assert_allclose(s(R_REAL), N2.lam(R_REAL), atol=1e-14)


def test_shell_term_adds_a_barrier():
    base = from_diatomic(N2)
    shell = SmoothR(f_inf=0.5, f_0=0.0, f_1=1.0, R_f=0.0)
    with_shell = FlexibleDiatomicModel(
        **{**base.__dict__, "shell": shell, "alpha_b": 2.0, "r_b": 3.0}
    )
    r = np.array([3.0])
    R = np.array([2.0])
    assert with_shell.v_int(r, R).real > base.v_int(r, R).real
    assert abs(with_shell.v_int(r, R) - base.v_int(r, R) - 0.5) < 1e-14


def test_complex_inputs_are_not_coerced_to_real():
    flex = from_diatomic(N2)
    out = flex.v0(R_CPLX)
    assert out.dtype == np.complex128 and np.any(out.imag != 0.0)


def test_params_round_trip_and_unknown_key():
    flex = from_diatomic(NO)
    p = params(flex)
    assert p["D_e"] == NO.D0 and p["beta0"] == NO.alpha0 and p["lam.f_inf"] == NO.lambda_inf
    assert "shell.f_inf" not in p
    back = with_params(flex, {"lam.f_inf": 7.0, "beta0": 1.5})
    assert back.lam.f_inf == 7.0 and back.betas == (1.5,) and back.D_e == NO.D0
    with pytest.raises(KeyError):
        with_params(flex, {"nope": 1.0})


def test_pack_unpack_are_inverse():
    flex = from_diatomic(F2)
    names = ["lam.f_inf", "lam.f_0", "alpha.f_inf"]
    x = pack(flex, names)
    assert x.shape == (3,)
    again = unpack(flex, names, x * 1.0)
    assert params(again) == params(flex)
    moved = unpack(flex, names, x + 0.1)
    assert abs(moved.lam.f_inf - flex.lam.f_inf - 0.1) < 1e-12
