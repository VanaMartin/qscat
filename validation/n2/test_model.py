import numpy as np

from validation.n2 import model


def test_morse_minimum_and_depth():
    p = model.PARAMS["potential"]
    R0, D0 = p["R_0"], p["D_0"]
    # V0(R0) == -D0 exactly
    assert model.v0(R0) == np.float64(-D0) or abs(model.v0(R0) + D0) < 1e-12
    # argmin over a fine grid is at R0
    R = np.linspace(1.0, 6.0, 200001)
    assert abs(R[np.argmin(model.v0(R))] - R0) < 1e-3
    # asymptote -> 0
    assert abs(model.v0(20.0)) < 1e-6


def test_lambda_at_Rc():
    p = model.PARAMS["potential"]
    assert abs(model.lam(p["R_c"]) - p["lambda_c"]) < 1e-12


def test_v_int_is_negative_decaying_well():
    R0 = model.PARAMS["potential"]["R_0"]
    assert model.v_int(1.0, R0) < 0.0
    assert abs(model.v_int(10.0, R0)) < abs(model.v_int(1.0, R0))  # decays in r


def test_v_eff_has_centrifugal_term():
    R0 = model.PARAMS["potential"]["R_0"]
    r = 2.0
    l = model.PARAMS["impulsemomentum"]
    assert abs(model.v_eff_el(r, R0) - (model.v_int(r, R0) + l * (l + 1) / (2 * r**2))) < 1e-14


def test_model_checks_all_pass():
    results = model.model_checks()
    assert results and all(ok for _n, ok, _d in results), results


def test_v_eff_el_is_complex_safe():
    # r may be an ECS-rotated (complex) tail point: v_eff_el must NOT coerce
    # to dtype=float internally (that would silently discard Im(r) and
    # corrupt the analytic continuation the ECS method relies on -- see the
    # docstring on model.v_eff_el and projects/n2_resonance/potential.py).
    R0 = model.PARAMS["potential"]["R_0"]
    r_complex = 2.0 * np.exp(1j * np.deg2rad(30.0))
    l = model.PARAMS["impulsemomentum"]
    expected = model.v_int(r_complex, R0) + l * (l + 1) / (2 * r_complex**2)
    out = model.v_eff_el(r_complex, R0)
    assert np.iscomplexobj(out)
    assert abs(out.imag) > 1e-6  # Im(r) must survive, not be discarded
    assert abs(out - expected) < 1e-12
