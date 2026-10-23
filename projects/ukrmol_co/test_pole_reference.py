"""Known-truth pole/residue, dense resolvent and hidden-state counterexamples."""

import unittest

import numpy as np

from projects.ukrmol_co.native_boundary import surface_matrix
from projects.ukrmol_co.outgoing_map import outgoing_map, physical_s
from projects.ukrmol_co.pole_reference import (
    bordered_search,
    contour_witness,
    derivative_residue,
    matching,
    projected_resolvent,
)


def level(z):
    poles = np.array([0.2])
    w = np.array([[0.3]])
    k = np.sqrt(2 * z)
    q = surface_matrix(poles, w, z)
    y = np.array([[1j * k]])
    a = np.array([[np.sqrt(k) * np.exp(-2j * k)]])
    return dict(q=q, y=y, a=a, poles=poles, w=w)


def exact_pole():
    return 0.2 - 0.3**4 - 1j * 0.3**2 * np.sqrt(0.4 - 0.3**4)


def circle(center, model=level):
    points = center + 0.001 * np.exp(2j * np.pi * np.arange(32) / 32)
    interior = []
    outgoing = []
    phases = []
    for z in points:
        data = model(z)
        x = np.eye(len(data["poles"]))
        r, t = projected_resolvent(data["poles"], data["w"], z, data["y"], data["a"], x, data["q"])
        interior.append(r)
        outgoing.append(t)
        phases.append(float(np.angle(np.linalg.slogdet(matching(data["q"], data["y"]))[0])))
    return contour_witness(center, points, interior, outgoing, phases)


def test_bordered_root_and_contour_against_exact_one_level_pole():
    root, report = bordered_search(
        level, [0.18 - 0.05j, 0.182 - 0.05j], [0.1, 0.25, -0.1, -0.01], 16
    )
    np.testing.assert_allclose(root, exact_pole(), atol=1e-10, rtol=0)
    assert report["evaluations"] <= 16
    witness, moments = circle(root)
    assert all(witness["checks"].values())
    np.testing.assert_allclose(moments[-1]["pole"], exact_pole(), atol=1e-10, rtol=0)


def test_schur_projection_against_small_dense_interior_inverse():
    rng = np.random.default_rng(72)
    poles = np.array([0.1, 0.2, 0.3])
    w = rng.normal(size=(3, 2)) * 0.2
    x, _ = np.linalg.qr(rng.normal(size=(3, 2)))
    z = 0.17 - 0.04j
    y = np.array([[-0.2 + 0.4j, 0.1], [0.1, -0.5 + 0.2j]])
    a = np.array([[0.5, 0.2j]])
    q = surface_matrix(poles, w, z)
    r, t = projected_resolvent(poles, w, z, y, a, x, q)
    inverse = np.linalg.inv(z * np.eye(3) - np.diag(poles) + w @ y @ w.T)
    np.testing.assert_allclose(r, x.T @ inverse @ x, atol=1e-12, rtol=0)
    np.testing.assert_allclose(t, -a @ w.T @ inverse @ w @ a.T, atol=1e-12, rtol=0)


def test_contour_and_derivative_scattering_residue_against_exact_result():
    root = exact_pole()
    data = level(root)
    k = np.sqrt(2 * root)
    expected = -2j * k * 0.3**2 * np.exp(-4j * k) / (1 + 1j * 0.3**2 / k)
    _, moments = circle(root)
    np.testing.assert_allclose(moments[-1]["residue"], [[expected]], atol=1e-10, rtol=0)
    for h in [1e-5, 5e-6]:
        plus, minus = level(root + h), level(root - h)
        derivative = (matching(plus["q"], plus["y"]) - matching(minus["q"], minus["y"])) / (2 * h)
        residue, _ = derivative_residue(data["q"], data["y"], data["a"], derivative)
        np.testing.assert_allclose(residue, [[expected]], atol=2e-8, rtol=0)


def test_wrong_width_and_invisible_interior_poles_do_not_pass():
    wrong = complex(exact_pole().real, exact_pole().imag * 1.1)
    witness, _ = circle(wrong)
    assert not all(witness["checks"].values())
    assert not witness["checks"]["interior_signal"]

    def invisible(z):
        poles = np.array([0.18, 0.4])
        w = np.array([[0.0], [0.3]])
        k = np.sqrt(2 * z)
        y = np.array([[1j * k]])
        a = np.array([[np.sqrt(k) * np.exp(-2j * k)]])
        return dict(poles=poles, w=w, y=y, a=a, q=surface_matrix(poles, w, z))

    witness, _ = circle(0.18 + 0j, invisible)
    assert witness["checks"]["interior_signal"]
    assert witness["open_residue_norm_hartree"] < 1e-10
    assert not witness["checks"]["open_signal"]


def test_complex_jost_rows_require_explicit_continuation_and_cannot_use_real_s():
    z = 0.3 - 0.05j
    k = np.sqrt(z)
    y, a, _ = outgoing_map(
        [0], [0], [[0]], [[0]], z, 2, 18, 30, amplitudes=True, continued_amplitudes=True, rtol=1e-12
    )
    np.testing.assert_allclose(y, [[1j * k]], atol=1e-11, rtol=0)
    np.testing.assert_allclose(a, [[np.sqrt(k) * np.exp(-2j * k)]], atol=1e-9, rtol=0)
    with unittest.TestCase().assertRaisesRegex(ValueError, "Complex boundary"):
        physical_s(np.array([[1 + 0.1j]]), y, a)


def test_search_box_and_cap_fail_without_automatic_reseeding():
    case = unittest.TestCase()
    with case.assertRaisesRegex(ValueError, "evaluation cap"):
        bordered_search(level, [0.18 - 0.05j, 0.182 - 0.05j], [0.1, 0.25, -0.1, -0.01], 1)
    with case.assertRaisesRegex(ValueError, "search box"):
        bordered_search(level, [0.18 - 0.05j, 0.182 - 0.05j], [0.1, 0.181, -0.1, -0.01], 16)


if __name__ == "__main__":
    suite = unittest.TestSuite(
        unittest.FunctionTestCase(fn)
        for name, fn in list(globals().items())
        if name.startswith("test_") and callable(fn)
    )
    raise SystemExit(0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1)
