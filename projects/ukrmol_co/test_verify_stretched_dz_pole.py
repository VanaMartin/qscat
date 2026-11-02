"""Independent analytic moment, unit and support controls; no CO reads or ODE."""

import json
from pathlib import Path

import numpy as np
import pytest
from qscat.units import hartree_to_ev

from projects.ukrmol_co.verify_stretched_dz_pole import character, difference, moments, schur


def test_position_and_full_width_differences_use_incident_hartree():
    delta = difference(0.03 - 0.005j, 0.03001 - 0.00502j)
    assert delta["position_ev"] == pytest.approx(hartree_to_ev(1e-5), abs=1e-12, rel=0)
    assert delta["full_width_ev"] == pytest.approx(hartree_to_ev(4e-5), abs=1e-12, rel=0)


def test_nested_analytic_simple_pole_moments_and_winding():
    center = 0.03 - 0.005j
    points = center + 0.001 * np.exp(2j * np.pi * np.arange(32) / 32)
    b = np.diag([1.0, 0.0])
    residue = np.diag([0.2j, 0.0])
    inside = np.asarray([b / (z - center) for z in points])
    outside = np.asarray([residue / (2j * (z - center)) for z in points])
    result = moments(center, points, inside, outside, np.angle(points - center))
    for row in result:
        np.testing.assert_allclose(row["c0"], b, atol=1e-13, rtol=0)
        np.testing.assert_allclose(row["residue"], residue, atol=1e-13, rtol=0)
        np.testing.assert_allclose(row["pole"], center, atol=1e-13, rtol=0)
        assert row["winding"] == pytest.approx(1, abs=1e-13, rel=0)


def test_scalar_schur_formula_matches_direct_interior_inverse():
    z = 0.03 - 0.005j
    poles = np.asarray([0.1])
    w = np.asarray([[0.2]])
    y = np.asarray([[0.3j]])
    a = np.asarray([[1.2j]])
    q = w.T @ (w / (poles - z)[:, None])
    interior, outgoing = schur(poles, w, z, q, y, a, np.ones((1, 1)))
    exact = np.linalg.inv(z * np.eye(1) - np.diag(poles) + w @ y @ w.T)
    np.testing.assert_allclose(interior, exact, atol=1e-13, rtol=0)
    np.testing.assert_allclose(outgoing, -a @ w.T @ exact @ w @ a.T, atol=1e-13, rtol=0)


def test_character_independently_distinguishes_pi_support_and_high_projection():
    c = json.loads(Path(__file__).with_name("stretched-dz-pole-proposal-20261102.json").read_text())
    labels = [[1, 1, 1, 1], [2, 1, 2, 1], [3, 1, 3, 1], [4, 1, 3, 3], [5, 1, 4, 1], [6, 1, 4, 3]]
    for index, expected in [(0, "pass"), (3, "reject")]:
        residue = np.zeros((6, 6), complex)
        residue[index, index] = -0.1j
        assert character(residue, labels, 0.01, c["character"])["disposition"] == expected
    assert (
        character(np.zeros((6, 6)), labels, 0.01, c["character"])["disposition"] == "inconclusive"
    )
