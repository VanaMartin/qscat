"""Analytic and independent norm/phase controls; no native fit or pole extraction."""

import numpy as np
import pytest

from projects.ukrmol_co.direct_scattering_screen import diagonal_bounds


def pole_factor(energy, center, width):
    return (energy - center - 0.5j * width) / (energy - center + 0.5j * width)


def test_diagonal_truth_and_wrong_widths_only_admit_representation():
    for width in (0.9, 1, 1.1):
        for energy in (0.5, 1, 1.5):
            s = np.diag([pole_factor(energy, 1, width), 1, 1, 1, 1, 1])
            result = diagonal_bounds(s)
            assert result["decision"] == "admit"
            assert result["width_identification_verdict"] is None
            assert result["common_potential_verdict"] is None
            assert result["lower_chord_bound"] == pytest.approx(0, abs=1e-14, rel=0)


def test_rotated_known_pole_rejects_at_resonance():
    q = np.eye(6)
    q[:2, :2] = np.array([[1, 1], [-1, 1]]) / np.sqrt(2)
    s = q @ np.diag([-1, 1, 1, 1, 1, 1]) @ q.T
    result = diagonal_bounds(s)
    assert result["decision"] == "reject"
    assert result["lower_chord_bound"] == pytest.approx(1, abs=1e-14, rel=0)


def test_chord_bound_and_independent_relative_eigenphase_identity():
    rng = np.random.default_rng(731)
    for _ in range(12):
        q, _ = np.linalg.qr(rng.normal(size=(6, 6)))
        s = q @ np.diag(np.exp(2j * rng.uniform(-1, 1, 6))) @ q.T
        result = diagonal_bounds(s)
        d = np.diag(
            np.asarray(result["candidate_diagonal_real"])
            + 1j * np.asarray(result["candidate_diagonal_imag"])
        )
        phase = np.max(abs(np.angle(np.linalg.eigvals(s.conj().T @ d)))) / 2
        assert result["constructive_upper_chord_bound"] == pytest.approx(
            2 * np.sin(phase), abs=1e-12, rel=0
        )
        arbitrary = np.diag(np.exp(1j * rng.uniform(-np.pi, np.pi, 6)))
        assert result["lower_chord_bound"] <= np.linalg.norm(s - arbitrary, ord=2) + 1e-12
        assert result["lower_chord_bound"] <= result["constructive_upper_chord_bound"] + 1e-12


def test_channel_rephasing_and_permutation_do_not_hide_mixing():
    q = np.eye(6)
    q[:2, :2] = np.array([[1, 1], [-1, 1]]) / np.sqrt(2)
    s = q @ np.diag([-1, 1, 1, 1, 1, 1]) @ q.T
    gauge = np.diag(np.exp(1j * np.arange(6)))
    permutation = [3, 1, 5, 0, 2, 4]
    transformed = (gauge.conj().T @ s @ gauge)[np.ix_(permutation, permutation)]
    a, b = diagonal_bounds(s), diagonal_bounds(transformed)
    for key in ("lower_chord_bound", "constructive_upper_chord_bound"):
        assert a[key] == pytest.approx(b[key], abs=1e-12, rel=0)
    assert a["decision"] == b["decision"] == "reject"


def test_bounds_can_be_inconclusive_without_an_automatic_optimizer():
    u = np.ones(6) / np.sqrt(6)
    s = np.eye(6) + (np.exp(0.16j) - 1) * np.outer(u, u)
    result = diagonal_bounds(s)
    assert result["lower_chord_bound"] < result["chord_budget"]
    assert result["constructive_upper_chord_bound"] > result["chord_budget"]
    assert result["decision"] == "inconclusive"


@pytest.mark.parametrize("s", [np.ones((2, 3)), np.diag([1, 0.9]), np.array([[complex("nan")]])])
def test_bad_input_is_not_a_scientific_rejection(s):
    with pytest.raises(ValueError):
        diagonal_bounds(s)
