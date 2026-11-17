"""Analytic and independent finite-radius controls; no molecular payloads."""

import json

import numpy as np
import pytest
from scipy.integrate import solve_bvp
from scipy.special import iv, ivp, kv, kvp

from projects.ukrmol_co.bound_box import (
    finite_box_dtn,
    nonnegative_exterior_witness,
    positive_exterior_witness,
)
from projects.ukrmol_co.bound_matching import BoundCaps, find_bound_state
from projects.ukrmol_co.test_bound_matching import free_map, well_spectrum


@pytest.mark.parametrize("boundary", ["dirichlet", "neumann"])
@pytest.mark.parametrize("energy", [-0.5, -0.002, -5e-9])
def test_free_box_tanh_and_coth(energy, boundary):
    thresholds = np.array([0.0, 0.4])
    kappa = np.sqrt(thresholds - 2 * energy)
    factor = np.tanh(kappa * 7)
    expected = -kappa / factor if boundary == "dirichlet" else -kappa * factor
    y, report = finite_box_dtn(
        energy,
        thresholds,
        np.zeros((2, 2)),
        np.zeros((2, 2)),
        2,
        endpoint_bohr=9,
        boundary=boundary,
        rtol=1e-12,
    )
    np.testing.assert_allclose(y, np.diag(expected), atol=3e-11, rtol=0)
    assert report["boundary"] == boundary and report["finite_box"]
    assert report["decaying"] is False and report["total_error_bound"] is None


def bessel_box_map(energy, angular, radius, endpoint, boundary):
    """Closed-form solutions of u''=(kappa**2+l(l+1)/r**2)u."""
    k = np.sqrt(-2 * energy)
    nu = angular + 0.5

    def basis(r):
        return np.sqrt(r) * np.array([iv(nu, k * r), kv(nu, k * r)])

    def derivative(r):
        return basis(r) / (2 * r) + np.sqrt(r) * k * np.array([ivp(nu, k * r), kvp(nu, k * r)])

    terminal = basis(endpoint) if boundary == "dirichlet" else derivative(endpoint)
    weights = np.array([terminal[1], -terminal[0]])
    return (derivative(radius) @ weights) / (basis(radius) @ weights)


@pytest.mark.parametrize("boundary", ["dirichlet", "neumann"])
@pytest.mark.parametrize("angular", [1, 2, 4])
def test_nonzero_angular_box_against_bessel_functions(boundary, angular):
    y, _ = finite_box_dtn(
        -0.02,
        [0],
        [[angular * (angular + 1)]],
        [[0]],
        2,
        endpoint_bohr=9,
        boundary=boundary,
        rtol=1e-12,
    )
    expected = bessel_box_map(-0.02, angular, 2, 9, boundary)
    np.testing.assert_allclose(y, [[expected]], atol=3e-11, rtol=0)


def coupled_bvp(boundary):
    """Second-order collocation, not a Riccati or inverse-Riccati solver."""
    thresholds = np.array([0.0, 0.3])
    c2 = np.array([[2.0, 0.25], [0.25, 6.0]])
    c3 = np.array([[0.1, 0.07], [0.07, -0.15]])
    x = np.linspace(2, 9, 100)

    def rhs(r, values):
        u = values[:4].reshape(2, 2, -1)
        potential = (
            np.diag(thresholds + 0.04)[:, :, None] + c2[:, :, None] / r**2 + c3[:, :, None] / r**3
        )
        return np.concatenate([values[4:], np.einsum("ijr,jkr->ikr", potential, u).reshape(4, -1)])

    def bc(left, right):
        terminal = right[:4] if boundary == "dirichlet" else right[4:]
        return np.r_[left[:4] - np.eye(2).ravel(), terminal]

    initial = np.zeros((8, len(x)))
    result = solve_bvp(rhs, bc, x, initial, tol=1e-10, max_nodes=12000)
    assert result.success, result.message
    return thresholds, c2, c3, result.y[4:, 0].reshape(2, 2)


@pytest.mark.parametrize("boundary", ["dirichlet", "neumann"])
def test_coupled_bvp_and_three_tolerances(boundary, record_property):
    thresholds, c2, c3, expected = coupled_bvp(boundary)
    errors = []
    for tolerance in [1e-4, 1e-7, 1e-10]:
        y, report = finite_box_dtn(
            -0.02,
            thresholds,
            c2,
            c3,
            2,
            endpoint_bohr=9,
            boundary=boundary,
            rtol=tolerance,
        )
        errors.append(float(np.max(abs(y - expected))))
        assert report["transpose_defect"] < 1e-8
    assert errors[-1] < errors[0] and errors[-1] < 2e-9
    record_property("tolerance_errors", json.dumps(errors))


def box_callback(endpoint, boundary):
    def callback(energy):
        return finite_box_dtn(
            energy,
            [0],
            [[0]],
            [[0]],
            1,
            endpoint_bohr=endpoint,
            boundary=boundary,
            rtol=1e-11,
        )

    return callback


def test_three_endpoint_shallow_well_brackets_same_finite_inner_root(record_property):
    poles, w = well_spectrum(1.25, 32)
    infinite = find_bound_state(poles, w, 0, (-0.01, -1e-8), free_map)["energy_hartree"]
    levels = []
    for endpoint in [80, 160, 320]:
        pair = []
        for boundary in ["neumann", "dirichlet"]:
            result = find_bound_state(
                poles,
                w,
                0,
                (-0.01, -1e-8),
                box_callback(endpoint, boundary),
                finite_box_boundary=boundary,
            )
            assert result["qualification"] == "finite-box-root-only"
            assert result["total_uncertainty_hartree"] is None
            pair.append(result["energy_hartree"])
        assert pair[0] < infinite < pair[1] < 0
        levels.append(pair)
    widths = np.diff(np.asarray(levels), axis=1).ravel()
    assert widths[2] < widths[1] < widths[0] and widths[-1] < 1e-7
    record_property("finite_inner_infinite_root_hartree", infinite)
    record_property("DN_levels_hartree", json.dumps(levels))


def test_short_box_and_subcritical_well_do_not_supply_a_bound_claim():
    poles, w = well_spectrum(1.25, 32)
    with pytest.raises(ValueError, match="bracket"):
        find_bound_state(
            poles,
            w,
            0,
            (-0.01, -1e-8),
            box_callback(10, "dirichlet"),
            finite_box_boundary="dirichlet",
        )
    poles, w = well_spectrum(1.2, 32)
    for boundary in ["neumann", "dirichlet"]:
        with pytest.raises(ValueError, match="bracket"):
            find_bound_state(
                poles,
                w,
                0,
                (-0.01, -1e-8),
                box_callback(80, boundary),
                finite_box_boundary=boundary,
            )


def test_finite_box_requires_matching_explicit_opt_in():
    with pytest.raises(ValueError, match="decaying"):
        find_bound_state([-0.5], [[0.4]], 0, (-0.49, -0.01), box_callback(10, "neumann"))
    with pytest.raises(ValueError, match="boundary"):
        find_bound_state(
            [-0.5],
            [[0.4]],
            0,
            (-0.49, -0.01),
            box_callback(10, "neumann"),
            finite_box_boundary="dirichlet",
        )
    with pytest.raises(ValueError, match="boundary"):
        find_bound_state(
            [-0.5], [[0.4]], 0, (-0.49, -0.01), free_map, finite_box_boundary="neumann"
        )
    with pytest.raises(ValueError, match="boundary"):
        find_bound_state([-0.5], [[0.4]], 0, (-0.49, -0.01), free_map, finite_box_boundary="robin")


def test_nonnegative_exterior_sufficient_witness_and_inconclusive_negative():
    report = nonnegative_exterior_witness(
        [0, 0.3], [[2, 0.2], [0.2, 0]], [[-0.1, 0.05], [0.05, 0.1]], 2
    )
    assert report["nonnegative_established"] and report["roundoff_certified"] is False
    assert report["radius_bohr"] == 2
    assert nonnegative_exterior_witness([0], [[0]], [[0]], 1)["nonnegative_established"]
    for c2, c3 in [([[-1]], [[0]]), ([[0]], [[-1]])]:
        with pytest.raises(ValueError, match="not established"):
            nonnegative_exterior_witness([0], c2, c3, 1)
    # This matrix is positive definite but not diagonally dominant: failure
    # of the sufficient witness is deliberately not a negativity diagnosis.
    with pytest.raises(ValueError, match="not established"):
        nonnegative_exterior_witness([0, 0], [[1, 2], [2, 5]], np.zeros((2, 2)), 1)


def test_positive_form_witness_without_diagonal_dominance():
    c2 = np.array([[1.0, 2.0], [2.0, 5.0]])
    report = positive_exterior_witness([0, 0], c2, np.zeros((2, 2)), 1)
    expected = 3 - np.sqrt(8)
    for endpoint in report["endpoints"]:
        np.testing.assert_allclose(endpoint["minimum_eigenvalue"], expected, atol=1e-14, rtol=0)
        assert 0 < endpoint["positive_lower_bound"] <= expected + 1e-14
    assert report["nonnegative_established"] and not report["roundoff_certified"]
    # The old API remains conservative; there is no automatic fallback.
    with pytest.raises(ValueError, match="not established"):
        nonnegative_exterior_witness([0, 0], c2, np.zeros((2, 2)), 1)


@pytest.mark.parametrize("scale", [1e-6, 1.0, 1e6])
def test_positive_form_rotation_and_whole_radial_domain(scale):
    rotation = np.array([[0.6, -0.8], [0.8, 0.6]])
    c2 = scale * rotation @ np.diag([0.2, 3.0]) @ rotation.T
    c2 = (c2 + c2.T) / 2  # Construct exactly symmetric synthetic input.
    c3 = -0.1 * c2
    report = positive_exterior_witness([0, 0.4 * scale], c2, c3, 2)
    assert all(row["positive_lower_bound"] > 0 for row in report["endpoints"])
    for radius in np.geomspace(2, 1e5, 31):
        potential = np.diag([0, 0.4 * scale]) + c2 / radius**2 + c3 / radius**3
        assert np.linalg.eigvalsh(potential)[0] > 0


@pytest.mark.parametrize("c2,c3", [([[0]], [[0]]), ([[-1]], [[0]]), ([[1]], [[-2]])])
def test_positive_form_rejects_unresolved_or_negative_endpoint(c2, c3):
    with pytest.raises(ValueError, match="not established"):
        positive_exterior_witness([0], c2, c3, 1)


@pytest.mark.parametrize("defect", ["residual", "orthogonality"])
def test_positive_form_does_not_trust_eigensolver_values_alone(monkeypatch, defect):
    from projects.ukrmol_co import bound_box

    def bad_eigh(matrix):
        return np.ones(2) * 2, np.eye(2) * (2 if defect == "orthogonality" else 1)

    monkeypatch.setattr(bound_box.np.linalg, "eigh", bad_eigh)
    with pytest.raises(ValueError, match="not established"):
        positive_exterior_witness([0, 0], np.diag([-1.0, 1.0]), np.zeros((2, 2)), 1)


def test_positive_form_caps_and_full_channel_synthetic_input():
    n = 298
    thresholds = np.linspace(0, 1, n)
    c2 = 2 * np.eye(n)
    c3 = np.zeros((n, n))
    with pytest.raises(ValueError, match="cap"):
        positive_exterior_witness(thresholds, c2, c3, 18)
    report = positive_exterior_witness(thresholds, c2, c3, 18, caps=BoundCaps(max_channels=n))
    assert min(row["positive_lower_bound"] for row in report["endpoints"]) > 1.99
    with pytest.raises(TimeoutError, match="cap"):
        positive_exterior_witness([0], [[1]], [[0]], 1, caps=BoundCaps(timeout_seconds=1e-12))
    with pytest.raises(ValueError, match="symmetric"):
        positive_exterior_witness([0, 0], [[1, 0.1], [0, 1]], np.zeros((2, 2)), 1)


@pytest.mark.parametrize("boundary", ["dirichlet", "neumann"])
def test_map_caps_dimensions_and_settings(boundary):
    for settings in [{"max_rhs": 1}, {"timeout_seconds": 1e-12}]:
        with pytest.raises(TimeoutError, match="cap"):
            finite_box_dtn(
                -0.02,
                [0],
                [[2]],
                [[0]],
                1,
                endpoint_bohr=10,
                boundary=boundary,
                caps=BoundCaps(**settings),
            )
    for energy in [0, 0.1, complex(-0.1, 0.01), np.nan]:
        with pytest.raises(ValueError):
            finite_box_dtn(energy, [0], [[0]], [[0]], 1, endpoint_bohr=10, boundary=boundary)
    with pytest.raises(ValueError, match="symmetric"):
        finite_box_dtn(
            -0.1,
            [0, 0.2],
            [[0, 1], [0, 0]],
            np.zeros((2, 2)),
            1,
            endpoint_bohr=10,
            boundary=boundary,
        )
    with pytest.raises(ValueError, match="cap"):
        finite_box_dtn(
            -0.1,
            [0, 0.2],
            np.zeros((2, 2)),
            np.zeros((2, 2)),
            1,
            endpoint_bohr=10,
            boundary=boundary,
            caps=BoundCaps(max_channels=1),
        )
    with pytest.raises(ValueError):
        finite_box_dtn(-0.1, [0], [[0]], [[0]], 2, endpoint_bohr=1, boundary=boundary)


def test_attractive_chart_failure_is_not_switched_to_other_boundary():
    with pytest.raises((ValueError, TimeoutError), match=r"chart|integration|cap"):
        finite_box_dtn(-0.001, [0], [[-100]], [[0]], 1, endpoint_bohr=10, boundary="dirichlet")


def test_finite_box_identical_input_replay_and_root_caps():
    seen = set()

    def changed(energy):
        y, report = box_callback(10, "neumann")(energy)
        if energy in seen:
            y[0, 0] += 1e-5
        seen.add(energy)
        return y, report

    with pytest.raises(ValueError, match="changed"):
        find_bound_state([-0.5], [[0.4]], 0, (-0.49, -0.01), changed, finite_box_boundary="neumann")
    with pytest.raises(TimeoutError, match="cap"):
        find_bound_state(
            [-0.5],
            [[0.4]],
            0,
            (-0.49, -0.01),
            box_callback(10, "neumann"),
            finite_box_boundary="neumann",
            caps=BoundCaps(max_evaluations=1),
        )


@pytest.mark.parametrize("failure", ["singular", "asymmetric", "changed-input"])
def test_unresolved_inverse_chart_symmetry_and_input_changes_stop(monkeypatch, failure):
    from types import SimpleNamespace

    from projects.ukrmol_co import bound_box

    second = np.eye(2)

    def fake_integrator(*args, **kwargs):
        if failure == "changed-input":
            second[0, 0] = 2
        chart = np.zeros((2, 2)) if failure == "singular" else np.array([[-1, 0.1], [0, -1]])
        return SimpleNamespace(success=True, y=chart.ravel()[:, None])

    monkeypatch.setattr(bound_box, "solve_ivp", fake_integrator)
    with pytest.raises(ValueError, match=r"chart|symmetric|changed"):
        finite_box_dtn(
            -0.1, [0, 0.2], second, np.zeros((2, 2)), 1, endpoint_bohr=10, boundary="dirichlet"
        )


@pytest.mark.parametrize("boundary", ["dirichlet", "neumann"])
def test_three_radii_free_map_approaches_infinite_tail(boundary, record_property):
    expected = -np.sqrt(0.004)
    errors = []
    for endpoint in [4, 8, 16]:
        y, _ = finite_box_dtn(
            -0.002, [0], [[0]], [[0]], 1, endpoint_bohr=endpoint, boundary=boundary, rtol=1e-12
        )
        errors.append(float(abs(y[0, 0] - expected)))
    assert errors[2] < errors[1] < errors[0]
    record_property("endpoint_errors", json.dumps(errors))


def test_witness_rejects_invalid_inputs_and_respects_time_cap():
    for thresholds, c2, c3, radius in [
        ([-0.1], [[0]], [[0]], 1),
        ([0], [[0]], [[0]], 0),
        ([0], [[np.nan]], [[0]], 1),
        ([0], [[1, 2]], [[0]], 1),
    ]:
        with pytest.raises(ValueError):
            nonnegative_exterior_witness(thresholds, c2, c3, radius)
    with pytest.raises(TimeoutError, match="cap"):
        nonnegative_exterior_witness([0], [[0]], [[0]], 1, caps=BoundCaps(timeout_seconds=1e-12))
