"""Analytic and independent small-matrix controls; no molecular inputs."""

import json

import numpy as np
import pytest
from qscat.units import hartree_to_ev
from scipy.integrate import quad, solve_bvp
from scipy.optimize import brentq

from projects.ukrmol_co.bound_matching import BoundCaps, decaying_dtn, find_bound_state
from projects.ukrmol_co.native_boundary import surface_matrix


def free_map(energy):
    return decaying_dtn(energy, [0], [[0]], [[0]], 1, endpoint_bohr=12)


def well_spectrum(depth, count):
    # Radial unit-radius square well, regular origin and Neumann inner boundary.
    # Normalized inner functions sqrt(2)*sin((n+1/2)*pi*r) give W=(-1)^n.
    n = np.arange(count)
    return 0.5 * ((n + 0.5) * np.pi) ** 2 - depth, ((-1.0) ** n)[:, None]


def exact_well_energy(depth):
    def condition(energy):
        q = np.sqrt(2 * (depth + energy))
        return q / np.tan(q) + np.sqrt(-2 * energy)

    return brentq(condition, -0.01, -1e-10, xtol=1e-15)


@pytest.mark.parametrize("energy", [-0.5, -0.02, -5e-9])
def test_free_closed_map_is_real_negative_kappa(energy):
    thresholds = np.array([0.0, 0.4])
    y, report = decaying_dtn(
        energy, thresholds, np.zeros((2, 2)), np.zeros((2, 2)), 2, endpoint_bohr=16
    )
    np.testing.assert_allclose(y, -np.diag(np.sqrt(thresholds - 2 * energy)), atol=2e-12, rtol=0)
    assert np.isrealobj(y)
    assert report["rhs_evaluations"] > 0
    assert report["transpose_defect"] < 1e-12


def test_native_exterior_norm_energy_derivative_and_radial_quadrature():
    energy, radius = -0.02, 2.0
    kappa = np.sqrt(-2 * energy)
    v = 0.3

    def normalized_radial(r):
        return np.exp(-kappa * (r - radius)) * (1 + 1 / (kappa * r)) / (1 + 1 / (kappa * radius))

    expected = 2 * v**2 * quad(lambda r: normalized_radial(r) ** 2, radius, np.inf)[0]
    errors = []
    for step in [1e-4, 5e-5, 2.5e-5]:
        values = [
            decaying_dtn(
                energy + sign * step,
                [0],
                [[2]],
                [[0]],
                radius,
                endpoint_bohr=20,
                tail_order=1,
                rtol=1e-13,
            )[0][0, 0]
            for sign in [-1, 1]
        ]
        exterior_norm = v**2 * (values[1] - values[0]) / (2 * step)
        errors.append(abs(exterior_norm - expected))
    assert errors[-1] < errors[1] < errors[0] and errors[-1] < 1e-6
    # Free-channel derivative is 1/kappa in Hartree, not 1/(2*kappa).
    np.testing.assert_allclose(v**2 / kappa, 2 * v**2 / (2 * kappa), atol=1e-15, rtol=0)


@pytest.mark.parametrize("energy", [-0.3, -0.002, -5e-9])
def test_p_wave_decaying_log_derivative(energy, record_property):
    kappa = np.sqrt(-2 * energy)
    radius = 2.0
    exact = -kappa - 1 / (radius * (kappa * radius + 1))
    errors = []
    for tolerance in [1e-6, 1e-9, 1e-12]:
        y, _ = decaying_dtn(
            energy, [0], [[2]], [[0]], radius, endpoint_bohr=20, tail_order=1, rtol=tolerance
        )
        errors.append(abs(y[0, 0] - exact))
    assert errors[-1] < errors[0]
    assert errors[-1] < 3e-11
    record_property("tolerance_errors", json.dumps(errors))


def test_direct_log_tail_three_orders_against_exact_p_wave():
    e, radius = -0.02, 2.0
    k = np.sqrt(-2 * e)
    exact = -k - 1 / (radius * (k * radius + 1))
    errors = []
    for order in [4, 8, 12]:
        y, d = decaying_dtn(
            e,
            [0],
            [[2]],
            [[0]],
            radius,
            endpoint_bohr=20,
            tail_order=order,
            tail_representation="log-derivative",
            rtol=1e-12,
        )
        errors.append(abs(y[0, 0] - exact))
        assert d["initial_transpose_defect"] == 0 and d["tail_representation"] == "log-derivative"
    assert errors[-1] < errors[1] < errors[0] and errors[-1] < 1e-10


def test_log_tail_rejects_unknown_representation_and_unresolved_branch():
    with pytest.raises(ValueError, match="representation"):
        decaying_dtn(-0.02, [0], [[2]], [[0]], 2, endpoint_bohr=20, tail_representation="fallback")
    with pytest.raises(ValueError, match="branch"):
        decaying_dtn(
            -0.02,
            [0],
            [[-100]],
            [[0]],
            2,
            endpoint_bohr=3,
            tail_order=2,
            tail_representation="log-derivative",
        )


def test_log_tail_coupled_independent_bvp_and_three_endpoints():
    thresholds = np.array([0.0, 0.3])
    c2 = np.array([[2.0, 0.25], [0.25, 6.0]])
    c3 = np.array([[0.1, 0.07], [0.07, -0.15]])
    expected = bvp_dtn(-0.0001, thresholds, c2, c3, 2.0, 2000.0)
    errors = []
    for endpoint, order, tolerance in [(60, 6, 1e-8), (120, 8, 1e-10), (240, 10, 1e-12)]:
        y, d = decaying_dtn(
            -0.0001,
            thresholds,
            c2,
            c3,
            2.0,
            endpoint_bohr=endpoint,
            tail_order=order,
            tail_representation="log-derivative",
            rtol=tolerance,
        )
        errors.append(float(np.max(abs(y - expected))))
        assert d["transpose_defect"] < 1e-10
    assert errors[-1] < errors[0] and errors[-1] < 1e-9


def bvp_dtn(energy, thresholds, c2, c3, radius, endpoint):
    """Independent second-order equation with leading Robin data far away."""
    n = len(thresholds)
    kappa = np.sqrt(thresholds - 2 * energy)
    x = np.linspace(radius, endpoint, 120)

    def rhs(r, values):
        u = values[: n * n].reshape(n, n, -1)
        derivative = values[n * n :].reshape(n, n, -1)
        potential = (
            np.diag(thresholds - 2 * energy)[:, :, None]
            + c2[:, :, None] / r**2
            + c3[:, :, None] / r**3
        )
        return np.concatenate(
            [
                derivative.reshape(n * n, -1),
                np.einsum("ijr,jkr->ikr", potential, u).reshape(n * n, -1),
            ]
        )

    def boundary(left, right):
        return np.concatenate(
            [
                left[: n * n] - np.eye(n).ravel(),
                (
                    right[n * n :].reshape(n, n) + kappa[:, None] * right[: n * n].reshape(n, n)
                ).ravel(),
            ]
        )

    u = np.zeros((n, n, len(x)))
    derivative = np.zeros_like(u)
    for i in range(n):
        u[i, i] = np.exp(-kappa[i] * (x - radius))
        derivative[i, i] = -kappa[i] * u[i, i]
    result = solve_bvp(
        rhs,
        boundary,
        x,
        np.concatenate([u.reshape(n * n, -1), derivative.reshape(n * n, -1)]),
        tol=1e-10,
        max_nodes=12000,
    )
    assert result.success, result.message
    return result.y[n * n :, 0].reshape(n, n)


def test_coupled_map_against_independent_real_bvp_and_tail_refinement(record_property):
    energy = -0.02
    thresholds = np.array([0.0, 0.3])
    c2 = np.array([[2.0, 0.25], [0.25, 6.0]])
    c3 = np.array([[0.1, 0.07], [0.07, -0.15]])
    reference = bvp_dtn(energy, thresholds, c2, c3, 2, 80)
    with pytest.raises(ValueError, match="symmetric"):
        decaying_dtn(energy, thresholds, c2, c3, 2, endpoint_bohr=12, tail_order=2, rtol=1e-12)
    errors = []
    for order in [0, 2, 4]:
        y, report = decaying_dtn(
            energy, thresholds, c2, c3, 2, endpoint_bohr=16, tail_order=order, rtol=1e-12
        )
        errors.append(np.max(abs(y - reference)))
        assert report["transpose_defect"] < 1e-8
    assert errors[-1] < errors[1] < errors[0], errors
    final, _ = decaying_dtn(
        energy, thresholds, c2, c3, 2, endpoint_bohr=40, tail_order=8, rtol=1e-12
    )
    np.testing.assert_allclose(final, reference, atol=2e-9, rtol=0)
    extent_errors = []
    for endpoint in [8, 16, 32]:
        y, _ = decaying_dtn(
            energy, thresholds, c2, c3, 2, endpoint_bohr=endpoint, tail_order=0, rtol=1e-12
        )
        extent_errors.append(np.max(abs(y - reference)))
    assert extent_errors[-1] < extent_errors[1] < extent_errors[0]
    record_property("tail_order_errors", json.dumps([float(e) for e in errors]))
    record_property("endpoint_errors", json.dumps([float(e) for e in extent_errors]))
    record_property("final_bvp_error", float(np.max(abs(final - reference))))


def test_well_root_converges_to_analytic_bound_energy(record_property):
    depth = 1.25  # Just above the first s-wave binding threshold pi**2/8.
    expected = exact_well_energy(depth)
    errors = []
    for count in [16, 32, 64]:
        poles, w = well_spectrum(depth, count)
        result = find_bound_state(poles, w, 0, (-0.01, -1e-8), free_map)
        errors.append(abs(result["energy_hartree"] - expected))
        assert result["residual_hartree"] < 1e-8
        assert result["isolation_gap_hartree"] > 1
        assert result["total_uncertainty_hartree"] is None
    assert errors[-1] < errors[1] < errors[0]
    assert hartree_to_ev(errors[-1]) < 1e-4
    record_property("analytic_energy_hartree", expected)
    record_property("spectrum_errors_hartree", json.dumps(errors))
    record_property("final_error_ev", float(hartree_to_ev(errors[-1])))


def test_subcritical_well_has_no_bound_root():
    poles, w = well_spectrum(1.2, 32)
    with pytest.raises(ValueError, match="bracket"):
        find_bound_state(poles, w, 0, (-0.1, -1e-8), free_map)


def test_q_and_inner_matching_agree_with_ground_total_shift():
    poles, w = well_spectrum(1.25, 32)
    ground = -10.0
    result = find_bound_state(poles + ground, w, ground, (-0.01, -1e-8), free_map)
    energy = result["energy_hartree"]
    q = surface_matrix(poles + ground, w, ground + energy)
    y, _ = free_map(energy)
    np.testing.assert_allclose(np.eye(1) - q @ y, 0, atol=1e-7, rtol=0)


def test_inner_pole_is_not_a_bound_root():
    with pytest.raises(ValueError, match="bracket"):
        find_bound_state([-0.2], [[0.4]], 0, (-0.21, -0.19), free_map)


@pytest.mark.parametrize("index", [0, 1])
def test_rotated_two_mode_roots_against_closed_form(index):
    angle = 0.37
    rotation = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    strengths = np.array([0.4, 0.7])
    w = rotation @ np.diag(strengths)

    def two_channels(energy):
        return -np.sqrt(-2 * energy) * np.eye(2), {}

    kappa = -(strengths[index] ** 2) + np.sqrt(strengths[index] ** 4 + 1)
    expected = -0.5 * kappa**2
    result = find_bound_state([-0.5, -0.5], w, 0, (-0.49, -0.01), two_channels, mode_index=index)
    np.testing.assert_allclose(result["energy_hartree"], expected, atol=2e-12, rtol=0)
    assert result["matching_transpose_defect_hartree"] < 1e-12


def test_coupled_root_against_opposite_boundary_matching():
    poles = np.array([-0.5, 0.1])
    w = np.array([[0.4, 0.15], [0.1, 0.3]])
    thresholds = np.array([0.0, 0.4])

    def coupled(energy):
        return -np.diag(np.sqrt(thresholds - 2 * energy)), {}

    def opposite_matching(energy):
        y, _ = coupled(energy)
        q = sum(np.outer(row, row) / (pole - energy) for pole, row in zip(poles, w, strict=True))
        return np.linalg.det(np.eye(2) - y @ q)

    expected = brentq(opposite_matching, -0.49, -0.01, xtol=1e-14)
    result = find_bound_state(poles, w, 0, (-0.49, -0.01), coupled)
    np.testing.assert_allclose(result["energy_hartree"], expected, atol=2e-12, rtol=0)


def test_fake_dtn_chart_pole_does_not_pass_as_a_root():
    def chart_pole(energy):
        return np.array([[1 / (energy + 0.2)]]), {}

    with pytest.raises((ValueError, ZeroDivisionError), match=r"residual|finite|zero"):
        find_bound_state([0.5], [[1]], 0, (-0.25, -0.15), chart_pole)


def test_degenerate_root_is_rejected():
    def diagonal_map(energy):
        return -np.sqrt(-2 * energy) * np.eye(2), {}

    with pytest.raises(ValueError, match="isolat"):
        find_bound_state([-0.5, -0.5], np.eye(2), 0, (-0.49, -0.01), diagonal_map)


def test_nondecaying_robin_branch_is_rejected():
    def wrong_branch(energy):
        return np.array([[np.sqrt(-2 * energy)]]), {"decaying": False}

    with pytest.raises(ValueError, match="decaying"):
        find_bound_state([-0.5], [[0.3]], 0, (-0.9, -0.01), wrong_branch)


def test_changed_inner_inputs_are_rejected():
    poles = np.array([-0.5])

    def changed(energy):
        poles[0] += 0.01
        return free_map(energy)

    with pytest.raises(ValueError, match="changed"):
        find_bound_state(poles, [[0.4]], 0, (-0.49, -0.01), changed)


@pytest.mark.parametrize("energy", [0.0, 0.1, -0.1 + 0.01j, np.nan])
def test_bound_map_rejects_nonnegative_complex_and_nonfinite_energy(energy):
    with pytest.raises(ValueError):
        decaying_dtn(energy, [0], [[0]], [[0]], 1, endpoint_bohr=10)


def test_invalid_dimensions_symmetry_and_caps():
    with pytest.raises(ValueError, match="symmetric"):
        decaying_dtn(-0.1, [0, 0.3], [[0, 1], [0, 0]], np.zeros((2, 2)), 1, endpoint_bohr=10)
    with pytest.raises(ValueError, match=r"[Cc]hannel"):
        decaying_dtn(
            -0.1,
            [0, 0],
            np.zeros((2, 2)),
            np.zeros((2, 2)),
            1,
            endpoint_bohr=10,
            caps=BoundCaps(max_channels=1),
        )
    with pytest.raises(ValueError, match="negative"):
        find_bound_state([-0.5], [[1]], 0, (-0.5, 0.0), free_map)
    with pytest.raises(ValueError, match="cap"):
        BoundCaps(max_rhs=0)
    with pytest.raises(TimeoutError, match="cap"):
        decaying_dtn(-0.1, [0], [[2]], [[0]], 1, endpoint_bohr=10, caps=BoundCaps(max_rhs=1))
    with pytest.raises(TimeoutError, match="cap"):
        find_bound_state(
            [-0.5], [[0.4]], 0, (-0.49, -0.01), free_map, caps=BoundCaps(max_evaluations=1)
        )


def test_identical_input_map_change_is_rejected():
    seen = set()

    def unstable(energy):
        y, report = free_map(energy)
        if energy in seen:
            y[0, 0] += 1e-5
        seen.add(energy)
        return y, report

    with pytest.raises(ValueError, match="changed"):
        find_bound_state([-0.5], [[0.4]], 0, (-0.49, -0.01), unstable)


def test_time_caps_are_checked(monkeypatch):
    from projects.ukrmol_co import bound_matching

    clock = iter([0.0, 2.0])
    monkeypatch.setattr(bound_matching, "monotonic", lambda: next(clock))
    with pytest.raises(TimeoutError, match="time cap"):
        decaying_dtn(
            -0.1, [0], [[0]], [[0]], 1, endpoint_bohr=10, caps=BoundCaps(timeout_seconds=1)
        )
    clock = iter([0.0, 2.0])
    with pytest.raises(TimeoutError, match="time cap"):
        find_bound_state(
            [-0.5], [[0.4]], 0, (-0.49, -0.01), free_map, caps=BoundCaps(timeout_seconds=1)
        )


@pytest.mark.parametrize("matrix", [[[np.nan]], [[1j]], [[np.inf]]])
def test_invalid_callback_data_stops(matrix):
    def invalid(energy):
        return np.asarray(matrix), {}

    with pytest.raises(ValueError, match="real and finite"):
        find_bound_state([-0.5], [[0.4]], 0, (-0.49, -0.01), invalid)


def test_finite_box_is_not_an_infinite_domain_oracle():
    poles, w = well_spectrum(1.25, 32)
    root = find_bound_state(poles, w, 0, (-0.01, -1e-8), free_map)
    assert root["energy_hartree"] < 0

    def short_box(energy):
        kappa = np.sqrt(-2 * energy)
        return np.array([[-kappa / np.tanh(9 * kappa)]]), {}

    # A short, Dirichlet-terminated outer box misses this shallow bound state.
    with pytest.raises(ValueError, match="bracket"):
        find_bound_state(poles, w, 0, (-0.01, -1e-8), short_box)
