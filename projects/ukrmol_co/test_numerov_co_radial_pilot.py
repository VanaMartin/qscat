"""Full-propagation instrumentation and independent finite-scattering gates."""

import json
from pathlib import Path

import numpy as np
import pytest
from scipy.special import spherical_jn, spherical_yn

from projects.ukrmol_co import numerov_co_radial_pilot as pilot


def analytic_arguments():
    energy = 0.3
    thresholds = np.array([0.0, 0.0, 0.8, 1.2])
    angular = np.ones(4, int)
    p = np.sqrt(energy)
    z = 18 * p
    f, g = z * spherical_jn(1, z), -z * spherical_yn(1, z)
    df = p * (spherical_jn(1, z) + z * spherical_jn(1, z, derivative=True))
    dg = -p * (spherical_yn(1, z) + z * spherical_yn(1, z, derivative=True))
    expected = np.diag([0.4, -0.3])
    q = np.diag([(f + g * k) / (df + dg * k) for k in np.diag(expected)] + [0.4, 0.5])
    return [q, thresholds, angular, 2 * np.eye(4), np.zeros((4, 4))], energy, expected


def test_full_instrumented_path_matches_analytic_inverse_square_and_restores(tmp_path):
    args, energy, expected = analytic_arguments()
    original = pilot.reference.refined_solve
    config = dict(id="test", endpoint_bohr=18.05, step_bohr=0.00625, backend="reference")
    row = pilot.propagate(args, energy, config, tmp_path / "reference", dtype=np.float64)
    np.testing.assert_allclose(row["k"], expected, rtol=0, atol=1e-10)
    assert row["solve_calls"] == 16
    assert row["samples"][-1]["role"] == "matching"
    assert row["samples"][-1]["radius_bohr"] == 18.05
    assert pilot.reference.refined_solve is original


def test_failure_restores_solver(tmp_path, monkeypatch):
    args, energy, _ = analytic_arguments()
    original = pilot.reference.refined_solve

    def fail(*args):
        raise np.linalg.LinAlgError("original failure")

    monkeypatch.setattr(pilot, "REFERENCE_SOLVE", fail)
    with pytest.raises(np.linalg.LinAlgError, match="original failure"):
        pilot.propagate(
            args,
            energy,
            dict(id="test", endpoint_bohr=18.05, step_bohr=0.00625, backend="reference"),
            tmp_path / "failed",
            dtype=np.float64,
        )
    assert pilot.reference.refined_solve is original


def test_native_k_gate_cannot_be_waived_by_small_s_and_phase_errors():
    row = dict(id="native-200-fine", endpoint_bohr=200, k=[[100.0001]])
    previous = dict(step_bohr=0.00625, endpoint_bohr=200, tail_order=8, k=[[100.0]])
    checks, _ = pilot.comparisons([row], np.array([[100.0]]), dict(rows=[previous]), 1e-6)
    assert not checks["native-200-fine_native64_k"]
    assert checks["native-200-fine_native64_s"]
    assert checks["native-200-fine_native64_phase_rad_mod_pi"]


def test_native_full_path_on_qualified_format(tmp_path):
    if np.finfo(np.longdouble).nmant != 63:
        pytest.skip("Requires qualified Linux x87 extended format")
    pytest.importorskip("qscat_extended_kernels")
    args, energy, expected = analytic_arguments()
    config = dict(id="test", endpoint_bohr=18.05, step_bohr=0.00625, backend="native")
    row = pilot.propagate(args, energy, config, tmp_path / "native")
    np.testing.assert_allclose(row["k"], expected, rtol=0, atol=1e-10)
    assert all(pilot.row_checks(row, 1e-6).values())
    assert len(row["captures"]) == 3
    assert all(all(r["checks"].values()) for r in pilot.oracle_checks(tmp_path / "native"))


def test_frozen_scope_and_long_range_sampling():
    path = Path(pilot.__file__).with_name("numerov-co-radial-pilot-contract.json")
    c = json.loads(path.read_text())
    ns = [round((r["endpoint_bohr"] - 18) / r["step_bohr"]) for r in c["configurations"]]
    assert sum(ns) == c["nominal_intervals"]
    assert 2 * sum(ns) == c["mixed_solves"]
    assert len(ns) == c["propagations"] == 5
    for config, n in zip(c["configurations"], ns, strict=True):
        indices = pilot.sample_indices(config["endpoint_bohr"], config["step_bohr"])
        assert indices <= set(range(2 * n))
        assert {0, 2 * n - 2, 2 * n - 1} <= indices
    assert len(pilot.sample_indices(200, 0.00625)) == 11
    assert len(pilot.sample_indices(400, 0.00625)) == 15
