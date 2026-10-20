"""Independent phase/K gates and exact scientific-record replay boundaries."""

import copy

import numpy as np
import pytest

from projects.ukrmol_co import numerov_co_radial_replay as replay


def test_right_solve_scattering_and_eigenphase_against_diagonal_analytic():
    k = np.diag([0.3, -7.0, 100.0])
    expected = np.diag((1 + 1j * np.diag(k)) / (1 - 1j * np.diag(k)))
    np.testing.assert_allclose(replay.scattering(k), expected, rtol=0, atol=1e-14)
    metrics = replay.compare(k, k + np.diag([0, 0, 1e-4]))
    assert metrics["k"] > 1e-6
    assert metrics["s"] < 1e-6
    assert metrics["phase_rad_mod_pi"] < 1e-6


def record():
    return dict(
        k=[[1.0]],
        diagnostics=dict(error=1e-20),
        samples=[dict(index=0, reciprocal_condition_1=0.25, x_sha256="exact")],
        captures=[dict(file="a.npz", sha256="storage", index=0)],
        oracles=[dict(file="a.npz", sha256="storage", references={"80": ["1e-20"]})],
        kernel_wall_seconds=1.0,
        role_solve_seconds={},
        unaudited_recurring_bins={},
    )


def test_only_timing_storage_and_auxiliary_condition_observations_are_nonexact():
    old, fresh = record(), record()
    fresh["kernel_wall_seconds"] = 10.0
    fresh["samples"][0]["reciprocal_condition_1"] += 2**-54
    fresh["captures"][0]["sha256"] = "new-storage"
    fresh["oracles"][0]["sha256"] = "new-storage"
    assert replay.scientific_record(old) == replay.scientific_record(fresh)
    assert replay.condition_difference(old, fresh) <= 32 * 2**-52
    fresh["samples"][0]["x_sha256"] = "different-scalar-values"
    assert replay.scientific_record(old) != replay.scientific_record(fresh)
    fresh = copy.deepcopy(old)
    fresh["oracles"][0]["references"]["80"] = ["2e-20"]
    assert replay.scientific_record(old) != replay.scientific_record(fresh)


def test_condition_observations_must_stay_positive_finite_and_within_budget():
    old, fresh = record(), record()
    for value in [0, -1, np.nan, np.inf]:
        fresh["samples"][0]["reciprocal_condition_1"] = value
        with pytest.raises(ValueError, match="condition observation"):
            replay.condition_difference(old, fresh)
    fresh["samples"][0]["reciprocal_condition_1"] = 0.251
    assert replay.condition_difference(old, fresh) > 32 * 2**-52


def test_missing_or_duplicate_radial_samples_reject_before_interpreting_counters(tmp_path):
    expected = replay.indices(200, 0.00625)
    assert len(expected) == 11
    assert len(replay.indices(400, 0.00625)) == 15
    for ids in [expected[:-1], [*expected, expected[-1]]]:
        row = dict(
            endpoint_bohr=200,
            step_bohr=0.00625,
            diagnostics={},
            samples=[dict(index=i) for i in ids],
        )
        with pytest.raises(ValueError, match="radial sample"):
            replay.saved_row(row, 1e-6, tmp_path)
