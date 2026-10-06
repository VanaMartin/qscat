"""Preserve failed preparations and compare refinement grids without interpolation."""

import json

import pytest

from projects.ukrmol_co.collect import collect


@pytest.mark.parametrize("scope", ["qc", "neutral"])
def test_qc_success_is_not_counted_as_independent_import_validation(tmp_path, scope):
    batch = tmp_path / "batches/qc"
    batch.mkdir(parents=True)
    (batch / "result.json").write_text(json.dumps({"jobs": [{"name": "qc", "exit_code": 0}]}))
    run = tmp_path / "runs/qc"
    run.mkdir(parents=True)
    (run / "config.json").write_text(json.dumps({f"{scope}_only": True}))
    (run / "resources.json").write_text('{"exit_code": 0}')
    (run / "result.json").write_text(json.dumps({f"{scope}_only": True}))
    (run / "stages.tsv").write_text("target\tpyscf\t\t\t10\t0\n")
    assert collect(tmp_path, [])["runs"]["qc"]["status"] == f"{scope}_validated"


@pytest.mark.parametrize("exit_code", [0, 2])
def test_job_without_a_run_directory_is_only_a_failure_if_the_batch_failed(tmp_path, exit_code):
    """A rejected ensemble can fail before config/resources/stage files exist."""
    batch = tmp_path / "batches/rejected-ensemble"
    batch.mkdir(parents=True)
    args = ["--target-roots", "1", "--sa-triplet-roots", "5", "5", "5", "5"]
    (batch / "jobs.json").write_text(json.dumps([{"name": "rejected", "args": args}]))
    (batch / "result.json").write_text(
        json.dumps({"jobs": [{"name": "rejected", "exit_code": exit_code}]})
    )
    if exit_code == 0:
        with pytest.raises(FileNotFoundError):
            collect(tmp_path, [])
    else:
        result = collect(tmp_path, [])
        rejected = result["runs"]["rejected"]
        assert rejected["status"] == "setup_failed"
        assert rejected["requested_args"] == args
        assert rejected["initial_batch_exit_code"] == 2
        assert "resources" not in rejected


@pytest.mark.parametrize("grid", ["shared", "disjoint", "equal-size-shifted"])
def test_phase_grid_refinement_uses_only_shared_native_energies(tmp_path, grid):
    """A phase branch jump and an unsampled sharp feature must not affect alignment."""
    import numpy as np

    batch = tmp_path / "batches/grid-refinement"
    batch.mkdir(parents=True)
    (batch / "result.json").write_text(
        json.dumps({"jobs": [{"name": name, "exit_code": 0} for name in ("coarse", "fine")]})
    )
    for name in ("coarse", "fine"):
        run = tmp_path / "runs" / name
        geometry = run / "output/CO/geom1"
        geometry.mkdir(parents=True)
        (run / "config.json").write_text("{}")
        (run / "resources.json").write_text('{"exit_code": 0}')
        (run / "result.json").write_text('{"native_resonances": {"B1": []}}')
        (run / "stages.tsv").write_text("")
    coarse = np.array([[0.01, 0.1], [0.03, 0.3], [0.05, 0.5]])
    fine = np.array([[0.01, 0.1 + np.pi], [0.02, 20], [0.03, 0.4 + np.pi], [0.04, -20]])
    if grid == "equal-size-shifted":
        fine = coarse.copy()
        fine[:, 0] += 0.0001
    elif grid == "disjoint":
        fine[:, 0] += 0.001
    for name, values in (("coarse", coarse), ("fine", fine)):
        np.savetxt(tmp_path / "runs" / name / "output/CO/geom1/eigenph.doublet.B1", values)
    if grid == "shared":
        comparison = collect(tmp_path, [["coarse", "fine"]])["comparisons"][0]
        assert comparison["energy_range_ev"] == [0.01, 0.03]
        assert comparison["compared_energy_points"] == 2
        assert comparison["left_energy_points"] == 3
        assert comparison["right_energy_points"] == 4
        assert comparison["max_phase_difference_rad_mod_pi"] == pytest.approx(0.1, abs=1e-14)
    else:
        with pytest.raises(ValueError, match="shared native energies"):
            collect(tmp_path, [["coarse", "fine"]])
