"""Keep pre-run CLI rejection evidence without hiding missing successful outputs."""

import json

import pytest

from projects.ukrmol_co.collect import collect


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
