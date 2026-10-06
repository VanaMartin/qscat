"""The shared profiler must preserve a non-UKRmol worker's command and failure."""

import json
import os
import sys

from projects.ukrmol_co.run import run_profiled


def test_explicit_worker_command_preserves_its_failure_and_workdir(tmp_path):
    command = [
        sys.executable,
        "-c",
        "from pathlib import Path; Path('worker.txt').write_text('done'); "
        "print('worker output'); raise SystemExit(3)",
    ]
    assert run_profiled(tmp_path, dict(os.environ), command) == 3
    assert (tmp_path / "worker.txt").read_text() == "done"
    assert "worker output" in (tmp_path / "run.log").read_text()
    resources = json.loads((tmp_path / "resources.json").read_text())
    assert resources["exit_code"] == 3
    assert resources["wall_seconds"] > 0
