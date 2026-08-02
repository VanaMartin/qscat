"""Architectural guard: `qscat_run` must NOT import `validation.*` or
`projects.*`. Preset deck values are copied in as literals (see
`qscat_run/presets.py`'s module docstring); this test fails loudly if that
boundary is ever crossed.

Checked in a FRESH interpreter (subprocess), mirroring
`libs/qscat/tests/test_core_no_model_import.py`.
"""

from __future__ import annotations

import subprocess
import sys


def test_qscat_run_does_not_import_validation_or_projects() -> None:
    code = (
        "import qscat_run.cli, qscat_run.config, qscat_run.presets, sys\n"
        "bad = [m for m in sys.modules "
        "if m == 'validation' or m.startswith('validation.') "
        "or m == 'projects' or m.startswith('projects.')]\n"
        "assert not bad, f'qscat_run pulled in forbidden modules at runtime: {bad}'\n"
        "print('OK')\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, (
        f"qscat_run violated the no-validation/projects-import boundary.\n"
        f"stdout: {result.stdout}\nstderr: {result.stderr}"
    )
    assert "OK" in result.stdout
