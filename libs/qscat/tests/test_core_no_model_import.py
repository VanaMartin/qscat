"""Architectural guard: `qscat.core` must NOT import `qscat.model` or `projects.*`.

The whole point of the `core` (model-independent engine) / `model` (per-molecule
form + parameters) split is that the solvers depend only on the
`qscat.model.ResonanceModel` *protocol* (a structural type, imported under
`TYPE_CHECKING` only), never on any concrete model or research project. If a
`qscat.core` module ever grows a runtime `import qscat.model` / `from projects
...`, adding a new molecule (or a different model form) would stop being a pure
`qscat.model` + validation change -- this test fails loudly if that boundary is
crossed.

Checked in a FRESH interpreter (subprocess): importing `qscat.core` must leave
`qscat.model` and every `projects.*` module absent from `sys.modules`.
"""

from __future__ import annotations

import subprocess
import sys


def test_qscat_core_does_not_import_model_or_projects() -> None:
    code = (
        "import qscat.core, sys\n"
        "bad = [m for m in sys.modules "
        "if m == 'qscat.model' or m.startswith('qscat.model.') "
        "or m == 'projects' or m.startswith('projects.')]\n"
        "assert not bad, f'qscat.core pulled in forbidden modules at runtime: {bad}'\n"
        "print('OK')\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, (
        f"qscat.core violated the core-never-imports-model/projects boundary.\n"
        f"stdout: {result.stdout}\nstderr: {result.stderr}"
    )
    assert "OK" in result.stdout
