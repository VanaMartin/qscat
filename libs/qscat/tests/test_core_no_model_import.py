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

`qscat.core.__init__` does NOT import every submodule (e.g. the `nrm`
subpackage is never pulled in by a plain `import qscat.core`), so checking
`sys.modules` right after that bare import would only ever see the modules
`__init__` happens to load and would pass trivially for anything it doesn't.
The subprocess therefore walks EVERY submodule of `qscat.core` with
`pkgutil.walk_packages` and imports each one explicitly before taking the
`sys.modules` snapshot -- so a future subpackage (like `nrm`) is covered
automatically, with no enumeration to keep in sync.
"""

from __future__ import annotations

import subprocess
import sys


def test_qscat_core_does_not_import_model_or_projects() -> None:
    code = (
        "import importlib, pkgutil, sys\n"
        "import qscat.core\n"
        "for _info in pkgutil.walk_packages(qscat.core.__path__, qscat.core.__name__ + '.'):\n"
        "    importlib.import_module(_info.name)\n"
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
