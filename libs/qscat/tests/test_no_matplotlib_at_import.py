"""Guard: importing `qscat.core`/`qscat.tuning` must NOT import matplotlib.

matplotlib is an OPTIONAL dependency (the `qscat[plot]` extra). A regression
where `qscat.core` eagerly imports it (e.g. a non-lazy `plot.py`) makes
`import qscat.core` -- the scattering engine -- and `import qscat.tuning`
crash `ModuleNotFoundError: matplotlib` on any clean install that did not
opt into the extra. This once shipped: `core/plot.py` imported matplotlib at
module scope and `core/__init__.py` re-exported it.

matplotlib IS installed in the dev venv, so absence can't be tested by
uninstalling it here. Instead we run each import in a fresh subprocess and
assert matplotlib never entered `sys.modules` -- i.e. the import graph does
not touch it, regardless of whether it happens to be installed. CI also runs
a true clean-venv install without the extra as the belt-and-suspenders check.
"""

from __future__ import annotations

import subprocess
import sys

import pytest


@pytest.mark.parametrize("module", ["qscat.core", "qscat.tuning", "qscat.viz", "qscat"])
def test_importing_module_does_not_import_matplotlib(module: str) -> None:
    code = (
        f"import sys; import {module}; "
        "assert 'matplotlib' not in sys.modules, "
        f"'{module} pulled in matplotlib at import time (it is an optional extra)'"
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
