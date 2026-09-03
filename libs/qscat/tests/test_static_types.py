"""The static-typing fixtures must type-check clean.

`libs/qscat/tests/static_typing/` holds modules that are never imported or
executed. They record, with `typing.assert_type` and with deliberately
suppressed errors, what a type checker infers at each public surface this
package promises to narrow: the `ScatteringProblem` detail flags and the
tuner's fixed-shape reports.

Those assertions are worth nothing unless a type checker actually reads them,
and no runtime test can stand in: `return_wavefunction=True` and a runtime
`bool` that happens to be `True` produce the same object at run time, so only
a checker can tell the two calls apart. This test is therefore the gate --
it runs mypy over the fixture directory under the repository's own
configuration and fails on the first diagnostic.

Skipped when mypy is not installed (it is a dev-group tool, not a runtime
dependency of `qscat`).
"""

from __future__ import annotations

import hashlib
import pathlib
import tempfile

import pytest

_HERE = pathlib.Path(__file__).resolve()
_FIXTURES = _HERE.parent / "static_typing"
_REPO_ROOT = _HERE.parents[3]
_PYPROJECT = _REPO_ROOT / "pyproject.toml"


def _cache_dir() -> pathlib.Path:
    """A private, incremental mypy cache for this checkout.

    Not `.mypy_cache`: a developer's own `uv run mypy` writes there, and a
    checker sharing a cache directory with a concurrent run can read a
    half-written entry. Keying the directory on the repository path keeps
    parallel worktrees off each other's cache while still letting repeat runs
    of this test reuse one (a cold run costs tens of seconds, a warm one ~1 s).
    """
    key = hashlib.sha256(str(_REPO_ROOT).encode()).hexdigest()[:16]
    return pathlib.Path(tempfile.gettempdir()) / f"qscat-static-typing-mypy-{key}"


def test_static_typing_fixtures_check_clean() -> None:
    api = pytest.importorskip("mypy.api")

    targets = sorted(str(p) for p in _FIXTURES.glob("*.py"))
    assert targets, f"no static-typing fixtures found in {_FIXTURES}"

    stdout, stderr, status = api.run(
        [
            "--config-file",
            str(_PYPROJECT),
            "--cache-dir",
            str(_cache_dir()),
            *targets,
        ]
    )
    assert status == 0, f"static-typing fixtures did not check clean:\n{stdout}\n{stderr}"


def test_fixtures_are_not_collected_as_tests() -> None:
    """The fixtures are type-check input, not runnable tests.

    They call solver methods on parameters that are never supplied, so
    importing one would fail. pytest only collects `test_*.py`, and this
    check keeps a future fixture from being named into collection by
    accident.
    """
    stray = [p.name for p in _FIXTURES.glob("*.py") if p.name.startswith("test_")]
    assert not stray, f"static-typing fixtures must not be named test_*: {stray}"
