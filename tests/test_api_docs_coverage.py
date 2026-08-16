"""Every public name in qscat's ``__all__`` must appear in the rendered API reference.

Guards the gap that motivated the split: the single-page ``docs/api.md``
documented eight submodules and silently omitted ``qscat.viz`` (9 public
names) and ``qscat.units`` (4), so 13 public names had no rendered
documentation at all. A layout can be re-prettified at any time; this test is
what stops the omission from recurring.

Lives at the repository root rather than under ``libs/qscat/tests/`` because it
reads ``docs/``, which ships in the repository but not in the qscat sdist.
"""

from __future__ import annotations

import importlib
import re
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
API_DOCS = REPO_ROOT / "docs" / "api"

# qscat._SUBMODULES plus `exceptions`, which is a module rather than a lazily
# exposed submodule but is equally part of the public surface (ADR 0004).
MODULES = (
    "core",
    "model",
    "dvr",
    "ecs",
    "linalg",
    "evolution",
    "special",
    "tuning",
    "viz",
    "units",
    "exceptions",
)

# Two shapes count as "documented", and both are whole-line matches:
#   * an autosummary entry -- a bare indented name on its own line
#   * an autodoc directive -- `.. autofunction:: qscat.core.ve_cross_section`
# Matching whole lines only is what keeps an incidental prose mention of a name
# from counting as documentation for it.
_DOCUMENTED = re.compile(r"^\s*(?:\.\.\s+auto\w+::\s+[\w.]*?\.)?(\w+)\s*$")


def _documented_names() -> set[str]:
    names: set[str] = set()
    for path in sorted(API_DOCS.rglob("*.md")):
        for line in path.read_text(encoding="utf-8").splitlines():
            match = _DOCUMENTED.match(line)
            if match:
                names.add(match.group(1))
    return names


def test_api_docs_directory_is_populated() -> None:
    # Without this, a mistyped path would make `_documented_names()` return an
    # empty set and every check below would fail loudly -- but a mistyped
    # *glob* would make it return nothing while the directory exists, so assert
    # both the directory and at least one page.
    assert API_DOCS.is_dir(), f"no API reference directory at {API_DOCS}"
    assert list(API_DOCS.rglob("*.md")), f"no .md pages under {API_DOCS}"


@pytest.mark.parametrize("module", MODULES)
def test_every_public_name_is_documented(module: str) -> None:
    mod = importlib.import_module(f"qscat.{module}")
    public = set(mod.__all__)
    assert public, f"qscat.{module} exports no public names -- is __all__ missing?"
    missing = sorted(public - _documented_names())
    assert not missing, (
        f"qscat.{module} exports {len(missing)} public name(s) with no entry "
        f"under docs/api/: {missing}"
    )


def test_documented_names_do_not_include_private_names() -> None:
    # The reference documents the public surface only (ADR 0004). A leading
    # underscore in an autodoc directive means a private name leaked into it.
    leaked = sorted(n for n in _documented_names() if n.startswith("_"))
    assert not leaked, f"private names documented under docs/api/: {leaked}"
