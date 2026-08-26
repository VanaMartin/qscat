"""Layering enforcement: projects/ must not depend on validation/.

The rule (CLAUDE.md, docs/adr): validation may import projects and qscat;
projects may import qscat only; qscat imports neither. The qscat and
qscat_run sides are already enforced (test_core_no_model_import.py,
apps/qscat-run/tests/test_no_validation_import.py); this test closes the
projects side, in BOTH forms the 2026-08-25 release review found in the
wild: a literal `import validation...`, and a filesystem traversal into
`validation/` via a path string (the config.json pattern). Docstrings may
mention validation/ in prose; string literals in CODE may not.
"""

from __future__ import annotations

import ast
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PROJECTS = REPO / "projects"


def _project_files() -> list[Path]:
    files = sorted(PROJECTS.rglob("*.py"))
    assert files, f"no Python files found under {PROJECTS}"
    return files


def _docstring_nodes(tree: ast.Module) -> set[int]:
    """id()s of every Constant node that is a module/class/function docstring."""
    out: set[int] = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            body = node.body
            if (
                body
                and isinstance(body[0], ast.Expr)
                and isinstance(body[0].value, ast.Constant)
                and isinstance(body[0].value.value, str)
            ):
                out.add(id(body[0].value))
    return out


def test_projects_never_import_validation() -> None:
    offenders: list[str] = []
    for path in _project_files():
        tree = ast.parse(path.read_text(), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                if any(
                    a.name == "validation" or a.name.startswith("validation.") for a in node.names
                ):
                    offenders.append(f"{path.relative_to(REPO)}:{node.lineno}")
            elif isinstance(node, ast.ImportFrom):
                mod = node.module or ""
                if node.level == 0 and (mod == "validation" or mod.startswith("validation.")):
                    offenders.append(f"{path.relative_to(REPO)}:{node.lineno}")
    assert not offenders, (
        "projects/ imports validation/ (forbidden direction; move the "
        f"dependency into validation/ or qscat): {offenders}"
    )


def test_projects_never_reference_validation_paths() -> None:
    offenders: list[str] = []
    for path in _project_files():
        tree = ast.parse(path.read_text(), filename=str(path))
        doc_ids = _docstring_nodes(tree)
        for node in ast.walk(tree):
            if (
                isinstance(node, ast.Constant)
                and isinstance(node.value, str)
                and id(node) not in doc_ids
                and ("validation/" in node.value or node.value == "validation")
            ):
                offenders.append(f"{path.relative_to(REPO)}:{node.lineno} {node.value!r}")
    assert not offenders, (
        "projects/ code contains a path string reaching into validation/ "
        "(the config.json-traversal pattern; read qscat.model instead): "
        f"{offenders}"
    )
