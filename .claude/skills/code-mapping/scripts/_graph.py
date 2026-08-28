"""The intra-repository import graph and an approximate caller table.

Python is dynamic: a symbol reached through a registry, a `getattr`, or a string in
a YAML config is invisible here. Every count in `collect_callers` is therefore a
LOWER BOUND, and a zero is a candidate for phase 1 to confirm, never a verdict.
"""

from __future__ import annotations

import ast
from collections import defaultdict
from pathlib import Path

from _ast_facts import module_qualname

_PACKAGE_PREFIXES = [
    ("libs/qscat/qscat", "qscat"),
    ("apps/", "apps"),
    ("projects/", "projects"),
    ("validation/", "validation"),
    ("benchmarks/", "benchmarks"),
    ("native/", "native"),
    ("tests/", "tests"),
]


def package_of(rel_path: str) -> str:
    """Coarse package label for a repo-relative path."""
    if "/tests/" in rel_path or Path(rel_path).name.startswith("test_"):
        return "tests"
    for prefix, label in _PACKAGE_PREFIXES:
        if rel_path.startswith(prefix):
            return label
    return "other"


def collect_imports(files: list[Path], repo_root: Path) -> list[dict]:
    """Every import edge between modules inside the repository, sorted."""
    known = {module_qualname(p, repo_root) for p in files}
    edges: list[dict] = []
    for path in files:
        importer = module_qualname(path, repo_root)
        is_init = path.name == "__init__.py"
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            for target in _import_targets(node, importer, is_init):
                if target in known:
                    edges.append({"importer": importer, "imported": target, "lineno": node.lineno})
    return sorted(edges, key=lambda e: (e["importer"], e["imported"], e["lineno"]))


def _import_targets(node: ast.AST, importer: str, is_init: bool) -> list[str]:
    """Candidate imported module names for one import statement."""
    if isinstance(node, ast.Import):
        return [alias.name for alias in node.names]
    if isinstance(node, ast.ImportFrom):
        base = _resolve_relative(node, importer, is_init)
        if base is None:
            return []
        return [base] + [f"{base}.{alias.name}" for alias in node.names]
    return []


def _resolve_relative(node: ast.ImportFrom, importer: str, is_init: bool) -> str | None:
    """Absolute module name for a possibly-relative `from ... import`.

    A relative import's level counts packages up from the MODULE'S OWN PACKAGE, not
    from its qualname — and for a package `__init__.py`, `module_qualname` already
    strips the `__init__` segment, so the importer's qualname IS its package (unlike a
    regular module, whose package is its qualname minus its own last component).
    """
    if not node.level:
        return node.module
    package_parts = importer.split(".") if is_init else importer.split(".")[:-1]
    drop = node.level - 1
    base_parts = package_parts[: len(package_parts) - drop] if drop else package_parts
    return ".".join(base_parts + ([node.module] if node.module else []))


def collect_callers(files: list[Path], repo_root: Path, symbols: list[dict]) -> list[dict]:
    """Approximate fan-in per symbol, with the packages the call sites live in."""
    by_name: dict[str, list[dict]] = defaultdict(list)
    own_by_file: dict[str, set[str]] = defaultdict(set)
    for record in symbols:
        by_name[record["qualname"].rsplit(".", 1)[-1]].append(record)
        own_by_file[record["file"]].add(record["qualname"])

    sites: dict[str, set[str]] = defaultdict(set)
    counts: dict[str, int] = defaultdict(int)
    for path in files:
        rel = str(path.resolve().relative_to(repo_root.resolve()))
        package = package_of(rel)
        own = own_by_file.get(rel, set())
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for name in _referenced_names(tree):
            for record in by_name.get(name, []):
                if record["qualname"] in own:
                    continue
                counts[record["qualname"]] += 1
                sites[record["qualname"]].add(package)

    return sorted(
        (
            {
                "qualname": r["qualname"],
                "sites": counts.get(r["qualname"], 0),
                "consumer_packages": sorted(sites.get(r["qualname"], set())),
                "ambiguous": len(by_name[r["qualname"].rsplit(".", 1)[-1]]) > 1,
            }
            for r in symbols
        ),
        key=lambda c: c["qualname"],
    )


def _referenced_names(tree: ast.Module) -> list[str]:
    """Every name used in a call or attribute position, plus bare name loads."""
    names: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            func = node.func
            if isinstance(func, ast.Name):
                names.append(func.id)
            elif isinstance(func, ast.Attribute):
                names.append(func.attr)
        elif isinstance(node, ast.Attribute):
            names.append(node.attr)
        elif isinstance(node, ast.Name) and isinstance(node.ctx, ast.Load):
            names.append(node.id)
    return names
