"""Per-symbol facts extracted from the AST: location, shape, and documentation."""

from __future__ import annotations

import ast
from pathlib import Path

_BRANCH_NODES = (
    ast.If,
    ast.For,
    ast.AsyncFor,
    ast.While,
    ast.ExceptHandler,
    ast.With,
    ast.AsyncWith,
    ast.BoolOp,
    ast.IfExp,
)

_NESTING_NODES = (
    ast.If,
    ast.For,
    ast.AsyncFor,
    ast.While,
    ast.Try,
    ast.With,
    ast.AsyncWith,
    ast.FunctionDef,
    ast.AsyncFunctionDef,
)


def module_qualname(path: Path, repo_root: Path) -> str:
    """Dotted module name for `path`, derived from its path below `repo_root`."""
    rel = path.resolve().relative_to(repo_root.resolve()).with_suffix("")
    parts = [p for p in rel.parts if p != "__init__"]
    return ".".join(parts)


def _branches(node: ast.AST) -> int:
    """Count branch-introducing nodes inside `node` — a cyclomatic proxy, not a score."""
    return sum(1 for child in ast.walk(node) if isinstance(child, _BRANCH_NODES))


def _max_nesting(node: ast.AST, depth: int = 0) -> int:
    """Deepest nesting of block-introducing statements inside `node`."""
    deepest = depth
    for child in ast.iter_child_nodes(node):
        step = 1 if isinstance(child, _NESTING_NODES) else 0
        deepest = max(deepest, _max_nesting(child, depth + step))
    return deepest


def _returns(node: ast.AST) -> int:
    """Count `return` statements with a value inside `node`."""
    return sum(1 for c in ast.walk(node) if isinstance(c, ast.Return) and c.value is not None)


def _docstring_lines(node: ast.AST) -> int:
    """Number of source lines the docstring occupies, or 0 when there is none."""
    doc = ast.get_docstring(node, clean=False)
    return 0 if doc is None else len(doc.splitlines())


def _signature(node: ast.FunctionDef | ast.AsyncFunctionDef) -> str:
    """Render the signature as source text."""
    return f"({ast.unparse(node.args)})" + (
        f" -> {ast.unparse(node.returns)}" if node.returns else ""
    )


def _param_count(node: ast.FunctionDef | ast.AsyncFunctionDef) -> int:
    """Count every declared parameter, including keyword-only and starred forms."""
    a = node.args
    return (
        len(a.posonlyargs)
        + len(a.args)
        + len(a.kwonlyargs)
        + (1 if a.vararg else 0)
        + (1 if a.kwarg else 0)
    )


def _comment_lines(source: str, lineno: int, end_lineno: int) -> int:
    """Count lines in [lineno, end_lineno] whose first non-space character is `#`."""
    lines = source.splitlines()[lineno - 1 : end_lineno]
    return sum(1 for line in lines if line.lstrip().startswith("#"))


def _exported_names(tree: ast.Module) -> set[str]:
    """Names listed in a module-level `__all__`."""
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
            isinstance(t, ast.Name) and t.id == "__all__" for t in node.targets
        ):
            if isinstance(node.value, (ast.List, ast.Tuple)):
                return {
                    e.value
                    for e in node.value.elts
                    if isinstance(e, ast.Constant) and isinstance(e.value, str)
                }
    return set()


def collect_symbols(files: list[Path], repo_root: Path) -> list[dict]:
    """Return one record per function, method and class across `files`, sorted by qualname."""
    records: list[dict] = []
    for path in files:
        source = path.read_text(encoding="utf-8")
        tree = ast.parse(source)
        module = module_qualname(path, repo_root)
        all_names = _exported_names(tree)
        rel = str(path.resolve().relative_to(repo_root.resolve()))
        is_test = path.name.startswith("test_") or "tests" in path.parts
        for node, owner in _iter_definitions(tree):
            name = node.name
            qual = f"{module}.{owner}.{name}" if owner else f"{module}.{name}"
            is_class = isinstance(node, ast.ClassDef)
            end = node.end_lineno or node.lineno
            doc_lines = _docstring_lines(node)
            records.append(
                {
                    "qualname": qual,
                    "kind": "class" if is_class else ("method" if owner else "function"),
                    "file": rel,
                    "lineno": node.lineno,
                    "end_lineno": end,
                    "signature": "" if is_class else _signature(node),
                    "decorators": sorted(ast.unparse(d) for d in node.decorator_list),
                    "code_lines": end - node.lineno + 1 - doc_lines,
                    "docstring_lines": doc_lines,
                    "comment_lines": _comment_lines(source, node.lineno, end),
                    "branches": _branches(node),
                    "max_nesting": _max_nesting(node),
                    "params": 0 if is_class else _param_count(node),
                    "returns": _returns(node),
                    "has_docstring": ast.get_docstring(node) is not None,
                    "exported": not name.startswith("_") and (not all_names or name in all_names),
                    "is_test": is_test,
                }
            )
    return sorted(records, key=lambda r: r["qualname"])


def _iter_definitions(tree: ast.Module) -> list[tuple[ast.AST, str]]:
    """Yield (node, owner_class_name) for every module-level def/class and its methods."""
    found: list[tuple[ast.AST, str]] = []
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            found.append((node, ""))
            if isinstance(node, ast.ClassDef):
                for member in node.body:
                    if isinstance(member, (ast.FunctionDef, ast.AsyncFunctionDef)):
                        found.append((member, node.name))
    return found
