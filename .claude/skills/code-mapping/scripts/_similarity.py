"""Structural clone detection and same-name clustering.

A function's body is normalized to a sequence of node-type tokens: docstrings and
comments vanish (they are not AST nodes, or are dropped explicitly), local names
become positional placeholders, and literal constants collapse to their type. Two
functions with identical sequences are clones; two whose sequences are 95% similar
are near-clones, which is where a copied function that drifted by a constant lands.
"""

from __future__ import annotations

import ast
import difflib
import hashlib
from collections import defaultdict
from pathlib import Path

MIN_NODES = 10
NEAR_RATIO = 0.95
BUCKET_TOLERANCE = 0.10


def _normalized_tokens(node: ast.AST) -> list[str]:
    """Node-type token sequence for a function body, with names and constants erased."""
    body = list(getattr(node, "body", []))
    if (
        body
        and isinstance(body[0], ast.Expr)
        and isinstance(body[0].value, ast.Constant)
        and isinstance(body[0].value.value, str)
    ):
        body = body[1:]  # drop the docstring (a leading STRING constant only, e.g. not `...`)
    tokens: list[str] = []
    for statement in body:
        for child in ast.walk(statement):
            if isinstance(child, ast.Name):
                tokens.append("NAME")
            elif isinstance(child, ast.Constant):
                tokens.append(f"CONST:{type(child.value).__name__}")
            elif isinstance(child, ast.Attribute):
                tokens.append(f"ATTR:{child.attr}")
            else:
                tokens.append(type(child).__name__)
    return tokens


def _functions(files: list[Path], repo_root: Path) -> list[tuple[str, list[str]]]:
    """(`path:lineno`, token sequence) for every function and method, sorted."""
    found: list[tuple[str, list[str]]] = []
    for path in files:
        rel = str(path.resolve().relative_to(repo_root.resolve()))
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                tokens = _normalized_tokens(node)
                if len(tokens) >= MIN_NODES:
                    found.append((f"{rel}:{node.lineno}", tokens))
    return sorted(found)


def _same_bucket(a: list[str], b: list[str]) -> bool:
    """True when two sequences are close enough in length to be worth comparing."""
    longer = max(len(a), len(b))
    return abs(len(a) - len(b)) <= longer * BUCKET_TOLERANCE


def collect_duplicates(files: list[Path], repo_root: Path) -> list[dict]:
    """Clone and near-clone clusters, sorted by first member.

    Near-clone grouping is ANCHOR-RELATIVE: each group is built from one anchor
    function compared against later candidates in sorted order, not from the
    transitive closure of "similar enough to any member". A chain of gradual
    drift (A~B and B~C both over NEAR_RATIO, but A~C under it) is therefore split
    across groups rather than gathered into one. A near-clone group is a LOWER
    BOUND on its family, not the complete family.
    """
    functions = _functions(files, repo_root)
    by_location = dict(functions)
    clusters: list[dict] = []
    exact: dict[str, list[str]] = defaultdict(list)
    for location, tokens in functions:
        digest = hashlib.sha256("|".join(tokens).encode()).hexdigest()
        exact[digest].append(location)

    clone_members: set[str] = set()
    for digest in sorted(exact):
        members = sorted(exact[digest])
        if len(members) > 1:
            clusters.append(
                {"kind": "clone", "nodes": len(by_location[members[0]]), "members": members}
            )
            clone_members.update(members)

    remaining = [(loc, tok) for loc, tok in functions if loc not in clone_members]
    seen: set[str] = set()
    for i, (loc_a, tok_a) in enumerate(remaining):
        if loc_a in seen:
            continue
        group = [loc_a]
        for loc_b, tok_b in remaining[i + 1 :]:
            if loc_b in seen or not _same_bucket(tok_a, tok_b):
                continue
            if difflib.SequenceMatcher(None, tok_a, tok_b).ratio() >= NEAR_RATIO:
                group.append(loc_b)
        if len(group) > 1:
            seen.update(group)
            clusters.append(
                {
                    "kind": "near-clone",
                    "nodes": len(by_location[loc_a]),
                    "members": sorted(group),
                }
            )

    clusters.sort(key=lambda c: (c["members"][0], c["kind"]))
    for index, cluster in enumerate(clusters):
        cluster["cluster"] = index
    return clusters


def collect_homonyms(symbols: list[dict]) -> list[dict]:
    """Clusters of symbols sharing a bare name, sorted by name."""
    by_name: dict[str, list[str]] = defaultdict(list)
    for record in symbols:
        by_name[record["qualname"].rsplit(".", 1)[-1]].append(
            f"{record['file']}:{record['lineno']}"
        )
    return sorted(
        (
            {"name": name, "count": len(locations), "members": sorted(locations)}
            for name, locations in by_name.items()
            if len(locations) > 1
        ),
        key=lambda h: h["name"],
    )
