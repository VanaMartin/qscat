"""Result-holder classes, their field sets, and the methods they repeat.

Two holders sharing most of their fields are a consolidation candidate; three
holders repeating the same method BODY are a stronger one. Two is coincidence, so
`method_echoes` starts at three.
"""

from __future__ import annotations

import ast
import hashlib
import itertools
from collections import defaultdict
from pathlib import Path

from _ast_facts import module_qualname
from _similarity import _normalized_tokens

JACCARD_THRESHOLD = 0.6
ECHO_THRESHOLD = 3

# A dataclass is marked by its DECORATOR; the other two flavours by their BASE class.
# Keeping them in one table would mean carrying an entry that every lookup must skip.
_BASE_FLAVOURS = {"NamedTuple": "namedtuple", "TypedDict": "typeddict"}


def _flavour(node: ast.ClassDef) -> str | None:
    """Holder flavour for a class, or None when it is an ordinary class."""
    if any("dataclass" in ast.unparse(d) for d in node.decorator_list):
        return "dataclass"
    for base in node.bases:
        text = ast.unparse(base)
        for needle, flavour in _BASE_FLAVOURS.items():
            if needle in text:
                return flavour
    return None


def _fields(node: ast.ClassDef) -> list[str]:
    """Annotated field names declared directly on the class."""
    return sorted(
        member.target.id
        for member in node.body
        if isinstance(member, ast.AnnAssign) and isinstance(member.target, ast.Name)
    )


def _methods(node: ast.ClassDef) -> dict[str, str]:
    """Method name to normalized-body hash."""
    out: dict[str, str] = {}
    for member in node.body:
        if isinstance(member, (ast.FunctionDef, ast.AsyncFunctionDef)):
            tokens = "|".join(_normalized_tokens(member))
            out[member.name] = hashlib.sha256(tokens.encode()).hexdigest()
    return out


def collect_holders(files: list[Path], repo_root: Path) -> dict:
    """Holders, pairwise field-overlap rows, and repeated methods."""
    holders: list[dict] = []
    method_index: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for path in files:
        rel = str(path.resolve().relative_to(repo_root.resolve()))
        module = module_qualname(path, repo_root)
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if not isinstance(node, ast.ClassDef):
                continue
            flavour = _flavour(node)
            if flavour is None:
                continue
            qual = f"{module}.{node.name}"
            holders.append(
                {
                    "qualname": qual,
                    "file": rel,
                    "lineno": node.lineno,
                    "flavour": flavour,
                    "fields": _fields(node),
                }
            )
            for name, digest in _methods(node).items():
                method_index[name].append((qual, digest))

    holders.sort(key=lambda h: h["qualname"])
    return {
        "holders": holders,
        "field_overlap_pairs": _field_overlap_pairs(holders),
        "method_echoes": _method_echoes(method_index),
    }


def _field_overlap_pairs(holders: list[dict]) -> list[dict]:
    """Pairs of holders whose field-name sets overlap above the Jaccard threshold."""
    pairs: list[dict] = []
    for a, b in itertools.combinations(holders, 2):
        fields_a, fields_b = set(a["fields"]), set(b["fields"])
        if not fields_a or not fields_b:
            continue
        overlap = len(fields_a & fields_b) / len(fields_a | fields_b)
        if overlap >= JACCARD_THRESHOLD:
            pairs.append(
                {
                    "members": sorted([a["qualname"], b["qualname"]]),
                    "shared": sorted(fields_a & fields_b),
                    "distinct": sorted(fields_a ^ fields_b),
                }
            )
    return sorted(pairs, key=lambda p: p["members"])


def _method_echoes(index: dict[str, list[tuple[str, str]]]) -> list[dict]:
    """Methods defined on at least ECHO_THRESHOLD holders, flagged when bodies match."""
    echoes: list[dict] = []
    for name in sorted(index):
        owners = sorted(index[name])
        if len(owners) < ECHO_THRESHOLD:
            continue
        echoes.append(
            {
                "method": name,
                "count": len(owners),
                "members": [qual for qual, _ in owners],
                "identical": len({digest for _, digest in owners}) == 1,
            }
        )
    return echoes
