#!/usr/bin/env python3
"""Emit the structure-audit map: deterministic JSON tables describing a source tree.

Run directly; the sibling `_`-prefixed modules are imported from this file's own
directory. Every table is sorted and written with sorted keys, so two runs over
identical source produce byte-identical output.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from _ast_facts import collect_symbols
from _graph import collect_callers, collect_imports
from _holders import collect_holders
from _ranking import collect_hotspots
from _similarity import collect_duplicates, collect_homonyms

# Directories named `fixtures` hold deliberate-defect test data — planted clones,
# missing docstrings, dead code. Scanning them reports the plants as findings.
SKIP_DIRS = {
    ".venv",
    ".git",
    ".hypothesis",
    "target",
    "node_modules",
    "__pycache__",
    "reference",
    ".superpowers",
    ".claude",
    "fixtures",
}


def python_files(roots: list[Path]) -> list[Path]:
    """Every `.py` file below `roots`, sorted, excluding vendored and generated trees."""
    found: list[Path] = []
    for root in roots:
        for path in sorted(root.rglob("*.py")):
            rel_parts = path.relative_to(root).parts
            if SKIP_DIRS.isdisjoint(rel_parts):
                found.append(path)
    return sorted(set(found))


def write_table(out_dir: Path, name: str, payload: object) -> None:
    """Write one table deterministically."""
    target = out_dir / f"{name}.json"
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> None:
    """Parse arguments, build every table, and write it to the output directory."""
    parser = argparse.ArgumentParser(description="Emit the structure-audit map.")
    parser.add_argument(
        "--root", action="append", required=True, type=Path, help="Source root to scan; repeatable."
    )
    parser.add_argument("--out", required=True, type=Path, help="Directory for the tables.")
    parser.add_argument(
        "--repo-root",
        required=True,
        type=Path,
        help="Repository root, used to make paths relative.",
    )
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    files = python_files(args.root)
    symbols = collect_symbols(files, args.repo_root)
    callers = collect_callers(files, args.repo_root, symbols)
    duplicates = collect_duplicates(files, args.repo_root)
    write_table(args.out, "symbols", symbols)
    write_table(args.out, "imports", collect_imports(files, args.repo_root))
    write_table(args.out, "callers", callers)
    write_table(args.out, "duplicates", duplicates)
    write_table(args.out, "homonyms", collect_homonyms(symbols))
    write_table(args.out, "holders", collect_holders(files, args.repo_root))
    write_table(
        args.out,
        "hotspots",
        collect_hotspots(symbols, callers, duplicates, files, args.repo_root),
    )


if __name__ == "__main__":
    main()
