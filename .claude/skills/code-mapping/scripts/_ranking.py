"""The review queue.

This scores REVIEW PRIORITY, not quality. A high score means "look here first",
never "this is bad"; a long, branchy function may be exactly right. The signals are
all measured, and each contributing signal is named in the record so a reader can
see why a symbol ranked where it did.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

_SIGNALS: list[tuple[str, int, Callable[[dict], bool]]] = [
    ("long", 2, lambda s: s["code_lines"] > 60),
    ("very-long", 2, lambda s: s["code_lines"] > 120),
    ("branchy", 2, lambda s: s["branches"] > 10),
    ("very-branchy", 2, lambda s: s["branches"] > 20),
    ("deep-nesting", 1, lambda s: s["max_nesting"] > 3),
    ("very-deep-nesting", 1, lambda s: s["max_nesting"] > 5),
    ("many-parameters", 1, lambda s: s["params"] > 5),
    ("very-many-parameters", 1, lambda s: s["params"] > 8),
]


def collect_hotspots(
    symbols: list[dict],
    callers: list[dict],
    duplicates: list[dict],
    files: list[Path],
    repo_root: Path,
) -> dict:
    """Rank symbols and files by how much they warrant a reviewer's attention."""
    fan_in = {c["qualname"]: c["sites"] for c in callers}
    cloned = {member for cluster in duplicates for member in cluster["members"]}

    ranked: list[dict] = []
    for record in symbols:
        signals: list[str] = []
        score = 0
        for name, weight, predicate in _SIGNALS:
            if predicate(record):
                signals.append(name)
                score += weight
        if record["exported"] and not record["has_docstring"]:
            signals.append("undocumented-export")
            score += 3
        if f"{record['file']}:{record['lineno']}" in cloned:
            signals.append("clone")
            score += 3
        if record["exported"] and not record["is_test"] and fan_in.get(record["qualname"], 0) == 0:
            signals.append("no-static-consumer")
            score += 4
        ranked.append(
            {
                "qualname": record["qualname"],
                "file": record["file"],
                "lineno": record["lineno"],
                "priority": score,
                "signals": sorted(signals),
            }
        )
    ranked.sort(key=lambda r: (-r["priority"], r["qualname"]))

    per_file: dict[str, int] = {}
    for record in ranked:
        per_file[record["file"]] = per_file.get(record["file"], 0) + record["priority"]
    file_rows: list[dict] = []
    for path in files:
        rel = str(path.resolve().relative_to(repo_root.resolve()))
        lines = len(path.read_text(encoding="utf-8").splitlines())
        score = per_file.get(rel, 0) + 2 * (lines > 400) + 2 * (lines > 800)
        file_rows.append({"file": rel, "lines": lines, "priority": score})
    file_rows.sort(key=lambda r: (-r["priority"], r["file"]))

    return {"symbols": ranked, "files": file_rows}
