"""Scaffolding tests do not reach `main`.

Developing a change means writing far more tests than the change should keep:
probes that pin behaviour long enough to understand it, duplicates of coverage
that turns out to exist already, checks of a property that cannot regress.
Writing them is right. Committing them is not — each one that stays has to be
paid for on every run, in every review, and by whoever later changes the code
it over-specifies.

The `wip` marker makes that decision explicit rather than tacit. A test marked
`@pytest.mark.wip` is declared scaffolding at the moment it is written, when
its author knows; this guard then refuses to let it merge, so the choice is
made deliberately rather than by forgetting. Removing the marker is the act of
claiming a test earns its keep.

Unmarked tests are not exempt from the judgement — they are just not caught by
a grep. The review checkpoint is `mastering-github`'s `/review-ready` Step 3.
"""

from __future__ import annotations

import re
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
_MARKER = re.compile(r"^\s*@pytest\.mark\.wip\b", re.MULTILINE)
_SEARCH = ("libs", "apps", "tests", "validation", "projects", "benchmarks")


def test_no_test_is_still_marked_wip() -> None:
    offenders = []
    for root in _SEARCH:
        for path in (_REPO / root).rglob("test_*.py"):
            if path.resolve() == Path(__file__).resolve():
                continue
            if _MARKER.search(path.read_text(encoding="utf-8")):
                offenders.append(str(path.relative_to(_REPO)))
    assert not offenders, (
        "these tests are still marked `wip`, i.e. declared scaffolding by their own author. "
        "Decide: drop the marker if the test earns its place on main, or delete it. "
        f"{sorted(offenders)}"
    )
