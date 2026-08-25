"""The `validate:*` suite map in .github/workflows/validation.yml must cover
every `slow` test.

The slow tier does not run in the default gate (docs/adr/0005), so the only
places it runs are a developer's machine, the Docker `test` image, and the
label-triggered Validation workflow. A slow test in a directory no suite
selects is therefore one that CI can never run, however many labels a reviewer
applies -- and nothing about it would look wrong: it collects, it passes
locally, and it silently never gates anything. (When this test was written it
found exactly that: `apps/qscat-run/tests` held 12 slow tests and no suite
named it.)

Collection runs ONCE, for the whole repo (`_all_slow_node_ids` is cached), and
suites are resolved by path prefix in-process. Measured on this tree: one
whole-repo `--collect-only -m slow` costs 1.7 s and finds every slow id, while
running the same collection once per suite costs 5.35 s in total (six
subprocesses, each re-paying interpreter start-up and conftest import) for the
same answer. Neither number is alarming on its own -- the point is that a guard
policing the tiering should be the cheapest thing in the tier, not something
that needs an exemption from it.
"""

from __future__ import annotations

import ast
import functools
import json
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
WORKFLOW = REPO / ".github" / "workflows" / "validation.yml"


def _suite_map() -> dict[str, str]:
    """The SUITES dict as the workflow's own selection script defines it."""
    match = re.search(r"SUITES = \{(.*?)\n\s*\}", WORKFLOW.read_text(), re.S)
    assert match, "could not find the SUITES map in the validation workflow"
    suites = ast.literal_eval("{" + match.group(1).replace("\n", " ") + "}")
    assert suites, "the SUITES map is empty"
    return suites


@functools.cache
def _all_slow_node_ids() -> tuple[str, ...]:
    """Every `slow` node id in the repo, collected once per test session.

    Cached because two tests below need the same answer and the collection is
    a subprocess: without this it ran twice for one identical result. A tuple
    rather than a list so a caller cannot mutate the cached value.
    """
    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "-m",
            "slow",
            "--collect-only",
            "-p",
            "no:cacheprovider",
        ],
        cwd=REPO,
        capture_output=True,
        text=True,
    )
    assert proc.returncode in (0, 5), f"collection failed:\n{proc.stdout}\n{proc.stderr}"
    ids = tuple(ln.strip() for ln in proc.stdout.splitlines() if "::" in ln)
    assert ids, "no slow tests collected at all -- the marker or the run is broken"
    return ids


def _prefixes(suites: dict[str, str]) -> dict[str, tuple[str, ...]]:
    return {name: tuple(paths.split()) for name, paths in suites.items()}


def test_every_slow_test_is_selected_by_some_suite() -> None:
    suites = _prefixes(_suite_map())
    every_prefix = tuple(p for paths in suites.values() for p in paths)
    orphans = [nid for nid in _all_slow_node_ids() if not nid.startswith(every_prefix)]
    assert not orphans, (
        f"{len(orphans)} slow test(s) are in no `validate:*` suite, so the Validation "
        "workflow can never run them:\n  "
        + "\n  ".join(sorted(orphans))
        + f"\n\nAdd their directory to SUITES in {WORKFLOW.relative_to(REPO)}."
    )


def test_no_suite_is_empty() -> None:
    """A suite selecting nothing is a typo'd path: the label appears to work,
    runs zero tests, and reports success."""
    node_ids = _all_slow_node_ids()
    for suite, prefixes in _prefixes(_suite_map()).items():
        assert any(nid.startswith(prefixes) for nid in node_ids), (
            f"suite 'validate:{suite}' selects no slow tests (paths: {prefixes}). "
            "Either the path is wrong or the suite should be removed."
        )


def test_suite_paths_exist() -> None:
    for suite, paths in _suite_map().items():
        for path in paths.split():
            assert (REPO / path).exists(), f"validate:{suite} points at missing path {path!r}"


def test_documented_labels_match_the_suite_map() -> None:
    """The header comment lists the labels a reviewer is told to apply. If it
    drifts from SUITES, reviewers apply labels that silently do nothing."""
    documented = set(re.findall(r"^#\s+validate:(\w+)\s{2,}", WORKFLOW.read_text(), re.M))
    documented.discard("all")  # the aggregate, not a suite
    assert documented == set(_suite_map()), (
        f"header comment documents {sorted(documented)} but SUITES defines {sorted(_suite_map())}"
    )


def _selection_script() -> str:
    match = re.search(
        r"python3 - <<'PY' >> \"\$GITHUB_OUTPUT\"\n(.*?)\n\s*PY\n", WORKFLOW.read_text(), re.S
    )
    assert match, "could not extract the selection script"
    return "\n".join(
        line[10:] if line.startswith(" " * 10) else line for line in match.group(1).split("\n")
    )


def test_selection_script_behaves() -> None:
    script = _selection_script()
    ast.parse(script)  # raises SyntaxError if the heredoc is malformed

    cases = [
        ({"LABELS": "validate:all", "DISPATCH": ""}, set(_suite_map())),
        ({"LABELS": "", "DISPATCH": ""}, set()),
        ({"LABELS": "bug docs", "DISPATCH": ""}, set()),
        ({"LABELS": "validate:core", "DISPATCH": ""}, {"core"}),
        ({"LABELS": "validate:nope", "DISPATCH": ""}, set()),
        ({"LABELS": "", "DISPATCH": "all"}, set(_suite_map())),
        ({"LABELS": "validate:all", "DISPATCH": "core"}, {"core"}),  # dispatch wins
    ]
    for env, expect in cases:
        proc = subprocess.run(
            [sys.executable, "-c", script],
            env={"PATH": "/usr/bin:/bin", **env},
            capture_output=True,
            text=True,
        )
        assert proc.returncode == 0, proc.stderr
        out = dict(ln.split("=", 1) for ln in proc.stdout.strip().splitlines())
        got = {e["suite"] for e in json.loads(out["matrix"])["include"]}
        assert got == expect, f"{env} -> {got}, expected {expect}"
        assert out["any"] == ("true" if expect else "false")
