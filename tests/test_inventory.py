"""The structure-audit mapper, driven as a subprocess over a fixture tree.

The mapper lives at .claude/skills/code-mapping/scripts/inventory.py, outside any
importable package, so these tests invoke its CLI rather than importing it.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / ".claude/skills/code-mapping/scripts/inventory.py"
FIXTURE = REPO_ROOT / "tests/fixtures/inventory_sample"


def run_inventory(out_dir: Path, root: Path = FIXTURE) -> Path:
    """Run the mapper over `root`, writing tables into `out_dir`."""
    out_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--root",
            str(root),
            "--out",
            str(out_dir),
            "--repo-root",
            str(REPO_ROOT),
        ],
        check=True,
        capture_output=True,
    )
    return out_dir


def load(out_dir: Path, name: str) -> object:
    """Read one emitted table."""
    return json.loads((out_dir / f"{name}.json").read_text())


def by_qualname(symbols: list[dict], suffix: str) -> dict:
    """Return the single symbol whose qualname ends with `suffix`."""
    hits = [s for s in symbols if s["qualname"].endswith(suffix)]
    assert len(hits) == 1, f"{suffix}: expected 1, got {[h['qualname'] for h in hits]}"
    return hits[0]


def test_symbols_records_location_and_shape(tmp_path: Path) -> None:
    out = run_inventory(tmp_path / "out")
    symbols = load(out, "symbols")
    clone = by_qualname(symbols, "alpha.core.clone_a")
    assert clone["kind"] == "function"
    assert clone["file"] == "tests/fixtures/inventory_sample/alpha/core.py"
    assert clone["params"] == 1
    assert clone["branches"] == 2  # one for, one if
    assert clone["max_nesting"] == 2
    assert clone["has_docstring"] is True
    assert clone["returns"] == 1


def test_symbols_flags_missing_docstring(tmp_path: Path) -> None:
    out = run_inventory(tmp_path / "out")
    symbols = load(out, "symbols")
    assert by_qualname(symbols, "alpha.core.undocumented_export")["has_docstring"] is False


def test_symbols_marks_all_exports(tmp_path: Path) -> None:
    out = run_inventory(tmp_path / "out")
    symbols = load(out, "symbols")
    assert by_qualname(symbols, "alpha.core.shared_helper")["exported"] is True
    assert by_qualname(symbols, "alpha.core._private_dead")["exported"] is False


def test_symbols_counts_parameters(tmp_path: Path) -> None:
    out = run_inventory(tmp_path / "out")
    assert by_qualname(load(out, "symbols"), "alpha.core.wide")["params"] == 9


def test_symbols_records_methods_with_their_class(tmp_path: Path) -> None:
    out = run_inventory(tmp_path / "out")
    symbols = load(out, "symbols")
    method = by_qualname(symbols, "alpha.holders.LevelsHolder.save")
    assert method["kind"] == "method"


def test_output_is_byte_identical_across_runs(tmp_path: Path) -> None:
    first = run_inventory(tmp_path / "a")
    second = run_inventory(tmp_path / "b")
    for table in ["symbols"]:
        assert (first / f"{table}.json").read_bytes() == (second / f"{table}.json").read_bytes()


def test_imports_records_intra_repo_edges(tmp_path: Path) -> None:
    out = run_inventory(tmp_path / "out")
    edges = load(out, "imports")
    pairs = {(e["importer"], e["imported"]) for e in edges}
    assert (
        "tests.fixtures.inventory_sample.beta.core",
        "tests.fixtures.inventory_sample.alpha.core",
    ) in pairs


def test_callers_reports_the_consuming_packages(tmp_path: Path) -> None:
    out = run_inventory(tmp_path / "out")
    callers = {c["qualname"]: c for c in load(out, "callers")}
    helper = callers["tests.fixtures.inventory_sample.alpha.core.shared_helper"]
    assert helper["sites"] >= 1
    assert "tests" in helper["consumer_packages"]


def test_callers_reports_zero_for_the_string_dispatched_symbol(tmp_path: Path) -> None:
    """The static scan cannot see a YAML-dispatched symbol. Phase 1 resolves it."""
    out = run_inventory(tmp_path / "out")
    callers = {c["qualname"]: c for c in load(out, "callers")}
    assert callers["tests.fixtures.inventory_sample.alpha.core.string_dispatched"]["sites"] == 0


def test_callers_marks_homonyms_ambiguous(tmp_path: Path) -> None:
    out = run_inventory(tmp_path / "out")
    callers = {c["qualname"]: c for c in load(out, "callers")}
    entry = callers["tests.fixtures.inventory_sample.alpha.core.documented_public"]
    assert entry["ambiguous"] is True


def test_imports_resolves_a_relative_import_from_a_package_init(tmp_path: Path) -> None:
    """A package __init__.py re-exporting from a sibling module must produce an edge."""
    out = run_inventory(tmp_path / "out")
    pairs = {(e["importer"], e["imported"]) for e in load(out, "imports")}
    assert (
        "tests.fixtures.inventory_sample.alpha",
        "tests.fixtures.inventory_sample.alpha.core",
    ) in pairs


def test_imports_include_the_library_package_reexports(tmp_path: Path) -> None:
    """qscat/core/__init__.py re-exports heavily; zero edges out of it means broken resolution."""
    out = run_inventory(tmp_path / "out", root=REPO_ROOT / "libs/qscat/qscat")
    out_edges = [e for e in load(out, "imports") if e["importer"].endswith("qscat.core")]
    assert len(out_edges) >= 10


def test_duplicates_finds_the_exact_clone_pair(tmp_path: Path) -> None:
    """clone_a and clone_b are byte-identical bodies in different packages."""
    out = run_inventory(tmp_path / "out")
    pairs = [c for c in load(out, "duplicates") if c["kind"] == "clone" and len(c["members"]) == 2]
    assert len(pairs) == 1
    assert {m.rsplit(":", 1)[0] for m in pairs[0]["members"]} == {
        "tests/fixtures/inventory_sample/alpha/core.py",
        "tests/fixtures/inventory_sample/beta/core.py",
    }


def test_duplicates_finds_the_repeated_holder_method(tmp_path: Path) -> None:
    """The three holders share one `save` body, so they cluster as a clone trio."""
    out = run_inventory(tmp_path / "out")
    trios = [c for c in load(out, "duplicates") if len(c["members"]) == 3]
    assert len(trios) == 1
    assert trios[0]["kind"] == "clone"


def test_duplicates_separates_near_clones_from_clones(tmp_path: Path) -> None:
    """near_clone_a and near_clone_b differ by one comparison operator, nothing else.

    Constants are erased by normalization, so a fixture pair differing only in a
    literal would cluster as an exact clone. The near-clone pair must differ
    STRUCTURALLY — here `<` against `<=` — for this table to have two kinds at all.
    """
    out = run_inventory(tmp_path / "out")
    near = [c for c in load(out, "duplicates") if c["kind"] == "near-clone"]
    assert len(near) == 1
    assert len(near[0]["members"]) == 2


def test_homonyms_lists_the_repeated_name(tmp_path: Path) -> None:
    out = run_inventory(tmp_path / "out")
    homonyms = {h["name"]: h for h in load(out, "homonyms")}
    assert homonyms["documented_public"]["count"] == 2
    assert homonyms["save"]["count"] == 3


def test_holders_lists_every_dataclass_with_its_fields(tmp_path: Path) -> None:
    out = run_inventory(tmp_path / "out")
    holders = {h["qualname"].rsplit(".", 1)[-1]: h for h in load(out, "holders")["holders"]}
    assert holders["LevelsHolder"]["flavour"] == "dataclass"
    assert holders["LevelsHolder"]["fields"] == ["energies", "grid", "label", "widths"]


def test_holders_pairs_overlapping_field_sets(tmp_path: Path) -> None:
    out = run_inventory(tmp_path / "out")
    pairs = load(out, "holders")["field_overlap_pairs"]
    assert pairs, "the three fixture holders share energies+grid+widths"
    shared = {tuple(p["shared"]) for p in pairs}
    assert ("energies", "grid", "widths") in shared


def test_holders_reports_identical_method_echoes(tmp_path: Path) -> None:
    out = run_inventory(tmp_path / "out")
    echoes = {e["method"]: e for e in load(out, "holders")["method_echoes"]}
    assert echoes["save"]["count"] == 3
    assert echoes["save"]["identical"] is True


def test_hotspots_scores_review_priority_not_quality(tmp_path: Path) -> None:
    out = run_inventory(tmp_path / "out")
    ranked = load(out, "hotspots")["symbols"]
    scores = {r["qualname"].rsplit(".", 1)[-1]: r for r in ranked}
    assert scores["wide"]["priority"] > scores["documented_public"]["priority"]
    assert "many-parameters" in scores["wide"]["signals"]


def test_hotspots_flags_an_undocumented_export(tmp_path: Path) -> None:
    out = run_inventory(tmp_path / "out")
    ranked = {r["qualname"].rsplit(".", 1)[-1]: r for r in load(out, "hotspots")["symbols"]}
    assert "undocumented-export" in ranked["undocumented_export"]["signals"]


def test_hotspots_flags_clone_membership(tmp_path: Path) -> None:
    out = run_inventory(tmp_path / "out")
    ranked = {r["qualname"].rsplit(".", 1)[-1]: r for r in load(out, "hotspots")["symbols"]}
    assert "clone" in ranked["clone_a"]["signals"]


def test_hotspots_is_sorted_by_descending_priority(tmp_path: Path) -> None:
    out = run_inventory(tmp_path / "out")
    priorities = [r["priority"] for r in load(out, "hotspots")["symbols"]]
    assert priorities == sorted(priorities, reverse=True)


def test_maps_the_real_library_deterministically_and_quickly(tmp_path: Path) -> None:
    """The whole-repo run is the working case; it must be fast and reproducible."""
    import time

    start = time.monotonic()
    first = run_inventory(tmp_path / "a", root=REPO_ROOT / "libs/qscat/qscat")
    elapsed = time.monotonic() - start
    assert elapsed < 60.0, f"library scan took {elapsed:.1f}s"
    second = run_inventory(tmp_path / "b", root=REPO_ROOT / "libs/qscat/qscat")
    for table in ["symbols", "imports", "callers", "duplicates", "homonyms", "holders", "hotspots"]:
        assert (first / f"{table}.json").read_bytes() == (second / f"{table}.json").read_bytes()


def test_the_real_library_has_the_known_homonyms(tmp_path: Path) -> None:
    """Guards the mapper against silently finding nothing."""
    out = run_inventory(tmp_path / "out", root=REPO_ROOT / "libs/qscat/qscat")
    names = {h["name"] for h in load(out, "homonyms")}
    assert {"da_cross_section", "ve_cross_section"} <= names


def test_real_scans_skip_the_deliberate_defect_fixture(tmp_path: Path) -> None:
    """A scan of tests/ must not report the fixture's planted clones as findings."""
    out = run_inventory(tmp_path / "out", root=REPO_ROOT / "tests")
    files = {s["file"] for s in load(out, "symbols")}
    assert not any(f.startswith("tests/fixtures/") for f in files)
    assert any(f.startswith("tests/") for f in files), "the rest of tests/ must still map"


CAPTURE = REPO_ROOT / ".claude/skills/code-mapping/scripts/capture_observables.py"


def run_capture(out_dir: Path) -> Path:
    """Run the observables capture into `out_dir`."""
    out_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [sys.executable, str(CAPTURE), "--out", str(out_dir)],
        check=True,
        capture_output=True,
        cwd=str(REPO_ROOT),
    )
    return out_dir


def test_capture_is_bit_identical_across_runs(tmp_path: Path) -> None:
    """The lane-B gate is only meaningful if an unchanged tree reproduces exactly."""
    first = run_capture(tmp_path / "a")
    second = run_capture(tmp_path / "b")
    manifest_a = json.loads((first / "manifest.json").read_text())
    manifest_b = json.loads((second / "manifest.json").read_text())
    assert manifest_a == manifest_b
    assert manifest_a, "the manifest must list at least one case"


def test_capture_manifest_covers_the_core_modules(tmp_path: Path) -> None:
    """The two shipped cases must cover the two hot-path modules the audit tracks."""
    out = run_capture(tmp_path / "out")
    covered = {
        m for entry in json.loads((out / "manifest.json").read_text()) for m in entry["modules"]
    }
    assert "qscat.linalg" in covered
    assert "qscat.dvr" in covered
