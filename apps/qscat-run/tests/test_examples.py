"""Coverage for the committed `apps/qscat-run/examples/*.yaml` configs.

Every example must schema-validate (fast, no solve) -- these are the
newcomer-facing sample configs the README points at, so a stale/broken
example would be a silent doc regression. `n2-ve.yaml` and `f2-da.yaml` are
additionally run end-to-end (`run_experiment` + `write_artifacts` for the
former, `run_experiment` for the latter): both use a deliberately tiny
explicit grid (see each file's own comments), so together they cost ~1 s and
stay in the fast gate. `h2p-dr.yaml` is the exception -- a multi-hundred-step
propagation on a real (if reduced) preset deck -- so it is schema-validated
here and solved only under `@slow`.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pytest
from qscat_run.artifacts import write_artifacts
from qscat_run.config import load_config, validate_config
from qscat_run.runner import run_experiment

EXAMPLES_DIR = Path(__file__).resolve().parent.parent / "examples"
EXAMPLE_PATHS = sorted(EXAMPLES_DIR.glob("*.yaml"))


def test_examples_directory_is_not_empty() -> None:
    # A guard against a typo'd glob silently collecting zero files and every
    # parametrized test below vacuously "passing".
    assert EXAMPLE_PATHS, f"no *.yaml files found under {EXAMPLES_DIR}"


@pytest.mark.parametrize("path", EXAMPLE_PATHS, ids=lambda p: p.name)
def test_example_validates_clean(path: Path) -> None:
    cfg = load_config(path)
    validate_config(cfg)  # raises ConfigError on any problem


def test_n2_ve_example_is_the_fast_end_to_end() -> None:
    """The one example small/fast enough to actually solve in the gate:
    `n2-ve.yaml`'s tiny explicit grid + two energies. Asserts the
    cross_section + manifest + resolved-config artifacts exist and are
    well-formed/finite."""
    cfg = load_config(EXAMPLES_DIR / "n2-ve.yaml")
    validate_config(cfg)

    result = run_experiment(cfg)
    assert result.cross_sections
    assert "ti:ve:v0->0" in result.cross_sections
    assert "ti:ve:v0->1" in result.cross_sections
    for series in result.cross_sections.values():
        assert np.all(np.isfinite(series))

    out_dir = Path(cfg.output_dir)
    assert not out_dir.exists(), (
        f"{out_dir} already exists on disk -- refusing to write into a "
        "pre-existing directory from a test run"
    )
    write_artifacts(result, cfg, out_dir, timestamp="2026-01-01T00:00:00")
    try:
        assert (out_dir / "cross_section.npz").exists()
        assert (out_dir / "cross_section.npz").exists()
        assert (out_dir / "cross_section.png").exists()
        assert (out_dir / "manifest.json").exists()
        assert (out_dir / "config.resolved.yaml").exists()

        arr = np.load(out_dir / "cross_section.npz")
        assert np.all(np.isfinite(arr["energy"]))
        for key in result.cross_sections:
            assert key in arr
            assert np.all(np.isfinite(arr[key]))
    finally:
        shutil.rmtree(out_dir, ignore_errors=True)


def test_f2_da_example_runs_end_to_end() -> None:
    """The mixed VE+DA, TI+TD example: the one config that exercises both
    observable kinds and both methods in a single run, including the
    per-kind `test_function` mapping a mixed run needs. Like `n2-ve.yaml`
    it is a plumbing smoke test on a deliberately tiny explicit grid, not a
    converged cross section, so it solves in ~1 s and belongs in the fast
    gate (docs/adr/0005 point 7)."""
    cfg = load_config(EXAMPLES_DIR / "f2-da.yaml")
    validate_config(cfg)
    result = run_experiment(cfg)
    assert result.cross_sections


@pytest.mark.slow
def test_h2p_dr_example_runs_end_to_end() -> None:
    """The H2+ proxy-grid DR example -- a multi-hundred-step TD propagation
    on a real (if reduced) preset deck; documented `@slow`, meant for Docker/
    a deliberate local run, not the fast gate."""
    cfg = load_config(EXAMPLES_DIR / "h2p-dr.yaml")
    validate_config(cfg)
    result = run_experiment(cfg)
    assert result.cross_sections


def test_the_o2_examples_specify_their_mesh_as_ranges_not_as_points() -> None:
    """The O2 sweeps are the reason `ranges` exists, so they are pinned to it.

    Point-by-point these three configs were 34 kB each and 3343 lines of
    energies -- a file nobody reads, in which the background step, the window
    width and the level each point belongs to are all invisible. Nothing stops
    a future regeneration from emitting `values` again, and the file would look
    plausible; this notices.
    """
    import yaml

    examples = Path(__file__).resolve().parents[1] / "examples"
    o2 = sorted(examples.glob("o2*-ve.yaml"))
    assert o2, "expected the O2 VE example configs"
    for path in o2:
        energies = yaml.safe_load(path.read_text())["energies"]
        assert "ranges" in energies, f"{path.name} should specify ranges"
        assert "values" not in energies, f"{path.name} should not list energies point by point"
        # the background sweep plus one window per level in the studied range
        assert len(energies["ranges"]) > 1, f"{path.name}: a level-aware mesh is not one segment"
