"""Coverage for the committed `apps/qscat-run/examples/*.yaml` configs.

Every example must schema-validate (fast, no solve) -- these are the
newcomer-facing sample configs the README points at, so a stale/broken
example would be a silent doc regression. `n2-ve.yaml` is additionally run
end-to-end (`run_experiment` + `write_artifacts`): it uses a deliberately
tiny explicit grid (see the file's own comments) so this costs a fraction of
a second, unlike `f2-da.yaml`/`h2p-dr.yaml` which are documented,
runnable-in-Docker examples on a real (if reduced) preset deck -- schema-
validated here but not solved, per the module docstring's `@slow` guidance.
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
        assert (out_dir / "cross_section.csv").exists()
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


@pytest.mark.slow
def test_f2_da_example_runs_end_to_end() -> None:
    """The heavier mixed VE+DA, TI+TD example -- documented `@slow` (not run
    in the fast gate; run explicitly with `pytest -m slow` or via
    `docker/run.sh`)."""
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
