"""The `resonance_levels` observable + artifact flag through the runner."""

from __future__ import annotations

import numpy as np
import pytest
from qscat_run import presets
from qscat_run.config import load_config, validate_config
from qscat_run.runner import run_experiment


def _resolved(path):
    cfg = load_config(path)
    validate_config(cfg)
    return presets.resolve_defaults(cfg)


LEVELS_ONLY = (
    "molecule: F2\n"
    "methods: [lcp]\n"
    "observables:\n"
    "  - kind: resonance_levels\n"
    "    channels: 4\n"
    "output_dir: out\n"
)


@pytest.mark.slow
def test_levels_only_run_produces_levels_and_no_cross_sections(tmp_path):
    cfg_path = tmp_path / "levels.yaml"
    cfg_path.write_text(LEVELS_ONLY)
    result = run_experiment(_resolved(cfg_path))

    assert len(result.resonance_levels) == 1
    run = result.resonance_levels[0]
    assert run.label == "lcp:resonance_levels"
    assert run.levels.energies.size == 4
    assert np.all(np.diff(run.levels.energies.real) > 0)
    assert run.levels.widths.size == 4
    assert run.R_axis.size == run.Vd.size == run.Gamma.size
    assert not result.cross_sections  # no sweep was requested or run


@pytest.mark.slow
def test_artifact_flag_adds_levels_to_a_da_run(tmp_path):
    cfg_path = tmp_path / "da.yaml"
    cfg_path.write_text(
        "molecule: F2\n"
        "methods: [lcp]\n"
        "observables: [{kind: da, channels: 1}]\n"
        "energies: {min: 0.02, max: 0.04, step: 0.02}\n"
        "artifacts: {resonance_levels: true}\n"
        "output_dir: out\n"
    )
    result = run_experiment(_resolved(cfg_path))
    assert "lcp:da:ch0" in result.cross_sections
    assert len(result.resonance_levels) == 1


@pytest.mark.slow
def test_omitted_channels_reports_every_selected_level(tmp_path):
    """The end-to-end consequence of `channels: None`: the run reports every
    angle-stable level the default window selects, not one and not four."""
    cfg_path = tmp_path / "all.yaml"
    cfg_path.write_text(
        "molecule: F2\nmethods: [lcp]\nobservables: [{kind: resonance_levels}]\noutput_dir: out\n"
    )
    result = run_experiment(_resolved(cfg_path))
    levels = result.resonance_levels[0].levels
    assert levels.energies.size > 4  # not the old `channels -> 1` default
    assert levels.energies.size == levels.widths.size == levels.residuals.size
