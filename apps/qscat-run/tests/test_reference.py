"""The `reference:` config block: overlay a published dataset on a computed
cross section.

The loader reads a data file BY PATH. It must never import `validation` --
`test_no_validation_import.py` enforces that layering, and the whole point of
naming the file in the config is that qscat_run stays independent of it.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from qscat_run.config import ConfigError, load_config, validate_config
from qscat_run.reference import load_reference

REPO_ROOT = Path(__file__).resolve().parents[3]
HOUFEK = REPO_ROOT / "validation" / "n2" / "data" / "CSVE.V00.J00"


def _write(tmp_path: Path, body: str) -> Path:
    p = tmp_path / "cfg.yaml"
    p.write_text(body)
    return p


BASE = """
molecule: N2
methods: [ti]
observables: [{kind: ve, channels: 2}]
output_dir: runs/x
"""


def test_config_without_reference_has_an_empty_tuple(tmp_path: Path) -> None:
    cfg = load_config(_write(tmp_path, BASE))
    assert cfg.reference == ()


def test_reference_block_is_parsed(tmp_path: Path) -> None:
    cfg = load_config(
        _write(
            tmp_path,
            BASE
            + f"""
reference:
  - path: {HOUFEK}
    format: houfek
    label: "Houfek (2006)"
    channels: [0, 1]
""",
        )
    )
    assert len(cfg.reference) == 1
    ref = cfg.reference[0]
    assert ref.format == "houfek"
    assert ref.label == "Houfek (2006)"
    assert ref.channels == (0, 1)


def test_unknown_format_is_an_actionable_config_error(tmp_path: Path) -> None:
    cfg = load_config(
        _write(tmp_path, BASE + f"\nreference:\n  - path: {HOUFEK}\n    format: nope\n")
    )
    with pytest.raises(ConfigError, match="nope"):
        validate_config(cfg)


def test_missing_file_fails_at_validate_time_not_at_plot_time(tmp_path: Path) -> None:
    cfg = load_config(
        _write(tmp_path, BASE + "\nreference:\n  - path: no/such/file.dat\n    format: houfek\n")
    )
    with pytest.raises(ConfigError, match="no/such/file.dat"):
        validate_config(cfg)


def test_bad_channel_index_fails_at_validate_time_not_after_the_solve(tmp_path: Path) -> None:
    """A typo'd `channels` index against a real file must be caught by
    `validate_config` itself -- not merely raise somewhere downstream (e.g.
    inside `load_reference`, which `write_artifacts` only reaches AFTER
    `run_experiment` has already solved)."""
    (tmp_path / "r.dat").write_text("0.1 1.0 2.0\n0.2 3.0 4.0\n")
    cfg = load_config(
        _write(
            tmp_path,
            BASE + "\nreference:\n  - path: r.dat\n    format: houfek\n    channels: [5]\n",
        )
    )
    with pytest.raises(ConfigError, match=r"\[5\]"):
        validate_config(cfg)


def test_relative_path_resolves_against_the_config_file(tmp_path: Path) -> None:
    (tmp_path / "data").mkdir()
    (tmp_path / "data" / "r.dat").write_text("0.1 1.0 2.0\n0.2 3.0 4.0\n")
    cfg = load_config(
        _write(tmp_path, BASE + "\nreference:\n  - path: data/r.dat\n    format: houfek\n")
    )
    validate_config(cfg)  # must not raise
    series = load_reference(cfg.reference[0], tmp_path)
    assert set(series) == {"ref:ve:ch0", "ref:ve:ch1"}
    energy, sigma = series["ref:ve:ch0"]
    assert energy.tolist() == [0.1, 0.2]
    assert sigma.tolist() == [1.0, 3.0]


def test_houfek_loader_reads_the_committed_dataset() -> None:
    from qscat_run.reference import ReferenceSpec

    spec = ReferenceSpec(path=str(HOUFEK), format="houfek", label=None, channels=(0, 1, 2))
    series = load_reference(spec, REPO_ROOT)
    assert set(series) == {"ref:ve:ch0", "ref:ve:ch1", "ref:ve:ch2"}
    energy, sigma = series["ref:ve:ch0"]
    assert energy.shape == (400,)
    assert sigma.shape == (400,)
    assert np.all(np.diff(energy) > 0.0)
    assert np.all(sigma >= 0.0)


def test_channels_omitted_loads_every_column(tmp_path: Path) -> None:
    (tmp_path / "r.dat").write_text("0.1 1.0 2.0 3.0\n0.2 4.0 5.0 6.0\n")
    from qscat_run.reference import ReferenceSpec

    spec = ReferenceSpec(path="r.dat", format="houfek", label=None, channels=None)
    series = load_reference(spec, tmp_path)
    assert set(series) == {"ref:ve:ch0", "ref:ve:ch1", "ref:ve:ch2"}
