"""Config schema: YAML -> `ExperimentConfig`, and `validate_config`'s
actionable rejections (unknown molecule, invalid observable-per-molecule,
missing `td` block, `dr` on a non-H2P molecule, explicit grid missing a
half, unknown extractor, empty observables)."""

from __future__ import annotations

import textwrap
from pathlib import Path

import pytest
from qscat_run.config import ConfigError, load_config, validate_config


def _write(tmp_path: Path, text: str) -> Path:
    p = tmp_path / "c.yaml"
    p.write_text(textwrap.dedent(text))
    return p


def test_minimal_config_loads_and_resolves(tmp_path: Path) -> None:
    cfg = load_config(
        _write(
            tmp_path,
            """
        molecule: F2
        methods: [ti]
        observables: [{kind: ve, channels: 2}, {kind: da, channels: 1}]
        output_dir: out
    """,
        )
    )
    validate_config(cfg)  # no raise
    assert cfg.molecule == "F2"
    assert [o.kind for o in cfg.observables] == ["ve", "da"]


def test_dr_on_neutral_rejected(tmp_path: Path) -> None:
    cfg = load_config(
        _write(
            tmp_path,
            """
        molecule: F2
        methods: [ti]
        observables: [{kind: dr, channels: 3}]
        output_dir: out
    """,
        )
    )
    with pytest.raises(ConfigError, match="dr.*H2P|not valid for F2"):
        validate_config(cfg)


def test_td_method_requires_td_block(tmp_path: Path) -> None:
    cfg = load_config(
        _write(
            tmp_path,
            """
        molecule: F2
        methods: [td]
        observables: [{kind: da, channels: 1}]
        output_dir: out
    """,
        )
    )
    with pytest.raises(ConfigError, match="td"):
        validate_config(cfg)


def test_unknown_molecule_rejected(tmp_path: Path) -> None:
    cfg = load_config(
        _write(
            tmp_path,
            """
        molecule: XX
        methods: [ti]
        observables: [{kind: ve, channels: 2}]
        output_dir: out
    """,
        )
    )
    with pytest.raises(ConfigError, match="XX"):
        validate_config(cfg)


def test_empty_observables_rejected(tmp_path: Path) -> None:
    cfg = load_config(
        _write(
            tmp_path,
            """
        molecule: F2
        methods: [ti]
        observables: []
        output_dir: out
    """,
        )
    )
    with pytest.raises(ConfigError, match="observable"):
        validate_config(cfg)


def test_n2_da_allowed_with_warning(tmp_path: Path) -> None:
    cfg = load_config(
        _write(
            tmp_path,
            """
        molecule: N2
        methods: [ti]
        observables: [{kind: da, channels: 1}]
        output_dir: out
    """,
        )
    )
    with pytest.warns(UserWarning, match="closed"):
        validate_config(cfg)


def test_ve_on_h2p_rejected(tmp_path: Path) -> None:
    cfg = load_config(
        _write(
            tmp_path,
            """
        molecule: H2P
        methods: [ti]
        observables: [{kind: ve, channels: 2}]
        output_dir: out
    """,
        )
    )
    with pytest.raises(ConfigError, match="ve"):
        validate_config(cfg)


def test_explicit_grid_missing_nuclear_rejected(tmp_path: Path) -> None:
    cfg = load_config(
        _write(
            tmp_path,
            """
        molecule: F2
        methods: [ti]
        observables: [{kind: ve, channels: 2}]
        grid:
          electronic: {real: [[9, 1.8]], ecs: {angle: 35, elements: 8, quadrature: 14}}
        output_dir: out
    """,
        )
    )
    with pytest.raises(ConfigError, match="nuclear"):
        validate_config(cfg)


def test_unknown_extractor_rejected(tmp_path: Path) -> None:
    cfg = load_config(
        _write(
            tmp_path,
            """
        molecule: F2
        methods: [td]
        observables: [{kind: da, channels: 1}]
        td: {dt: 1.0, n_steps: 10, order: 3, extractors: [bogus]}
        output_dir: out
    """,
        )
    )
    with pytest.raises(ConfigError, match="bogus"):
        validate_config(cfg)


def test_dr_on_h2p_allowed(tmp_path: Path) -> None:
    cfg = load_config(
        _write(
            tmp_path,
            """
        molecule: H2P
        methods: [ti]
        observables: [{kind: dr, channels: 3}]
        output_dir: out
    """,
        )
    )
    validate_config(cfg)


def test_energies_values_form(tmp_path: Path) -> None:
    cfg = load_config(
        _write(
            tmp_path,
            """
        molecule: F2
        methods: [ti]
        observables: [{kind: ve, channels: 2}]
        energies: {values: [0.03, 0.04]}
        output_dir: out
    """,
        )
    )
    validate_config(cfg)
    assert cfg.energies is not None
    assert cfg.energies.values == (0.03, 0.04)
