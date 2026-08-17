"""Config schema: YAML -> `ExperimentConfig`, and `validate_config`'s
actionable rejections (unknown molecule, invalid observable-per-molecule,
missing `td` block, `dr` on a non-H2P molecule, explicit grid missing a
half, unknown extractor, empty observables)."""

from __future__ import annotations

import textwrap
from pathlib import Path

import numpy as np
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


def test_lcp_method_accepted_for_f2_with_da(tmp_path: Path) -> None:
    cfg = load_config(
        _write(
            tmp_path,
            """
        molecule: F2
        methods: [ti, lcp]
        observables: [{kind: da, channels: 1}]
        output_dir: out
    """,
        )
    )
    validate_config(cfg)  # no raise
    assert "lcp" in cfg.methods


def test_lcp_rejected_for_n2(tmp_path: Path) -> None:
    cfg = load_config(
        _write(
            tmp_path,
            """
        molecule: N2
        methods: [lcp]
        observables: [{kind: da, channels: 1}]
        output_dir: out
    """,
        )
    )
    with pytest.raises(ConfigError, match="lcp.*not available"):
        validate_config(cfg)


def test_lcp_without_da_observable_rejected(tmp_path: Path) -> None:
    cfg = load_config(
        _write(
            tmp_path,
            """
        molecule: F2
        methods: [lcp]
        observables: [{kind: ve, channels: 2}]
        output_dir: out
    """,
        )
    )
    with pytest.raises(ConfigError, match="no 'da' or 'resonance_levels' observable"):
        validate_config(cfg)


def test_lcp_with_explicit_grid_rejected(tmp_path: Path) -> None:
    cfg = load_config(
        _write(
            tmp_path,
            """
        molecule: F2
        methods: [lcp]
        observables: [{kind: da, channels: 1}]
        output_dir: out
        grid:
          electronic: {real: [[4, 2.0], [2, 8.0]], ecs: {angle: 30, elements: 3, quadrature: 5}}
          nuclear: {real: [[3, 1.8], [3, 10.0]], ecs: {angle: 30, elements: 3, quadrature: 6}}
    """,
        )
    )
    with pytest.raises(ConfigError, match="lcp.*does not support an explicit grid"):
        validate_config(cfg)


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


def test_resonance_levels_observable_needs_no_energies(tmp_path):
    from qscat_run import presets
    from qscat_run.config import load_config, validate_config

    cfg_path = tmp_path / "levels.yaml"
    cfg_path.write_text(
        "molecule: F2\n"
        "methods: [lcp]\n"
        "observables:\n"
        "  - kind: resonance_levels\n"
        "    channels: 6\n"
        "output_dir: out\n"
    )
    cfg = load_config(cfg_path)
    validate_config(cfg)
    resolved = presets.resolve_defaults(cfg)
    assert resolved.observables[0].kind == "resonance_levels"
    assert resolved.energies is None  # levels-only run: no sweep to resolve


def test_resonance_levels_is_rejected_for_n2(tmp_path):
    from qscat_run.config import ConfigError, load_config, validate_config

    cfg_path = tmp_path / "bad.yaml"
    cfg_path.write_text(
        "molecule: N2\nmethods: [lcp]\nobservables: [{kind: resonance_levels}]\noutput_dir: out\n"
    )
    with pytest.raises(ConfigError, match="not valid for N2"):
        validate_config(load_config(cfg_path))


def test_lcp_without_da_is_still_rejected_when_no_levels_requested(tmp_path):
    from qscat_run.config import ConfigError, load_config, validate_config

    cfg_path = tmp_path / "bad2.yaml"
    cfg_path.write_text(
        "molecule: F2\nmethods: [lcp]\nobservables: [{kind: ve, channels: 2}]\n"
        "energies: {min: 0.01, max: 0.05, step: 0.02}\n"
        "output_dir: out\n"
    )
    with pytest.raises(ConfigError, match="no 'da' or 'resonance_levels' observable"):
        validate_config(load_config(cfg_path))


def test_artifacts_resonance_levels_flag_parses(tmp_path):
    from qscat_run.config import load_config

    cfg_path = tmp_path / "art.yaml"
    cfg_path.write_text(
        "molecule: F2\nmethods: [lcp]\n"
        "observables: [{kind: da, channels: 1}]\n"
        "energies: {min: 0.01, max: 0.05, step: 0.02}\n"
        "artifacts: {resonance_levels: true}\n"
        "output_dir: out\n"
    )
    assert load_config(cfg_path).artifacts.resonance_levels is True


def test_nuclear_angle_b_defaults_to_ten_degrees_below(tmp_path):
    from qscat_run import presets
    from qscat_run.config import load_config

    cfg_path = tmp_path / "ang.yaml"
    cfg_path.write_text(
        "molecule: F2\nmethods: [lcp]\nobservables: [{kind: resonance_levels}]\noutput_dir: out\n"
    )
    cfg = presets.resolve_defaults(load_config(cfg_path))
    g_a, _ea, _eb = presets.resolve_lcp_grids(cfg)
    g_b = presets.nuclear_grid_at_angle(cfg, presets.nuclear_angle_b(cfg))
    ang_a = max(el.angle_deg for el in g_a.spec.elements)
    ang_b = max(el.angle_deg for el in g_b.spec.elements)
    assert ang_b == pytest.approx(ang_a - 10.0)
    # The real nodes must be shared -- that is what makes the two spectra
    # comparable in `lcp_resonance_levels`.
    ra, rb = g_a.points, g_b.points
    assert np.array_equal(ra[ra.imag == 0.0], rb[rb.imag == 0.0])


def test_omitted_channels_on_resonance_levels_means_all_levels(tmp_path):
    """`channels` is the number of levels to report; omitting it must mean
    "every angle-stable level in the default window" (`n_levels=None`), which
    is what the design spec and `examples/f2-resonance-levels.yaml`'s own
    comment promise. It used to resolve to the `da`/`dr` default of 1, so an
    omitted `channels` silently reported a single level."""
    from qscat_run import presets
    from qscat_run.config import load_config, validate_config

    cfg_path = tmp_path / "levels.yaml"
    cfg_path.write_text(
        "molecule: F2\nmethods: [lcp]\nobservables: [{kind: resonance_levels}]\noutput_dir: out\n"
    )
    cfg = load_config(cfg_path)
    validate_config(cfg)
    resolved = presets.resolve_defaults(cfg)
    assert resolved.observables[0].channels is None


def test_omitted_channels_still_defaults_for_da(tmp_path):
    """The `None` default is `resonance_levels`-specific: `da` keeps its
    single-channel default."""
    from qscat_run import presets
    from qscat_run.config import load_config, validate_config

    cfg_path = tmp_path / "da.yaml"
    cfg_path.write_text(
        "molecule: F2\n"
        "methods: [lcp]\n"
        "observables: [{kind: da}]\n"
        "energies: {min: 0.02, max: 0.04, step: 0.02}\n"
        "output_dir: out\n"
    )
    cfg = load_config(cfg_path)
    validate_config(cfg)
    resolved = presets.resolve_defaults(cfg)
    assert resolved.observables[0].channels == 1


def test_bare_NO_molecule_gives_an_actionable_yaml_boolean_error(tmp_path: Path) -> None:
    """PyYAML is YAML 1.1, where bare `NO` is the boolean false -- so a
    hand-written `molecule: NO` used to fail with "unknown molecule 'False'",
    which says nothing about the real cause. `qscat-run init` quotes the name,
    so only hand-written configs hit this.
    """
    p = tmp_path / "no.yaml"
    p.write_text(
        "molecule: NO\nmethods: [ti]\nobservables: [{kind: da, channels: 1}]\noutput_dir: runs/x\n"
    )
    with pytest.raises(ConfigError, match="YAML boolean"):
        load_config(p)


def test_quoted_NO_molecule_loads(tmp_path: Path) -> None:
    p = tmp_path / "no.yaml"
    p.write_text(
        "molecule: 'NO'\nmethods: [ti]\n"
        "observables: [{kind: da, channels: 1}]\noutput_dir: runs/x\n"
    )
    assert load_config(p).molecule == "NO"
