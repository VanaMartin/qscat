"""Config schema: YAML -> `ExperimentConfig`, and `validate_config`'s
actionable rejections (unknown molecule, invalid observable-per-molecule,
missing `td` block, `dr` on a non-H2P molecule, explicit grid missing a
half, unknown extractor, empty observables)."""

from __future__ import annotations

import textwrap
from pathlib import Path

import numpy as np
import pytest
from qscat_run.config import ConfigError, _load_energies, load_config, validate_config


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


def test_lcp_ve_accepted_for_n2(tmp_path: Path) -> None:
    cfg = load_config(
        _write(
            tmp_path,
            """
        molecule: N2
        methods: [lcp]
        observables: [{kind: ve, channels: 2}]
        output_dir: out
    """,
        )
    )
    validate_config(cfg)  # no raise


def test_lcp_rejected_for_h2p(tmp_path: Path) -> None:
    cfg = load_config(
        _write(
            tmp_path,
            """
        molecule: H2P
        methods: [lcp]
        observables: [{kind: dr, channels: 1}]
        output_dir: out
    """,
        )
    )
    with pytest.raises(ConfigError, match=r"lcp.*not available"):
        validate_config(cfg)


def test_lcp_with_ve_observable_accepted_for_f2(tmp_path: Path) -> None:
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
    validate_config(cfg)  # no raise


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
    with pytest.raises(ConfigError, match=r"lcp.*does not support an explicit grid"):
        validate_config(cfg)


def test_nrm_block_loads_with_both_choices(tmp_path: Path) -> None:
    cfg = load_config(
        _write(
            tmp_path,
            """
        molecule: F2
        methods: [ti, lcp, nrm]
        observables: [{kind: da, channels: 1}]
        nrm: {choices: [a, b], n_states: 80}
        output_dir: out
    """,
        )
    )
    validate_config(cfg)  # no raise
    assert cfg.nrm is not None
    assert cfg.nrm.choices == ("a", "b")
    assert cfg.nrm.n_states == 80


def test_nrm_defaults_are_materialized_by_resolve_defaults(tmp_path: Path) -> None:
    """An omitted `nrm:` block still records what actually ran: `resolve_defaults`
    fills the measured defaults so `config.resolved.yaml` is not `nrm: null`."""
    from qscat_run import presets

    cfg = load_config(
        _write(
            tmp_path,
            """
        molecule: F2
        methods: [nrm]
        observables: [{kind: da, channels: 1}]
        output_dir: out
    """,
        )
    )
    validate_config(cfg)
    assert cfg.nrm is None
    resolved = presets.resolve_defaults(cfg)
    assert resolved.nrm is not None
    assert resolved.nrm.choices == ("b",)
    assert resolved.nrm.n_states == 100


def test_nrm_rejected_for_the_ion(tmp_path: Path) -> None:
    """H2+ has no NRM deck: the model's discrete state and Eq. (60) state sum
    are set up for the neutral diatomics, so the config must say so rather
    than fail deep inside a solve."""
    cfg = load_config(
        _write(
            tmp_path,
            """
        molecule: H2P
        methods: [nrm]
        observables: [{kind: dr, channels: 1}]
        output_dir: out
    """,
        )
    )
    with pytest.raises(ConfigError, match=r"nrm.*not available"):
        validate_config(cfg)


def test_nrm_accepts_a_ve_only_config(tmp_path: Path) -> None:
    """The NRM approximates vibrational excitation as well as DA (PRA 77
    Eq. 28/31/37), including for N2 -- the molecule the committed figure
    compares against Figs. 4 and 8."""
    cfg = load_config(
        _write(
            tmp_path,
            """
        molecule: N2
        methods: [ti, nrm]
        observables: [{kind: ve, channels: 2}]
        nrm: {choices: [a, b], n_states: 100, include_background: false}
        output_dir: out
    """,
        )
    )
    validate_config(cfg)
    assert cfg.nrm is not None
    assert cfg.nrm.choices == ("a", "b")
    assert cfg.nrm.include_background is False


def test_nrm_without_ve_or_da_observable_rejected(tmp_path: Path) -> None:
    """The NRM approximates `ve` and `da` only, so a config asking for neither
    must be rejected rather than silently produce no NRM series."""
    cfg = load_config(
        _write(
            tmp_path,
            """
        molecule: F2
        methods: [lcp, nrm]
        observables: [{kind: resonance_levels, channels: 2}]
        output_dir: out
    """,
        )
    )
    with pytest.raises(ConfigError, match="no 've' or 'da' observable"):
        validate_config(cfg)


def test_nrm_with_explicit_grid_rejected(tmp_path: Path) -> None:
    cfg = load_config(
        _write(
            tmp_path,
            """
        molecule: F2
        methods: [nrm]
        observables: [{kind: da, channels: 1}]
        output_dir: out
        grid:
          electronic: {real: [[4, 2.0], [2, 8.0]], ecs: {angle: 30, elements: 3, quadrature: 5}}
          nuclear: {real: [[3, 1.8], [3, 10.0]], ecs: {angle: 30, elements: 3, quadrature: 6}}
    """,
        )
    )
    with pytest.raises(ConfigError, match=r"nrm.*does not support an explicit grid"):
        validate_config(cfg)


def test_nrm_unknown_choice_rejected(tmp_path: Path) -> None:
    """PRA 77's third ('compact') discrete state is not implemented, so a `c`
    must fail loudly rather than resolve to something else."""
    cfg = load_config(
        _write(
            tmp_path,
            """
        molecule: F2
        methods: [nrm]
        observables: [{kind: da, channels: 1}]
        nrm: {choices: [c]}
        output_dir: out
    """,
        )
    )
    with pytest.raises(ConfigError, match="unknown nrm discrete-state choice"):
        validate_config(cfg)


def test_nrm_nonpositive_n_states_rejected(tmp_path: Path) -> None:
    cfg = load_config(
        _write(
            tmp_path,
            """
        molecule: F2
        methods: [nrm]
        observables: [{kind: da, channels: 1}]
        nrm: {n_states: 0}
        output_dir: out
    """,
        )
    )
    with pytest.raises(ConfigError, match="n_states"):
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
    with pytest.raises(ConfigError, match=r"dr.*H2P|not valid for F2"):
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


def test_lcp_with_ve_only_and_energies_accepted(tmp_path):
    from qscat_run.config import load_config, validate_config

    cfg_path = tmp_path / "ok2.yaml"
    cfg_path.write_text(
        "molecule: F2\nmethods: [lcp]\nobservables: [{kind: ve, channels: 2}]\n"
        "energies: {min: 0.01, max: 0.05, step: 0.02}\n"
        "output_dir: out\n"
    )
    validate_config(load_config(cfg_path))  # no raise


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


class TestMissingSubkeys:
    """A typo'd or missing sub-key must produce an actionable ConfigError
    naming the block and the key -- never a bare KeyError (exp-M4)."""

    def _load(self, tmp_path: Path, body: str):
        return load_config(_write(tmp_path, body))

    def test_observable_missing_kind_is_actionable(self, tmp_path: Path) -> None:
        # the reproduced `{kine: ve}` typo
        with pytest.raises(ConfigError, match=r"observables\[0\].*'kind'"):
            self._load(
                tmp_path,
                """
            molecule: F2
            methods: [ti]
            observables: [{kine: ve}]
            output_dir: out
        """,
            )

    def test_td_missing_n_steps_is_actionable(self, tmp_path: Path) -> None:
        with pytest.raises(ConfigError, match=r"'td'.*'n_steps'"):
            self._load(
                tmp_path,
                """
            molecule: F2
            methods: [ti, td]
            observables: [{kind: ve, channels: 2}]
            td: {dt: 0.5}
            output_dir: out
        """,
            )

    def test_energies_missing_step_is_actionable(self, tmp_path: Path) -> None:
        with pytest.raises(ConfigError, match=r"'energies'.*'step'"):
            self._load(
                tmp_path,
                """
            molecule: F2
            methods: [ti]
            observables: [{kind: ve, channels: 2}]
            energies: {min: 0.01, max: 0.05}
            output_dir: out
        """,
            )

    def test_td_incident_missing_sigma_is_actionable(self, tmp_path: Path) -> None:
        with pytest.raises(ConfigError, match=r"'td.incident'.*'sigma'"):
            self._load(
                tmp_path,
                """
            molecule: F2
            methods: [ti, td]
            observables: [{kind: ve, channels: 2}]
            td: {dt: 0.5, n_steps: 10, incident: {r0: 20.0, p0: -0.4}}
            output_dir: out
        """,
            )


# --- energy ranges -----------------------------------------------------------


def _ranges_cfg(tmp_path: Path, ranges_yaml: str) -> object:
    return load_config(
        _write(
            tmp_path,
            f"""
        molecule: F2
        methods: [ti]
        observables: [{{kind: ve, channels: 2}}]
        energies:
{ranges_yaml}
        output_dir: out
    """,
        )
    )


class TestEnergyRanges:
    """`energies.ranges` -- a list of `np.arange(start, stop, step)` segments.

    A level-aware mesh is a union of uniform segments: one coarse background
    sweep plus a dense window around each resonance level. Written out point by
    point that is a 3343-line config built from 27 numbers, and a file that long
    is not read, only scrolled past. The segments ARE the specification.
    """

    def test_a_single_range_expands_like_arange(self, tmp_path: Path) -> None:
        cfg = _ranges_cfg(tmp_path, "          ranges: [{start: 0.01, stop: 0.05, step: 0.01}]")
        np.testing.assert_allclose(cfg.energies.as_array(), [0.01, 0.02, 0.03, 0.04])

    def test_the_upper_bound_is_exclusive_as_arange_documents(self, tmp_path: Path) -> None:
        """`stop` is not a sample. That is the one thing about `arange` that
        surprises people, so it is pinned here rather than left to be
        discovered: 0.05 is absent above, and a caller who wants it pads
        `stop` itself -- exactly as they would when calling numpy."""
        cfg = _ranges_cfg(tmp_path, "          ranges: [{start: 0.01, stop: 0.055, step: 0.01}]")
        np.testing.assert_allclose(cfg.energies.as_array(), [0.01, 0.02, 0.03, 0.04, 0.05])

    def test_overlapping_ranges_are_merged_sorted_and_deduplicated(self, tmp_path: Path) -> None:
        """A background sweep and a level window overlap by construction, so
        the mesh is their UNION, not their concatenation. A duplicated energy
        is solved twice and drawn as a vertical segment."""
        cfg = _ranges_cfg(
            tmp_path,
            "          ranges:\n"
            "            - {start: 0.01, stop: 0.045, step: 0.01}\n"
            "            - {start: 0.02, stop: 0.065, step: 0.01}",
        )
        got = cfg.energies.as_array()
        np.testing.assert_allclose(got, [0.01, 0.02, 0.03, 0.04, 0.05, 0.06])
        assert got.size == np.unique(got).size

    def test_ranges_survive_the_yaml_round_trip_config_resolved_performs(
        self, tmp_path: Path
    ) -> None:
        """`config.resolved.yaml` is written by dumping the dataclass, so the
        ranges must survive that dump and reload to the same mesh -- it is the
        offline reproduction path, and a recipe that does not round-trip is
        worse than no recipe.

        Scoped to the `energies` block on purpose. Reloading a WHOLE resolved
        config is separately broken (`_load_segment` rejects the `ecs: null`
        that `asdict` emits for a segment without an ECS tail), for the
        `values` form just as much as this one, so asserting on the whole file
        here would fail for a reason that has nothing to do with ranges.
        """
        import dataclasses

        import yaml

        cfg = _ranges_cfg(
            tmp_path, "          ranges: [{start: 0.002, stop: 0.0105, step: 0.0005}]"
        )
        dumped = yaml.safe_dump(dataclasses.asdict(cfg.energies), sort_keys=False)
        reloaded = _load_energies(yaml.safe_load(dumped))
        np.testing.assert_array_equal(reloaded.as_array(), cfg.energies.as_array())

    def test_a_non_positive_step_is_rejected_by_name(self, tmp_path: Path) -> None:
        """`np.arange` with step 0 raises ZeroDivisionError from inside numpy,
        and with a negative step returns an empty sweep in silence. Neither
        names the config line that is wrong."""
        with pytest.raises(ConfigError, match=r"step.*positive"):
            _ranges_cfg(tmp_path, "          ranges: [{start: 0.01, stop: 0.05, step: 0}]")

    def test_a_stop_below_start_is_rejected_by_name(self, tmp_path: Path) -> None:
        with pytest.raises(ConfigError, match=r"stop.*start"):
            _ranges_cfg(tmp_path, "          ranges: [{start: 0.05, stop: 0.01, step: 0.01}]")

    def test_mixing_ranges_with_values_is_rejected(self, tmp_path: Path) -> None:
        """Two meshes in one block have no defined answer; picking one would
        silently run a sweep the author did not write."""
        with pytest.raises(ConfigError, match="exactly one"):
            _ranges_cfg(
                tmp_path,
                "          values: [0.01]\n"
                "          ranges: [{start: 0.01, stop: 0.05, step: 0.01}]",
            )

    def test_an_empty_ranges_list_is_rejected(self, tmp_path: Path) -> None:
        with pytest.raises(ConfigError, match="at least one"):
            _ranges_cfg(tmp_path, "          ranges: []")
