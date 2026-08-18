"""End-to-end NRM-method coverage through `run_experiment` + `write_artifacts`:
F2 dissociative attachment and N2 vibrational excitation.

The NRM solve runs on the preset grids (the molecule's electronic deck at two
ECS angles + its nuclear deck), so the solving tests are `@slow` -- the fast
gate is the config-level NRM validation in `test_config.py` plus the
grid-sharing checks here. They are nonetheless real NUMERICAL checks, not smoke
tests: each pins sigma against the value `validation/diatomic/{nrm,ve_nrm}.py`
measured through the library directly, so a future rewiring of the app path
that changed the grids, the vibrational basis, the state-sum truncation,
`v_init` or the background setting would move it.
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import numpy as np
import pytest
from qscat_run.artifacts import write_artifacts
from qscat_run.config import load_config, validate_config
from qscat_run.presets import resolve_nrm_grids
from qscat_run.runner import run_experiment


def _write(tmp_path: Path, text: str) -> Path:
    p = tmp_path / "c.yaml"
    p.write_text(textwrap.dedent(text))
    return p


def test_nrm_grids_share_the_ti_electronic_deck(tmp_path: Path) -> None:
    """The NRM must run on the SAME electronic discretisation as the exact `ti`
    solve, or a `methods: [ti, nrm]` ratio compares two different problems. Fast
    -- grid construction only, no solve."""
    cfg = load_config(
        _write(
            tmp_path,
            """
        molecule: F2
        methods: [ti, nrm]
        observables: [{kind: da, channels: 1}]
        energies: {values: [0.03]}
        grid: {preset: emoscat}
        output_dir: out
    """,
        )
    )
    validate_config(cfg)
    from qscat_run.presets import resolve_grid

    tg = resolve_grid(cfg, "ti")
    nuc, elec_a, elec_b = resolve_nrm_grids(cfg)

    np.testing.assert_allclose(elec_a.points, tg.grids[0].points, rtol=0, atol=0)
    np.testing.assert_allclose(nuc.points, tg.grids[1].points, rtol=0, atol=0)
    # The second angle must genuinely differ, or choice A's two-angle pole walk
    # has nothing to match against.
    assert not np.allclose(elec_b.points, elec_a.points)


@pytest.mark.slow
def test_nrm_choice_b_reproduces_the_recorded_f2_cross_section(tmp_path: Path) -> None:
    """F2, choice B, E = 0.030 Ha on the eMoScat production deck.

    1.65514 bohr^2 is what `validation/diatomic/nrm.py` measured driving
    `qscat.core.nrm` directly (against an exact-2D oracle of 1.65611 -- a 0.059 %
    deviation). The app path must reproduce it, so the tolerance is tight: 1e-4
    relative, ~100x the run-to-run round-off and ~6x under the physical
    deviation from the oracle it is quoted against.
    """
    out_dir = tmp_path / "out"
    cfg = load_config(
        _write(
            tmp_path,
            f"""
        molecule: F2
        methods: [nrm]
        observables: [{{kind: da, channels: 1}}]
        energies: {{values: [0.030]}}
        grid: {{preset: emoscat}}
        nrm: {{choices: [b], n_states: 100}}
        backend: auto
        output_dir: {out_dir}
    """,
        )
    )
    validate_config(cfg)
    result = run_experiment(cfg)

    assert "nrm-b:da:ch0" in result.cross_sections
    sigma = result.cross_sections["nrm-b:da:ch0"]
    assert sigma.shape == (1,)
    np.testing.assert_allclose(sigma[0], 1.65514, rtol=1e-4)

    write_artifacts(result, cfg, out_dir, timestamp="2026-01-01T00:00:00")
    arr = np.load(out_dir / "cross_section.npz")
    assert "nrm-b:da:ch0" in arr


@pytest.mark.slow
def test_both_choices_produce_disjoint_series(tmp_path: Path) -> None:
    """`choices: [a, b]` must give TWO series under distinct keys -- the whole
    point of the comparison figure -- and they must not be numerically equal
    (PRA 77's finding is that the discrete-state choice matters; on F2 at
    E = 0.030 the recorded values are 1.39751 for A against 1.65514 for B)."""
    out_dir = tmp_path / "out"
    cfg = load_config(
        _write(
            tmp_path,
            f"""
        molecule: F2
        methods: [nrm]
        observables: [{{kind: da, channels: 1}}]
        energies: {{values: [0.030]}}
        grid: {{preset: emoscat}}
        nrm: {{choices: [a, b], n_states: 100}}
        backend: auto
        output_dir: {out_dir}
    """,
        )
    )
    validate_config(cfg)
    result = run_experiment(cfg)

    assert set(result.cross_sections) == {"nrm-a:da:ch0", "nrm-b:da:ch0"}
    sigma_a = result.cross_sections["nrm-a:da:ch0"][0]
    sigma_b = result.cross_sections["nrm-b:da:ch0"][0]
    assert np.isfinite(sigma_a) and sigma_a > 0.0
    np.testing.assert_allclose(sigma_a, 1.39751, rtol=1e-4)
    np.testing.assert_allclose(sigma_b, 1.65514, rtol=1e-4)


def test_nrm_ve_config_validates_and_keys_match_the_ti_shape(tmp_path: Path) -> None:
    """A `ve` + `nrm` config is accepted and its knobs survive loading. Fast --
    no solve. The `include_background` flag is the paper's own "nonlocal" vs
    "nonlocal + background" distinction, so it has to reach `NrmSpec` intact."""
    cfg = load_config(
        _write(
            tmp_path,
            """
        molecule: N2
        methods: [ti, nrm]
        observables: [{kind: ve, channels: 2}]
        energies: {values: [0.10]}
        grid: {preset: emoscat}
        nrm: {choices: [a, b], n_states: 100, include_background: false}
        output_dir: out
    """,
        )
    )
    validate_config(cfg)
    assert cfg.nrm is not None
    assert cfg.nrm.choices == ("a", "b")
    assert cfg.nrm.n_states == 100
    assert cfg.nrm.include_background is False

    # N2 has no `lcp_grids`; the NRM resolves its decks off `ti_grid()` instead,
    # so `ti` and `nrm` in this run share one discretisation.
    from qscat_run.presets import resolve_grid

    tg = resolve_grid(cfg, "ti")
    nuc, elec_a, elec_b = resolve_nrm_grids(cfg)
    np.testing.assert_allclose(elec_a.points, tg.grids[0].points, rtol=0, atol=0)
    np.testing.assert_allclose(nuc.points, tg.grids[1].points, rtol=0, atol=0)
    assert not np.allclose(elec_b.points, elec_a.points)


@pytest.mark.slow
def test_nrm_ve_reproduces_the_recorded_n2_cross_sections(tmp_path: Path) -> None:
    """N2 vibrational excitation, both discrete-state choices, E = 0.100 Ha on
    the eMoScat production deck.

    The four values pinned here are what `validation/diatomic/ve_nrm.py`
    measured driving `qscat.core.nrm` directly -- reproduced through the app
    path digit for digit, which is the point: `qscat-run` must not be a second
    implementation. They are the same numbers behind
    `docs/physics/figures/n2-ve-nrm-vs-exact.png` at that energy, where the
    exact 2-D oracle gives 23.67447212 and 6.12299519 bohr^2. Choice B lands
    within 0.11 % / 0.08 % of it; choice A is 0.01 % / 10.0 % off, the
    Born-Oppenheimer breakdown of PRA 77 Sec. VI A showing up in the inelastic
    channel.

    rtol = 1e-4 -- tight enough that a changed grid, vibrational basis,
    state-sum truncation, `v_init` or background setting moves it, and far
    under the physical spread it separates.
    """
    out_dir = tmp_path / "out"
    cfg = load_config(
        _write(
            tmp_path,
            f"""
        molecule: N2
        methods: [nrm]
        observables: [{{kind: ve, channels: 2}}]
        energies: {{values: [0.100]}}
        grid: {{preset: emoscat}}
        nrm: {{choices: [a, b], n_states: 100, include_background: true}}
        backend: auto
        output_dir: {out_dir}
    """,
        )
    )
    validate_config(cfg)
    result = run_experiment(cfg)

    expected = {
        "nrm-a:ve:v0->0": 23.67134387,
        "nrm-a:ve:v0->1": 5.51229923,
        "nrm-b:ve:v0->0": 23.64751979,
        "nrm-b:ve:v0->1": 6.11807041,
    }
    assert set(result.cross_sections) == set(expected)
    for key, value in expected.items():
        series = result.cross_sections[key]
        assert series.shape == (1,)
        np.testing.assert_allclose(series[0], value, rtol=1e-4)

    write_artifacts(result, cfg, out_dir, timestamp="2026-01-01T00:00:00")
    arr = np.load(out_dir / "cross_section.npz")
    for key in expected:
        assert key in arr


@pytest.mark.slow
def test_nrm_ve_background_flag_is_load_bearing(tmp_path: Path) -> None:
    """`include_background: false` must give PRA 77's bare "nonlocal" curve,
    not the same numbers.

    Pinned rather than merely "different": at N2 E = 0.100 Ha, choice B, the
    elastic cross section drops 23.64751979 -> 22.59611344 bohr^2 when the
    Eq. (37) background is removed (4.4 %), and the exact oracle is
    23.67447212 -- so dropping it moves the answer AWAY from exact. A flag that
    silently did nothing would pass an inequality check against a wrongly
    recorded value; it cannot pass this one.
    """
    out_dir = tmp_path / "out"
    cfg = load_config(
        _write(
            tmp_path,
            f"""
        molecule: N2
        methods: [nrm]
        observables: [{{kind: ve, channels: 2}}]
        energies: {{values: [0.100]}}
        grid: {{preset: emoscat}}
        nrm: {{choices: [b], n_states: 100, include_background: false}}
        backend: auto
        output_dir: {out_dir}
    """,
        )
    )
    validate_config(cfg)
    result = run_experiment(cfg)

    np.testing.assert_allclose(result.cross_sections["nrm-b:ve:v0->0"][0], 22.59611344, rtol=1e-4)
    np.testing.assert_allclose(result.cross_sections["nrm-b:ve:v0->1"][0], 6.1455567, rtol=1e-4)
