"""End-to-end NRM-method coverage: F2 `methods: [nrm]` through
`run_experiment` + `write_artifacts`.

The NRM solve runs on the F2 preset grids (the electronic deck at two ECS
angles + the fine nuclear deck), so it is `@slow` -- the fast gate is the
config-level NRM validation in `test_config.py`. It is nonetheless a real
NUMERICAL check, not a smoke test: it pins sigma_DA against the value
`validation/diatomic/test_nrm.py` measured through the library directly, so a
future rewiring of the app path that changed the grids, the vibrational basis,
the state-sum truncation or `v_init` would move it.
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
