"""End-to-end TI-runner coverage: a TINY explicit-grid F2 config (fast, no
preset-sized deck), run through `run_experiment` + `write_artifacts`, and a
`CliRunner` `--dry-run` smoke test."""

from __future__ import annotations

import textwrap
from pathlib import Path

import numpy as np
from click.testing import CliRunner
from qscat_run.artifacts import write_artifacts
from qscat_run.cli import main
from qscat_run.config import (
    EcsSpec,
    ExperimentConfig,
    GridSpec,
    Observable,
    SegmentSpec,
    load_config,
    validate_config,
)
from qscat_run.runner import _n_vib, run_experiment


def _tiny_f2_ti_yaml(output_dir: str) -> str:
    # A deliberately tiny explicit grid: small electronic r_max (~8) and
    # nuclear r_max (~10), low quadrature order, so the whole run costs
    # seconds -- this is a plumbing smoke test, not a convergence study.
    return textwrap.dedent(
        f"""
        molecule: F2
        methods: [ti]
        observables:
          - {{kind: ve, channels: 2}}
          - {{kind: da, channels: 1}}
        energies: {{values: [0.03, 0.05]}}
        v_init: 0
        grid:
          electronic:
            real: [[4, 2.0], [2, 4.0], [2, 8.0]]
            ecs: {{angle: 30, elements: 3, quadrature: 5}}
          nuclear:
            real: [[3, 1.8], [1, 2.0], [2, 2.5], [2, 2.6], [2, 2.7], [3, 10.0]]
            ecs: {{angle: 30, elements: 3, quadrature: 6}}
        artifacts:
          cross_section: true
          wavefunction_snapshots:
            ti_energies: [0.05]
        backend: auto
        output_dir: {output_dir}
        """
    )


def _write(tmp_path: Path, text: str) -> Path:
    p = tmp_path / "c.yaml"
    p.write_text(text)
    return p


def test_ti_run_writes_cross_section_and_manifest(tmp_path: Path) -> None:
    out_dir = tmp_path / "out"
    cfg = load_config(_write(tmp_path, _tiny_f2_ti_yaml(str(out_dir))))
    validate_config(cfg)

    result = run_experiment(cfg)
    assert result.cross_sections  # ve + da series present
    write_artifacts(result, cfg, out_dir, timestamp="2026-01-01T00:00:00")

    assert (out_dir / "cross_section.csv").exists()
    assert (out_dir / "cross_section.png").exists()
    assert (out_dir / "manifest.json").exists()
    assert (out_dir / "config.resolved.yaml").exists()

    arr = np.load(out_dir / "cross_section.npz")
    assert np.all(np.isfinite(arr["energy"]))
    for key in result.cross_sections:
        assert key in arr
        assert np.all(np.isfinite(arr[key]))

    # ve[2] -> v0->0, v0->1; da[1] -> ch0.
    assert "ti:ve:v0->0" in result.cross_sections
    assert "ti:ve:v0->1" in result.cross_sections
    assert "ti:da:ch0" in result.cross_sections

    # the requested TI wavefunction snapshot at E=0.05.
    assert len(result.wavefunctions) == 1
    wf = result.wavefunctions[0]
    assert np.all(np.isfinite(wf.rho_r))
    assert np.all(np.isfinite(wf.rho_R))
    assert (out_dir / "wavefunction" / f"psi_{wf.label}.npz").exists()
    assert (out_dir / "wavefunction" / f"psi_{wf.label}.png").exists()


def test_run_cli_dry_run_prints_plan_without_solving(tmp_path: Path) -> None:
    out_dir = tmp_path / "out"
    cfg_path = _write(tmp_path, _tiny_f2_ti_yaml(str(out_dir)))

    r = CliRunner().invoke(main, ["run", str(cfg_path), "--dry-run"])
    assert r.exit_code == 0, r.output
    assert "molecule: F2" in r.output
    assert "grid[ti]" in r.output
    assert not out_dir.exists()


def _tiny_explicit_segment() -> SegmentSpec:
    return SegmentSpec(
        real=((4, 2.0), (2, 4.0), (2, 8.0)),
        ecs=EcsSpec(angle=30.0, elements=3, quadrature=5),
    )


def test_n_vib_explicit_grid_uses_required_not_preset_floor() -> None:
    """An EXPLICIT grid must not be floored at the molecule's preset `n_vib`
    (N2:emoscat's is 6): a config asking only for the v_init channel
    (`required=1`) should get `n_vib=1`, not 6 -- forcing 6 bound states on a
    tiny/coarse custom grid that may not support them is exactly the
    spurious `vibrational_states` ValueError this fixes."""
    cfg = ExperimentConfig(
        molecule="N2",
        methods=("ti",),
        observables=(Observable(kind="ve", channels=1),),
        output_dir="unused",
        v_init=0,
        grid=GridSpec(electronic=_tiny_explicit_segment(), nuclear=_tiny_explicit_segment()),
    )
    assert _n_vib(cfg, required=1) == 1


def test_n_vib_preset_grid_still_floors_at_preset_n_vib() -> None:
    """No explicit grid (preset only): the old floor-at-preset behavior is
    unchanged -- N2:emoscat's `n_vib=6` still wins over a smaller `required`."""
    cfg = ExperimentConfig(
        molecule="N2",
        methods=("ti",),
        observables=(Observable(kind="ve", channels=1),),
        output_dir="unused",
        v_init=0,
        grid=GridSpec(),
    )
    assert _n_vib(cfg, required=1) == 6
