"""End-to-end LCP-method coverage (R1): F2 `methods: [lcp]` (and `[ti, lcp]`)
through `run_experiment` + `write_artifacts`. The LCP solve runs on the F2
preset grids (two ECS-angle electronic decks + the fine nuclear deck), so it
is `@slow` -- the fast gate is the config-level LCP validation in
`test_config.py`."""

from __future__ import annotations

import textwrap
from pathlib import Path

import numpy as np
import pytest
from qscat_run.artifacts import write_artifacts
from qscat_run.config import load_config, validate_config
from qscat_run.runner import run_experiment


def _write(tmp_path: Path, text: str) -> Path:
    p = tmp_path / "c.yaml"
    p.write_text(textwrap.dedent(text))
    return p


@pytest.mark.slow
def test_lcp_run_produces_da_cross_section(tmp_path: Path) -> None:
    out_dir = tmp_path / "out"
    cfg = load_config(
        _write(
            tmp_path,
            f"""
        molecule: F2
        methods: [lcp]
        observables: [{{kind: da, channels: 1}}]
        energies: {{values: [0.03, 0.04]}}
        grid: {{preset: emoscat}}
        backend: auto
        output_dir: {out_dir}
    """,
        )
    )
    validate_config(cfg)
    result = run_experiment(cfg)

    assert "lcp:da:ch0" in result.cross_sections
    sigma = result.cross_sections["lcp:da:ch0"]
    assert sigma.shape == (2,)
    assert np.all(np.isfinite(sigma)) and np.all(sigma >= 0.0)

    write_artifacts(result, cfg, out_dir, timestamp="2026-01-01T00:00:00")
    assert (out_dir / "cross_section.csv").exists()
    arr = np.load(out_dir / "cross_section.npz")
    assert "lcp:da:ch0" in arr


@pytest.mark.slow
def test_lcp_emits_resonance_state_and_scattering_wavefunctions(tmp_path: Path) -> None:
    # #1 + #2: with eigenstates + full_field wavefunction snapshots, the LCP run
    # emits the resonance electronic eigenstate (complex pole + width) AND the
    # 1-D nuclear scattering states psi_sc(R) at the snapshot energies.
    out_dir = tmp_path / "out"
    cfg = load_config(
        _write(
            tmp_path,
            f"""
        molecule: F2
        methods: [lcp]
        observables: [{{kind: da, channels: 1}}]
        energies: {{values: [0.02, 0.03, 0.04]}}
        grid: {{preset: emoscat}}
        artifacts:
          cross_section: true
          eigenstates: true
          wavefunction_snapshots: {{ti_energies: [0.03, 0.04], full_field: true}}
        backend: auto
        output_dir: {out_dir}
    """,
        )
    )
    validate_config(cfg)
    result = run_experiment(cfg)

    # #1 resonance state
    assert len(result.resonance_states) == 1
    rs = result.resonance_states[0]
    assert rs.label == "lcp:resonance"
    assert rs.width > 0.0 and abs(-2.0 * rs.energy.imag - rs.width) < 1e-12
    assert rs.state.shape == (rs.axis.size,) and rs.state.dtype == np.complex128
    assert 1.0 < rs.R < 3.0

    # #2 LCP scattering states (F2 exothermic -> both snapshot energies open)
    scat = [es for es in result.eigenstates if es.kind == "lcp_scattering"]
    assert len(scat) == 1
    assert scat[0].states.shape == (2, scat[0].axis.size)
    assert np.all(np.isfinite(scat[0].states))

    write_artifacts(result, cfg, out_dir, timestamp="2026-01-01T00:00:00")
    assert (out_dir / "resonance" / "resonance_lcp_resonance.npz").exists()
    assert (out_dir / "resonance" / "resonance_lcp_resonance.png").exists()
    assert (out_dir / "eigenstates" / "eigenstates_lcp_scattering.npz").exists()


@pytest.mark.slow
def test_ti_and_lcp_overlay_disjoint_keys(tmp_path: Path) -> None:
    # methods: [ti, lcp] must produce BOTH the exact (ti:da) and approximate
    # (lcp:da) DA cross sections under disjoint keys -> one overlaid figure.
    out_dir = tmp_path / "out"
    cfg = load_config(
        _write(
            tmp_path,
            f"""
        molecule: F2
        methods: [ti, lcp]
        observables: [{{kind: da, channels: 1}}]
        energies: {{values: [0.03, 0.04]}}
        grid: {{preset: emoscat}}
        backend: auto
        output_dir: {out_dir}
    """,
        )
    )
    validate_config(cfg)
    result = run_experiment(cfg)
    assert "ti:da:ch0" in result.cross_sections
    assert "lcp:da:ch0" in result.cross_sections


@pytest.mark.slow
def test_lcp_run_produces_ve_cross_section_for_n2(tmp_path: Path) -> None:
    out_dir = tmp_path / "out"
    cfg = load_config(
        _write(
            tmp_path,
            f"""
        molecule: N2
        methods: [lcp]
        observables: [{{kind: ve, channels: [0, 1]}}]
        energies: {{values: [0.05, 0.1]}}
        output_dir: {out_dir}
    """,
        )
    )
    validate_config(cfg)
    result = run_experiment(cfg)
    for key in ("lcp:ve:v0->0", "lcp:ve:v0->1"):
        assert key in result.cross_sections
        sigma = result.cross_sections[key]
        assert sigma.shape == (2,)
        assert np.all(np.isfinite(sigma)) and np.all(sigma >= 0.0)
    # the Pi_g resonance region beats near-threshold on 0->1
    assert result.cross_sections["lcp:ve:v0->1"][1] > result.cross_sections["lcp:ve:v0->1"][0]
