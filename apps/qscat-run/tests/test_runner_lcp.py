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
