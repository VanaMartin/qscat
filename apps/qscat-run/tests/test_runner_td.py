"""End-to-end TD-runner coverage: a TINY explicit-grid F2 config (fast, few
propagation steps), mixing a VE and a DA observable in ONE propagation, run
through `run_experiment` + `write_artifacts`.

Grid/incident/test-function parameters mirror
`libs/qscat/tests/test_td_extractors.py`'s deliberately tiny, fast N2
fixture (seconds, not a converged cross section) -- this is a plumbing smoke
test, not a convergence study. The explicit grid is the SAME one
`test_runner_ti.py` already proved supports a bound anion electronic state
(its `da` observable solves successfully there), so `da` is safe to include
here too.

`td.test_function` uses the per-observable-kind mapping (`{ve: {...}, da:
{...}}`, an explicit custom grid with no matching preset -- see
`presets.resolve_test_function`): the `ve` packet is electronic (in `r`,
launched near the SAME grid's electronic real region), the `da` packet is
nuclear (in `R`) -- two deliberately DIFFERENT small numbers, covering the
`presets.py`/`runner.py` fix that a single shared packet used to conflate.
`v_init=0` is one of `ve`'s requested channels, so this run also exercises
the elastic VE free-reference propagation (Fix 2, `test_runner_td_fixes.py`
asserts on it directly; this file just checks the resulting artifacts are
all present/finite)."""

from __future__ import annotations

import textwrap
from pathlib import Path

import numpy as np
from qscat_run.artifacts import write_artifacts
from qscat_run.config import load_config, validate_config
from qscat_run.runner import run_experiment


def _tiny_f2_td_yaml(output_dir: str, *, correlations: bool = False) -> str:
    return textwrap.dedent(
        f"""
        molecule: F2
        methods: [td]
        observables:
          - {{kind: ve, channels: 2}}
          - {{kind: da, channels: 1}}
        v_init: 0
        energies: {{values: [0.03, 0.05]}}
        grid:
          electronic:
            real: [[4, 2.0], [2, 4.0], [2, 8.0]]
            ecs: {{angle: 30, elements: 3, quadrature: 5}}
          nuclear:
            real: [[3, 1.8], [1, 2.0], [2, 2.5], [2, 2.6], [2, 2.7], [3, 10.0]]
            ecs: {{angle: 30, elements: 3, quadrature: 6}}
        td:
          dt: 0.2
          n_steps: 10
          order: 3
          extractors: [tw, delta]
          incident: {{r0: 4.0, p0: -0.5, sigma: 1.0}}
          test_function:
            ve: {{r0_out: 6.0, p0_out: 0.5, sigma_out: 1.0}}
            da: {{r0_out: 5.0, p0_out: 3.0, sigma_out: 0.5}}
        artifacts:
          cross_section: true
          cross_section_vs_time:
            moments: [1.0, 2.0]
          wavefunction_snapshots:
            td_times: [0.0, 1.0]
          correlations: {correlations}
        backend: auto
        output_dir: {output_dir}
        """
    )


def _write(tmp_path: Path, text: str) -> Path:
    p = tmp_path / "c.yaml"
    p.write_text(text)
    return p


def test_td_run_writes_all_td_artifacts(tmp_path: Path) -> None:
    out_dir = tmp_path / "out"
    cfg = load_config(_write(tmp_path, _tiny_f2_td_yaml(str(out_dir))))
    validate_config(cfg)

    result = run_experiment(cfg)
    write_artifacts(result, cfg, out_dir, timestamp="2026-01-01T00:00:00")

    # Cross sections: ve[2] x {tw, delta} -> v0->0, v0->1; da[1] x {tw, delta} -> ch0.
    assert result.cross_sections
    for extractor in ("tw", "delta"):
        assert f"td:ve:{extractor}:v0->0" in result.cross_sections
        assert f"td:ve:{extractor}:v0->1" in result.cross_sections
        assert f"td:da:{extractor}:ch0" in result.cross_sections
    for key, arr in result.cross_sections.items():
        assert np.all(np.isfinite(arr)), key

    assert (out_dir / "cross_section.npz").exists()
    assert (out_dir / "cross_section.png").exists()

    # Moment-resolved cross_section_vs_time: one entry per (cross-section
    # key, moment).
    assert (out_dir / "cross_section_vs_time.npz").exists()
    assert (out_dir / "cross_section_vs_time.png").exists()
    assert result.cross_section_vs_time
    npz = np.load(out_dir / "cross_section_vs_time.npz")
    for t_i in (1.0, 2.0):
        key = f"td:ve:tw:v0->0@t{t_i:g}"
        assert key in result.cross_section_vs_time
        assert key in npz
        assert np.all(np.isfinite(npz[key]))
    # A later moment differs from an earlier one (a genuinely truncated read,
    # not an accidental full-series repeat). With the elastic free-reference
    # fix (Fix 2) the v0->0 channel's absolute magnitude is now genuinely
    # tiny (the spurious ~500x literal-1-fallback background is gone), so
    # `atol=0` forces a RELATIVE comparison -- otherwise `allclose`'s default
    # absolute tolerance would swallow the (still large, ~7x) relative
    # difference between the two moments.
    early = result.cross_section_vs_time["td:ve:tw:v0->0@t1"]
    late = result.cross_section_vs_time["td:ve:tw:v0->0@t2"]
    assert not np.allclose(early, late, rtol=1e-3, atol=0.0)
    # The last moment (t=2.0 == n_steps*dt) matches the untruncated sigma.
    np.testing.assert_allclose(late, result.cross_sections["td:ve:tw:v0->0"])

    # TD wavefunction density snapshots at the requested td_times.
    assert len(result.wavefunctions) == 2
    wf_dir = out_dir / "wavefunction"
    npz_files = sorted(wf_dir.glob("psi_t*.npz"))
    assert npz_files
    for wf in result.wavefunctions:
        assert wf.kind == "td"
        assert np.all(np.isfinite(wf.rho_r))
        assert np.all(np.isfinite(wf.rho_R))
        assert (wf_dir / f"psi_{wf.label}.npz").exists()
        assert (wf_dir / f"psi_{wf.label}.png").exists()

    # correlations opt-in (default False here): no correlations.npz.
    assert not (out_dir / "correlations.npz").exists()
    assert not result.correlations


def test_td_run_with_correlations_writes_correlations_npz(tmp_path: Path) -> None:
    out_dir = tmp_path / "out"
    cfg = load_config(_write(tmp_path, _tiny_f2_td_yaml(str(out_dir), correlations=True)))
    validate_config(cfg)

    result = run_experiment(cfg)
    write_artifacts(result, cfg, out_dir, timestamp="2026-01-01T00:00:00")

    assert result.correlations
    corr_path = out_dir / "correlations.npz"
    assert corr_path.exists()
    npz = np.load(corr_path)
    # tw -> "{label}:t"/"{label}:c"; delta (Dirac) -> the same "t"/"c" shape.
    assert "td:ve:tw:t" in npz
    assert "td:ve:tw:c" in npz
    assert "td:ve:delta:t" in npz
    assert "td:ve:delta:c" in npz
    assert "td:da:tw:t" in npz
    assert "td:da:delta:t" in npz
    for key in npz.files:
        assert np.all(np.isfinite(npz[key])), key
