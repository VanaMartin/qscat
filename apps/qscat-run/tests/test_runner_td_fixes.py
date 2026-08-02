"""Covering tests for the two TD-runner correctness fixes:

  Fix 1 -- per-kind test-function/surface resolution: `da`/`dr` must get a
  NUCLEAR outgoing packet + analysis surface, never the electronic `ve`
  default (F2/NO's `r0_out=24`, which is off the F2 nuclear grid's real
  region, R0~10.7).
  Fix 2 -- the elastic VE channel (`v_init in vprimes`) must get a SECOND
  `V_int=0` free-reference propagation feeding `ext.sigma(E, free=...)`,
  instead of silently reading `free=None` (the less-accurate literal-`S_ref=1`
  fallback).

All configs here are deliberately tiny (a handful of propagation steps) --
plumbing/wiring checks, not convergence studies.
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import numpy as np
import pytest
from qscat.core import TannorWeeks
from qscat.core import propagate as _original_propagate
from qscat_run import presets
from qscat_run import runner as runner_module
from qscat_run.config import (
    EnergySpec,
    ExperimentConfig,
    Observable,
    TdSpec,
    load_config,
    validate_config,
)
from qscat_run.runner import _nuclear_index_near, run_experiment

# --- Fix 1: per-kind test-function/surface resolution -----------------------


def _da_only_td_cfg() -> ExperimentConfig:
    """A minimal F2 TD config requesting only `da` (preset grid, no explicit
    grid) -- no propagation is run by these tests, only the resolution
    helpers, so `td.dt`/`n_steps` are arbitrary placeholders."""
    return ExperimentConfig(
        molecule="F2",
        methods=("td",),
        observables=(Observable(kind="da", channels=1),),
        output_dir="unused",
        energies=EnergySpec(values=(0.03,)),
        td=TdSpec(dt=0.2, n_steps=2, extractors=("flow",)),
    )


def test_da_only_resolves_nuclear_test_function_not_electronic_default() -> None:
    """F2's VALIDATED DA nuclear packet (r0_out=8, p0_out=72, sigma_out=0.07)
    must be what a `da` observable resolves to -- NOT the electronic `ve`
    default (r0_out=24) the pre-fix bug silently reused as the nuclear
    packet too."""
    resolved = presets.resolve_defaults(_da_only_td_cfg())
    tf = presets.resolve_test_function(resolved, "da")
    assert tf.r0_out == pytest.approx(8.0)
    assert tf.p0_out == pytest.approx(72.0)
    assert tf.sigma_out == pytest.approx(0.07)
    assert tf.r0_out != pytest.approx(24.0)


def test_da_only_resolves_surface_r_distinct_from_packet_r0_out() -> None:
    """The DA analysis SURFACE (R=6.0) is a genuinely different point from
    the nuclear packet's own center (R=8.0) -- see `MoleculePreset`'s
    docstring / `libs/qscat/tests/test_td_extractors.py`'s SP2 validation."""
    resolved = presets.resolve_defaults(_da_only_td_cfg())
    surface = presets.resolve_surface_r(resolved, "da")
    assert surface == pytest.approx(6.0)
    assert surface != presets.resolve_test_function(resolved, "da").r0_out


def test_da_only_nuclear_index_lands_inside_real_region_not_off_grid() -> None:
    """Wiring check on the ACTUAL F2 TD grid: the resolved surface (R=6.0)
    lands well inside the nuclear real region (R0~10.7 for F2's deck); the
    OLD buggy value (`r0_out=24.0`, the electronic default) falls past the
    real region entirely -- `_nuclear_index_near` masks anything past `R0`
    to `inf`, so `argmin` silently degenerates to the real-region EDGE point
    instead of a chosen physical analysis point. This test pins that the fix
    avoids that degenerate case."""
    resolved = presets.resolve_defaults(_da_only_td_cfg())
    tg = presets.resolve_grid(resolved, "td")
    nuc = tg.grids[1]
    assert nuc.R0 < 24.0  # confirms r0_out=24 would indeed be off the real region

    good_idx = _nuclear_index_near(tg, presets.resolve_surface_r(resolved, "da"))
    bad_idx = _nuclear_index_near(tg, 24.0)  # the pre-fix (buggy) value
    edge_idx = _nuclear_index_near(tg, nuc.R0)  # what "off-grid" degenerates to

    assert nuc.real_points[good_idx] == pytest.approx(6.0, abs=0.5)
    assert bad_idx == edge_idx  # the bug's failure mode: r0_out=24 hits the edge
    assert good_idx != bad_idx  # the fix picks a genuinely different point


# --- Fix 1 + 2 combined: a tiny mixed VE+DA end-to-end run -------------------


def _tiny_mixed_td_yaml(output_dir: str, *, ve_channels: str) -> str:
    """A TINY explicit-grid F2 TD config (seconds, not a convergence study),
    per-kind `test_function` (electronic `ve` / nuclear `da`, deliberately
    different small numbers -- the two-scale distinction the fix restores).
    `ve_channels` is a raw YAML fragment (`"2"` for `range(2)` -- includes
    the elastic v0->0 channel; `"[1]"` for just the inelastic v0->1 channel).
    """
    return textwrap.dedent(
        f"""
        molecule: F2
        methods: [td]
        observables:
          - {{kind: ve, channels: {ve_channels}}}
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
          n_steps: 6
          order: 3
          extractors: [tw]
          incident: {{r0: 4.0, p0: -0.5, sigma: 1.0}}
          test_function:
            ve: {{r0_out: 6.0, p0_out: 0.5, sigma_out: 1.0}}
            da: {{r0_out: 5.0, p0_out: 3.0, sigma_out: 0.5}}
        backend: auto
        output_dir: {output_dir}
        """
    )


def _write(tmp_path: Path, text: str) -> Path:
    p = tmp_path / "c.yaml"
    p.write_text(text)
    return p


def test_mixed_ve_da_run_uses_per_kind_packets_one_propagation_plus_free(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A mixed ve[2] (includes elastic v0->0) + da[1] run: ONE main
    propagation drives both observables' extractors, plus a SECOND
    free-reference propagation for the elastic VE channel (Fix 2) -- exactly
    2 `propagate()` calls total, never re-propagating per observable/kind."""
    calls: list[object] = []

    def _counting_propagate(*args: object, **kwargs: object) -> object:
        calls.append(kwargs.get("hamiltonian"))
        return _original_propagate(*args, **kwargs)

    monkeypatch.setattr(runner_module, "propagate", _counting_propagate)

    out_dir = tmp_path / "out"
    cfg = load_config(_write(tmp_path, _tiny_mixed_td_yaml(str(out_dir), ve_channels="2")))
    validate_config(cfg)
    result = run_experiment(cfg)

    assert len(calls) == 2  # main + free-reference
    assert "td:propagate" in result.timings
    assert "td:propagate_free" in result.timings

    for key in ("td:ve:tw:v0->0", "td:ve:tw:v0->1", "td:da:tw:ch0"):
        assert key in result.cross_sections
        assert np.all(np.isfinite(result.cross_sections[key])), key


def test_ve_without_elastic_channel_skips_free_propagation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The free-reference propagation is a no-op skip when no `ve`
    observable requests the elastic (`v_init`) channel -- no extra
    propagation cost when it isn't needed."""
    calls: list[object] = []

    def _counting_propagate(*args: object, **kwargs: object) -> object:
        calls.append(kwargs.get("hamiltonian"))
        return _original_propagate(*args, **kwargs)

    monkeypatch.setattr(runner_module, "propagate", _counting_propagate)

    out_dir = tmp_path / "out"
    cfg = load_config(_write(tmp_path, _tiny_mixed_td_yaml(str(out_dir), ve_channels="[1]")))
    validate_config(cfg)
    result = run_experiment(cfg)

    assert len(calls) == 1  # no free-reference propagation
    assert "td:propagate_free" not in result.timings
    assert np.all(np.isfinite(result.cross_sections["td:ve:tw:v0->1"]))


def test_elastic_ve_extractor_sigma_called_with_free_reference(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Direct evidence the elastic VE channel is NOT read via the literal-1
    fallback: `TannorWeeks.sigma` is spied on, and the ELECTRONIC-axis (`ve`)
    call must have a non-`None` `free`, while the nuclear-axis (`da`) call
    -- which has no elastic diagonal to subtract a reference from -- never
    does."""
    recorded: list[tuple[str, object]] = []
    original_sigma = TannorWeeks.sigma

    def _spy_sigma(
        self: TannorWeeks, E: object, *, free: object = None, n_steps: object = None
    ) -> object:
        recorded.append((self._axis, free))
        return original_sigma(self, E, free=free, n_steps=n_steps)

    monkeypatch.setattr(TannorWeeks, "sigma", _spy_sigma)

    out_dir = tmp_path / "out"
    cfg = load_config(_write(tmp_path, _tiny_mixed_td_yaml(str(out_dir), ve_channels="2")))
    validate_config(cfg)
    run_experiment(cfg)

    axes_with_free = {axis for axis, free in recorded if free is not None}
    axes_without_free = {axis for axis, free in recorded if free is None}
    assert axes_with_free == {"electronic"}
    assert "nuclear" in axes_without_free
