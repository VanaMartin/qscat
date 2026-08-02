"""The TI (time-independent) experiment runner: resolve a config, run ONE
driven-Psi+ sweep per requested observable on the shared TI grid, and hand
back an in-memory `ExperimentResult` for `artifacts.write_artifacts`.

TD (time-dependent) is not implemented here yet (Task 3) -- `run_experiment`
raises `NotImplementedError` if `"td"` is requested.

Observable dispatch mirrors the design spec exactly:
  - `ve`: `Observable.channels` is either an explicit tuple of final
    vibrational levels or an `int` count, expanded to `list(range(count))`;
    `qscat.core.ve_cross_section` returns one column per requested `v'`.
  - `da`/`dr`: `Observable.channels` is a plain `int` channel count, passed
    straight through as `n_channels`.

Every cross-section series is stored under a provenance-carrying key
`"{method}:{kind}:{channel_label}"` (e.g. `"ti:ve:v0->1"`, `"ti:da:ch0"`) --
`artifacts.py` writes these keys verbatim as CSV columns / npz names.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt
from qscat.core import da_cross_section, dr_cross_section, ve_cross_section, vibrational_states
from qscat.dvr import TensorGrid
from qscat.linalg import set_default_backend

from qscat_run import presets
from qscat_run.config import ConfigError, ExperimentConfig, Observable

__all__ = ["WavefunctionSnapshot", "ExperimentResult", "run_experiment"]


@dataclass(frozen=True)
class WavefunctionSnapshot:
    """One TI Psi+ density snapshot, projected onto each axis.

    `rho_r`/`rho_R` are `|Psi+|^2` summed over the OTHER axis (masked to the
    unscaled/real region first, via `TensorGrid.real_mask()`); `r`/`R` are
    the matching real (unscaled) coordinate axes (`FemDvrEcsGrid.real_points`),
    same length as `rho_r`/`rho_R` respectively.
    """

    kind: str
    label: str
    rho_r: npt.NDArray[np.float64]
    rho_R: npt.NDArray[np.float64]
    r: npt.NDArray[np.float64]
    R: npt.NDArray[np.float64]


@dataclass
class ExperimentResult:
    """Everything a config run produced, ready for `artifacts.write_artifacts`."""

    energies: npt.NDArray[np.float64]
    cross_sections: dict[str, npt.NDArray[np.float64]]
    wavefunctions: list[WavefunctionSnapshot]
    resolved_cfg: ExperimentConfig
    timings: dict[str, float] = field(default_factory=dict)
    grids: dict[str, TensorGrid] = field(default_factory=dict)


def _vprimes(obs: Observable) -> list[int]:
    """`ve`'s final-state list: an explicit tuple used verbatim, or an `int`
    count expanded to `range(count)`."""
    if obs.channels is None:
        raise ConfigError(f"observable {obs.kind!r} has no 'channels' (not resolved?)")
    if isinstance(obs.channels, tuple):
        return list(obs.channels)
    return list(range(obs.channels))


def _n_channels(obs: Observable) -> int:
    """`da`/`dr`'s channel count: a plain `int`, or the length of an
    explicit tuple (accepted for symmetry with `ve`)."""
    if obs.channels is None:
        raise ConfigError(f"observable {obs.kind!r} has no 'channels' (not resolved?)")
    if isinstance(obs.channels, tuple):
        return len(obs.channels)
    return obs.channels


def _n_vib(cfg: ExperimentConfig, required: int) -> int:
    """The vibrational basis size to diagonalize: the molecule's preset
    `n_vib`, widened (never narrowed) to cover every level `run_experiment`
    actually needs (`v_init` plus every requested `ve` final state)."""
    variant = cfg.grid.preset or presets.DEFAULT_PRESET
    preset = presets.PRESETS.get(f"{cfg.molecule}:{variant}")
    base = preset.n_vib if preset is not None else required
    return max(base, required)


def _project_density(
    tg: TensorGrid, psi: npt.NDArray[np.complex128]
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """`|psi|^2`, masked to the real region, summed onto each axis."""
    masked = psi.copy()
    masked[~tg.real_mask()] = 0.0
    dens = np.abs(masked.reshape(tg.shape)) ** 2
    rho_r = np.asarray(dens.sum(axis=1), dtype=np.float64)
    rho_R = np.asarray(dens.sum(axis=0), dtype=np.float64)
    return rho_r, rho_R


def _run_ti(
    cfg: ExperimentConfig,
    timings: dict[str, float],
) -> tuple[dict[str, npt.NDArray[np.float64]], list[WavefunctionSnapshot], TensorGrid]:
    if cfg.energies is None:
        raise ConfigError("no energies resolved for this config (missing 'energies' block?)")
    energies = cfg.energies.as_array()

    t0 = time.time()
    tg = presets.resolve_grid(cfg, "ti")
    timings["ti:grid"] = time.time() - t0

    required = cfg.v_init + 1
    for obs in cfg.observables:
        if obs.kind == "ve":
            required = max(required, max(_vprimes(obs), default=-1) + 1)
    model = presets.MODELS[cfg.molecule]
    n_vib = _n_vib(cfg, required)

    t0 = time.time()
    eps, chi = vibrational_states(tg.grids[1], model.mu, n_vib, model.v0)
    timings["ti:vibrational_states"] = time.time() - t0

    cross_sections: dict[str, npt.NDArray[np.float64]] = {}
    for obs in cfg.observables:
        t0 = time.time()
        if obs.kind == "ve":
            vprimes = _vprimes(obs)
            sigma = ve_cross_section(tg, model, eps, chi, cfg.v_init, vprimes, energies)
            for j, vp in enumerate(vprimes):
                cross_sections[f"ti:ve:v{cfg.v_init}->{vp}"] = sigma[:, j]
        elif obs.kind == "da":
            n_channels = _n_channels(obs)
            sigma_da = da_cross_section(
                tg, model, eps, chi, cfg.v_init, energies, n_channels=n_channels
            )
            for c in range(n_channels):
                cross_sections[f"ti:da:ch{c}"] = sigma_da[:, c]
        elif obs.kind == "dr":
            n_channels = _n_channels(obs)
            sigma_dr = dr_cross_section(
                tg, model, eps, chi, cfg.v_init, energies, n_channels=n_channels
            )
            for c in range(n_channels):
                cross_sections[f"ti:dr:ch{c}"] = sigma_dr[:, c]
        else:  # pragma: no cover -- validate_config already rejects unknown kinds
            raise ConfigError(f"unknown observable kind {obs.kind!r}")
        timings[f"ti:{obs.kind}"] = timings.get(f"ti:{obs.kind}", 0.0) + (time.time() - t0)

    wavefunctions: list[WavefunctionSnapshot] = []
    wf_spec = cfg.artifacts.wavefunction_snapshots
    if wf_spec is not None and wf_spec.ti_energies:
        t0 = time.time()
        for e in wf_spec.ti_energies:
            _, psis = ve_cross_section(
                tg, model, eps, chi, cfg.v_init, [cfg.v_init], np.array([e]),
                return_wavefunction=True,
            )
            psi_plus = psis[0] if isinstance(psis, list) else psis
            if psi_plus is None:  # E <= 0: no driven solve, nothing to snapshot
                continue
            rho_r, rho_R = _project_density(tg, psi_plus)
            wavefunctions.append(
                WavefunctionSnapshot(
                    kind="ti",
                    label=f"E{e:g}",
                    rho_r=rho_r,
                    rho_R=rho_R,
                    r=tg.grids[0].real_points,
                    R=tg.grids[1].real_points,
                )
            )
        timings["ti:wavefunction_snapshots"] = time.time() - t0

    return cross_sections, wavefunctions, tg


def run_experiment(cfg: ExperimentConfig) -> ExperimentResult:
    """Resolve `cfg` fully (defaults, backend), run every requested method's
    observables on its own grid, and return the collected result.

    Only `"ti"` is implemented here; `"td" in cfg.methods` raises
    `NotImplementedError` (Task 3's job).
    """
    t_start = time.time()
    resolved = presets.resolve_defaults(cfg)

    if resolved.backend != "auto":
        set_default_backend(resolved.backend)

    if "td" in resolved.methods:
        raise NotImplementedError(
            "run_experiment: the TD (time-dependent) path is not implemented yet "
            "(Task 3); remove 'td' from 'methods' or run only 'ti' for now."
        )

    timings: dict[str, float] = {}
    cross_sections: dict[str, npt.NDArray[np.float64]] = {}
    wavefunctions: list[WavefunctionSnapshot] = []
    grids: dict[str, TensorGrid] = {}
    energies = np.zeros(0, dtype=np.float64)

    if "ti" in resolved.methods:
        cross_sections, wavefunctions, tg = _run_ti(resolved, timings)
        grids["ti"] = tg
        if resolved.energies is not None:
            energies = resolved.energies.as_array()

    timings["total"] = time.time() - t_start

    return ExperimentResult(
        energies=energies,
        cross_sections=cross_sections,
        wavefunctions=wavefunctions,
        timings=timings,
        resolved_cfg=resolved,
        grids=grids,
    )
