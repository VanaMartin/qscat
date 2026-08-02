"""The TI (time-independent) and TD (time-dependent) experiment runners:
resolve a config, run the requested method(s)' observables on their own
grid, and hand back an in-memory `ExperimentResult` for
`artifacts.write_artifacts`.

Observable dispatch mirrors the design spec exactly:
  - `ve`: `Observable.channels` is either an explicit tuple of final
    vibrational levels or an `int` count, expanded to `list(range(count))`;
    `qscat.core.ve_cross_section` (TI) / the electronic-axis TD extractors
    return one column per requested `v'`.
  - `da`/`dr`: `Observable.channels` is a plain `int` channel count, passed
    straight through as `n_channels`.

Every cross-section series is stored under a provenance-carrying key
`"{method}:{kind}:{channel_label}"` (TI: e.g. `"ti:ve:v0->1"`, `"ti:da:ch0"`;
TD: e.g. `"td:ve:tw:v0->1"`, `"td:da:flow:ch0"` -- the TD key also carries the
extractor name, since `td.extractors` may request several) --
`artifacts.py` writes these keys verbatim as CSV columns / npz names.
`methods: [ti, td]` runs both and merges their cross sections into ONE dict
(disjoint key prefixes, so nothing collides) -- `cross_section.png` then
overlays both automatically.

## The TD path (`_run_td`)

ONE Pade propagation (`qscat.core.propagate`) drives every requested
observable's extractor(s) at once: electronic-axis `TannorWeeks`/`Dirac`/
`Flux` for `ve` observables, nuclear-axis siblings (`axis="nuclear"`) for
`da`/`dr`. `td.extractors` (a subset of `{tw, delta, flow}`) selects which
extractor(s) to build per observable; `_build_extractor` is the one factory
both axes share.

`td.test_function` is the ONE test-packet config the schema carries; it
doubles as the electronic outgoing packet (VE's `wp_out`, in `r`) AND, for
`da`/`dr`, the nuclear outgoing packet (in `R`) -- a real approximation
across two physically different scales, but the only choice available
without a second schema field (see the design brief). The SAME point,
`test_function.r0_out`, also gives the fixed electronic/nuclear DVR index
`Dirac`/`Flux` analyze at (`_electronic_index_near`/`_nuclear_index_near`,
following `validation/h2plus/td_dr.py`'s helper of the same name) --
mirroring the pattern `libs/qscat/tests/test_td_extractors.py` uses
throughout ("reuses `POSITION`... colocated with `WP_OUT`").

LIMITATION (documented, not fixed here): unlike `qscat.core.
td_ve_cross_section`/`td_ve_cross_sections_all`, this mixed-observable path
does NOT run a second `V_int=0` free-reference propagation -- the design
brief calls for exactly ONE propagation feeding every extractor. An elastic
VE channel (`v_init in vprimes`) therefore reads `ext.sigma(E)` with
`free=None`, the literal-`S_ref=1` fallback (see `qscat.core.
time_dependent._sigma_one_energy`'s own docstring for why that fallback is
the less-accurate one) rather than the free-particle reference -- a real,
bounded inaccuracy on the elastic channel only; inelastic VE and DA/DR
channels have no diagonal to subtract a reference from in the first place
and are unaffected.

Moment-resolved `cross_section_vs_time` reads `ext.sigma(E, n_steps=n_i)`
against the SAME already-completed propagation (`n_i = round(t_i / dt) + 1`
samples, i.e. through time `t_i`) -- the `n_steps=` truncation kwarg added to
`TannorWeeks`/`Dirac`/`Flux.sigma` in `qscat.core.td_extractors` for this
task (`n_steps=None`, the default used everywhere else, is byte-identical to
the pre-truncation behavior).
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import numpy as np
import numpy.typing as npt
from qscat.core import (
    Dirac,
    Flux,
    TannorWeeks,
    da_cross_section,
    dr_cross_section,
    initial_state,
    propagate,
    ve_cross_section,
    vibrational_states,
)
from qscat.dvr import FemDvrEcsGrid, TensorGrid
from qscat.linalg import set_default_backend

from qscat_run import presets
from qscat_run.config import VALID_EXTRACTORS, ConfigError, ExperimentConfig, Observable

if TYPE_CHECKING:
    from qscat.model import ResonanceModel

__all__ = ["WavefunctionSnapshot", "ExperimentResult", "run_experiment"]

# The three extractor classes `td.extractors` can select, sharing one
# `Extractor`-protocol-conformant interface (`record`/`sigma`) plus the
# `n_steps=` truncated read added for this task's moment-resolved artifact.
TdExtractor = TannorWeeks | Dirac | Flux


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
    """Everything a config run produced, ready for `artifacts.write_artifacts`.

    `cross_section_vs_time` (TD only, opt-in): keyed
    `"{cross_sections key}@t{t_i:g}"`, one entry per (extractor, channel,
    requested moment) -- the moment-resolved sigma(E) series for the
    `cross_section_vs_time` artifact. `correlations` (TD only, opt-in):
    keyed `"{label}:t"`/`"{label}:c"` (`TannorWeeks`/`Dirac`) or
    `"{label}:t"`/`"{label}:b"`/`"{label}:d"` (`Flux`) -- the raw recorded
    per-step series behind each extractor's transform.
    """

    energies: npt.NDArray[np.float64]
    cross_sections: dict[str, npt.NDArray[np.float64]]
    wavefunctions: list[WavefunctionSnapshot]
    resolved_cfg: ExperimentConfig
    timings: dict[str, float] = field(default_factory=dict)
    grids: dict[str, TensorGrid] = field(default_factory=dict)
    cross_section_vs_time: dict[str, npt.NDArray[np.float64]] = field(default_factory=dict)
    correlations: dict[str, npt.NDArray[Any]] = field(default_factory=dict)


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
    """The vibrational basis size to diagonalize.

    Mirrors `presets.resolve_grid`'s explicit-grid precedence exactly: an
    explicit grid (`cfg.grid.electronic`/`nuclear` both given) bypasses the
    preset entirely, so `n_vib` is just `required` (`v_init` plus every
    requested `ve` final state) -- an explicit/custom coarse grid must not
    be forced to support more bound states than the request actually needs.
    Only when resolving to a NAMED PRESET grid is the count widened (never
    narrowed) to the preset's `n_vib`."""
    if cfg.grid.electronic is not None and cfg.grid.nuclear is not None:
        return required
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


# --- TD (time-dependent) runner ---------------------------------------------


def _index_near(grid: FemDvrEcsGrid, r_value: float) -> int:
    """Nearest REAL-region (unscaled) DVR index to `r_value` (bohr) on a
    single `FemDvrEcsGrid`, the shared body of `_electronic_index_near`/
    `_nuclear_index_near`, following `validation/h2plus/td_dr.py`'s
    `_nuclear_index_near` helper of the same idea (nearest real-region
    index, complex-tail points masked to `inf` first so they never win the
    `argmin`)."""
    real = grid.real_points
    masked = np.where(real <= grid.R0, real, np.inf)
    return int(np.argmin(np.abs(masked - r_value)))


def _electronic_index_near(tg: TensorGrid, r_value: float) -> int:
    """The fixed electronic DVR index `Dirac`/`Flux` analyze at, nearest
    `r_value` (bohr) -- `_run_td` passes `td.test_function.r0_out` (the same
    point the VE `TannorWeeks` outgoing test packet is centered on; see
    module docstring)."""
    return _index_near(tg.grids[0], r_value)


def _nuclear_index_near(tg: TensorGrid, r_value: float) -> int:
    """The nuclear-axis twin of `_electronic_index_near` (DA/DR's fixed
    analysis point)."""
    return _index_near(tg.grids[1], r_value)


def _build_extractor(
    name: str,
    tg: TensorGrid,
    model: ResonanceModel,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    vprimes: list[int],
    wp_out: dict[str, float],
    *,
    wp_in: dict[str, float],
    dt: float,
    position: int,
    surface: int,
    axis: str,
    n_channels: int,
) -> TdExtractor:
    """The one factory shared by every (observable kind, extractor name,
    axis) combination `_run_td` builds: `name` selects the class
    (`"tw"` -> `TannorWeeks`, `"delta"` -> `Dirac`, `"flow"` -> `Flux`);
    `axis`/`n_channels`/`vprimes` are passed straight through (electronic VE:
    `axis="electronic"`, `vprimes` the requested final states, `n_channels`
    unused; nuclear DA/DR: `axis="nuclear"`, `vprimes=[]`, `n_channels` the
    anion-channel count)."""
    if name == "tw":
        return TannorWeeks(
            tg, model, eps, chi, v_init, vprimes, wp_out,
            wp_in=wp_in, dt=dt, axis=axis, n_channels=n_channels,
        )
    if name == "delta":
        return Dirac(
            tg, model, eps, chi, v_init, vprimes, position,
            wp_in=wp_in, dt=dt, axis=axis, n_channels=n_channels,
        )
    if name == "flow":
        return Flux(
            tg, model, eps, chi, v_init, vprimes, surface,
            wp_in=wp_in, dt=dt, axis=axis, n_channels=n_channels,
        )
    raise ConfigError(  # pragma: no cover -- validate_config already rejects unknown names
        f"unknown td extractor {name!r}; choose from {sorted(VALID_EXTRACTORS)}"
    )


def _run_td(
    cfg: ExperimentConfig,
    timings: dict[str, float],
) -> tuple[
    dict[str, npt.NDArray[np.float64]],
    dict[str, npt.NDArray[np.float64]],
    dict[str, npt.NDArray[Any]],
    list[WavefunctionSnapshot],
    TensorGrid,
]:
    """The TD path: ONE Pade propagation drives every requested observable's
    extractor(s); returns `(cross_sections, cross_section_vs_time,
    correlations, wavefunctions, grid)` -- see module docstring for the full
    design (extractor construction, the free-reference limitation, the
    moment-truncation mechanism).
    """
    if cfg.td is None:
        raise ConfigError("no 'td' block resolved for this config")
    if cfg.energies is None:
        raise ConfigError("no energies resolved for this config (missing 'energies' block?)")
    td = cfg.td
    if td.incident is None:
        raise ConfigError("td.incident not resolved (missing preset defaults?)")
    if td.test_function is None:
        raise ConfigError("td.test_function not resolved (missing preset defaults?)")
    energies = cfg.energies.as_array()

    t0 = time.time()
    tg = presets.resolve_grid(cfg, "td")
    timings["td:grid"] = time.time() - t0

    required = cfg.v_init + 1
    for obs in cfg.observables:
        if obs.kind == "ve":
            required = max(required, max(_vprimes(obs), default=-1) + 1)
    model = presets.MODELS[cfg.molecule]
    n_vib = _n_vib(cfg, required)

    t0 = time.time()
    eps, chi = vibrational_states(tg.grids[1], model.mu, n_vib, model.v0)
    timings["td:vibrational_states"] = time.time() - t0

    wp_in = {"r0": td.incident.r0, "p0": td.incident.p0, "sigma": td.incident.sigma}
    wp_out = {
        "r0_out": td.test_function.r0_out,
        "p0_out": td.test_function.p0_out,
        "sigma_out": td.test_function.sigma_out,
    }
    elec_index = _electronic_index_near(tg, td.test_function.r0_out)
    nuc_index = _nuclear_index_near(tg, td.test_function.r0_out)

    psi0 = initial_state(tg, chi[cfg.v_init], **wp_in)

    # One (label, extractor, channel-label-list) entry per (observable,
    # requested extractor name); `channel_labels[k]` is the suffix appended
    # to `label` for `ext.sigma(...)[:, k]`.
    entries: list[tuple[str, TdExtractor, list[str]]] = []
    for obs in cfg.observables:
        if obs.kind == "ve":
            vprimes = _vprimes(obs)
            channel_labels = [f"v{cfg.v_init}->{vp}" for vp in vprimes]
            for name in td.extractors:
                ext = _build_extractor(
                    name, tg, model, eps, chi, cfg.v_init, vprimes, wp_out,
                    wp_in=wp_in, dt=td.dt, position=elec_index, surface=elec_index,
                    axis="electronic", n_channels=1,
                )
                entries.append((f"td:ve:{name}", ext, channel_labels))
        elif obs.kind in ("da", "dr"):
            n_channels = _n_channels(obs)
            channel_labels = [f"ch{c}" for c in range(n_channels)]
            for name in td.extractors:
                ext = _build_extractor(
                    name, tg, model, eps, chi, cfg.v_init, [], wp_out,
                    wp_in=wp_in, dt=td.dt, position=nuc_index, surface=nuc_index,
                    axis="nuclear", n_channels=n_channels,
                )
                entries.append((f"td:{obs.kind}:{name}", ext, channel_labels))
        else:  # pragma: no cover -- validate_config already rejects unknown kinds
            raise ConfigError(f"unknown observable kind {obs.kind!r}")

    wf_spec = cfg.artifacts.wavefunction_snapshots
    snapshot_times = list(wf_spec.td_times) if wf_spec is not None and wf_spec.td_times else None

    t0 = time.time()
    result = propagate(
        tg,
        psi0,
        [],
        dt=td.dt,
        n_steps=td.n_steps,
        hamiltonian=model.hamiltonian(tg),
        order=td.order,
        extractors=[ext for _, ext, _ in entries],
        snapshot_times=snapshot_times,
    )
    timings["td:propagate"] = time.time() - t0

    cross_sections: dict[str, npt.NDArray[np.float64]] = {}
    for label, ext, channel_labels in entries:
        sigma = ext.sigma(energies)
        for k, suffix in enumerate(channel_labels):
            cross_sections[f"{label}:{suffix}"] = sigma[:, k]

    cross_section_vs_time: dict[str, npt.NDArray[np.float64]] = {}
    cvt_spec = cfg.artifacts.cross_section_vs_time
    if cvt_spec is not None and cvt_spec.moments:
        t0 = time.time()
        for t_i in cvt_spec.moments:
            n_i = round(t_i / td.dt) + 1
            for label, ext, channel_labels in entries:
                sigma_t = ext.sigma(energies, n_steps=n_i)
                for k, suffix in enumerate(channel_labels):
                    cross_section_vs_time[f"{label}:{suffix}@t{t_i:g}"] = sigma_t[:, k]
        timings["td:cross_section_vs_time"] = time.time() - t0

    correlations: dict[str, npt.NDArray[Any]] = {}
    if cfg.artifacts.correlations:
        for label, ext, _channel_labels in entries:
            if isinstance(ext, Flux):
                t_arr, b, d = ext._arrays
                correlations[f"{label}:t"] = t_arr
                correlations[f"{label}:b"] = b
                correlations[f"{label}:d"] = d
            else:  # TannorWeeks | Dirac
                res = ext.result
                correlations[f"{label}:t"] = res.t
                correlations[f"{label}:c"] = res.c

    wavefunctions: list[WavefunctionSnapshot] = []
    if snapshot_times is not None:
        for snap in result.snapshots:
            wavefunctions.append(
                WavefunctionSnapshot(
                    kind="td",
                    label=f"t{snap.time:g}",
                    rho_r=snap.rho_r,
                    rho_R=snap.rho_R,
                    r=tg.grids[0].real_points,
                    R=tg.grids[1].real_points,
                )
            )

    return cross_sections, cross_section_vs_time, correlations, wavefunctions, tg


def run_experiment(cfg: ExperimentConfig) -> ExperimentResult:
    """Resolve `cfg` fully (defaults, backend), run every requested method's
    observables on its own grid, and return the collected result.

    `methods: [ti, td]` runs both (`_run_ti` then `_run_td`) and merges their
    `cross_sections`/`wavefunctions` into one result -- their key prefixes
    (`"ti:"`/`"td:"`) never collide, so the merge is a plain dict update.
    """
    t_start = time.time()
    resolved = presets.resolve_defaults(cfg)

    if resolved.backend != "auto":
        set_default_backend(resolved.backend)

    timings: dict[str, float] = {}
    cross_sections: dict[str, npt.NDArray[np.float64]] = {}
    cross_section_vs_time: dict[str, npt.NDArray[np.float64]] = {}
    correlations: dict[str, npt.NDArray[Any]] = {}
    wavefunctions: list[WavefunctionSnapshot] = []
    grids: dict[str, TensorGrid] = {}
    energies = np.zeros(0, dtype=np.float64)

    if "ti" in resolved.methods:
        ti_cross_sections, ti_wavefunctions, tg_ti = _run_ti(resolved, timings)
        cross_sections.update(ti_cross_sections)
        wavefunctions.extend(ti_wavefunctions)
        grids["ti"] = tg_ti
        if resolved.energies is not None:
            energies = resolved.energies.as_array()

    if "td" in resolved.methods:
        (
            td_cross_sections,
            td_cross_section_vs_time,
            td_correlations,
            td_wavefunctions,
            tg_td,
        ) = _run_td(resolved, timings)
        cross_sections.update(td_cross_sections)
        cross_section_vs_time.update(td_cross_section_vs_time)
        correlations.update(td_correlations)
        wavefunctions.extend(td_wavefunctions)
        grids["td"] = tg_td
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
        cross_section_vs_time=cross_section_vs_time,
        correlations=correlations,
    )
