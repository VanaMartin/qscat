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

`td.test_function` carries EITHER a flat, back-compat single test-packet
block, or a per-observable-kind mapping (`{ve: {...}, da: {...}, dr: {...}}`)
-- `presets.resolve_test_function(cfg, kind)` is the ONE resolution entry
point, called separately per observable KIND: `ve`'s packet is electronic
(`wp_out`, in `r`); `da`/`dr`'s is nuclear (in `R`) -- physically different
scales, no longer conflated (the pre-fix bug: a single electronic-scale
packet silently doubling as the nuclear one too). The fixed electronic/
nuclear DVR analysis index `Dirac`/`Flux` read at
(`_electronic_index_near`/`_nuclear_index_near`) is likewise resolved
per-kind: `ve` reuses its own packet's `r0_out` (unchanged); `da`/`dr` use
`presets.resolve_surface_r(cfg, kind)`, which is generally a DIFFERENT point
from that kind's own nuclear `r0_out` (see `MoleculePreset`'s docstring --
e.g. F2's validated DA deck: packet centered at R=8, surface read at R=6).

An elastic VE channel (`cfg.v_init in vprimes`, requested by a `ve`
observable) gets a SECOND `V_int=0` free-reference propagation
(`qscat.core.time_dependent.free_hamiltonian`), mirroring `qscat.core.
td_ve_cross_sections_all`: a matching free extractor is built per (ve
observable, extractor name) and driven by that second propagation, then
`ext.sigma(E, free=free_ext)` subtracts the free-particle `S_free(E)` on the
diagonal channel instead of the less-accurate literal-`S_ref=1` fallback
(see `qscat.core.time_dependent._sigma_one_energy`'s docstring). This only
runs when an elastic VE channel is actually requested (a no-op skip
otherwise -- no extra propagation cost). DA/DR extractors take no free
reference (a pure rearrangement channel, no elastic diagonal to subtract
from -- `axis="nuclear"` extractors reject a non-`None` `free`).

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
from qscat.core.lcp import (
    lcp_da_cross_section,
    local_complex_potential,
    resonance_eigenstate_at_peak_width,
)
from qscat.core.time_dependent import free_hamiltonian  # same helper td_ve_cross_sections_all uses
from qscat.dvr import TensorGrid
from qscat.linalg import set_default_backend

from qscat_run import presets
from qscat_run.config import (
    VALID_EXTRACTORS,
    ConfigError,
    ExperimentConfig,
    Observable,
    TestFunctionSpec,
)

if TYPE_CHECKING:
    from qscat.model import ResonanceModel

__all__ = [
    "WavefunctionSnapshot",
    "EigenStates",
    "ResonanceState",
    "ExperimentResult",
    "run_experiment",
]

# The three extractor classes `td.extractors` can select, sharing one
# `Extractor`-protocol-conformant interface (`record`/`sigma`) plus the
# `n_steps=` truncated read added for this task's moment-resolved artifact.
type TdExtractor = TannorWeeks | Dirac | Flux


@dataclass(frozen=True)
class WavefunctionSnapshot:
    """One TI Psi+ (or TD Psi(t)) snapshot: per-axis density + optional field.

    `rho_r`/`rho_R` are `|Psi|^2` summed over the OTHER axis (masked to the
    unscaled/real region first, via `TensorGrid.real_mask()`); `r`/`R` are
    the matching real (unscaled) coordinate axes (`FemDvrEcsGrid.real_points`),
    same length as `rho_r`/`rho_R` respectively.

    `psi` is the FULL complex field reshaped to `(len(r), len(R))` and masked to
    the real region (zeros in the ECS tail), present only when the config asks
    for `full_field` -- the phase-carrying field `qscat.viz` domain-colours
    (the marginals throw the phase away). `None` otherwise.
    """

    kind: str
    label: str
    rho_r: npt.NDArray[np.float64]
    rho_R: npt.NDArray[np.float64]
    r: npt.NDArray[np.float64]
    R: npt.NDArray[np.float64]
    psi: npt.NDArray[np.complex128] | None = None


@dataclass(frozen=True)
class EigenStates:
    """A set of labelled complex fields on one axis, at labelled (real) energies.

    Two producers today: the target's vibrational levels + eigenfunctions
    (`kind="vibrational"`, from the `eps`/`chi` already diagonalized for the
    cross section; `label` tags method+grid, e.g. `"ti:vibrational"`), and the
    LCP nuclear SCATTERING states (`kind="lcp_scattering"`, the 1-D resolvent
    `psi_sc(R)` at each open collision energy). `energies` are the level/collision
    energies (Hartree), `states` shape `(n, len(axis))` the complex fields, `axis`
    the real (unscaled) coordinate the fields live on (`FemDvrEcsGrid.real_points`).
    """

    kind: str
    label: str
    energies: npt.NDArray[np.float64]
    states: npt.NDArray[np.complex128]
    axis: npt.NDArray[np.float64]


@dataclass(frozen=True)
class ResonanceState:
    """The resonance electronic eigenstate at the nuclear width peak (opt-in via
    `artifacts.eigenstates`, LCP-capable molecules).

    `energy` is the complex resonance pole `E_r - i*Gamma/2` (Hartree); `width`
    is `Gamma` (Hartree); `R` the nuclear geometry (bohr) it was solved at (the
    peak of `Gamma(R)`); `state` shape `(len(axis),)` the c-product-normalized
    electronic eigenfunction; `axis` the real electronic coordinate.
    """

    label: str
    R: float
    energy: complex
    width: float
    state: npt.NDArray[np.complex128]
    axis: npt.NDArray[np.float64]


@dataclass
class ExperimentResult:
    """Everything a config run produced, ready for `artifacts.write_artifacts`.

    `cross_section_vs_time` (TD only, opt-in): keyed
    `"{cross_sections key}@t{t_i:g}"`, one entry per (extractor, channel,
    requested moment) -- the moment-resolved sigma(E) series for the
    `cross_section_vs_time` artifact. `correlations` (TD only, opt-in):
    keyed `"{label}:t"`/`"{label}:c"` (`TannorWeeks`/`Dirac`) or
    `"{label}:t"`/`"{label}:b"`/`"{label}:d"` (`Flux`) -- the raw recorded
    per-step series behind each extractor's transform. `eigenstates` (opt-in):
    one entry per method run, the vibrational levels + eigenfunctions.
    """

    energies: npt.NDArray[np.float64]
    cross_sections: dict[str, npt.NDArray[np.float64]]
    wavefunctions: list[WavefunctionSnapshot]
    resolved_cfg: ExperimentConfig
    timings: dict[str, float] = field(default_factory=dict)
    grids: dict[str, TensorGrid] = field(default_factory=dict)
    cross_section_vs_time: dict[str, npt.NDArray[np.float64]] = field(default_factory=dict)
    correlations: dict[str, npt.NDArray[Any]] = field(default_factory=dict)
    eigenstates: list[EigenStates] = field(default_factory=list)
    resonance_states: list[ResonanceState] = field(default_factory=list)


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


def _masked_field(tg: TensorGrid, psi: npt.NDArray[np.complex128]) -> npt.NDArray[np.complex128]:
    """`psi` reshaped to `tg.shape` with the ECS tail zeroed (real region only)."""
    masked = psi.copy()
    masked[~tg.real_mask()] = 0.0
    return np.asarray(masked.reshape(tg.shape), dtype=np.complex128)


def _project_density(
    tg: TensorGrid, psi: npt.NDArray[np.complex128]
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """`|psi|^2`, masked to the real region, summed onto each axis."""
    dens = np.abs(_masked_field(tg, psi)) ** 2
    rho_r = np.asarray(dens.sum(axis=1), dtype=np.float64)
    rho_R = np.asarray(dens.sum(axis=0), dtype=np.float64)
    return rho_r, rho_R


def _vibrational_eigenstates(
    label: str,
    tg: TensorGrid,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
) -> EigenStates:
    """Package the already-computed nuclear `eps`/`chi` as an `EigenStates`."""
    return EigenStates(
        kind="vibrational",
        label=label,
        energies=np.asarray(eps, dtype=np.float64),
        states=np.asarray(chi, dtype=np.complex128),
        axis=tg.grids[1].real_points,
    )


def _run_ti(
    cfg: ExperimentConfig,
    timings: dict[str, float],
) -> tuple[
    dict[str, npt.NDArray[np.float64]], list[WavefunctionSnapshot], list[EigenStates], TensorGrid
]:
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

    eigenstates: list[EigenStates] = []
    if cfg.artifacts.eigenstates:
        eigenstates.append(_vibrational_eigenstates("ti:vibrational", tg, eps, chi))

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
                tg,
                model,
                eps,
                chi,
                cfg.v_init,
                [cfg.v_init],
                np.array([e]),
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
                    psi=_masked_field(tg, psi_plus) if wf_spec.full_field else None,
                )
            )
        timings["ti:wavefunction_snapshots"] = time.time() - t0

    return cross_sections, wavefunctions, eigenstates, tg


# --- TD (time-dependent) runner ---------------------------------------------


def _electronic_index_near(tg: TensorGrid, r_value: float) -> int:
    """The fixed electronic DVR index `Dirac`/`Flux` analyze at, nearest
    `r_value` (bohr) -- `_run_td` passes the resolved `ve` test function's
    own `r0_out` (the same point the VE `TannorWeeks` outgoing test packet
    is centered on; see module docstring)."""
    return tg.grids[0].real_index_near(r_value)


def _nuclear_index_near(tg: TensorGrid, r_value: float) -> int:
    """The nuclear-axis twin of `_electronic_index_near` (DA/DR's fixed
    analysis point) -- `_run_td` passes `presets.resolve_surface_r(cfg,
    kind)`, generally a DIFFERENT point from that kind's own nuclear test
    function's `r0_out` (see module docstring / `MoleculePreset`)."""
    return tg.grids[1].real_index_near(r_value)


def _wp_out_dict(tf: TestFunctionSpec) -> dict[str, float]:
    """A resolved per-kind `TestFunctionSpec` -> the `wp_out` dict the
    `qscat.core` wavepacket/extractor constructors take (`r0_out`/`p0_out`/
    `sigma_out`)."""
    return {"r0_out": tf.r0_out, "p0_out": tf.p0_out, "sigma_out": tf.sigma_out}


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
            tg,
            model,
            eps,
            chi,
            v_init,
            vprimes,
            wp_out,
            wp_in=wp_in,
            dt=dt,
            axis=axis,
            n_channels=n_channels,
        )
    if name == "delta":
        return Dirac(
            tg,
            model,
            eps,
            chi,
            v_init,
            vprimes,
            position,
            wp_in=wp_in,
            dt=dt,
            axis=axis,
            n_channels=n_channels,
        )
    if name == "flow":
        return Flux(
            tg,
            model,
            eps,
            chi,
            v_init,
            vprimes,
            surface,
            wp_in=wp_in,
            dt=dt,
            axis=axis,
            n_channels=n_channels,
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
    list[EigenStates],
    TensorGrid,
]:
    """The TD path: ONE Pade propagation drives every requested observable's
    extractor(s), plus a SECOND `V_int=0` free-reference propagation when an
    elastic VE channel is requested; returns `(cross_sections,
    cross_section_vs_time, correlations, wavefunctions, eigenstates, grid)` --
    see module docstring for the full design (per-kind test-function
    resolution, the free-reference propagation, the moment-truncation
    mechanism).
    """
    if cfg.td is None:
        raise ConfigError("no 'td' block resolved for this config")
    if cfg.energies is None:
        raise ConfigError("no energies resolved for this config (missing 'energies' block?)")
    td = cfg.td
    if td.incident is None:
        raise ConfigError("td.incident not resolved (missing preset defaults?)")
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

    eigenstates: list[EigenStates] = []
    if cfg.artifacts.eigenstates:
        eigenstates.append(_vibrational_eigenstates("td:vibrational", tg, eps, chi))

    wp_in = {"r0": td.incident.r0, "p0": td.incident.p0, "sigma": td.incident.sigma}
    psi0 = initial_state(tg, chi[cfg.v_init], **wp_in)

    # One (label, extractor, channel-label-list, kind) entry per (observable,
    # requested extractor name); `channel_labels[k]` is the suffix appended
    # to `label` for `ext.sigma(...)[:, k]`. `free_entries[label]` (VE
    # observables requesting the elastic channel only) is the companion
    # `V_int=0`-propagation extractor `ext.sigma`'s `free=` reads (see module
    # docstring) -- DA/DR entries never get one (no elastic diagonal).
    entries: list[tuple[str, TdExtractor, list[str], str]] = []
    free_entries: dict[str, TdExtractor] = {}
    for obs in cfg.observables:
        if obs.kind == "ve":
            vprimes = _vprimes(obs)
            channel_labels = [f"v{cfg.v_init}->{vp}" for vp in vprimes]
            tf_ve = presets.resolve_test_function(cfg, "ve")
            wp_out_ve = _wp_out_dict(tf_ve)
            elec_index = _electronic_index_near(tg, tf_ve.r0_out)
            wants_free = cfg.v_init in vprimes
            for name in td.extractors:
                ext = _build_extractor(
                    name,
                    tg,
                    model,
                    eps,
                    chi,
                    cfg.v_init,
                    vprimes,
                    wp_out_ve,
                    wp_in=wp_in,
                    dt=td.dt,
                    position=elec_index,
                    surface=elec_index,
                    axis="electronic",
                    n_channels=1,
                )
                label = f"td:ve:{name}"
                entries.append((label, ext, channel_labels, "ve"))
                if wants_free:
                    free_entries[label] = _build_extractor(
                        name,
                        tg,
                        model,
                        eps,
                        chi,
                        cfg.v_init,
                        vprimes,
                        wp_out_ve,
                        wp_in=wp_in,
                        dt=td.dt,
                        position=elec_index,
                        surface=elec_index,
                        axis="electronic",
                        n_channels=1,
                    )
        elif obs.kind in ("da", "dr"):
            n_channels = _n_channels(obs)
            channel_labels = [f"ch{c}" for c in range(n_channels)]
            tf_nuc = presets.resolve_test_function(cfg, obs.kind)
            wp_out_nuc = _wp_out_dict(tf_nuc)
            surface_r = presets.resolve_surface_r(cfg, obs.kind)
            nuc_index = _nuclear_index_near(tg, surface_r)
            for name in td.extractors:
                ext = _build_extractor(
                    name,
                    tg,
                    model,
                    eps,
                    chi,
                    cfg.v_init,
                    [],
                    wp_out_nuc,
                    wp_in=wp_in,
                    dt=td.dt,
                    position=nuc_index,
                    surface=nuc_index,
                    axis="nuclear",
                    n_channels=n_channels,
                )
                entries.append((f"td:{obs.kind}:{name}", ext, channel_labels, obs.kind))
        else:  # pragma: no cover -- validate_config already rejects unknown kinds
            raise ConfigError(f"unknown observable kind {obs.kind!r}")

    wf_spec = cfg.artifacts.wavefunction_snapshots
    snapshot_times = list(wf_spec.td_times) if wf_spec is not None and wf_spec.td_times else None
    # Keep the full Psi(t) field at the snapshot times only when full_field is
    # asked for (each kept field is a tg.size copy -- opt-in for memory).
    keep_psi_at = snapshot_times if (wf_spec is not None and wf_spec.full_field) else None

    t0 = time.time()
    result = propagate(
        tg,
        psi0,
        [],
        dt=td.dt,
        n_steps=td.n_steps,
        hamiltonian=model.hamiltonian(tg),
        order=td.order,
        extractors=[ext for _, ext, _, _ in entries],
        snapshot_times=snapshot_times,
        keep_psi_at=keep_psi_at,
    )
    timings["td:propagate"] = time.time() - t0

    if free_entries:
        # The elastic-VE free-reference propagation (Fix 2): SAME incident
        # packet/grid, `V_int=0` Hamiltonian -- no snapshots/other channels
        # needed, only the free extractors' recorded `c_{v'}(t)` series.
        t0 = time.time()
        propagate(
            tg,
            psi0,
            [],
            dt=td.dt,
            n_steps=td.n_steps,
            hamiltonian=free_hamiltonian(model, tg),
            order=td.order,
            extractors=list(free_entries.values()),
        )
        timings["td:propagate_free"] = time.time() - t0

    cross_sections: dict[str, npt.NDArray[np.float64]] = {}
    for label, ext, channel_labels, kind in entries:
        if kind == "ve":
            sigma = ext.sigma(energies, free=free_entries.get(label))
        else:
            sigma = ext.sigma(energies)
        for k, suffix in enumerate(channel_labels):
            cross_sections[f"{label}:{suffix}"] = sigma[:, k]

    cross_section_vs_time: dict[str, npt.NDArray[np.float64]] = {}
    cvt_spec = cfg.artifacts.cross_section_vs_time
    if cvt_spec is not None and cvt_spec.moments:
        t0 = time.time()
        for t_i in cvt_spec.moments:
            n_i = round(t_i / td.dt) + 1
            for label, ext, channel_labels, kind in entries:
                if kind == "ve":
                    sigma_t = ext.sigma(energies, free=free_entries.get(label), n_steps=n_i)
                else:
                    sigma_t = ext.sigma(energies, n_steps=n_i)
                for k, suffix in enumerate(channel_labels):
                    cross_section_vs_time[f"{label}:{suffix}@t{t_i:g}"] = sigma_t[:, k]
        timings["td:cross_section_vs_time"] = time.time() - t0

    correlations: dict[str, npt.NDArray[Any]] = {}
    if cfg.artifacts.correlations:
        for label, ext, _channel_labels, _kind in entries:
            if isinstance(ext, Flux):
                t_arr, b, d = ext.series
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
                    psi=_masked_field(tg, snap.psi) if snap.psi is not None else None,
                )
            )

    return cross_sections, cross_section_vs_time, correlations, wavefunctions, eigenstates, tg


def _run_lcp(
    cfg: ExperimentConfig,
    timings: dict[str, float],
) -> tuple[dict[str, npt.NDArray[np.float64]], list[EigenStates], list[ResonanceState]]:
    """The LCP path: reduce the exact 2-D DA resonance to a 1-D local complex
    potential `(V_d(R), Gamma(R))` (`qscat.core.lcp.local_complex_potential`,
    two ECS-angle electronic decks + the fine nuclear deck), then solve the 1-D
    driven equation (`lcp_da_cross_section`) per requested `da` observable.
    Keyed `"lcp:da:ch0"` so a `methods: [ti, lcp]` run overlays the exact and
    approximate DA cross sections on the SAME figure (disjoint key prefixes).

    Also emits (opt-in) the vibrational levels + the resonance electronic
    eigenstate at the width peak (`artifacts.eigenstates`) and the 1-D nuclear
    scattering states `psi_sc(R)` at the requested snapshot energies
    (`artifacts.wavefunction_snapshots.full_field` + `ti_energies`).
    """
    if cfg.energies is None:
        raise ConfigError("no energies resolved for this config (missing 'energies' block?)")
    energies = cfg.energies.as_array()

    t0 = time.time()
    g_R, elec_a, elec_b = presets.resolve_lcp_grids(cfg)
    timings["lcp:grid"] = time.time() - t0

    model = presets.MODELS[cfg.molecule]
    n_vib = _n_vib(cfg, cfg.v_init + 1)

    t0 = time.time()
    eps, chi = vibrational_states(g_R, model.mu, n_vib, model.v0)
    timings["lcp:vibrational_states"] = time.time() - t0

    t0 = time.time()
    v_d, gamma = local_complex_potential(model, g_R, elec_a, elec_b)
    timings["lcp:local_complex_potential"] = time.time() - t0

    cross_sections: dict[str, npt.NDArray[np.float64]] = {}
    for obs in cfg.observables:
        if obs.kind != "da":
            continue  # LCP only produces the DA cross section
        t0 = time.time()
        sigma = lcp_da_cross_section(g_R, model.mu, v_d, gamma, eps, chi, cfg.v_init, energies)
        cross_sections["lcp:da:ch0"] = np.asarray(sigma, dtype=np.float64)
        timings["lcp:da"] = timings.get("lcp:da", 0.0) + (time.time() - t0)

    eigenstates: list[EigenStates] = []
    resonance_states: list[ResonanceState] = []
    if cfg.artifacts.eigenstates:
        eigenstates.append(
            EigenStates(
                kind="vibrational",
                label="lcp:vibrational",
                energies=np.asarray(eps, dtype=np.float64),
                states=np.asarray(chi, dtype=np.complex128),
                axis=g_R.real_points,
            )
        )
        # The resonance electronic eigenstate at the width peak (#1).
        t0 = time.time()
        try:
            R_star, e_pole, phi_res = resonance_eigenstate_at_peak_width(
                model, g_R, elec_a, elec_b
            )
            resonance_states.append(
                ResonanceState(
                    label="lcp:resonance",
                    R=R_star,
                    energy=e_pole,
                    width=float(-2.0 * e_pole.imag),
                    state=phi_res,
                    axis=elec_a.real_points,
                )
            )
        except RuntimeError:
            pass  # no resolvable resonance width for this molecule/grid -- skip
        timings["lcp:resonance_state"] = time.time() - t0

    # The 1-D nuclear scattering states psi_sc(R) at the snapshot energies (#2).
    wf_spec = cfg.artifacts.wavefunction_snapshots
    if wf_spec is not None and wf_spec.full_field and wf_spec.ti_energies:
        t0 = time.time()
        snap_E = np.asarray(wf_spec.ti_energies, dtype=np.float64)
        _, psis = lcp_da_cross_section(
            g_R, model.mu, v_d, gamma, eps, chi, cfg.v_init, snap_E, return_wavefunction=True
        )
        psi_list = psis if isinstance(psis, list) else [psis]
        open_E = [float(e) for e, p in zip(snap_E, psi_list, strict=True) if p is not None]
        open_psi = [p for p in psi_list if p is not None]
        if open_psi:
            eigenstates.append(
                EigenStates(
                    kind="lcp_scattering",
                    label="lcp:scattering",
                    energies=np.asarray(open_E, dtype=np.float64),
                    states=np.asarray(open_psi, dtype=np.complex128),
                    axis=g_R.real_points,
                )
            )
        timings["lcp:scattering_states"] = time.time() - t0

    return cross_sections, eigenstates, resonance_states


def run_experiment(cfg: ExperimentConfig) -> ExperimentResult:
    """Resolve `cfg` fully (defaults, backend), run every requested method's
    observables on its own grid, and return the collected result.

    `methods: [ti, td, lcp]` runs each and merges their `cross_sections`/
    `wavefunctions`/`eigenstates` into one result -- their key prefixes
    (`"ti:"`/`"td:"`/`"lcp:"`) never collide, so the merge is a plain dict
    update and `cross_section.png` overlays exact vs approximation.
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
    eigenstates: list[EigenStates] = []
    resonance_states: list[ResonanceState] = []
    grids: dict[str, TensorGrid] = {}
    energies = np.zeros(0, dtype=np.float64)

    if "ti" in resolved.methods:
        ti_cross_sections, ti_wavefunctions, ti_eigenstates, tg_ti = _run_ti(resolved, timings)
        cross_sections.update(ti_cross_sections)
        wavefunctions.extend(ti_wavefunctions)
        eigenstates.extend(ti_eigenstates)
        grids["ti"] = tg_ti
        if resolved.energies is not None:
            energies = resolved.energies.as_array()

    if "td" in resolved.methods:
        (
            td_cross_sections,
            td_cross_section_vs_time,
            td_correlations,
            td_wavefunctions,
            td_eigenstates,
            tg_td,
        ) = _run_td(resolved, timings)
        cross_sections.update(td_cross_sections)
        cross_section_vs_time.update(td_cross_section_vs_time)
        correlations.update(td_correlations)
        wavefunctions.extend(td_wavefunctions)
        eigenstates.extend(td_eigenstates)
        grids["td"] = tg_td
        if resolved.energies is not None:
            energies = resolved.energies.as_array()

    if "lcp" in resolved.methods:
        lcp_cross_sections, lcp_eigenstates, lcp_resonance_states = _run_lcp(resolved, timings)
        cross_sections.update(lcp_cross_sections)
        eigenstates.extend(lcp_eigenstates)
        resonance_states.extend(lcp_resonance_states)
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
        eigenstates=eigenstates,
        resonance_states=resonance_states,
    )
