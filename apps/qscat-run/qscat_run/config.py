"""The `ExperimentConfig` schema: YAML -> frozen dataclasses, + validation.

`load_config` is a tolerant YAML -> dataclass parser: every block beyond the
four required top-level keys (`molecule`, `methods`, `observables`,
`output_dir`) is optional and, if omitted, is left as `None` (or a field
default) rather than raising -- `presets.resolve_defaults` fills the gaps
from the molecule's preset later. `validate_config` is the actionable-error
gate: it checks the molecule is known, the (molecule, observable.kind) combo
is allowed per the validity matrix, a `td` block is present whenever `"td"`
is requested, `td.extractors` names only known extractors, and an explicit
grid supplies both halves.

`Observable.channels` (and `energies`, `td.incident`/`test_function`) are
typed as optional (`None` when the YAML omits them) precisely so
`presets.resolve_defaults` has something to fill -- a deliberate widening of
the "always present" reading of the schema's YAML example.

This module must NOT import `qscat_run.presets` at module scope: `presets`
imports this module's dataclasses (`EnergySpec`, `IncidentSpec`, ...), so a
top-level `config -> presets` edge would close a cycle. `validate_config`
instead imports `presets` lazily, inside the function body, once both
modules have finished loading.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import click
import numpy as np
import numpy.typing as npt
import yaml

__all__ = [
    "ConfigError",
    "Observable",
    "EnergySpec",
    "EcsSpec",
    "SegmentSpec",
    "GridSpec",
    "IncidentSpec",
    "TestFunctionSpec",
    "TdSpec",
    "CrossSectionVsTimeSpec",
    "WavefunctionSnapshotsSpec",
    "ArtifactSpec",
    "ExperimentConfig",
    "load_config",
    "validate_config",
    "VALID_METHODS",
    "VALID_EXTRACTORS",
]


class ConfigError(click.ClickException):
    """An actionable, user-facing config problem (unknown molecule, an
    invalid (molecule, observable) combo, a missing `td` block, ...).

    A `click.ClickException` subclass so `cli.py` commands can simply let it
    propagate: click's own error handling prints `Error: <message>` and sets
    a non-zero exit code without any command-level try/except.
    """


VALID_METHODS = frozenset({"ti", "td"})
VALID_EXTRACTORS = frozenset({"flow", "delta", "tw"})


# --- observables -------------------------------------------------------------


@dataclass(frozen=True)
class Observable:
    """One requested observable: `kind` (`"ve"` | `"da"` | `"dr"`) + the
    final-state channel count or explicit list. `channels is None` means
    "omitted -- fill from the preset" (`presets.resolve_defaults`)."""

    kind: str
    channels: int | tuple[int, ...] | None = None


def _parse_channels(raw: Any) -> int | tuple[int, ...] | None:
    if raw is None:
        return None
    if isinstance(raw, list):
        return tuple(int(c) for c in raw)
    return int(raw)


def _load_observables(raw: list[Any]) -> tuple[Observable, ...]:
    return tuple(
        Observable(kind=str(item["kind"]), channels=_parse_channels(item.get("channels")))
        for item in raw
    )


# --- energies ----------------------------------------------------------------


@dataclass(frozen=True)
class EnergySpec:
    """Either a `min/max/step` sweep or an explicit `values` list."""

    min: float | None = None
    max: float | None = None
    step: float | None = None
    values: tuple[float, ...] | None = None

    def as_array(self) -> npt.NDArray[np.float64]:
        """The concrete energy sweep: `values` verbatim if given, else
        `np.arange(min, max + step/2, step)` (the half-step pad makes the
        sweep inclusive of `max` despite `arange`'s exclusive upper bound
        and float round-off), rounded to 10 decimals to kill the last-ULP
        noise `arange` accumulates over many steps."""
        if self.values is not None:
            return np.asarray(self.values, dtype=np.float64)
        if self.min is None or self.max is None or self.step is None:
            raise ConfigError(
                "EnergySpec.as_array() needs either 'values' or all of 'min'/'max'/'step'"
            )
        arr = np.arange(self.min, self.max + 0.5 * self.step, self.step, dtype=np.float64)
        return np.round(arr, decimals=10)


def _load_energies(raw: dict[str, Any] | None) -> EnergySpec | None:
    if raw is None:
        return None
    if "values" in raw:
        return EnergySpec(values=tuple(float(v) for v in raw["values"]))
    return EnergySpec(min=float(raw["min"]), max=float(raw["max"]), step=float(raw["step"]))


# --- grid ----------------------------------------------------------------


@dataclass(frozen=True)
class EcsSpec:
    """The ECS tail of an explicit grid segment: `angle` (degrees),
    `elements` (the tail element count, i.e. `fem_grid_exp_tail`'s
    `tail_n`), `quadrature`."""

    angle: float
    elements: int
    quadrature: int


@dataclass(frozen=True)
class SegmentSpec:
    """One explicit-grid axis: `real` is `(n_elements, endpoint)` segment
    pairs (eMoScat `grids.txt` format), `ecs` the exponentially-growing
    complex tail beyond the last real endpoint."""

    real: tuple[tuple[int, float], ...]
    ecs: EcsSpec


def _load_segment(raw: dict[str, Any]) -> SegmentSpec:
    real = tuple((int(n), float(end)) for n, end in raw["real"])
    ecs_raw = raw["ecs"]
    ecs = EcsSpec(
        angle=float(ecs_raw["angle"]),
        elements=int(ecs_raw["elements"]),
        quadrature=int(ecs_raw["quadrature"]),
    )
    return SegmentSpec(real=real, ecs=ecs)


@dataclass(frozen=True)
class GridSpec:
    """`preset` selects a named per-molecule deck (`None` -> the default
    preset); an explicit `electronic`/`nuclear` pair overrides the preset
    for BOTH methods (the user then owns the TD launch-box requirement)."""

    preset: str | None = None
    electronic: SegmentSpec | None = None
    nuclear: SegmentSpec | None = None


def _load_grid(raw: dict[str, Any] | None) -> GridSpec:
    raw = raw or {}
    electronic = _load_segment(raw["electronic"]) if "electronic" in raw else None
    nuclear = _load_segment(raw["nuclear"]) if "nuclear" in raw else None
    preset = raw.get("preset")
    return GridSpec(
        preset=str(preset) if preset is not None else None,
        electronic=electronic,
        nuclear=nuclear,
    )


# --- td ----------------------------------------------------------------


@dataclass(frozen=True)
class IncidentSpec:
    """The TD incident electronic Gaussian wavepacket: position/momentum/width."""

    r0: float
    p0: float
    sigma: float


@dataclass(frozen=True)
class TestFunctionSpec:
    """One TD outgoing test packet (electronic `wp_out`, in `r`, for `ve`; or
    nuclear, in `R`, for `da`/`dr` -- the SAME shape, different physical
    scale/axis depending which observable kind it is resolved for; see
    `TdSpec.test_function`/`test_functions`): position/impulse/width."""

    r0_out: float
    p0_out: float
    sigma_out: float


@dataclass(frozen=True)
class TdSpec:
    """The `td` block, required iff `"td" in methods`.

    `test_function`/`test_functions` are the two shapes `td.test_function`'s
    YAML may take (`_load_td` picks exactly one, based on whether the block
    has an `r0_out` key): a flat `TestFunctionSpec` (back-compat -- the
    historical single packet, unambiguous for a single-observable-kind run,
    see `presets.resolve_test_function`), or a per-observable-kind mapping
    (`{"ve": ..., "da": ..., "dr": ...}`) -- required to disambiguate a
    MIXED `ve`+`da`/`dr` run, since `ve`'s outgoing packet is electronic (in
    `r`) and `da`/`dr`'s is nuclear (in `R`), physically different scales.
    `presets.resolve_test_function`/`resolve_surface_r` are the resolution
    entry points; `presets.resolve_defaults` also fills `test_functions` from
    the preset for convenience, but callers should always resolve through
    those functions rather than reading either field directly.
    """

    dt: float
    n_steps: int
    order: int = 3
    extractors: tuple[str, ...] = ()
    incident: IncidentSpec | None = None
    test_function: TestFunctionSpec | None = None
    test_functions: dict[str, TestFunctionSpec] | None = None


def _load_td(raw: dict[str, Any] | None) -> TdSpec | None:
    if raw is None:
        return None
    incident = None
    if "incident" in raw:
        ir = raw["incident"]
        incident = IncidentSpec(r0=float(ir["r0"]), p0=float(ir["p0"]), sigma=float(ir["sigma"]))
    test_function = None
    test_functions = None
    if "test_function" in raw:
        tr = raw["test_function"]
        if "r0_out" in tr:
            test_function = TestFunctionSpec(
                r0_out=float(tr["r0_out"]),
                p0_out=float(tr["p0_out"]),
                sigma_out=float(tr["sigma_out"]),
            )
        else:
            test_functions = {
                str(kind): TestFunctionSpec(
                    r0_out=float(block["r0_out"]),
                    p0_out=float(block["p0_out"]),
                    sigma_out=float(block["sigma_out"]),
                )
                for kind, block in tr.items()
            }
    return TdSpec(
        dt=float(raw["dt"]),
        n_steps=int(raw["n_steps"]),
        order=int(raw.get("order", 3)),
        extractors=tuple(str(e) for e in raw.get("extractors", ())),
        incident=incident,
        test_function=test_function,
        test_functions=test_functions,
    )


# --- artifacts ----------------------------------------------------------------


@dataclass(frozen=True)
class CrossSectionVsTimeSpec:
    moments: tuple[float, ...] = ()


@dataclass(frozen=True)
class WavefunctionSnapshotsSpec:
    td_times: tuple[float, ...] = ()
    ti_energies: tuple[float, ...] = ()
    # When true, snapshots also carry the FULL complex Psi field (masked to the
    # real region) -- emitted as `psi` in the npz + a domain-coloured png, for
    # `qscat.viz`. Default false keeps only the cheap per-axis density marginals.
    full_field: bool = False


@dataclass(frozen=True)
class ArtifactSpec:
    cross_section: bool = True
    cross_section_vs_time: CrossSectionVsTimeSpec | None = None
    correlations: bool = False
    wavefunction_snapshots: WavefunctionSnapshotsSpec | None = None


def _load_artifacts(raw: dict[str, Any] | None) -> ArtifactSpec:
    raw = raw or {}
    cvt_raw = raw.get("cross_section_vs_time")
    cvt = (
        CrossSectionVsTimeSpec(moments=tuple(float(m) for m in cvt_raw.get("moments", ())))
        if cvt_raw
        else None
    )
    wf_raw = raw.get("wavefunction_snapshots")
    wf = (
        WavefunctionSnapshotsSpec(
            td_times=tuple(float(t) for t in wf_raw.get("td_times", ())),
            ti_energies=tuple(float(e) for e in wf_raw.get("ti_energies", ())),
            full_field=bool(wf_raw.get("full_field", False)),
        )
        if wf_raw
        else None
    )
    return ArtifactSpec(
        cross_section=bool(raw.get("cross_section", True)),
        cross_section_vs_time=cvt,
        correlations=bool(raw.get("correlations", False)),
        wavefunction_snapshots=wf,
    )


# --- the top-level config ----------------------------------------------------


@dataclass(frozen=True)
class ExperimentConfig:
    """A parsed (not necessarily yet default-resolved) experiment config.

    `energies`/`grid.preset`/`td.incident`/`td.test_function`/`test_functions`/
    `Observable.channels` may all be `None` (omitted) straight out of
    `load_config` -- `qscat_run.presets.resolve_defaults` fills them from
    the molecule's preset. `validate_config` accepts either state (resolved
    or not) since it only checks *validity*, not completeness.
    """

    molecule: str
    methods: tuple[str, ...]
    observables: tuple[Observable, ...]
    output_dir: str
    energies: EnergySpec | None = None
    grid: GridSpec = field(default_factory=GridSpec)
    v_init: int = 0
    td: TdSpec | None = None
    artifacts: ArtifactSpec = field(default_factory=ArtifactSpec)
    backend: str = "auto"


_REQUIRED_KEYS = ("molecule", "methods", "observables", "output_dir")


def load_config(path: str | Path) -> ExperimentConfig:
    """Parse a YAML config into an `ExperimentConfig`. Tolerant of every
    optional block being omitted; raises `ConfigError` only if one of the
    four required top-level keys is missing."""
    raw = yaml.safe_load(Path(path).read_text()) or {}
    missing = [k for k in _REQUIRED_KEYS if k not in raw]
    if missing:
        raise ConfigError(f"config is missing required key(s): {missing}")

    return ExperimentConfig(
        molecule=str(raw["molecule"]),
        methods=tuple(str(m) for m in raw["methods"]),
        observables=_load_observables(raw["observables"]),
        output_dir=str(raw["output_dir"]),
        energies=_load_energies(raw.get("energies")),
        grid=_load_grid(raw.get("grid")),
        v_init=int(raw.get("v_init", 0)),
        td=_load_td(raw.get("td")),
        artifacts=_load_artifacts(raw.get("artifacts")),
        backend=str(raw.get("backend", "auto")),
    )


def validate_config(cfg: ExperimentConfig) -> None:
    """Raise `ConfigError` with an actionable message for the first problem
    found; emit a `UserWarning` (not an error) for an allowed-but-suspect
    combo (N2 `da`, closed-in-range).

    Checks, in order: molecule known -> methods are a non-empty subset of
    `{ti, td}` -> observables non-empty and each valid for the molecule ->
    `td` block present iff `"td" in methods` -> `td.extractors` all known ->
    an explicit grid supplies both `electronic` and `nuclear` -> a named
    preset (if given, with no explicit grid) exists for the molecule.
    """
    from qscat_run import presets  # local import: breaks the config<->presets cycle

    if cfg.molecule not in presets.MODELS:
        raise ConfigError(
            f"unknown molecule {cfg.molecule!r}; choose one of {sorted(presets.MODELS)}"
        )

    if not cfg.methods:
        raise ConfigError("'methods' must list at least one of {'ti', 'td'}")
    unknown_methods = sorted(set(cfg.methods) - VALID_METHODS)
    if unknown_methods:
        raise ConfigError(
            f"unknown method(s) {unknown_methods}; choose from {sorted(VALID_METHODS)}"
        )

    if not cfg.observables:
        raise ConfigError(
            "'observables' must list at least one observable, "
            "e.g. `observables: [{kind: ve, channels: 2}]`"
        )

    valid_kinds = presets.VALIDITY.get(cfg.molecule, frozenset())
    warn_kinds = presets.WARN_OBSERVABLES.get(cfg.molecule, frozenset())
    for obs in cfg.observables:
        if obs.kind not in valid_kinds:
            hint = " ('dr' is only defined for H2P)" if obs.kind == "dr" else ""
            raise ConfigError(
                f"observable {obs.kind!r} is not valid for {cfg.molecule}; "
                f"{cfg.molecule} supports {sorted(valid_kinds)}{hint}"
            )
        if obs.kind in warn_kinds:
            warnings.warn(
                f"observable {obs.kind!r} is closed-in-range for {cfg.molecule}: that "
                "channel only opens near the top of the studied energy range; treating "
                "it as valid, but expect strongly energy-dependent (and possibly zero) "
                "results away from there.",
                UserWarning,
                stacklevel=2,
            )

    if "td" in cfg.methods and cfg.td is None:
        raise ConfigError(
            "methods includes 'td' but no 'td' block is present; add a "
            "`td: {dt: ..., n_steps: ..., order: 3}` block (see the design spec's "
            "config schema for the full set of td keys)"
        )

    if cfg.td is not None:
        bad_extractors = sorted(e for e in cfg.td.extractors if e not in VALID_EXTRACTORS)
        if bad_extractors:
            raise ConfigError(
                f"unknown extractor(s) {bad_extractors} in 'td.extractors'; "
                f"choose from {sorted(VALID_EXTRACTORS)}"
            )

    grid = cfg.grid
    if (grid.electronic is None) != (grid.nuclear is None):
        missing_half = "nuclear" if grid.electronic is not None else "electronic"
        raise ConfigError(
            f"explicit grid is missing '{missing_half}': an explicit grid must provide "
            "both 'electronic' and 'nuclear' segments (or omit both and use "
            "`grid: {preset: ...}`)"
        )
    if grid.preset is not None and grid.electronic is None and grid.nuclear is None:
        available = presets.available_presets(cfg.molecule)
        if grid.preset not in available:
            raise ConfigError(
                f"unknown preset {grid.preset!r} for {cfg.molecule}; available: {sorted(available)}"
            )
