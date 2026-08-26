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

from qscat_run import reference as _reference
from qscat_run.reference import ReferenceSpec

__all__ = [
    "VALID_EXTRACTORS",
    "VALID_METHODS",
    "VALID_NRM_CHOICES",
    "ArtifactSpec",
    "ConfigError",
    "CrossSectionVsTimeSpec",
    "EcsSpec",
    "EnergySpec",
    "ExperimentConfig",
    "GridSpec",
    "IncidentSpec",
    "NrmSpec",
    "Observable",
    "ReferenceSpec",
    "SegmentSpec",
    "TdSpec",
    "TestFunctionSpec",
    "WavefunctionSnapshotsSpec",
    "load_config",
    "validate_config",
]


class ConfigError(click.ClickException):
    """An actionable, user-facing config problem (unknown molecule, an
    invalid (molecule, observable) combo, a missing `td` block, ...).

    A `click.ClickException` subclass so `cli.py` commands can simply let it
    propagate: click's own error handling prints `Error: <message>` and sets
    a non-zero exit code without any command-level try/except.
    """


VALID_METHODS = frozenset({"ti", "td", "lcp", "nrm"})
VALID_EXTRACTORS = frozenset({"flow", "delta", "tw"})
# PRA 77's two implemented discrete-state choices (Sec. VI A / VI B); the
# paper's third, "compact" choice C is not implemented -- see
# docs/physics/nonlocal-resonance-model.md.
VALID_NRM_CHOICES = frozenset({"a", "b"})


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
    # Tail ECS angle of the SECOND nuclear grid used for two-angle level
    # selection. Defaults to `angle_a - 10` degrees (eMoScat's decks pair
    # 44/35 and 40/30). Real segments and quadrature are always copied from
    # grid a -- the shared real region is what makes the comparison valid.
    nuclear_angle_b: float | None = None


def _load_grid(raw: dict[str, Any] | None) -> GridSpec:
    raw = raw or {}
    electronic = _load_segment(raw["electronic"]) if "electronic" in raw else None
    nuclear = _load_segment(raw["nuclear"]) if "nuclear" in raw else None
    preset = raw.get("preset")
    return GridSpec(
        preset=str(preset) if preset is not None else None,
        electronic=electronic,
        nuclear=nuclear,
        nuclear_angle_b=(
            float(raw["nuclear_angle_b"]) if raw.get("nuclear_angle_b") is not None else None
        ),
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


# --- nrm ----------------------------------------------------------------------


@dataclass(frozen=True)
class NrmSpec:
    """The `nrm` block -- optional even when `"nrm" in methods` (these
    defaults then apply, and `presets.resolve_defaults` writes them into
    `config.resolved.yaml`).

    `choices` names PRA 77's discrete-state choices to run: `"a"` the
    R-dependent "physical" state (Sec. VI A), `"b"` the R-independent
    asymptotic bound state (Sec. VI B). Each requested choice gets its own
    cross-section key (`nrm-a:da:ch0`, `nrm-b:ve:v0->1`, ...), so one run
    overlays both against `ti`/`lcp`. `"b"` alone is the default because it is
    the choice PRA 77 shows reproducing the exact F2 DA cross section, and the
    one this repo measured at 0.06-1.9 % of the exact oracle for DA and within
    0.7 % of it for vibrational excitation.

    `n_states` truncates the Eq. (60) sum over projected electronic states.
    100 is the measured value: the F2/NO x A/B ladders reproduce the
    untruncated sum to numerical identity there, and the sum is NOT
    front-loaded (n=50 is still 33 % off) -- see
    docs/physics/nonlocal-resonance-model.md.

    `include_background` adds the Eq. (37) background T-matrix to the resonant
    one before squaring -- PRA 77's "nonlocal + background" curve as against
    its bare "nonlocal" one (Figs. 4-6 and 8 plot both, and the difference is
    the paper's own argument for why a bare LCP curve is missing something).
    It applies to `ve` only: `da` has no background term in this model, and a
    `da` observable ignores the flag.
    """

    choices: tuple[str, ...] = ("b",)
    n_states: int = 100
    include_background: bool = True


def _load_nrm(raw: dict[str, Any] | None) -> NrmSpec | None:
    if raw is None:
        return None
    choices = raw.get("choices")
    return NrmSpec(
        choices=("b",) if choices is None else tuple(str(c).lower() for c in choices),
        n_states=int(raw.get("n_states", 100)),
        include_background=bool(raw.get("include_background", True)),
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
    # Emit the target's vibrational energy levels + their eigenstate wavefunctions
    # (the eps/chi already diagonalized for the cross section) as an artifact.
    eigenstates: bool = False
    # Emit the quasi-bound vibrational levels of the anion in the LCP complex
    # potential (`qscat.core.lcp.resonance_levels`) -- the BO approximation to
    # the resonance energies. LCP path only.
    resonance_levels: bool = False


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
        eigenstates=bool(raw.get("eigenstates", False)),
        resonance_levels=bool(raw.get("resonance_levels", False)),
    )


# --- reference -----------------------------------------------------------


def _load_reference(raw: list[Any] | None) -> tuple[ReferenceSpec, ...]:
    if not raw:
        return ()
    out: list[ReferenceSpec] = []
    for i, item in enumerate(raw):
        if not isinstance(item, dict) or "path" not in item:
            raise ConfigError(f"reference[{i}] must be a mapping with a 'path' key")
        chans = item.get("channels")
        out.append(
            ReferenceSpec(
                path=str(item["path"]),
                format=str(item.get("format", "houfek")),
                label=None if item.get("label") is None else str(item["label"]),
                channels=None if chans is None else tuple(int(c) for c in chans),
            )
        )
    return tuple(out)


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
    nrm: NrmSpec | None = None
    artifacts: ArtifactSpec = field(default_factory=ArtifactSpec)
    reference: tuple[ReferenceSpec, ...] = ()
    backend: str = "auto"
    # The directory the config YAML itself lives in, so relative `reference`
    # paths resolve against the config file rather than the CWD `qscat-run`
    # happens to be invoked from. Kept LAST so existing positional
    # `ExperimentConfig(...)` construction (there is none left, but future
    # callers) is unaffected by this field's addition.
    config_dir: str | None = None


_REQUIRED_KEYS = ("molecule", "methods", "observables", "output_dir")


def _molecule_name(raw: Any) -> str:
    """`raw` as a molecule name, with a specific error for the YAML-1.1 trap.

    PyYAML implements YAML 1.1, where a bare ``NO`` is the BOOLEAN false, not
    the string "NO". A hand-written ``molecule: NO`` therefore arrives here as
    ``False`` and, coerced with ``str()``, produced "unknown molecule 'False'"
    -- an error that says nothing about the real cause. ``qscat-run init``
    already quotes the name, so this only bites configs written by hand.
    """
    if isinstance(raw, bool):
        return _raise_yaml_bool_molecule(raw)
    return str(raw)


def _raise_yaml_bool_molecule(raw: bool) -> str:
    got = "true" if raw else "false"
    raise ConfigError(
        f"'molecule' parsed as the YAML boolean {got}, not a name. YAML 1.1 reads "
        f"bare NO/No/no (and ON/OFF/YES) as booleans -- quote it: molecule: 'NO'"
    )


def load_config(path: str | Path) -> ExperimentConfig:
    """Parse a YAML config into an `ExperimentConfig`. Tolerant of every
    optional block being omitted; raises `ConfigError` only if one of the
    four required top-level keys is missing."""
    raw = yaml.safe_load(Path(path).read_text()) or {}
    missing = [k for k in _REQUIRED_KEYS if k not in raw]
    if missing:
        raise ConfigError(f"config is missing required key(s): {missing}")

    return ExperimentConfig(
        molecule=_molecule_name(raw["molecule"]),
        methods=tuple(str(m) for m in raw["methods"]),
        observables=_load_observables(raw["observables"]),
        output_dir=str(raw["output_dir"]),
        energies=_load_energies(raw.get("energies")),
        grid=_load_grid(raw.get("grid")),
        v_init=int(raw.get("v_init", 0)),
        td=_load_td(raw.get("td")),
        nrm=_load_nrm(raw.get("nrm")),
        artifacts=_load_artifacts(raw.get("artifacts")),
        reference=_load_reference(raw.get("reference")),
        backend=str(raw.get("backend", "auto")),
        config_dir=str(Path(path).resolve().parent),
    )


def validate_config(cfg: ExperimentConfig) -> None:
    """Raise `ConfigError` with an actionable message for the first problem
    found; emit a `UserWarning` (not an error) for an allowed-but-suspect
    combo (N2 `da`, closed-in-range).

    Checks, in order: molecule known -> methods are a non-empty subset of
    `VALID_METHODS` -> observables non-empty and each valid for the molecule
    -> `td` block present iff `"td" in methods` -> `lcp`/`nrm` each get the
    grid form and molecule they need -> `td.extractors` all known
    -> an explicit grid supplies both `electronic` and `nuclear` -> a named
    preset (if given, with no explicit grid) exists for the molecule ->
    each `reference` entry names a known `format`, resolves to a file that
    actually exists, and (if `channels` is given) requests only channels
    the file actually has.
    """
    from qscat_run import presets  # local import: breaks the config<->presets cycle

    if cfg.molecule not in presets.MODELS:
        raise ConfigError(
            f"unknown molecule {cfg.molecule!r}; choose one of {sorted(presets.MODELS)}"
        )

    if not cfg.methods:
        raise ConfigError(f"'methods' must list at least one of {sorted(VALID_METHODS)}")
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

    if "lcp" in cfg.methods:
        # LCP is the local-complex-potential APPROXIMATION of the exact 2-D
        # solve, so it needs a molecule with an LCP path (N2/F2/NO -- H2P is
        # DR, not DA/VE-LCP) and the preset grids (no explicit-grid schema for
        # the two ECS-angle electronic decks + fine nuclear deck). No
        # per-observable-kind check is needed here: every kind that survives
        # the per-molecule validity check above (`ve`/`da` for N2, plus
        # `resonance_levels` for F2/NO) is one LCP can serve, so a run that
        # reaches this point already requests only LCP-servable observables.
        if cfg.grid.electronic is not None or cfg.grid.nuclear is not None:
            raise ConfigError(
                "the 'lcp' method does not support an explicit grid (it needs the "
                "preset's paired two-ECS-angle electronic + fine nuclear decks); "
                "use `grid: {preset: ...}` (or omit grid) with methods including 'lcp'"
            )
        variant = cfg.grid.preset or presets.DEFAULT_PRESET
        lcp_preset = presets.PRESETS.get(f"{cfg.molecule}:{variant}")
        if lcp_preset is None or lcp_preset.lcp_grids is None:
            raise ConfigError(
                f"the 'lcp' method is not available for {cfg.molecule}; LCP is "
                "wired for the neutral diatomics (N2, F2, NO); H2P's observable is DR"
            )

    if "nrm" in cfg.methods:
        # The nonlocal resonance model approximates the SAME two observables the
        # exact `ti` solve gives: vibrational excitation (PRA 77 Eq. 28/31/37)
        # and dissociative attachment (Eq. 52-54). Like `lcp` it needs the
        # preset's electronic deck at two ECS angles plus the molecule's own
        # nuclear deck, so it has no explicit-grid form.
        variant = cfg.grid.preset or presets.DEFAULT_PRESET
        nrm_preset = presets.PRESETS.get(f"{cfg.molecule}:{variant}")
        if nrm_preset is None or nrm_preset.nrm_elec_b is None:
            # Checked first: "not available for this molecule" is the more
            # actionable message than "no ve/da observable" for a molecule
            # (H2P) whose only observable the NRM could never serve anyway.
            raise ConfigError(
                f"the 'nrm' method is not available for {cfg.molecule}; the nonlocal "
                "resonance model is wired for the neutral diatomics (N2, NO, F2)"
            )
        kinds = {obs.kind for obs in cfg.observables}
        if not kinds & {"ve", "da"}:
            raise ConfigError(
                "methods includes 'nrm' but no 've' or 'da' observable is requested; "
                "the nonlocal resonance model approximates those two cross sections "
                "-- add `{kind: ve, channels: 2}` or `{kind: da, channels: 1}`"
            )
        if cfg.grid.electronic is not None or cfg.grid.nuclear is not None:
            raise ConfigError(
                "the 'nrm' method does not support an explicit grid (it needs the "
                "preset's electronic deck at two ECS angles + the fine nuclear "
                "deck); use `grid: {preset: ...}` (or omit grid) with methods "
                "including 'nrm'"
            )
        if cfg.nrm is not None:
            if not cfg.nrm.choices:
                raise ConfigError(
                    "'nrm.choices' must name at least one discrete-state choice; "
                    f"choose from {sorted(VALID_NRM_CHOICES)}"
                )
            bad_choices = sorted(set(cfg.nrm.choices) - VALID_NRM_CHOICES)
            if bad_choices:
                raise ConfigError(
                    f"unknown nrm discrete-state choice(s) {bad_choices}; choose from "
                    f"{sorted(VALID_NRM_CHOICES)} ('a' = the R-dependent physical "
                    "state, 'b' = the R-independent asymptotic bound state; PRA 77's "
                    "third 'compact' choice is not implemented)"
                )
            if cfg.nrm.n_states < 1:
                raise ConfigError(
                    f"'nrm.n_states' must be >= 1, got {cfg.nrm.n_states}; 100 is the "
                    "measured converged value (see docs/physics/"
                    "nonlocal-resonance-model.md)"
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

    base = _reference.config_base_dir(cfg.config_dir)
    for i, ref in enumerate(cfg.reference):
        if ref.format not in _reference.REFERENCE_FORMATS:
            raise ConfigError(
                f"reference[{i}]: unknown format {ref.format!r}; "
                f"choose one of {sorted(_reference.REFERENCE_FORMATS)}"
            )
        resolved = _reference.resolve_path(ref, base)
        if not resolved.is_file():
            raise ConfigError(f"reference[{i}]: no such file {ref.path!r} (looked at {resolved})")
        # Bounds-check `channels` here too (not just at `load_reference` time,
        # which `write_artifacts` only reaches AFTER `run_experiment` has
        # already solved) -- a typo'd channel index should fail fast, before
        # the solve, like every other config problem. `peek_n_channels` reads
        # only the file's first line, so this stays cheap even on a large
        # production dataset.
        try:
            n_channels = _reference.peek_n_channels(resolved)
        except (OSError, ValueError) as exc:
            raise ConfigError(f"reference[{i}]: could not read {resolved}: {exc}") from exc
        bad = _reference.bad_channels(ref.channels, n_channels)
        if bad:
            raise ConfigError(
                f"reference[{i}]: requested channel(s) {bad} but {resolved} has "
                f"{n_channels} column(s) (valid 0..{n_channels - 1})"
            )
