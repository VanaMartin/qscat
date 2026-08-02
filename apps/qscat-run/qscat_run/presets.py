"""The single source of the validated per-molecule PRODUCTION decks.

Consolidates what today lives in `validation/diatomic/config.py` (F2/NO) and
`validation/h2plus/config.py` (H2P) -- as literal DATA, copied by value, not
imported: `qscat_run` must never depend on `validation.*`/`projects.*` (see
`tests/test_no_validation_import.py`). N2's deck has no `validation/`
counterpart of its own (its numbers live in `projects.n2_2d_cross_section
.convergence.WORKING_GRID`/`projects.n2_2d_td_cross_section.convergence
.TD_WORKING_GRID`) and are transcribed here the same way.

Two grid roles per molecule, per the CLI design spec:
  - `ti_grid()` -- the driven-solve deck (small electronic box; F2/NO reuse
    their `da_grid()`-style fine per-molecule nuclear deck so DA is resolved
    too).
  - `td_grid()` -- the LAUNCH-BOX deck: a larger electronic `r_max` so an
    incident wavepacket launches, interacts, and (for DA/DR) its outgoing
    test function sits, entirely inside the real region. For F2/NO this is
    documented as a REASONABLE, not yet independently convergence-tested,
    choice (r_max=30, within the spec's suggested 25-40 bohr range) --
    unlike N2 (r_max=50) and H2+ (r_max=1300/60), which come from an actual
    Task-4-style convergence study or an eMoScat deck. See the module's
    `CONCERNS` note in the Task 1 report.

`PRESETS` keys are `"{molecule}:{variant}"` (e.g. `"F2:emoscat"`,
`"H2P:proxy"`) -- a molecule maps to one variant (N2/NO/F2, laptop-feasible
already) or two (H2P: `emoscat` the ~1.15M-unknown Docker/MUMPS deck,
`proxy` a laptop-sized reduction, both from `validation/h2plus/config.py`).
`DEFAULT_PRESET` (`"emoscat"`) is what an omitted `grid.preset` resolves to.

Kept pure: data + builder functions only, no CLI/YAML concerns.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, replace

from qscat.core.grids import electronic_grid, fem_grid_exp_tail, nuclear_grid, segmented_grid
from qscat.dvr import FemDvrEcsGrid, TensorGrid
from qscat.model import F2, H2P, N2, NO, ResonanceModel

from qscat_run.config import (
    ConfigError,
    EnergySpec,
    ExperimentConfig,
    IncidentSpec,
    SegmentSpec,
    TestFunctionSpec,
)

__all__ = [
    "MODELS",
    "VALIDITY",
    "WARN_OBSERVABLES",
    "DEFAULT_PRESET",
    "MoleculePreset",
    "PRESETS",
    "available_presets",
    "resolve_grid",
    "resolve_defaults",
]

MODELS: dict[str, ResonanceModel] = {"N2": N2, "NO": NO, "F2": F2, "H2P": H2P}

# The (molecule, observable.kind) validity matrix from the design spec.
# N2 "da" is CLOSED-IN-RANGE -- allowed, not rejected, but flagged
# (`WARN_OBSERVABLES`) since the channel is only weakly open near the top of
# the studied energy range.
VALIDITY: dict[str, frozenset[str]] = {
    "N2": frozenset({"ve", "da"}),
    "NO": frozenset({"ve", "da"}),
    "F2": frozenset({"ve", "da"}),
    "H2P": frozenset({"dr"}),
}
WARN_OBSERVABLES: dict[str, frozenset[str]] = {"N2": frozenset({"da"})}

DEFAULT_PRESET = "emoscat"


@dataclass(frozen=True)
class MoleculePreset:
    """One molecule/variant's numerical deck: the TI and TD grid builders,
    the default energy sweep, the default TD incident/test-function
    parameters, the observables this molecule supports, and the vibrational
    basis size to diagonalize."""

    molecule: str
    variant: str
    ti_grid: Callable[[], TensorGrid]
    td_grid: Callable[[], TensorGrid]
    default_energies: EnergySpec
    default_incident: IncidentSpec
    default_test_function: TestFunctionSpec
    valid_observables: frozenset[str]
    n_vib: int


# --- N2 -----------------------------------------------------------------
# `projects.n2_2d_cross_section.convergence.WORKING_GRID` (TI) /
# `projects.n2_2d_td_cross_section.convergence.TD_WORKING_GRID` (TD).


def _n2_ti_grid() -> TensorGrid:
    return TensorGrid(
        [
            electronic_grid(r_max=16.0, angle_deg=35.0, order=7, n_complex=5),
            nuclear_grid(angle_deg=35.0, r_max=20.0, n_complex=5, quadrature=10),
        ]
    )


def _n2_td_grid() -> TensorGrid:
    return TensorGrid(
        [
            electronic_grid(r_max=50.0, angle_deg=35.0, order=8, n_complex=6),
            nuclear_grid(angle_deg=35.0, r_max=22.0, n_complex=5, quadrature=10),
        ]
    )


# --- NO / F2 --------------------------------------------------------------
# `validation.diatomic.config.CONFIGS["NO"/"F2"]`: electronic_grid (VE deck)
# x the eMoScat per-molecule nuclear deck (`da_grid()`, via `segmented_grid`).

_NO_NUC_REAL = ((1, 1.0), (1, 1.6), (37, 9.0))
_NO_NUC_COMPLEX = ((1, 9.25), (1, 10.0), (1, 12.0), (4, 42.0))
_NO_NUC_ANGLE, _NO_NUC_QUAD = 45.0, 14

_F2_NUC_REAL = ((9, 1.8), (1, 2.0), (5, 2.5), (4, 2.596908), (4, 2.7), (40, 10.7))
_F2_NUC_COMPLEX = (
    (1, 10.8),
    (1, 11.0),
    (1, 11.5),
    (1, 12.5),
    (1, 14.0),
    (1, 18.0),
    (4, 30.0),
    (2, 101.0),
)
_F2_NUC_ANGLE, _F2_NUC_QUAD = 35.0, 14


def _no_nuc_grid() -> FemDvrEcsGrid:
    return segmented_grid(
        _NO_NUC_REAL, _NO_NUC_COMPLEX, angle_deg=_NO_NUC_ANGLE, quadrature=_NO_NUC_QUAD
    )


def _no_ti_grid() -> TensorGrid:
    return TensorGrid([electronic_grid(r_max=16.0, order=8, n_complex=6), _no_nuc_grid()])


def _no_td_grid() -> TensorGrid:
    # Launch-box: same fine nuclear deck, a larger electronic r_max (see
    # module docstring's caveat -- not independently convergence-tested).
    return TensorGrid([electronic_grid(r_max=30.0, order=8, n_complex=6), _no_nuc_grid()])


def _f2_nuc_grid() -> FemDvrEcsGrid:
    return segmented_grid(
        _F2_NUC_REAL, _F2_NUC_COMPLEX, angle_deg=_F2_NUC_ANGLE, quadrature=_F2_NUC_QUAD
    )


def _f2_ti_grid() -> TensorGrid:
    return TensorGrid([electronic_grid(r_max=16.0, order=8, n_complex=6), _f2_nuc_grid()])


def _f2_td_grid() -> TensorGrid:
    return TensorGrid([electronic_grid(r_max=30.0, order=8, n_complex=6), _f2_nuc_grid()])


# --- H2+ -------------------------------------------------------------------
# `validation.h2plus.config.full_grid`/`proxy_grid`, transcribed verbatim.
# TI and TD share the SAME grid here (no separate launch box): the full
# deck's electronic real region already runs to 1300 bohr, comfortably
# holding both the driven-solve deck and eMoScat's r0=800 incident.


def _h2p_full_grid() -> TensorGrid:
    electronic = fem_grid_exp_tail(
        ((10, 1.0), (10, 4.0), (16, 20.0), (20, 100.0), (120, 1300.0)),
        angle_deg=5.0,
        quadrature=8,
        tail_n=25,
    )
    nuclear = fem_grid_exp_tail(
        ((5, 1.0), (20, 4.0), (67, 14.0)),
        angle_deg=22.0,
        quadrature=8,
        tail_n=25,
    )
    return TensorGrid([electronic, nuclear])


def _h2p_proxy_grid() -> TensorGrid:
    electronic = fem_grid_exp_tail(
        ((10, 1.0), (10, 4.0), (16, 20.0), (10, 60.0)),
        angle_deg=5.0,
        quadrature=8,
        tail_n=8,
    )
    nuclear = fem_grid_exp_tail(
        ((5, 1.0), (20, 4.0), (40, 14.0)),
        angle_deg=22.0,
        quadrature=8,
        tail_n=8,
    )
    return TensorGrid([electronic, nuclear])


PRESETS: dict[str, MoleculePreset] = {
    "N2:emoscat": MoleculePreset(
        molecule="N2",
        variant="emoscat",
        ti_grid=_n2_ti_grid,
        td_grid=_n2_td_grid,
        # Approximates the dense curve's np.linspace(0.005, 0.2, 60)
        # (validation/n2/ti_curve.py) as a min/max/step sweep.
        default_energies=EnergySpec(min=0.005, max=0.20, step=0.005),
        default_incident=IncidentSpec(r0=25.0, p0=-0.5, sigma=5.0),
        default_test_function=TestFunctionSpec(r0_out=35.0, p0_out=0.5, sigma_out=4.0),
        valid_observables=VALIDITY["N2"],
        n_vib=6,
    ),
    "NO:emoscat": MoleculePreset(
        molecule="NO",
        variant="emoscat",
        ti_grid=_no_ti_grid,
        td_grid=_no_td_grid,
        default_energies=EnergySpec(min=0.004, max=0.120, step=0.004),
        # NO has no validated TD experiment yet (see CLAUDE.md's diatomic
        # note) -- these incident/test-function defaults are a reasonable,
        # UNVALIDATED scaling of N2's (see module docstring's caveat).
        default_incident=IncidentSpec(r0=20.0, p0=-0.4, sigma=4.0),
        default_test_function=TestFunctionSpec(r0_out=24.0, p0_out=0.4, sigma_out=3.0),
        valid_observables=VALIDITY["NO"],
        n_vib=4,
    ),
    "F2:emoscat": MoleculePreset(
        molecule="F2",
        variant="emoscat",
        ti_grid=_f2_ti_grid,
        td_grid=_f2_td_grid,
        default_energies=EnergySpec(min=0.004, max=0.100, step=0.004),
        # Same caveat as NO: unvalidated TD defaults, scaled from N2's.
        default_incident=IncidentSpec(r0=20.0, p0=-0.4, sigma=4.0),
        default_test_function=TestFunctionSpec(r0_out=24.0, p0_out=0.4, sigma_out=3.0),
        valid_observables=VALIDITY["F2"],
        n_vib=4,
    ),
    "H2P:emoscat": MoleculePreset(
        molecule="H2P",
        variant="emoscat",
        ti_grid=_h2p_full_grid,
        td_grid=_h2p_full_grid,
        default_energies=EnergySpec(min=0.001, max=0.050, step=0.001),
        default_incident=IncidentSpec(r0=800.0, p0=-0.25, sigma=8.0),
        default_test_function=TestFunctionSpec(r0_out=12.0, p0_out=15.0, sigma_out=0.4),
        valid_observables=VALIDITY["H2P"],
        n_vib=4,
    ),
    "H2P:proxy": MoleculePreset(
        molecule="H2P",
        variant="proxy",
        ti_grid=_h2p_proxy_grid,
        td_grid=_h2p_proxy_grid,
        default_energies=EnergySpec(min=0.001, max=0.050, step=0.001),
        # r0 scaled DOWN from the full deck's 800 to fit inside the proxy
        # grid's ~60-bohr electronic real region (an off-box incident lands
        # in the ECS tail and diverges -- see validation/h2plus/td_dr.py).
        default_incident=IncidentSpec(r0=40.0, p0=-0.25, sigma=8.0),
        default_test_function=TestFunctionSpec(r0_out=12.0, p0_out=15.0, sigma_out=0.4),
        valid_observables=VALIDITY["H2P"],
        n_vib=4,
    ),
}


def available_presets(molecule: str) -> frozenset[str]:
    """The preset variant names (e.g. `{"emoscat", "proxy"}`) defined for `molecule`."""
    prefix = f"{molecule}:"
    return frozenset(key[len(prefix) :] for key in PRESETS if key.startswith(prefix))


def _build_explicit_segment(seg: SegmentSpec) -> FemDvrEcsGrid:
    """One explicit-grid axis -> a `FemDvrEcsGrid`, via `fem_grid_exp_tail`
    (the ECS spec gives `angle`/`elements`/`quadrature` but no explicit tail
    endpoint, so the exponential-tail builder -- not `segmented_grid`, which
    needs an endpoint for every segment -- is the one that fits)."""
    return fem_grid_exp_tail(
        seg.real,
        angle_deg=seg.ecs.angle,
        quadrature=seg.ecs.quadrature,
        tail_n=seg.ecs.elements,
    )


def resolve_grid(cfg: ExperimentConfig, method: str) -> TensorGrid:
    """The `TensorGrid` for `method` (`"ti"` or `"td"`): an explicit grid
    (if `cfg.grid.electronic`/`nuclear` are both given) used for BOTH
    methods, else the named preset's `ti_grid()`/`td_grid()`."""
    grid = cfg.grid
    if grid.electronic is not None and grid.nuclear is not None:
        return TensorGrid(
            [_build_explicit_segment(grid.electronic), _build_explicit_segment(grid.nuclear)]
        )

    variant = grid.preset or DEFAULT_PRESET
    key = f"{cfg.molecule}:{variant}"
    preset = PRESETS.get(key)
    if preset is None:
        raise ConfigError(
            f"no preset {variant!r} for molecule {cfg.molecule!r}; "
            f"available: {sorted(available_presets(cfg.molecule))}"
        )
    if method == "ti":
        return preset.ti_grid()
    if method == "td":
        return preset.td_grid()
    raise ConfigError(f"unknown method {method!r}; choose 'ti' or 'td'")


def _default_channels(kind: str, preset: MoleculePreset) -> int | tuple[int, ...]:
    """A reasonable channel default when the config omits `channels`: for
    `ve`, the first few excited vibrational levels (bounded by `n_vib`);
    for `da`/`dr`, a single channel (the usual case)."""
    if kind == "ve":
        return tuple(range(1, min(preset.n_vib, 4)))
    return 1


def resolve_defaults(cfg: ExperimentConfig) -> ExperimentConfig:
    """Fill omitted `energies`/`td.incident`/`td.test_function`/
    `Observable.channels` from the resolved preset. A no-op (returns `cfg`
    unchanged) if no preset matches (an explicit grid with no matching
    preset key, or an unknown molecule/variant) -- `validate_config` is
    responsible for rejecting that case with an actionable message."""
    variant = cfg.grid.preset or DEFAULT_PRESET
    preset = PRESETS.get(f"{cfg.molecule}:{variant}")
    if preset is None:
        return cfg

    energies = cfg.energies if cfg.energies is not None else preset.default_energies
    observables = tuple(
        obs
        if obs.channels is not None
        else replace(obs, channels=_default_channels(obs.kind, preset))
        for obs in cfg.observables
    )
    td = cfg.td
    if td is not None:
        incident = td.incident if td.incident is not None else preset.default_incident
        test_function = (
            td.test_function if td.test_function is not None else preset.default_test_function
        )
        td = replace(td, incident=incident, test_function=test_function)

    return replace(cfg, energies=energies, observables=observables, td=td)
