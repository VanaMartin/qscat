"""The one-shot a-priori grid assembler: `analyze` -> `mesh` -> `ecs`, wired
into a complete `FemDvrEcsGrid`.

`propose_grid` is the a-priori HALF of the hybrid tuner design (see
`docs/superpowers/specs/2026-07-28-discretisation-tuner-design.md`): given a
`ResonanceModel`, a coordinate ("nuclear" or "electronic"), and a target
energy range, it builds the potential-adaptive grid directly from the
potential curve -- no eigensolve, no probing. The `discretisation-tuner`
skill runs the probe/refine loop ON TOP of this a-priori starting point; this
module owns only the a-priori half.

The MODEL ADAPTER (`_nuclear_adapter`/`_electronic_adapter`) is the only
model-aware code here: it picks the per-coordinate 1-D potential, mass, real
extent, and channel wavenumber from the model + energy range. Everything
downstream (`analyze_potential` -> `optimal_real_mesh` -> `max_stable_angle`
+ `tune_ecs_tail`) is the same model-independent pipeline already validated
in Tasks 1-3.
"""

from __future__ import annotations

import dataclasses
import math
from collections.abc import Callable
from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

import numpy as np
import numpy.typing as npt

from qscat.core.grids import electronic_grid
from qscat.dvr import ElementSpec, FemDvrEcsGrid, GridSpec

from .analyze import PotentialProfile, analyze_potential
from .ecs import max_stable_angle, tune_ecs_tail
from .mesh import combined_profile, optimal_real_mesh
from .resonance import resonance_curve

if TYPE_CHECKING:
    from qscat.model import ResonanceModel

__all__ = ["propose_grid"]

Coordinate = Literal["nuclear", "electronic"]

# Sane real-region cutoff defaults (bohr) -- see the design spec's per-
# coordinate guidance ("14-22 bohr" nuclear, "~16-30 bohr" electronic); a
# fixed default here, refined by the h/p mesh sweep and (Task 6) extended by
# an `incident` requirement. Task 8 calibrated the mesh's phase constant `C`
# against the eMoScat decks but did NOT retune these fixed extents: they are
# a per-molecule-INDEPENDENT default, which Task 8's calibration found
# exceeds N2's (12 bohr) / NO's (9 bohr) committed nuclear real regions --
# the reason those two molecules' proposed grids cost more DVR points than
# their decks even at the calibrated C (see docs/physics/
# discretisation-tuning.md). Deriving x_max from the potential profile itself
# (rather than this fixed constant) is a follow-on, not done here.
_NUCLEAR_X_MAX_DEFAULT = 18.0
_ELECTRONIC_X_MAX_DEFAULT = 20.0

# Element-length bounds (bohr) fed to `optimal_real_mesh`'s equidistribution
# sweep -- sane, not independently calibrated (Task 8 calibrated only the
# phase constant `C`). The nuclear floor is deliberately not tiny: a heavy
# reduced mass mu inflates the local
# wavenumber k = sqrt(2*mu*(e_max-V)) deep in the classically-forbidden well
# wall, which (via the kappa-decay-length branch) would otherwise carve that
# region into far more sub-min_len elements than the physics needs -- 0.15
# bohr (close to the eMoScat N2 deck's own finest equilibrium-region
# spacing, see qscat.core.grids.nuclear_grid's _REAL_SEGMENTS) keeps the
# proposed grid's point count comparable to that deck.
_NUCLEAR_MIN_LEN = 0.15
_NUCLEAR_MAX_LEN = 2.0
_ELECTRONIC_MIN_LEN = 0.05
_ELECTRONIC_MAX_LEN = 3.0

# How far out (bohr) `max_stable_angle` probes the rotated tail for growth --
# generous enough to catch a diverging continuation before the actual tail
# length (from `tune_ecs_tail`) is known; a fixed, molecule-independent
# constant (matches the scale used in `test_tuning_ecs.py`).
_ANGLE_PROBE_TAIL_EXTENT = 40.0

PotentialFn = Callable[[npt.ArrayLike], npt.NDArray[np.complexfloating]]

# Default two-angle electronic grids for the resonance-pole match
# (`qscat.tuning.resonance.resonance_curve`) driving the resonant nuclear
# path (`channel="dissociation"`) -- two different ECS rotation angles are
# needed so `find_resonance_pole` can triangulate the pole across them.
# Overridable via `propose_grid`'s `elec_grids` so tests can inject small
# grids instead of paying for these full-size eigensolves.
_RESONANCE_ELEC_R_MAX = 16.0
_RESONANCE_ELEC_ORDER = 7
_RESONANCE_ELEC_N_COMPLEX = 6
_RESONANCE_ELEC_ANGLE_A = 35.0
_RESONANCE_ELEC_ANGLE_B = 44.0


@dataclass(frozen=True)
class _CoordinateSpec:
    """What the model adapter hands the coordinate-independent pipeline.

    `charge` is not consumed by this pipeline (the potential callable `V`
    already embeds any Coulomb tail for a charged residual channel -- see
    `qscat.model.IonicResonanceModel`); it is exposed here only because the
    design spec's per-coordinate adapter names it as one of the quantities
    the adapter picks, for a future consumer (e.g. a channel-representation
    probe run on top of this grid) to read off the same model without
    re-deriving it.
    """

    V: PotentialFn
    mass: float
    charge: int
    x_min: float
    x_max: float
    channel_k: float
    min_len: float
    max_len: float


def _nuclear_adapter(model: ResonanceModel, e_max: float) -> _CoordinateSpec:
    """Nuclear coordinate: `V = model.v0`, `m = model.mu`, real region [0,
    x_max], outgoing dissociation wavenumber `K = sqrt(2*mu*e_max)` (heavy
    reduced mass -> large K, per the design spec's nuclear-adapter note).
    """
    return _CoordinateSpec(
        V=model.v0,
        mass=model.mu,
        charge=0,  # the nuclear channel itself carries no Coulomb charge
        x_min=0.0,
        x_max=_NUCLEAR_X_MAX_DEFAULT,
        channel_k=math.sqrt(2.0 * model.mu * e_max),
        min_len=_NUCLEAR_MIN_LEN,
        max_len=_NUCLEAR_MAX_LEN,
    )


def _well_minimum(model: ResonanceModel) -> float:
    """The bond length `R` minimizing `Re(model.v0(R))` -- the Morse-like
    neutral-curve equilibrium -- found by a plain scan over a generous
    diatomic bond-length range. Uses only `v0` (part of the `ResonanceModel`
    protocol), so this works for any model rather than reaching into a
    concrete class's private equilibrium-distance field.
    """
    r_scan = np.linspace(0.5, 8.0, 4000)
    v_scan = np.real(model.v0(r_scan))
    return float(r_scan[np.argmin(v_scan)])


def _electronic_adapter(model: ResonanceModel, e_max: float) -> _CoordinateSpec:
    """Electronic coordinate: the effective electronic potential
    `V(r) = model.surface(r, R_eq) - model.v0(R_eq)` at the neutral-curve
    well minimum `R_eq` (subtracting off the R-dependent neutral-curve
    constant so the profile reflects only the centrifugal + interaction
    terms an electron at fixed R_eq feels), mass 1, incident wavenumber
    `k = sqrt(2*e_max)`, and the model's residual-channel `charge`.
    """
    r_eq = _well_minimum(model)

    def v_elec(r: npt.ArrayLike) -> npt.NDArray[np.complexfloating]:
        return np.asarray(model.surface(r, r_eq) - model.v0(r_eq), dtype=np.complex128)

    return _CoordinateSpec(
        V=v_elec,
        mass=1.0,
        charge=model.charge,
        x_min=0.0,
        x_max=_ELECTRONIC_X_MAX_DEFAULT,
        channel_k=math.sqrt(2.0 * e_max),
        min_len=_ELECTRONIC_MIN_LEN,
        max_len=_ELECTRONIC_MAX_LEN,
    )


def _resonant_nuclear_profile(
    model: ResonanceModel,
    spec: _CoordinateSpec,
    x_max: float,
    e_max_mesh: float,
    elec_grids: tuple[FemDvrEcsGrid, FemDvrEcsGrid] | None,
    resonance_n_dense: int,
) -> PotentialProfile:
    """The `channel="dissociation"` nuclear profile: the worst-case merge of
    `v0`-alone with the adiabatic resonance curve `V_d(R)`, plus a turning-
    point injected at the `Gamma(R)` peak so `equidistribution_elements`'s
    existing feature-refinement halves the elements straddling it.

    Runs the (expensive) two-angle pole match `resonance_curve` -- callers
    needing this cheap may inject small `elec_grids` (see `propose_grid`).
    """
    if elec_grids is None:
        ga = electronic_grid(
            r_max=_RESONANCE_ELEC_R_MAX,
            order=_RESONANCE_ELEC_ORDER,
            n_complex=_RESONANCE_ELEC_N_COMPLEX,
            angle_deg=_RESONANCE_ELEC_ANGLE_A,
        )
        gb = electronic_grid(
            r_max=_RESONANCE_ELEC_R_MAX,
            order=_RESONANCE_ELEC_ORDER,
            n_complex=_RESONANCE_ELEC_N_COMPLEX,
            angle_deg=_RESONANCE_ELEC_ANGLE_B,
        )
    else:
        ga, gb = elec_grids

    R, Vd, Gamma = resonance_curve(model, ga, gb, R_max=x_max, n_dense=resonance_n_dense)
    Vd_real = np.real(Vd)
    vd_left = float(Vd_real[0])
    vd_right = float(Vd_real[-1])

    def Vd_of_R(rr: npt.ArrayLike) -> npt.NDArray[np.float64]:
        # Constant-extrapolate both ends -- outside the sampled range V_d
        # has saturated to its asymptote (see `resonance_curve`'s docstring:
        # a single far point standing in for the whole outer region).
        return np.asarray(
            np.interp(np.real(np.asarray(rr)), R, Vd_real, left=vd_left, right=vd_right),
            dtype=np.float64,
        )

    profile_v0 = analyze_potential(spec.V, spec.x_min, x_max, spec.mass, e_max_mesh)
    profile_vd = analyze_potential(Vd_of_R, spec.x_min, x_max, spec.mass, e_max_mesh)
    combined = combined_profile([profile_v0, profile_vd])

    if Gamma.max() > 0.0:
        r_peak = float(R[int(np.argmax(Gamma))])
        turning_points = np.unique(np.concatenate([combined.turning_points, [r_peak]]))
        combined = dataclasses.replace(combined, turning_points=turning_points)

    return combined


def propose_grid(
    model: ResonanceModel,
    coordinate: Coordinate,
    energy_range: tuple[float, float],
    *,
    rtol: float = 1e-3,
    incident: object | None = None,
    phase_coeff: float | None = None,
    channel: str = "ve",
    elec_grids: tuple[FemDvrEcsGrid, FemDvrEcsGrid] | None = None,
    resonance_n_dense: int = 25,
) -> FemDvrEcsGrid:
    """The one-shot a-priori `FemDvrEcsGrid` for `model`/`coordinate` over
    `energy_range = (E_min, E_max)`.

    Pipeline: a small per-coordinate MODEL ADAPTER (`_nuclear_adapter` /
    `_electronic_adapter`) picks `V`/mass/extent/channel-`k` -> `analyze_
    potential` profiles it at `E_max` -> `optimal_real_mesh` sweeps the h/p
    equidistribution mesh -> `max_stable_angle` finds the ECS rotation angle
    -> `tune_ecs_tail` sizes the absorbing tail for the channel wavenumber ->
    the resulting real + tail `ElementSpec` list becomes a `GridSpec` /
    `FemDvrEcsGrid`.

    `rtol` is accepted for interface parity with the probe/refine loop this
    feeds (the `discretisation-tuner` skill); this a-priori assembler itself
    does not probe or refine -- it has no eigensolve to converge.

    `phase_coeff`, if given, overrides `optimal_real_mesh`'s calibrated
    default de-Broglie phase-per-`(order-1)` constant `C` -- the knob
    `validation.tuning.calibrate` sweeps to find that calibrated value
    against the eMoScat decks. `None` (the default) leaves
    `optimal_real_mesh` at its own calibrated default; ordinary callers never
    need to pass this.

    `incident` (`qscat.tuning.incident.IncidentSpec`, Task 6) is accepted
    here as BOTH an extent floor AND a resolution floor:

    - EXTENT: `getattr(incident, "required_extent", lambda: 0.0)()` extends
      the real-region cutoff to at least that value before the mesh/ECS
      steps run.
    - RESOLUTION: `getattr(incident, "incident_energy", lambda: 0.0)()`
      raises the effective `E_max` fed to `analyze_potential`/
      `optimal_real_mesh` to `max(E_max, incident_energy)`. Without this, a
      HAND-BUILT `IncidentSpec` whose `impulse` implies an energy above
      `energy_range`'s own `E_max` would still get a mesh sized only for the
      (lower) `energy_range`'s local wavenumber -- silently under-resolving
      the incident's actual wave. (A `tw_analysis`-produced `IncidentSpec`
      never triggers this: its energy is bounded by `energy_range` by
      construction, so `max(...)` is then a no-op.) The ECS-tail
      `channel_k` is deliberately left keyed to `energy_range`'s own
      `E_max` alone, not this raised value: the incident wavepacket is a
      TD object CONTAINED in the real region by the extent floor above, not
      an outgoing wave the tail need absorb -- widening `channel_k` too is
      a plausible future refinement, not required by this baseline.

    Both getattrs default to `0.0` (a no-op) for any duck-typed `incident`
    that does not define the corresponding method; `IncidentSpec` defines
    both, precisely so these duck-typed calls work against it unchanged --
    see `qscat.tuning.incident`'s docstring for the reconciliation. The
    placement logic itself (impulse/width/observation boundary,
    `tw_analysis`) lives there; this is only the extent/resolution floor.

    `channel` selects which physical channel the mesh targets:

    - `"ve"` (the default): the VE (vibrational-excitation) path, v0-alone,
      exactly as before this parameter existed -- BYTE-IDENTICAL to the
      pre-`channel` behavior; nothing below reads `elec_grids` or
      `resonance_n_dense` on this path.
    - `"dissociation"`, `coordinate="nuclear"` only: the resonance-aware
      nuclear path. Builds the adiabatic resonance curve
      `(R, V_d(R), Gamma(R))` (`qscat.tuning.resonance.resonance_curve`,
      via a two-angle ECS pole match -- `elec_grids`, if given, overrides
      the default electronic grids used for that match, and
      `resonance_n_dense` overrides its dense-sampling point count; both
      exist so tests can inject small/cheap grids), then feeds
      `optimal_real_mesh` the WORST CASE of `v0`-alone and `V_d(R)`
      (`qscat.tuning.mesh.combined_profile`: elementwise max-`k`/min-
      `kappa`, unioned turning points), with an extra turning point
      injected at the `Gamma(R)` peak so the existing near-feature
      refinement halves the elements straddling the resonance crossing.
      `"dissociation"` with `coordinate="electronic"` raises `ValueError`
      (this resonant path is nuclear-only in this sub-project).
    - any other value raises `ValueError`.
    """
    del rtol  # interface parity with the probe/refine loop; unused here

    if channel not in ("ve", "dissociation"):
        raise ValueError(f"channel must be 've' or 'dissociation', got {channel!r}")
    if channel == "dissociation" and coordinate != "nuclear":
        raise ValueError(f"channel='dissociation' is nuclear-only, got coordinate={coordinate!r}")

    e_min, e_max = energy_range
    if e_max <= e_min:
        raise ValueError(f"energy_range must have E_max > E_min, got {energy_range}")

    adapter = _nuclear_adapter if coordinate == "nuclear" else _electronic_adapter
    spec = adapter(model, e_max)

    x_max = spec.x_max
    e_max_mesh = e_max
    if incident is not None:
        required_extent = getattr(incident, "required_extent", lambda: 0.0)()
        x_max = max(x_max, float(required_extent))
        incident_energy = getattr(incident, "incident_energy", lambda: 0.0)()
        e_max_mesh = max(e_max_mesh, float(incident_energy))

    if channel == "dissociation":
        profile = _resonant_nuclear_profile(
            model, spec, x_max, e_max_mesh, elec_grids, resonance_n_dense
        )
    else:
        profile = analyze_potential(spec.V, spec.x_min, x_max, spec.mass, e_max_mesh)
    real_lengths, order = (
        optimal_real_mesh(profile, min_len=spec.min_len, max_len=spec.max_len)
        if phase_coeff is None
        else optimal_real_mesh(
            profile, phase_coeff=phase_coeff, min_len=spec.min_len, max_len=spec.max_len
        )
    )

    R0 = spec.x_min + sum(real_lengths)
    angle = max_stable_angle(spec.V, R0, _ANGLE_PROBE_TAIL_EXTENT)
    tail_lengths = tune_ecs_tail(spec.channel_k, R0, angle=angle, order=order)

    elements = [ElementSpec(h) for h in real_lengths] + [
        ElementSpec(h, angle) for h in tail_lengths
    ]
    grid_spec = GridSpec(quadrature=order, elements=elements, x_min=spec.x_min)
    return FemDvrEcsGrid(grid_spec)
