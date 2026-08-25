"""The one-shot a-priori grid assembler: `analyze` -> `mesh` -> `ecs`, wired
into a complete `FemDvrEcsGrid`.

`propose_grid` is the a-priori HALF of the hybrid tuner (see
`docs/physics/discretisation-tuning.md`): given a `ResonanceModel`, a
coordinate ("nuclear" or "electronic"), and a target energy range, it builds
the potential-adaptive grid directly from the potential curve -- no
eigensolve, no probing. The `discretisation-tuner` skill runs the
probe/refine loop ON TOP of this a-priori starting point; this module owns
only the a-priori half.

The MODEL ADAPTER (`_nuclear_adapter`/`_electronic_adapter`) is the only
model-aware code here: it picks the per-coordinate 1-D potential, mass, real
extent, and channel wavenumber from the model + energy range. Everything
downstream (`analyze_potential` -> `optimal_real_mesh` -> `max_stable_angle`
+ `tune_ecs_tail`) is the same model-independent pipeline the `analyze`,
`mesh`, and `ecs` modules validate on their own.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

import numpy as np
import numpy.typing as npt

from qscat.core.grids import electronic_grid
from qscat.dvr import ElementSpec, FemDvrEcsGrid, GridSpec

from .analyze import analyze_potential
from .ecs import max_stable_angle, tune_ecs_tail
from .mesh import optimal_real_mesh, order_for_wavenumber, refine_elements_in_window
from .resonance import interaction_region, resonance_curve

if TYPE_CHECKING:
    from qscat.model import ResonanceModel

__all__ = ["propose_grid"]

Coordinate = Literal["nuclear", "electronic"]

# Real-region cutoff defaults (bohr), in the ranges the eMoScat decks use
# (14-22 bohr nuclear, ~16-30 bohr electronic). Each is a fixed starting
# value, refined by the h/p mesh sweep and extended when an `incident`
# requires more room. The mesh's phase constant `C` is calibrated against the
# eMoScat decks (`validation.tuning.calibrate`); these extents are NOT -- they
# are one molecule-independent default, and it exceeds N2's (12 bohr) / NO's
# (9 bohr) committed nuclear real regions, which is why those two molecules'
# proposed grids cost more DVR points than their decks even at the calibrated
# C (see docs/physics/discretisation-tuning.md). Deriving x_max from the
# potential profile itself is a follow-on, not done here.
_NUCLEAR_X_MAX_DEFAULT = 18.0
_ELECTRONIC_X_MAX_DEFAULT = 20.0

# Element-length bounds (bohr) fed to `optimal_real_mesh`'s equidistribution
# sweep -- reasonable, but not independently calibrated (the calibration
# targets only the phase constant `C`). The nuclear floor is deliberately not
# tiny: a heavy reduced mass mu inflates the local
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

# Real-region extent (bohr) for the RESONANT (`channel="dissociation"`)
# nuclear path -- overrides `_NUCLEAR_X_MAX_DEFAULT` for that path only. The
# resonant mesh's ECS tail absorbs the outgoing dissociation wave, so the
# real region need only host it a few de Broglie wavelengths past the
# interaction region, not the blunt VE-path default (18 bohr) sized for a
# bound-state-like vibrational problem. 10.5 bohr ~ the eMoScat F2 DA deck's
# own real-region extent (~10.7 bohr) -- verified to keep F2 sigma_DA 2-D-
# converged (1.656 at E=0.03) while bringing the grid to deck-parity (n=1000
# vs the ~974-pt deck). A per-molecule derivation from the potential profile
# (where v0/V_d flatten + an outgoing-wave margin) is the documented follow-on
# -- 10.5 is an F2-anchored default; a longer-range channel may need more (an
# `incident`/extent override still applies). See docs/physics/
# discretisation-tuning.md.
_RESONANT_NUCLEAR_X_MAX_DEFAULT = 10.5

# Points-per-wavelength target for `order_for_wavenumber` when sizing the
# resonant path's DVR order against the exit-wave wavenumber `K_exit` -- see
# `qscat.tuning.mesh.order_for_wavenumber`'s docstring for what this counts.
_EXIT_TARGET_PPW = 6.0

# Target element length (bohr) the resonance-crossing window is super-
# refined to -- overrides `_NUCLEAR_MIN_LEN` LOCALLY inside that window
# (`qscat.tuning.mesh.refine_elements_in_window`); this is the sub-0.1-bohr
# spacing the eMoScat F2 deck hand-places around its own R~2.5-2.7 bohr
# interaction feature (see docs/physics/discretisation-tuning.md's 2-D
# spot-check finding).
_CROSSING_TARGET_LEN = 0.03

# Half-width (bohr) of the crossing super-refinement window, clamped from
# the Gamma-closing width computed around the crossing -- a floor so a
# vanishingly narrow closing width still refines a physically sane region,
# and a ceiling so a contaminated (frozen-plateau) closing-width estimate
# cannot blow the window out to the whole domain.
_CROSSING_DELTA_MIN = 0.15
# Ceiling on the crossing half-width. The Gamma-closing-width estimate is
# contaminated by the frozen-plateau artifact (Gamma held constant across the
# inner range where the pole finder broke down), which would otherwise push
# the window out to ~0.6 bohr and over-refine. The genuine crossing feature
# is narrow (the eMoScat deck super-refines only ~[2.5,2.7], a ~0.1-bohr
# half-width); cap here so the contaminated estimate cannot bloat the grid.
_CROSSING_DELTA_MAX = 0.18

# Fraction of the peak (region-restricted) Gamma a point must clear to count
# as part of the "Gamma-significant" closing region bracketing the crossing.
_GAMMA_SIGNIFICANT_FRAC = 0.1

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
    `qscat.model.IonicResonanceModel`). It is carried here so a future
    consumer (e.g. a channel-representation probe run on top of this grid)
    can read it off the same adapter output rather than re-deriving it from
    the model.
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
    x_max], outgoing dissociation wavenumber `K = sqrt(2*mu*e_max)` -- the
    heavy reduced mass makes this K large, which is what forces the fine
    exit-wave resolution the nuclear grids need.
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


def _outermost_crossing(R: npt.NDArray[np.float64], diff: npt.NDArray[np.float64]) -> float:
    """The largest-`R` (interpolated) zero of `diff(R)`, or `R[argmin(diff)]`
    if `diff` never changes sign -- the resonance-crossing localizer (`R* =`
    where `Re(V_d(R))` crosses `v0(R)`; see `_resonant_nuclear_mesh`).
    """
    sign = np.sign(diff)
    crossings = np.nonzero(np.diff(sign))[0]
    if crossings.size == 0:
        return float(R[int(np.argmin(diff))])
    i = int(crossings[-1])  # outermost: R ascending, so the LAST crossing
    x0, x1 = R[i], R[i + 1]
    f0, f1 = diff[i], diff[i + 1]
    if f1 == f0:
        return float(x0)
    t = -f0 / (f1 - f0)
    return float(x0 + t * (x1 - x0))


def _crossing_half_width(
    model: ResonanceModel, R: npt.NDArray[np.float64], Gamma: npt.NDArray[np.float64], r_star: float
) -> float:
    """Half-width `delta` of the `Gamma`-closing region bracketing `r_star`,
    clamped to `[_CROSSING_DELTA_MIN, _CROSSING_DELTA_MAX]`.

    Restricts to `R >= R_lo` (`interaction_region`) to exclude the walk's
    frozen inner plateau (see `resonance_curve`/`resonance_pole_walk`'s
    docstrings: on breakdown the LAST accepted `Gamma` freezes for all
    remaining -- smaller -- `R`, an artifact, not physics), then takes the
    `R`-span where `Gamma` clears `_GAMMA_SIGNIFICANT_FRAC` of its
    (region-restricted) peak.
    """
    R_lo, _R_hi = interaction_region(model)
    mask = R >= R_lo
    if not np.any(mask):
        return _CROSSING_DELTA_MIN
    gamma_peak = float(np.max(Gamma[mask]))
    if gamma_peak <= 0.0:
        return _CROSSING_DELTA_MIN
    sig_R = R[mask & (Gamma >= _GAMMA_SIGNIFICANT_FRAC * gamma_peak)]
    if sig_R.size == 0:
        return _CROSSING_DELTA_MIN
    delta = max(r_star - float(np.min(sig_R)), float(np.max(sig_R)) - r_star)
    return float(np.clip(delta, _CROSSING_DELTA_MIN, _CROSSING_DELTA_MAX))


def _resonant_nuclear_mesh(
    model: ResonanceModel,
    spec: _CoordinateSpec,
    x_max: float,
    e_max: float,
    e_max_mesh: float,
    elec_grids: tuple[FemDvrEcsGrid, FemDvrEcsGrid] | None,
    resonance_n_dense: int,
    phase_coeff: float | None,
) -> tuple[list[float], int, float]:
    """The `channel="dissociation"` real-region mesh: `(real_lengths, order,
    channel_k)`.

    Builds the adiabatic resonance curve `(R, V_d(R), Gamma(R))`
    (`resonance_curve`; the (expensive) two-angle pole match -- callers
    needing this cheap may inject small `elec_grids`, as `propose_grid`
    does), then:

    1. **Exit wave**: `K_exit = sqrt(2*mass*max(e_max - V_d_asym, e_max))`,
       `V_d_asym` the asymptotic `Re(V_d)` at the largest sampled `R` -- the
       fast outgoing dissociation wavenumber the ECS tail must absorb
       (`channel_k`, returned). The DVR `order` is sized to resolve THIS
       wave at the base `min_len` via `order_for_wavenumber`, not swept for
       point-count as the ve path does.
    2. **Base mesh**: `v0`-alone (`analyze_potential` + `optimal_real_mesh`
       restricted to the single chosen `order`).
    3. **Crossing**: `R* =` the outermost `Re(V_d) - v0` sign change
       (`_outermost_crossing`); `delta` the Gamma-closing half-width
       (`_crossing_half_width`). `[R* - delta, R* + delta]` is super-refined
       to `_CROSSING_TARGET_LEN` via `refine_elements_in_window` -- this
       LOCALLY overrides `min_len`, which is the point (see module
       docstring / the design note this supersedes: a floored `min_len` is
       exactly why the prior worst-case-`k`-merge design was inert).
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

    vd_asym = float(Vd_real[-1])  # R ascending -> the largest sampled R
    channel_k = math.sqrt(2.0 * spec.mass * max(e_max - vd_asym, e_max))
    order = order_for_wavenumber(channel_k, spec.min_len, target_ppw=_EXIT_TARGET_PPW)

    # Exit-region element floor: the COARSEST element that still resolves the
    # fast exit wave `channel_k` at `target_ppw` points per wavelength for the
    # chosen `order` (`order * lambda / ppw`, lambda = 2*pi/channel_k). Since
    # `order` is picked high enough to resolve the exit wave, the base
    # `spec.min_len` (0.15) OVER-resolves it -- floor at this coarser length
    # instead (the eMoScat deck's ~0.2-bohr exit elements at q=14), then let
    # the crossing window super-refine locally below it. This is the size
    # lever that keeps the resonant grid competitive with the hand deck.
    exit_min_len = max(spec.min_len, order * (2.0 * math.pi / channel_k) / _EXIT_TARGET_PPW)

    profile_v0 = analyze_potential(spec.V, spec.x_min, x_max, spec.mass, e_max_mesh)
    real_lengths, order = (
        optimal_real_mesh(profile_v0, orders=(order,), min_len=exit_min_len, max_len=spec.max_len)
        if phase_coeff is None
        else optimal_real_mesh(
            profile_v0,
            orders=(order,),
            phase_coeff=phase_coeff,
            min_len=exit_min_len,
            max_len=spec.max_len,
        )
    )

    v0_at_R = np.real(np.asarray(spec.V(R), dtype=np.complex128))
    r_star = _outermost_crossing(R, Vd_real - v0_at_R)
    delta = _crossing_half_width(model, R, Gamma, r_star)
    real_lengths = refine_elements_in_window(
        real_lengths, spec.x_min, r_star - delta, r_star + delta, _CROSSING_TARGET_LEN
    )

    return real_lengths, order, channel_k


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

    `incident` (`qscat.tuning.incident.IncidentSpec`) is accepted here as
    BOTH an extent floor AND a resolution floor:

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
      nuclear path (see `_resonant_nuclear_mesh`). Builds the adiabatic
      resonance curve `(R, V_d(R), Gamma(R))`
      (`qscat.tuning.resonance.resonance_curve`, via a two-angle ECS pole
      match -- `elec_grids`, if given, overrides the default electronic
      grids used for that match, and `resonance_n_dense` overrides its
      dense-sampling point count; both exist so tests can inject
      small/cheap grids). A REDUCED real-region extent
      (`_RESONANT_NUCLEAR_X_MAX_DEFAULT`, not the VE path's
      `_NUCLEAR_X_MAX_DEFAULT`) is used, since the ECS tail absorbs the
      outgoing wave. The DVR order is sized (`order_for_wavenumber`) to
      resolve the fast dissociation EXIT wave `K_exit` at the base
      `min_len`, and that same `K_exit` (not `spec.channel_k`) drives the
      ECS tail. The narrow resonance CROSSING `R*` (the outermost
      `Re(V_d) - v0` sign change) is then LOCALLY super-refined
      (`refine_elements_in_window`, overriding `min_len` only inside a
      Gamma-closing-width window around `R*`). The local override is the
      point: merging the resonance's wavenumber into the global profile
      instead would be inert, because the merged profile's finer elements
      get floored right back up by the shared global `min_len`.
      `"dissociation"` with `coordinate="electronic"` raises `ValueError` --
      the resonance-aware path is nuclear-only.
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

    # The resonant path uses a REDUCED default extent (the ECS tail absorbs
    # the exit wave, so the real region need only host it a few wavelengths
    # past the interaction region) -- not `spec.x_max` (the VE-path default,
    # sized for a bound-state-like vibrational problem). `incident`, if
    # given, may still widen either default below.
    x_max = _RESONANT_NUCLEAR_X_MAX_DEFAULT if channel == "dissociation" else spec.x_max
    e_max_mesh = e_max
    if incident is not None:
        required_extent = getattr(incident, "required_extent", lambda: 0.0)()
        x_max = max(x_max, float(required_extent))
        incident_energy = getattr(incident, "incident_energy", lambda: 0.0)()
        e_max_mesh = max(e_max_mesh, float(incident_energy))

    if channel == "dissociation":
        real_lengths, order, channel_k = _resonant_nuclear_mesh(
            model, spec, x_max, e_max, e_max_mesh, elec_grids, resonance_n_dense, phase_coeff
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
        channel_k = spec.channel_k

    R0 = spec.x_min + sum(real_lengths)
    angle = max_stable_angle(spec.V, R0, _ANGLE_PROBE_TAIL_EXTENT)
    tail_lengths = tune_ecs_tail(channel_k, R0, angle=angle, order=order)

    elements = [ElementSpec(h) for h in real_lengths] + [
        ElementSpec(h, angle) for h in tail_lengths
    ]
    grid_spec = GridSpec(quadrature=order, elements=elements, x_min=spec.x_min)
    return FemDvrEcsGrid(grid_spec)
