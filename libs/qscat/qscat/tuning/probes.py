"""Decoupled 1-D convergence probes: the empirical validators that tell the
tuner whether a candidate FEM-DVR-ECS grid resolves the physics it is meant
to carry.

Each probe compares a quantity computed on the CANDIDATE grid against the
same quantity on a `refine`-d (one h-refinement step) grid, or -- for the
cheapest probe -- against a fine-quadrature reference, and reports whether
the change is within `rtol`. All three probes return a `ProbeResult`
(`value`, `converged`, `cost`, `detail`), so the tuner can treat them
uniformly regardless of what each one actually computes.

- `refine` -- ONE h-refinement step: every REAL `ElementSpec` is split into
  two elements of half its length (so the real-region span is unchanged, its
  resolution doubles); the ECS tail, quadrature order, and `x_min` are left
  untouched. `FemDvrEcsGrid` already carries the `GridSpec` it was built from
  (`grid.spec`), so `refine` needs only the grid -- no separate spec plumbing.

- `probe_nuclear` -- the `n_vib` lowest vibrational eigenvalues
  (`qscat.core.vibrational.vibrational_states`) on the grid vs. its refinement;
  `converged` iff every eigenvalue's relative shift is `< rtol`. Needs an
  eigensolve -- moderate cost.

- `probe_electronic` -- the lowest anion bound electronic-state energy at a
  fixed internuclear distance `R` (`qscat.core.dissociation.
  anion_electronic_states`, `n_states=1`) on the grid vs. its refinement.
  This is a CHEAPER PROXY for the full two-angle resonance-pole convergence
  check (`qscat.ecs.find_resonance_pole` needs two ECS-angle spectra, which
  is a follow-on): the bound-state energy is angle-independent on a
  converged grid exactly like the nuclear vibrational levels are (see
  `vibrational_states`'s module docstring), so it is a legitimate, much
  cheaper stand-in for electronic-grid convergence. `window` is accepted for
  interface parity with a future two-angle pole probe but is UNUSED by this
  proxy (there is no pole search here) -- kept so callers/tuner code can be
  written once against a stable signature and swapped later without a
  breaking change.

- `probe_channel_representation` -- THE cheapest and most diagnostic probe:
  no eigensolve at all. It checks whether the free (`charge=0`, mass-`mass`
  `riccati_bessel_en_mass`) or Coulomb (`coulomb_f_en`) scattering channel
  function at wavenumber `k` and partial wave `l` is REPRESENTABLE on the
  grid's real region, by comparing the DVR quadrature estimate of `integral
  |F(r)|^2 dr` over `[x_min, R0)` against a reference computed on a uniform
  grid fine enough to resolve the fastest oscillation present (`>= 40`
  samples per de Broglie wavelength `2*pi/k`) over that SAME `[x_min, R0]`
  span, via Simpson's rule. A grid whose element lengths are large compared
  to `1/k` aliases this integral badly -- this is exactly the failure mode
  of a K~58 wave on ~1.0-bohr elements (the coarse-grid dissociative-
  attachment failure this probe is designed to catch), and it costs nothing
  beyond evaluating a special function. `mass` matters for a non-electron
  channel (e.g. a mass-`mu` nuclear dissociation wave); it is honored in
  BOTH the neutral (`charge=0`) and Coulomb branches, and in the reference.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, NamedTuple

import numpy as np
from scipy.integrate import simpson

from qscat.core.dissociation import anion_electronic_states
from qscat.core.vibrational import vibrational_states
from qscat.dvr import ElementSpec, FemDvrEcsGrid, GridSpec
from qscat.special import coulomb_f_en, riccati_bessel_en_mass

if TYPE_CHECKING:
    from qscat.model import ResonanceModel

__all__ = [
    "ProbeResult",
    "probe_channel_representation",
    "probe_electronic",
    "probe_nuclear",
    "refine",
]


class ProbeResult(NamedTuple):
    """A convergence probe's verdict: the computed `value`, whether it is
    `converged` (change under one refinement / vs. a fine reference is
    `< rtol`), the discretisation `cost` (grid point count) it was computed
    at, and a `detail` dict with the raw numbers behind the verdict (refined
    value, relative error, reference value, ...) for diagnostics/reporting.
    """

    value: Any
    converged: bool
    cost: int
    detail: dict[str, Any]


# Number of uniform samples per de Broglie wavelength (2*pi/k) used to build
# the fine-reference integral in `probe_channel_representation` -- generous
# enough that the reference itself is never the accuracy bottleneck.
_SAMPLES_PER_WAVELENGTH = 40
_MIN_FINE_SAMPLES = 4001

# Below this magnitude a "relative" error would divide by ~0; fall back to
# an absolute comparison against this floor instead (mirrors the pattern
# `qscat.core.vibrational`/`dissociation` use for near-zero bound-state
# guards).
_REL_FLOOR = 1e-12


def refine(grid: FemDvrEcsGrid) -> FemDvrEcsGrid:
    """One h-refinement step: split every REAL element into two half-length
    elements; leave the ECS tail, quadrature order, and `x_min` untouched.

    `grid.spec` is the `GridSpec` the grid was built from (`FemDvrEcsGrid`
    stores it directly), so no separate spec-plumbing is needed -- this is
    the "SIMPLEST" option the design brief calls out: reconstruct from the
    grid's own stored spec rather than threading a `GridSpec` through every
    probe signature.
    """
    spec = grid.spec
    new_elements: list[ElementSpec] = []
    for el in spec.elements:
        if el.angle_deg == 0.0:
            half = el.length / 2.0
            new_elements.append(ElementSpec(half, 0.0))
            new_elements.append(ElementSpec(half, 0.0))
        else:
            new_elements.append(ElementSpec(el.length, el.angle_deg))
    new_spec = GridSpec(quadrature=spec.quadrature, elements=new_elements, x_min=spec.x_min)
    return FemDvrEcsGrid(new_spec)


def probe_nuclear(
    model: ResonanceModel,
    nuclear_grid: FemDvrEcsGrid,
    n_vib: int,
    *,
    rtol: float = 1e-3,
) -> ProbeResult:
    """Convergence of the `n_vib` lowest vibrational eigenvalues under one
    h-refinement of `nuclear_grid`.

    `converged = max_v |eps_v(refined) - eps_v| / |eps_v| < rtol` (an
    absolute floor is used for any `eps_v` too close to zero). `cost` is the
    (unrefined) grid's point count `nuclear_grid.n`.
    """
    eps, _ = vibrational_states(nuclear_grid, model.mu, n_vib, model.v0)
    refined = refine(nuclear_grid)
    eps_ref, _ = vibrational_states(refined, model.mu, n_vib, model.v0)

    denom = np.where(np.abs(eps) > _REL_FLOOR, np.abs(eps), 1.0)
    rel = np.abs(eps_ref - eps) / denom
    max_rel = float(np.max(rel))
    return ProbeResult(
        value=eps,
        converged=bool(max_rel < rtol),
        cost=nuclear_grid.n,
        detail={
            "eps_refined": eps_ref,
            "max_rel_delta": max_rel,
            "refined_cost": refined.n,
        },
    )


def probe_electronic(
    model: ResonanceModel,
    elec_grid: FemDvrEcsGrid,
    R: float,
    *,
    window: tuple[float, float] | None,
    rtol: float = 1e-3,
) -> ProbeResult:
    """Convergence of the lowest anion bound electronic-state energy at
    internuclear distance `R`, under one h-refinement of `elec_grid`.

    See the module docstring for why this bound-state energy (rather than a
    full two-angle resonance-pole match) is used as the electronic-grid
    convergence proxy; `window` is accepted for signature parity with a
    future pole-based probe but is UNUSED here.
    """
    eps_e, _ = anion_electronic_states(elec_grid, model, R, n_states=1)
    refined = refine(elec_grid)
    eps_e_ref, _ = anion_electronic_states(refined, model, R, n_states=1)

    e0 = float(eps_e[0])
    e0_ref = float(eps_e_ref[0])
    denom = abs(e0) if abs(e0) > _REL_FLOOR else 1.0
    rel = abs(e0_ref - e0) / denom
    return ProbeResult(
        value=e0,
        converged=bool(rel < rtol),
        cost=elec_grid.n,
        detail={
            "e0_refined": e0_ref,
            "rel_delta": rel,
            "refined_cost": refined.n,
        },
    )


def _fine_reference_norm(
    x_min: float, R0: float, k: float, l: int, *, charge: int, mass: float
) -> float:
    """`integral_{x_min}^R0 |F(r)|^2 dr`, via Simpson's rule on a uniform grid
    fine enough (`_SAMPLES_PER_WAVELENGTH` samples per `2*pi/k`) to resolve
    the fastest oscillation `F` exhibits over `[x_min, R0]` -- the same real
    domain `probe_channel_representation`'s DVR-quadrature sum covers.
    """
    wavelength = 2.0 * np.pi / k
    span = R0 - x_min
    n_pts = max(_MIN_FINE_SAMPLES, int(_SAMPLES_PER_WAVELENGTH * span / wavelength) + 1)
    r_fine = np.linspace(x_min, R0, n_pts)
    f_fine = (
        riccati_bessel_en_mass(r_fine, k, l, mass)
        if charge == 0
        else coulomb_f_en(r_fine, k, float(charge), mass, l)
    )
    return float(simpson(np.abs(f_fine) ** 2, x=r_fine))


def probe_channel_representation(
    grid: FemDvrEcsGrid,
    k: float,
    l: int,
    *,
    charge: int = 0,
    mass: float = 1.0,
    rtol: float = 1e-3,
) -> ProbeResult:
    """Is the energy-normalized channel function `F_{k,l}` REPRESENTABLE on
    `grid`'s real region?

    Compares the grid's own composite Gauss-Lobatto quadrature estimate of
    `integral |F|^2 dr` over the real region `[x_min, R0)` -- `sum(Re(w_j) *
    F(r_j)^2)` at the grid's own real-region nodes/weights -- against a
    fine-uniform-grid Simpson reference (`_fine_reference_norm`, integrated
    over that SAME `[x_min, R0]` span -- `grid.spec.x_min`, not hardcoded to
    0.0; every grid built by `qscat.core.grids` uses `x_min=0.0` today, but
    `GridSpec.x_min` is a supported nonzero field). The node exactly at `R0`
    (shared with the first ECS-tail element) is EXCLUDED: its global
    bridge-summed weight mixes in that neighbor's complex Jacobian (see the
    inline comment at the mask below), which would corrupt a real-region-only
    quadrature at any nonzero ECS angle. `charge=0` uses `riccati_bessel_en_
    mass(..., mass)` (mass-`mass` free wave; reduces to the mass-1 electron
    channel at the default `mass=1.0`) and `charge != 0` uses `coulomb_f_en`
    (which already takes `mass`) -- both the DVR sum and the reference honor
    `mass` identically. No eigensolve is needed, so this is by far the
    cheapest of the three probes -- and the most diagnostic: a grid whose
    element lengths are not small compared to the channel's wavelength
    `2*pi/k` cannot represent `F` and this probe catches it directly (this is
    what would have caught the K~58-on-1.0-bohr-elements coarse-grid DA
    failure).
    """
    # Strict `<`, not `<=`: the node exactly at R0 is the bridge point shared
    # with the first ECS-tail element, and its GLOBAL bridge-summed weight
    # (grid.weights) mixes in that neighbor's complex Jacobian (the "classic
    # assembly trap", qscat.dvr.kinetic's term for using the global weight
    # where a purely real, real-region-local one is needed) -- e.g. its real
    # part alone still carries a spurious `Re(hz_tail) * w[0]` term at a
    # nonzero ECS angle. Dropping that single shared node is harmless: on a
    # PURE real grid (no ECS tail) R0 == x_max, which is never itself a
    # retained node (the Dirichlet endpoint is already dropped by
    # `FemDvrEcsGrid`), so `<` and `<=` coincide there and this never excises
    # a genuine interior real-region point.
    real_mask = grid.real_points < grid.R0
    r_real = grid.real_points[real_mask]
    w_real = grid.weights[real_mask].real

    f_vals = (
        riccati_bessel_en_mass(r_real, k, l, mass)
        if charge == 0
        else coulomb_f_en(r_real, k, float(charge), mass, l)
    )
    f2 = np.abs(np.asarray(f_vals, dtype=np.complex128)) ** 2
    quad_val = float(np.sum(w_real * f2))

    ref_val = _fine_reference_norm(
        float(grid.spec.x_min), float(grid.R0), k, l, charge=charge, mass=mass
    )
    denom = abs(ref_val) if abs(ref_val) > _REL_FLOOR else 1.0
    rel_err = abs(quad_val - ref_val) / denom
    return ProbeResult(
        value=quad_val,
        converged=bool(rel_err < rtol),
        cost=int(np.count_nonzero(real_mask)),
        detail={
            "reference": ref_val,
            "rel_error": rel_err,
        },
    )
