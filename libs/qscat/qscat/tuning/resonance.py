"""Interaction-region localisation for the FEM-DVR-ECS discretisation
tuner (the resonance-aware nuclear mesh follow-on).

`interaction_region` finds the R-window over which the electron-molecule
interaction `V_int(r, R)` is still TRANSITIONING between its low-R and
high-R regimes -- the window a later, denser resonance-curve sampling step
confines itself to, freezing at a single asymptotic value outside it --
only inside `[R_lo, R_hi]` does `V_d` differ from its asymptote (see
docs/physics/discretisation-tuning.md). It uses only
`model.v_int` (the `ResonanceModel` protocol), so it works for any
registered model (diatomic or ionic), not just the shared
Morse+sigmoid+Gaussian-in-r form.

**Why a percentile-of-range crossing, not a raw magnitude threshold**: for
every registered `DiatomicResonanceModel` (N2/NO/F2), `lambda(R)` (and
hence `max_r |v_int(r, R)|`) is a *sigmoid* in `R` that saturates to a
substantial, generally NON-ZERO plateau on both the low-R and high-R sides
(verified numerically: F2's plateaus are ~9.7 and ~18.3, both large
compared to their difference; N2's low-R plateau is even the far larger of
the two once you account for its sign flip). A literal "`s(R) >=
frac * s.max()`" threshold (`frac` a small fraction like the 0.02 default)
is then satisfied at BOTH ends of any scanned domain, collapsing the
bracket to the scan's own edges rather than localising anything. What
actually varies over a bounded region is where `s(R)` is still ascending
between its own floor and ceiling -- so `s.min()`/`s.max()` over the scan
set the reference RANGE, and `R_lo`/`R_hi` bracket the middle
`(1 - 2*frac)` fraction of that rise (excluding the outer `frac` at each
end), which is exactly the transition zone the downstream sampler needs to
resolve densely.
"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from qscat.core.dissociation import anion_electronic_states
from qscat.core.lcp import resonance_pole_walk

if TYPE_CHECKING:
    from qscat.dvr import FemDvrEcsGrid
    from qscat.model import ResonanceModel

__all__ = ["interaction_region", "resonance_curve_arrays"]

FloatArray = NDArray[np.float64]

# Default electronic probe set: a modest grid spanning the r-range where
# every registered model's Gaussian-in-r interaction has non-negligible
# support (alpha_c in [0.4, 3.0] for N2/NO/F2, decaying well inside 15 bohr).
# The r-shape factors out of the R-profile entirely (v_int(r, R) =
# -lambda(R) * f(r) for the shared diatomic form), so the exact r_probe
# choice only needs to be consistent, not tuned per model.
_DEFAULT_R_PROBE_MIN = 0.1
_DEFAULT_R_PROBE_MAX = 15.0
_DEFAULT_R_PROBE_N = 200


def interaction_region(
    model: ResonanceModel,
    *,
    r_probe: FloatArray | None = None,
    R_max: float = 8.0,
    frac: float = 0.02,
    n: int = 400,
) -> tuple[float, float]:
    """The `(R_lo, R_hi)` bracket over which `V_int(r, R)` is still
    transitioning between its low-R and high-R regimes.

    Defines `s(R) = max_r |Re(model.v_int(r_probe, R))|` over a fixed
    electronic probe set `r_probe`, scans `R` on
    `linspace(1e-3, R_max, n)`, and returns the first/last `R` whose
    `s(R)` falls within the middle `(1 - 2*frac)` fraction of `s`'s own
    `[min, max]` range on that scan -- i.e. `R_lo` is the first `R` with
    `s(R) - s.min() >= frac * (s.max() - s.min())` and `R_hi` is the last
    `R` with `s(R) - s.min() <= (1 - frac) * (s.max() - s.min())`. See the
    module docstring for why this is the range that localises (a raw
    threshold against `s.max()` alone does not, for these sigmoid-shaped
    interaction profiles).

    `R` starts at `1e-3` rather than 0 to stay clear of any `R=0`
    singularity in a model's `v_int` (e.g. a Coulomb-tail ionic model).
    """
    if r_probe is None:
        r_probe = np.linspace(_DEFAULT_R_PROBE_MIN, _DEFAULT_R_PROBE_MAX, _DEFAULT_R_PROBE_N)

    R_scan = np.linspace(1e-3, R_max, n)
    v = model.v_int(r_probe[:, None], R_scan[None, :])
    s = np.max(np.abs(np.real(v)), axis=0)

    s_min, s_max = float(np.min(s)), float(np.max(s))
    span = s_max - s_min
    if span <= 0.0:
        # s(R) is (numerically) constant over the scan -- no transition to
        # localise; report the full scanned domain.
        return float(R_scan[0]), float(R_scan[-1])

    rel = (s - s_min) / span
    lo_above = np.nonzero(rel >= frac)[0]
    hi_below = np.nonzero(rel <= 1.0 - frac)[0]
    R_lo = R_scan[int(lo_above[0])] if lo_above.size else float(R_scan[0])
    R_hi = R_scan[int(hi_below[-1])] if hi_below.size else float(R_scan[-1])
    return float(R_lo), float(R_hi)


def resonance_curve_arrays(
    model: ResonanceModel,
    elec_grid_a: FemDvrEcsGrid,
    elec_grid_b: FemDvrEcsGrid,
    *,
    R_max: float = 22.0,
    n_dense: int = 25,
    region: tuple[float, float] | None = None,
) -> tuple[FloatArray, FloatArray, FloatArray]:
    """Efficient adiabatic resonance curve `(R, V_d(R), Gamma(R))`.

    Samples the pole continuation (`qscat.core.lcp.resonance_pole_walk`)
    DENSELY (`n_dense` points) only inside the interaction window
    `region = region or interaction_region(model)`, plus a few inner points
    below `R_lo` and a SINGLE far point at `R_max` standing in for the
    saturated asymptote -- the efficiency constraint that keeps this to
    `O(n_dense)` electronic solves rather than a scan over hundreds of `R`
    points. The walk runs over the samples DESCENDING (outer -> inner) and
    its freeze (see `resonance_pole_walk`) correctly holds `V_d` constant
    across the big `R_max -> R_hi` gap, since `Gamma` has already saturated
    to ~0 out there.

    Renamed from `resonance_curve` (2026-08-25 API surface pass) to end the
    collision with `qscat.core.bo.resonance_curve`, which shares the same
    underlying pole walk but returns an `ElectronicCurves` carrying the
    eigenVECTORS, for building Born-Oppenheimer basis states. This one
    returns plain `(R, V_d, Gamma)` arrays -- the name says so -- and exists
    only to size a grid, so it discards the states and samples as sparsely
    as it can.

    Seeded from the bound anion state at `R_max`
    (`qscat.core.dissociation.anion_electronic_states`). Returns arrays
    ascending in `R`; `V_d(R) = Re(model.v0(R)) + shift`.
    """
    R_lo, R_hi = region or interaction_region(model)

    inner = np.linspace(max(1e-3, R_lo - 0.6), R_lo, 3)[:-1]
    dense = np.linspace(R_lo, R_hi, n_dense)
    R_samples = np.unique(np.concatenate([[R_max], dense, inner]))
    R_descending = R_samples[::-1]

    eps_e, _ = anion_electronic_states(elec_grid_a, model, R_max, 1)
    seed_window = (eps_e[0] - 0.05, eps_e[0] + 0.05, -0.05, 0.05)

    shift, gamma = resonance_pole_walk(
        model,
        R_descending,
        elec_grid_a,
        elec_grid_b,
        seed_window,
    )

    order = np.argsort(R_descending)
    R = R_descending[order]
    Vd = np.real(model.v0(R)) + shift[order]
    Gamma = gamma[order]
    return R, Vd, Gamma


# --- Deprecated aliases (2026-08-25 API surface pass) ------------------------
# One release cycle per ADR 0004, then delete this block. Not in `__all__`:
# the public surface is the new name; the alias only keeps old imports alive.

_DEPRECATED = {"resonance_curve": "resonance_curve_arrays"}


def __getattr__(name: str) -> object:
    if name in _DEPRECATED:
        new = _DEPRECATED[name]
        warnings.warn(
            f"{__name__}.{name} was renamed to {new} in the 2026-08-25 API "
            "surface pass; the old name is a deprecated alias for one release "
            "cycle (docs/adr/0004-public-api-stability-policy.md)",
            DeprecationWarning,
            stacklevel=2,
        )
        return globals()[new]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
