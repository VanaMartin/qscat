"""Interaction-region localisation for the FEM-DVR-ECS discretisation
tuner (the resonance-aware nuclear mesh follow-on).

`interaction_region` finds the R-window over which the electron-molecule
interaction `V_int(r, R)` is still TRANSITIONING between its low-R and
high-R regimes -- the window a later, denser resonance-curve sampling step
confines itself to, freezing at a single asymptotic value outside it (see
`docs/superpowers/specs/2026-07-28-resonance-aware-mesh-design.md`: "only
in `[R_lo, R_hi]` does `V_d` differ from the asymptote"). It uses only
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

from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    from qscat.model import ResonanceModel

__all__ = ["interaction_region"]

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
