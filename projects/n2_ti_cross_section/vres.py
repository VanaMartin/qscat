"""V_d(R) / Gamma(R) recomputed per nuclear-R via the sub-project #2
electronic resonance pole finder (sub-project #3, Task 2).

the eMoScat TI extraction's "Efficiency note":
V_d(R)/Gamma(R) is recomputed at *every* nuclear grid point R by an
independent electronic-structure pole search (the user's choice, not
interpolation from a coarse R-scan). At each R: `E_res(R) = Re(E_pole)`,
`Gamma(R) = max(0, -2*Im(E_pole))`, `V_d(R) = v0(R) + E_res(R)`, where
`E_pole` is the two-angle-matched electronic resonance pole from
`projects/n2_resonance/pole.find_pole(R, grid_a, grid_b, window)` (imported
package-absolute, the pattern already used by `vibrational.py`). The two
electronic ECS grids (35 deg / 44 deg) are built once and reused across all
R -- only `V_eff_el`'s R argument changes, per the efficiency note.

Continuation over the REAL nuclear grid points:
  Sort the real grid points by R, seed a window near the R0-equilibrium
  pole region (`_WINDOW0`, matching `projects/n2_resonance/pole.py`'s
  `resonance_curve`), and walk outward in both directions from the point
  nearest `_SEED_R`, each step recentering the window on the previous step's
  matched pole (+/- `_RE_HALF_WIDTH`/`_IM_HALF_WIDTH`). This is the same
  continuation strategy sub-project #2's `pole.resonance_curve` already uses
  and validates (`test_curve.py`, R in [1.6, 3.0]); this module extends the
  walk to the FULL real nuclear grid (R roughly 0.01 to 12 bohr).

Breakdown region (small R, deep repulsive wall):
  Walking the continuation down to very small R (below roughly 0.9-1.2
  bohr), the two-angle match degrades and eventually fails outright
  (`find_pole` raises `ValueError`: the window empties in one of the two
  spectra; the per-step call also defensively catches
  `np.linalg.LinAlgError`, in case `qscat.dvr.eigen`'s diagonalization ever
  fails to converge at a pathological R -- treated identically to a failed
  match). This is NOT a solver bug -- `v0(R)` there is already several
  Hartree above the neutral dissociation limit (e.g. v0(1.0) ~= 3.0 Ha,
  v0(0.5) ~= 16.3 Ha; the interaction strength `lambda(R)` itself swings
  from ~6.2 at large R to negative at very small R), so the model's shape
  resonance stops being a well-defined, angle-stable pole there. Physically
  this region is irrelevant: every neutral vibrational level sits between
  -D_0 = -0.751 Ha and 0 Ha, so the nuclear wavefunctions have utterly
  negligible amplitude wherever v0(R) is several Hartree above zero -- but
  Task 2 still requires Vd/Gamma to be *defined* (finite) at every nuclear
  grid point. Policy (documented, matches the complex-tail treatment below):
  a step is accepted only if `find_pole` succeeds AND its angle-matching
  residual is below `_RESID_TOL_HA` (1e-3 Ha, the same threshold
  sub-project #2's own `test_pole.py::test_V2_pole_is_stable` uses to call a
  pole "angle-stable"). The first time a step in a given walk direction is
  rejected (raises, or residual too large), the walk stops advancing in
  that direction; every remaining, more-extreme real grid point in that
  direction reuses the last-accepted `E_res`/`Gamma` (with `V_d` still
  evaluated at that point's own, correct, rapidly-rising `v0(R)` -- only the
  resonance-shift part is frozen, not the whole curve). On the validated
  default nuclear grid (`nuclear_grid.n2_nuclear_grid()`), this affects only
  the innermost ~22-30 of 299 real grid points (R below ~1.1-1.2 bohr).

Complex tail (R > R0_nuclear_ECS_pivot = 12 bohr):
  `v_eff_el(r, R)` is analytic in R, so `find_pole` *can* in principle be
  called at complex nuclear R too, but that stacks the electronic ECS
  continuation on top of a second, independent nuclear-ECS continuation --
  "delicate" per the task brief, and unnecessary here because by R ~ 12 bohr
  the resonance has long since closed into a real, angle-independent bound
  anion curve (Gamma is already ~1e-14 Ha at the largest real grid point,
  see `test_gamma_closes_beyond_crossing_and_on_complex_tail`). We therefore
  use the recommended simpler, robust path: on the complex tail, clamp
  `Gamma = 0` and continue `E_res` as the constant asymptote taken from the
  largest real grid point (`E_res_asymptote`), while still evaluating
  `v0(R)` at the ECS-rotated complex `R` itself (matching how
  `vibrational.py` evaluates `v0(grid.points)` directly on the same nuclear
  grid -- the local potential is evaluated at the true, possibly complex,
  ECS coordinate; only the resonance-shift correction is frozen at its
  converged real-R asymptote). This keeps `Gamma >= 0` and `V_d` finite and
  smooth everywhere, and is exact at the real/complex-tail boundary by
  construction (the boundary point itself supplies `E_res_asymptote`).
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from qscat.dvr import FemDvrEcsGrid

# Reuse the sub-project #2 electronic pole finder.
from projects.n2_resonance.grid_n2 import n2_electronic_grid
from projects.n2_resonance.pole import find_pole
from projects.n2_resonance.potential import v0

__all__ = ["vres_on_grid"]

_ANGLE_A_DEG = 35.0
_ANGLE_B_DEG = 44.0

# Continuation seed/window, matching `projects/n2_resonance/pole.py`'s
# `resonance_curve` defaults (validated there against the R0 equilibrium
# pole and the R in [1.6, 3.0] smooth-curve test).
_SEED_R = 2.01943
_WINDOW0: tuple[float, float, float, float] = (0.04, 0.16, -0.05, 0.0)
_RE_HALF_WIDTH = 0.05
_IM_HALF_WIDTH = 0.05

# Angle-stability threshold: same bar sub-project #2's own
# `test_pole.py::test_V2_pole_is_stable` uses to call a matched pole
# genuinely angle-independent (residual << Gamma) rather than a coincidental
# near-miss between unrelated eigenvalues.
_RESID_TOL_HA = 1e-3


def vres_on_grid(
    grid: FemDvrEcsGrid,
) -> tuple[npt.NDArray[np.complex128], npt.NDArray[np.float64]]:
    """Recompute `(V_d(R), Gamma(R))` at every nuclear grid point.

    `V_d(R) = v0(R) + E_res(R)`, `Gamma(R) = max(0, -2*Im(E_pole(R)))`,
    where `E_pole(R)` is the two-angle-matched electronic resonance pole
    (`projects/n2_resonance/pole.find_pole`) evaluated at the (real or
    complex) nuclear grid point `R = grid.points[i]`. See the module
    docstring for the continuation walk over real grid points, the
    small-R breakdown-region fallback, and the complex-tail treatment.

    Returns `(Vd, Gamma)`, each shape `(grid.n,)`: `Vd` is complex128
    (matching how `v0(grid.points)` is evaluated at the true, possibly
    ECS-rotated complex coordinate elsewhere in this project, e.g.
    `vibrational.py`); `Gamma` is real float64, `>= 0` everywhere.
    """
    points = grid.points
    n = grid.n

    real_mask = points.imag == 0.0
    real_idx = np.flatnonzero(real_mask)
    tail_idx = np.flatnonzero(~real_mask)

    real_R_unsorted = points[real_idx].real
    order = np.argsort(real_R_unsorted)
    sorted_local_idx = real_idx[order]  # global grid index, ascending in R
    sorted_R = real_R_unsorted[order]
    m = sorted_R.size

    E_res = np.empty(m, dtype=np.float64)
    Gamma_sorted = np.empty(m, dtype=np.float64)

    ga = n2_electronic_grid(_ANGLE_A_DEG)
    gb = n2_electronic_grid(_ANGLE_B_DEG)

    seed_pos = int(np.argmin(np.abs(sorted_R - _SEED_R)))

    def _walk(positions: range, seed_from: int | None) -> None:
        """Walk `positions` (a contiguous run away from the seed), tracking
        the pole by window continuation. `seed_from`, if given, is the
        position whose already-computed (E_res, Gamma) anchors the walk's
        initial window and breakdown fallback (used for the decreasing-R
        walk, which continues from the increasing-R walk's seed point).
        """
        window = _WINDOW0
        last_good: int | None = seed_from
        broken = False
        for pos in positions:
            if broken:
                assert last_good is not None
                E_res[pos] = E_res[last_good]
                Gamma_sorted[pos] = Gamma_sorted[last_good]
                continue

            R = float(sorted_R[pos])
            try:
                E_pole, residual = find_pole(R, ga, gb, window)
            except (ValueError, np.linalg.LinAlgError):
                # ValueError: the two-angle window match found no candidate pole
                # (see module docstring's "Breakdown region"). LinAlgError:
                # defensive -- `find_pole` diagonalizes H_el(R) via `qscat.dvr.eigen`,
                # which could in principle raise on a pathological (non-converging)
                # eigendecomposition at some R; treated the same as a failed match
                # (per the Task-2 review Minor).
                residual = np.inf
            else:
                if residual < _RESID_TOL_HA:
                    E_res[pos] = E_pole.real
                    Gamma_sorted[pos] = max(0.0, -2.0 * E_pole.imag)
                    window = (
                        E_pole.real - _RE_HALF_WIDTH,
                        E_pole.real + _RE_HALF_WIDTH,
                        E_pole.imag - _IM_HALF_WIDTH,
                        E_pole.imag + _IM_HALF_WIDTH,
                    )
                    last_good = pos
                    continue

            # Rejected step (raised, or residual >= tolerance): the
            # two-angle match is no longer trustworthy here (see module
            # docstring's "Breakdown region"). Freeze at the last accepted
            # point and stop advancing further in this direction.
            if last_good is None:
                raise RuntimeError(
                    f"vres_on_grid: pole finder failed at the very first "
                    f"continuation step (R={R}, window={window}, "
                    f"residual={residual}) -- no prior accepted point to "
                    "fall back on."
                )
            broken = True
            E_res[pos] = E_res[last_good]
            Gamma_sorted[pos] = Gamma_sorted[last_good]

    _walk(range(seed_pos, m), seed_from=None)
    _walk(range(seed_pos - 1, -1, -1), seed_from=seed_pos)

    Vd = np.empty(n, dtype=np.complex128)
    Gamma = np.zeros(n, dtype=np.float64)

    Vd[sorted_local_idx] = v0(sorted_R) + E_res
    Gamma[sorted_local_idx] = Gamma_sorted

    if tail_idx.size:
        e_res_asymptote = E_res[-1]  # E_res at the largest real R
        Vd[tail_idx] = v0(points[tail_idx]) + e_res_asymptote
        Gamma[tail_idx] = 0.0

    return Vd, Gamma
