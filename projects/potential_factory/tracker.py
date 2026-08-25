# projects/potential_factory/tracker.py
"""Fixed-R resonance-pole solving with a spurious-pole gate, the c-product
Hellmann-Feynman gradient, Newton in (lam, alpha), and continuation over R."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from qscat.core.grids import electronic_grid
from qscat.dvr import FemDvrEcsGrid, eigen, kinetic
from qscat.ecs import find_resonance_pole
from qscat.exceptions import ConvergenceError
from qscat.linalg import c_product

__all__ = [
    "DEFAULT_BOUND_WINDOW",
    "DEFAULT_RES_WINDOW",
    "MAX_STEP",
    "ElectronicPair",
    "Pole",
    "TrackResult",
    "WellParams",
    "Window",
    "pole_sensitivity",
    "pole_vector",
    "solve_pole_params",
    "track_curve",
    "well_potential",
]

Window = tuple[float, float, float, float]
VFn = Callable[[npt.NDArray[np.complex128]], npt.NDArray[np.complex128]]

DEFAULT_RES_WINDOW: Window = (0.006, 0.8, -0.6, 0.0)
DEFAULT_BOUND_WINDOW: Window = (-2.0, -1e-6, -1e-6, 1e-6)
_BOUND_IM_TOL = 1e-6

# BRIDGING: `DEFAULT_BOUND_WINDOW`'s real range is hard
# capped at -1e-6, so recentering candidate search on it while `pole` is
# still bound-classified can NEVER see a candidate that has crossed into the
# resonant regime (real >= -1e-6) -- the walk asymptotes toward the branch
# point but can never cross it (measured: geometric convergence of the real
# part toward 0 from below, every larger step rejected as "no gated pole").
# While classes are mismatched, search a window wide enough to see BOTH
# regimes (the union of DEFAULT_BOUND_WINDOW's and DEFAULT_RES_WINDOW's real/
# imag extents); `ElectronicPair.pole`'s own bound/e_floor/residual gates
# still apply unchanged, so this only widens what the SEARCH can see, not
# what is ultimately accepted.
_BRIDGE_WINDOW: Window = (-2.0, 0.8, -0.6, 1e-6)

# `_BRIDGE_WINDOW`'s union search has no tie-break toward the pole
# being TRACKED -- `find_resonance_pole` returns the global residual-argmin
# over the whole window, so a second, more angle-stable state anywhere in
# (-2.0, 0.8) (e.g. another bound state deeper in the well) could be silently
# selected instead of the state Newton is actually walking. `_MAX_STEP` is
# the continuity guard's radius (also the recentred candidate window's
# half-width, unifying the two uses of "how far can one Newton step move the
# tracked pole").
_MAX_STEP = 0.15
# Public alias so `extract.py`'s per-node T1
# gated walk can build the same recentred-search half-width as `track_curve`
# without reaching into a private name.
MAX_STEP = _MAX_STEP


@dataclass(frozen=True)
class Pole:
    energy: complex
    residual: float

    @property
    def gamma(self) -> float:
        return max(0.0, -2.0 * self.energy.imag)

    @property
    def shift(self) -> float:
        return float(self.energy.real)


class ElectronicPair:
    """Two electronic FEM-DVR-ECS grids differing only in the ECS angle, with
    their kinetic matrices precomputed once."""

    def __init__(
        self,
        angles: tuple[float, float] = (35.0, 44.0),
        r_max: float = 16.0,
        order: int = 7,
        n_complex: int = 6,
    ) -> None:
        self.grid_a: FemDvrEcsGrid = electronic_grid(
            r_max=r_max, order=order, n_complex=n_complex, angle_deg=angles[0]
        )
        self.grid_b: FemDvrEcsGrid = electronic_grid(
            r_max=r_max, order=order, n_complex=n_complex, angle_deg=angles[1]
        )
        self._Ta = kinetic(self.grid_a, 1.0)
        self._Tb = kinetic(self.grid_b, 1.0)

    def hamiltonians(
        self, v_fn: VFn
    ) -> tuple[npt.NDArray[np.complex128], npt.NDArray[np.complex128]]:
        return (
            self._Ta + np.diag(v_fn(self.grid_a.points)),
            self._Tb + np.diag(v_fn(self.grid_b.points)),
        )

    def pole(
        self,
        v_fn: VFn,
        window: Window,
        *,
        resid_tol: float = 1e-3,
        gate_frac: float = 0.05,
        e_floor: float = 0.006,
    ) -> Pole | None:
        """The angle-stable pole in `window`, or None if no candidate passes the gate.

        `window` is ABSOLUTE: pass `v_fn` on whatever energy scale your window
        is on. The factory's own callers use the bare well (`v_int +
        ell(ell+1)/2r^2`, i.e. WITHOUT `v0(R)`), so a returned `Pole.shift` is
        an electronic-shift energy relative to `v0(R)`, not an absolute one.
        """
        Ha, Hb = self.hamiltonians(v_fn)
        try:
            E, resid = find_resonance_pole(eigen(Ha)[0], eigen(Hb)[0], window)
        except ValueError:
            return None
        p = Pole(complex(E), float(resid))
        if p.residual >= resid_tol:
            return None
        bound = p.energy.real < 0.0 and abs(p.energy.imag) < _BOUND_IM_TOL
        if bound:
            return p
        if p.gamma <= 0.0 or p.residual >= gate_frac * p.gamma:
            return None
        if p.energy.real < e_floor:
            return None
        return p


def pole_vector(H: npt.NDArray[np.complex128], E: complex) -> npt.NDArray[np.complex128]:
    """Right eigenvector of the eigenvalue of `H` nearest `E`, c-normalised (psi^T psi = 1)."""
    vals, vecs = eigen(H)
    j = int(np.argmin(np.abs(vals - E)))
    psi = np.asarray(vecs[:, j], dtype=np.complex128)
    norm = c_product(psi, psi)
    if abs(norm) < 1e-12:
        raise ValueError("eigenvector is c-product self-orthogonal; cannot normalise")
    return psi / np.sqrt(norm)


def pole_sensitivity(H: npt.NDArray[np.complex128], E: complex) -> npt.NDArray[np.complex128]:
    """`dE/dV_i = psi_i^2` for complex-symmetric `H` (Hellmann-Feynman under the c-product)."""
    psi = pole_vector(H, E)
    return psi * psi


@dataclass(frozen=True)
class WellParams:
    lam: float
    alpha: float


def well_potential(ell: int, lam: float, alpha: float, extra: VFn | None) -> VFn:
    def v(r: npt.NDArray[np.complex128]) -> npt.NDArray[np.complex128]:
        rr = np.asarray(r, dtype=np.complex128)
        out = -lam * np.exp(-alpha * rr**2) + ell * (ell + 1) / (2.0 * rr**2)
        if extra is not None:
            out = out + extra(rr)
        return np.asarray(out, dtype=np.complex128)

    return v


def _is_bound_target(target: complex) -> bool:
    # Tolerate up to `_BOUND_IM_TOL` (the
    # SAME threshold `ElectronicPair.pole` already uses to call a FOUND pole
    # "bound", line ~135) rather than requiring `target.imag == 0.0` exactly.
    # `fit.py::pole_target` always constructs a bound target with EXACTLY
    # `0.0`, so this widening changes nothing for the package's real callers
    # (every resonant `gamma_target` used there clears `tol.gamma_floor >=
    # 2e-3`, three orders of magnitude above `_BOUND_IM_TOL`) -- it only
    # matters for a target that itself carries the pole finder's own noise
    # floor (~1e-7, see `_BOUND_IM_TOL`; measured on a genuinely real
    # state), which must still be routed to the bound branch below
    # for its real-part-only convergence check to apply.
    return abs(target.imag) < _BOUND_IM_TOL and target.real < 0.0


def _is_bound_pole(p: Pole) -> bool:
    return p.energy.real < 0.0 and abs(p.energy.imag) < 1e-6


def _window_about(E: complex, bound: bool) -> Window:
    # 0.15 Ha half-width (not 0.05) -- these windows recentre on a
    # Newton CANDIDATE pole close to the current iterate, not on the target,
    # which may still be several tenths of an eV away mid-iteration.
    if bound:
        return DEFAULT_BOUND_WINDOW
    # im_hi = 1e-6 (not 0.0) so a pole that is momentarily bound,
    # or that sits exactly on the real axis mid-bridging, stays visible to
    # the window rather than being excluded by a hard Im <= 0 edge.
    return (E.real - _MAX_STEP, E.real + _MAX_STEP, min(E.imag - _MAX_STEP, -_MAX_STEP), 1e-6)


def solve_pole_params(
    pair: ElectronicPair,
    ell: int,
    target: complex,
    seed: WellParams,
    *,
    extra: VFn | None = None,
    window: Window | None = None,
    max_iter: int = 30,
    tol: float = 1e-8,
) -> tuple[WellParams, Pole]:
    """Newton on (lam, alpha) so that the gated pole equals `target`.

    For a bound target (`target.imag == 0` and `target.real < 0`) only `lam`
    is solved; `alpha` is held fixed at `seed.alpha`.
    """
    bound = _is_bound_target(target)
    lam, alpha = seed.lam, seed.alpha
    r = pair.grid_a.points
    # The seed is deliberately far from the target, so the FIRST
    # search uses the default (or caller-supplied) window, not one recentred
    # on the (possibly distant) target -- a narrow window about the target
    # can miss the seed's own pole entirely.
    win = window if window is not None else (DEFAULT_BOUND_WINDOW if bound else DEFAULT_RES_WINDOW)
    pole = pair.pole(well_potential(ell, lam, alpha, extra), win)
    if pole is None and window is None:
        # Bridging must also cover the BOOTSTRAP, not just later
        # iterations: a seed frozen from a different (lam, alpha) regime can
        # have its own true class differ from the target's -- e.g. a target
        # is resonant but the seed itself is still a genuinely bound state.
        # `win` above, chosen from the TARGET's class, then searches the
        # wrong energy range and finds nothing even though the seed's own
        # pole exists (just under the other window). Try the other default
        # window before giving up, so the first Newton iteration has a real
        # (mismatched-class) pole to bridge from.
        other = DEFAULT_RES_WINDOW if bound else DEFAULT_BOUND_WINDOW
        pole = pair.pole(well_potential(ell, lam, alpha, extra), other)
    if pole is None:
        raise ConvergenceError("no gated pole at the seed parameters")
    for _ in range(max_iter):
        f = pole.energy - target
        # A BOUND target only ever moves `lam`
        # (the step below is real-part-only, `alpha` frozen) -- `f.imag` for a
        # bound target is entirely the ECS pole-finder's own irreducible
        # numerical noise floor (measured ~1e-7 on a genuinely real
        # state), not a residual the
        # lam-only step can or should try to drive down. Requiring
        # `abs(f) < tol` on the FULL complex residual made that noise floor
        # unreachable and stalled Newton even after `f.real` converged to
        # round-off. `bound` targets converge on `abs(f.real) < tol`; every
        # other target (resonant, or still bridging into resonant) keeps the
        # full complex test, since there `alpha` genuinely participates (or,
        # for bridging, is expected to once the pole's own class catches up).
        if (abs(f.real) < tol) if bound else (abs(f) < tol):
            return WellParams(lam, alpha), pole
        Ha, _ = pair.hamiltonians(well_potential(ell, lam, alpha, extra))
        s = pole_sensitivity(Ha, pole.energy)
        g = np.exp(-alpha * r**2)
        dlam = complex(np.sum(s * (-g)))
        dalpha = complex(np.sum(s * (lam * r**2 * g)))
        # CLASSIFICATION BRIDGING. A true bound state has Gamma == 0
        # identically over a whole neighbourhood of (lam, alpha) -- not just
        # approximately -- so the Im-energy row of the full 2x2 Jacobian is
        # genuinely degenerate (cond ~1e12 measured) whenever the pole's OWN
        # current class doesn't match the target's class (e.g. continuing a
        # resonant target from a seed frozen at a bound state). While that
        # mismatch persists, fall back to a lam-only step on the real part
        # (the same reduced step the bound-target branch below always uses)
        # and accept candidates by real-part-only improvement; once the
        # pole's own class catches up with the target's, resume the full 2x2
        # step with the ordinary |f|-decrease acceptance.
        bridging = (not bound) and _is_bound_pole(pole)
        if bound or bridging:
            if dlam.real == 0.0:
                raise ConvergenceError(
                    f"zero real dlam sensitivity at lam={lam:.6g}, alpha={alpha:.6g}; "
                    "cannot take a lam-only step"
                )
            step = np.array([-f.real / dlam.real, 0.0])
        else:
            J = np.array([[dlam.real, dalpha.real], [dlam.imag, dalpha.imag]])
            try:
                step = np.linalg.solve(J, -np.array([f.real, f.imag]))
            except np.linalg.LinAlgError as e:
                raise ConvergenceError(
                    f"singular Jacobian at lam={lam:.6g}, alpha={alpha:.6g}"
                ) from e
        damp = 1.0
        while damp > 1e-4:
            lam_n, alpha_n = lam + damp * step[0], alpha + damp * step[1]
            if lam_n <= 0.0 or alpha_n <= 0.0:
                damp *= 0.5
                continue
            v_n = well_potential(ell, lam_n, alpha_n, extra)
            win_n = _BRIDGE_WINDOW if bridging else _window_about(pole.energy, _is_bound_pole(pole))
            cand = pair.pole(v_n, win_n)
            if cand is not None:
                if bound or bridging:
                    # Real-part-only metric: applies BOTH to bridging (where
                    # the Im row is genuinely degenerate at a bound-classified
                    # pole, so |f| is meaningless) and to the pre-existing,
                    # unconditional bound-TARGET path, which used the
                    # full complex `abs(cand.energy - target) < abs(f)` -- this
                    # DOES change that path's acceptance metric, but not its
                    # outcome: for a bound target `target.imag == 0` and
                    # `f.imag`/`cand.energy.imag` stay near round-off, so the
                    # two metrics agree on which candidate is closer there.
                    improved = abs(cand.energy.real - target.real) < abs(
                        pole.energy.real - target.real
                    )
                else:
                    improved = abs(cand.energy - target) < abs(f)
                # CONTINUITY GUARD. `_BRIDGE_WINDOW` (and, in
                # principle, any window) can contain more than one gated
                # pole -- `find_resonance_pole` returns the global residual-
                # argmin over the whole window, with no preference for the
                # state actually being tracked. Without this guard a second,
                # more angle-stable state elsewhere in the window (e.g. a
                # deeper bound state) could be silently substituted for the
                # one Newton is walking. Require the candidate to also be
                # within one Newton-step radius of the CURRENT pole; a
                # closer-to-target but far-away candidate is rejected exactly
                # like a non-improving one (halve the damping and retry).
                continuous = abs(cand.energy - pole.energy) <= _MAX_STEP
                if improved and continuous:
                    lam, alpha, pole = lam_n, alpha_n, cand
                    break
            damp *= 0.5
        else:
            msg = f"Newton stalled at lam={lam:.6g}, alpha={alpha:.6g}, |f|={abs(f):.3g}"
            raise ConvergenceError(msg)
    raise ConvergenceError(f"Newton did not converge in {max_iter} iterations")


@dataclass(frozen=True)
class TrackResult:
    R: npt.NDArray[np.float64]
    lam: npt.NDArray[np.float64]
    alpha: npt.NDArray[np.float64]
    converged: npt.NDArray[np.bool_]
    poles: list[Pole | None]


def track_curve(
    pair: ElectronicPair,
    ell: int,
    R_desc: npt.ArrayLike,
    target_of_R: Callable[[float], complex],
    seed: WellParams,
    *,
    extra_of_R: Callable[[float], VFn | None] | None = None,
    alpha_of_R: Callable[[float], float] | None = None,
) -> TrackResult:
    """Continuation of `solve_pole_params` over a strictly-descending `R` grid.

    Each node is seeded from the previous accepted node's `(lam, alpha)`. On a
    `ConvergenceError`, or a `target_of_R` that reports "no target" via a NaN
    energy (`target.real` is NaN -- e.g. the pole exists but is gated out at
    that R), the node is flagged `converged=False`, its `(lam, alpha)` are
    copied from the last accepted node (freeze-and-flag), and tracking
    continues to the next node.

    `alpha_of_R`, when given, overrides the
    `alpha` component of every node's seed with `alpha_of_R(R_j)` -- the
    node's `lam` still comes from the previous accepted node (`seed.lam` at
    the first node). Without this, a BOUND target's Newton step never moves
    `alpha` (`solve_pole_params` solves `lam` only there), so the whole bound
    branch would otherwise run at whatever `alpha` the very first node's seed
    carried, not at a caller-supplied `alpha(R)` curve -- invisible when that
    curve happens to be constant, wrong whenever it is not.
    """
    R = np.asarray(R_desc, dtype=np.float64)
    if R.size > 1 and np.any(np.diff(R) >= 0.0):
        raise ValueError("R_desc must be strictly descending")
    lam = np.empty(R.size)
    alpha = np.empty(R.size)
    ok = np.zeros(R.size, dtype=bool)
    poles: list[Pole | None] = []
    cur = seed
    for j, Rj in enumerate(R):
        extra = extra_of_R(float(Rj)) if extra_of_R is not None else None
        seed_j = cur if alpha_of_R is None else WellParams(cur.lam, float(alpha_of_R(float(Rj))))
        target = target_of_R(float(Rj))
        if np.isnan(target.real):
            poles.append(None)
        else:
            try:
                cur, pole = solve_pole_params(pair, ell, target, seed_j, extra=extra)
                ok[j] = True
                poles.append(pole)
            except ConvergenceError:
                poles.append(None)
        lam[j], alpha[j] = cur.lam, cur.alpha
    return TrackResult(R=R, lam=lam, alpha=alpha, converged=ok, poles=poles)
