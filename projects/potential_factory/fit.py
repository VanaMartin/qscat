"""Staged fitting: T0 (neutral) -> T1 (pole curves) -> T3 (energy-dependent width).
Each stage seeds the next; a stage that misses its tolerance stops and reports."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import replace

import numpy as np
import numpy.typing as npt
from qscat.core.dissociation import anion_electronic_states
from qscat.core.grids import nuclear_grid
from qscat.core.nrm.coupling import gamma_from_coupling, v_dk_plus
from qscat.core.nrm.discrete_state import AsymptoticDiscreteState
from qscat.core.vibrational import vibrational_states
from qscat.exceptions import ConvergenceError, GridError
from scipy.optimize import least_squares
from scipy.special import expit

from projects.potential_factory.ansatz import (
    FlexibleDiatomicModel,
    SmoothR,
    pack,
    params,
    unpack,
    with_params,
    y_p,
)
from projects.potential_factory.extract import walk_t1
from projects.potential_factory.report import FitReport, TierResult, Tolerances, ecs_bounded
from projects.potential_factory.target import CouplingTarget, NeutralTarget, ResonanceTarget, Target
from projects.potential_factory.tracker import (
    DEFAULT_BOUND_WINDOW,
    DEFAULT_RES_WINDOW,
    MAX_STEP,
    ElectronicPair,
    Pole,
    WellParams,
    Window,
    pole_sensitivity,
    track_curve,
    well_potential,
)

__all__ = ["fit", "fit_coupling", "fit_neutral", "fit_resonance", "model_gamma_tilde"]

_R_INF = 10.0
# What one node of `_joint_polish` costs when it yields no gated pole
# at all. ~1000x a converged per-node residual, so a
# trial that loses nodes is rejected outright -- see `_joint_polish`.
_POLISH_MISS_PENALTY = 1.0
# How much of the v=0/v=1 probability may sit OUTSIDE a target table
# before that table is judged too narrow to define a vibrational
# spacing at all -- see `_omega_e_reference`.
_TAIL_WEIGHT_TOL = 1e-8


def _with_n_beta(model: FlexibleDiatomicModel, n_beta: int) -> FlexibleDiatomicModel:
    """Grow (with zeros) or truncate the EMO beta expansion to `n_beta` terms."""
    betas = list(model.betas) + [0.0] * max(0, n_beta - len(model.betas))
    return replace(model, betas=tuple(betas[:n_beta]))


def _omega_e_10(
    model_mu: float, v0_of_R: Callable[[np.ndarray], np.ndarray]
) -> tuple[float, np.ndarray, np.ndarray]:
    """The `v=0 -> v=1` vibrational spacing of `v0_of_R` on the package's
    standard nuclear grid, plus the two eigenvectors and the grid points.

    This is the ANHARMONIC spacing `G(1) - G(0)`, not the harmonic constant
    `omega_e` a spectroscopic table quotes -- see `_omega_e_reference`.
    """
    g = nuclear_grid()
    eps, chi = vibrational_states(g, model_mu, 2, v0_of_R)
    return float(eps[1] - eps[0]), np.asarray(chi), np.asarray(g.points)


def _omega_e_reference(
    target: NeutralTarget, fitted: FlexibleDiatomicModel
) -> tuple[float | None, str]:
    """The `v=0 -> v=1` spacing the TARGET implies, and where it came from.

    Two routes, in order:

    1. `constants["omega_e"]`: a spectroscopic table quotes the HARMONIC
       constant, which is NOT the 0->1 spacing. For the Morse curve this
       branch of `fit_neutral` builds from those same constants, the exact
       spacing is `G(1) - G(0) = omega_e - 2 omega_e x_e` with
       `omega_e x_e = omega_e^2 / (4 D_e)`. Comparing against the bare
       `omega_e` instead would report a fixed anharmonic offset as fit error:
       measured 0.83% (N2), 1.96% (NO), 3.33% (F2) -- two of the three over
       the 1% `omega_e_rel` default, for a curve reproduced to 1e-13.
    2. The target CURVE, solved on the same grid. The curve is a table with
       NaN outside its span and is real-argument only, so the grid's ECS tail
       and any real point outside the table are filled with the FITTED curve.
       That fill is only admissible if the two levels do not reach it, which
       is CHECKED (`_TAIL_WEIGHT_TOL`) rather than assumed; a target table too
       narrow to confine `v=0,1` reports "not checked" instead.

    `None` means neither route applies.
    """
    c = target.constants
    if "omega_e" in c and "D_e" in c:
        w = float(c["omega_e"])
        return w - 2.0 * w * w / (4.0 * float(c["D_e"])), "constants (Morse G(1)-G(0))"
    if target.curve is None or target.curve.x is None:
        return None, "no table"
    x_lo, x_hi = float(target.curve.x.min()), float(target.curve.x.max())

    def v0_ref(Rq: npt.ArrayLike) -> np.ndarray:
        Rc = np.asarray(Rq, dtype=np.complex128)
        out = np.asarray(fitted.v0(Rc), dtype=np.complex128).copy()
        real_pt = np.flatnonzero(Rc.imag == 0.0)
        if real_pt.size:
            xr = Rc[real_pt].real
            inside = (xr >= x_lo) & (xr <= x_hi)
            if inside.any():
                out[real_pt[inside]] = target.curve(xr[inside])
        return out

    spacing, chi, pts = _omega_e_10(fitted.mu, v0_ref)
    outside = (pts.imag != 0.0) | (pts.real < x_lo) | (pts.real > x_hi)
    leak = float(np.max(np.sum(np.abs(chi[:, outside]) ** 2, axis=1)))
    if leak > _TAIL_WEIGHT_TOL:
        return None, f"table [{x_lo:.3g}, {x_hi:.3g}] does not confine v=0,1 (leak={leak:.1e})"
    return spacing, f"table [{x_lo:.3g}, {x_hi:.3g}] (leak={leak:.1e})"


def _check_omega_e(
    target: NeutralTarget, fitted: FlexibleDiatomicModel, tol: Tolerances
) -> tuple[bool, str]:
    """Spec step 1's ladder check: the fitted curve's own `v=0 -> v=1` spacing
    against the target's, on `qscat.core.grids.nuclear_grid()`.

    Returns `(ok, detail)`. When no reference is available `ok` is True and the
    detail says `omega_e: not checked` -- a missing check must not fail a tier
    that met its curve tolerance, but it must be VISIBLE in the report.
    """
    ref, src = _omega_e_reference(target, fitted)
    if ref is None:
        return True, f"omega_e: not checked ({src})"
    try:
        fit_w, _, _ = _omega_e_10(fitted.mu, fitted.v0)
    except (GridError, ValueError) as err:
        return False, f"omega_e: fitted curve binds no v=0,1 on the nuclear grid ({str(err)[:80]})"
    rel = abs(fit_w - ref) / abs(ref)
    return rel <= tol.omega_e_rel, (
        f"omega_e (v=0->1) fit={fit_w:.6e} target={ref:.6e} Ha, "
        f"omega_e_rel={rel:.2e} vs tol {tol.omega_e_rel:.2e} [{src}]"
    )


def fit_neutral(
    target: NeutralTarget,
    seed: FlexibleDiatomicModel,
    *,
    n_beta: int = 1,
    tol: Tolerances,
) -> tuple[FlexibleDiatomicModel, TierResult]:
    model = _with_n_beta(seed, n_beta)
    if target.curve is None:
        c = target.constants
        beta0 = c["omega_e"] * np.sqrt(model.mu / (2.0 * c["D_e"]))
        model = with_params(model, {"D_e": c["D_e"], "R_e": c["R_e"], "beta0": float(beta0)})
        ok_w, w_detail = _check_omega_e(target, model, tol)
        return model, TierResult(
            "T0",
            "met" if ok_w else "not met",
            0.0,
            0.0,
            f"constants only; Morse relation for beta0; {w_detail}",
        )
    lo, hi = target.R_range
    if target.curve.x is not None:
        # Fit AT the curve's own table nodes: a `CubicSpline` interpolates
        # exactly through its input data, so the residual is exactly zero
        # there at the true parameters -- an off-node probe grid inherits the
        # spline's O(h^4) interpolation error, which is systematic (not
        # noise least squares averages away) and caps recovery at ~1e-6/1e-7
        # instead of float64 precision.
        nodes = target.curve.x
        R = nodes[(nodes >= lo) & (nodes <= hi)]
    else:
        R = np.linspace(lo, hi, 200)
    y = target.curve(R)
    names = ["D_e", "R_e"] + [f"beta{i}" for i in range(n_beta)]

    # The EMO's beta(R) must stay positive everywhere the nuclear grid reaches,
    # or the curve blows up at small R and binds nothing: bounds on the
    # constants plus a soft penalty on beta(R) < 0 over [0.3, R_e].
    R_pen = np.linspace(0.3, max(lo, 0.31), 12)
    lower = [1e-4, 0.5, 0.05] + [-2.0] * (n_beta - 1)
    upper = [2.0, 10.0, 10.0] + [2.0] * (n_beta - 1)

    def resid(x):
        m = unpack(model, names, x)
        pen = 10.0 * np.maximum(0.0, -m.beta_R(R_pen).real)
        return np.concatenate([m.v0(R).real - y, pen])

    x0 = np.clip(pack(model, names), np.array(lower) + 1e-9, np.array(upper) - 1e-9)
    sol = least_squares(
        resid, x0, bounds=(lower, upper), xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=2000
    )
    fitted = unpack(model, names, sol.x)
    r = resid(sol.x)[: R.size]
    rms, mx = float(np.sqrt(np.mean(r**2))), float(np.max(np.abs(r)))
    ok_w, w_detail = _check_omega_e(target, fitted, tol)
    status = "met" if (rms <= tol.v0_rms and ok_w) else "not met"
    return fitted, TierResult(
        "T0", status, rms, mx, f"EMO with {n_beta} beta(s) on {R.size} points; {w_detail}"
    )


def _grow_coeffs(s: SmoothR, n: int) -> SmoothR:
    c = list(s.coeffs) + [0.0] * max(0, n - len(s.coeffs))
    return replace(s, coeffs=tuple(c[:n]))


def _smooth_rel_rms(resid_arr: np.ndarray, y: np.ndarray) -> float:
    """`rms(|f(R_j) - y_j| / |y_j|)` -- the relative node residual `_fit_smooth`
    reports and `fit_resonance` gates on. `y` floored at 1e-12 so a
    node whose target happens to be exactly 0 doesn't divide by zero; every
    physical lam(R)/alpha(R) node in this package is O(1), so the floor never
    engages in practice."""
    rel = resid_arr / np.where(np.abs(y) > 1e-12, y, 1e-12)
    return float(np.sqrt(np.mean(rel**2)))


def _smooth_reparam(
    s: SmoothR, R_min: float, R_max: float
) -> tuple[np.ndarray, list[float], list[float], Callable[[np.ndarray], SmoothR]]:
    """The log-amplitude reparametrization of ONE `SmoothR`:
    `(x0, lo, hi, decode)` where `decode(x) -> SmoothR` maps a flat vector
    back to the dataclass. `f_0 == 0` (an R-INDEPENDENT term -- every
    published `alpha(R)` in this package) is degenerate (`sig(R) == 0`
    regardless of `f_1`/`R_f`) and held fixed; the flat vector then carries
    only `f_inf` + `coeffs`. Shared by `_fit_smooth` and
    `_joint_polish` so both optimize the exact SAME
    parametrization of the exact same model -- see `_fit_smooth`'s docstring
    for the conditioning story this reparametrization fixes.

    `_smooth_grad` differentiates a decoded `SmoothR` with respect to THIS
    exact layout (the analytic Jacobian `_joint_polish` passes to
    `least_squares`); keep the two in step.
    """
    n_c = len(s.coeffs)
    if s.f_0 != 0.0:
        sign0 = 1.0 if s.f_0 >= 0.0 else -1.0
        log_amp0 = float(np.log(abs(s.f_0)))
        # See `_fit_smooth`'s docstring: these bounds are generous safety
        # rails against float64 `exp` overflow during a multi-start, never
        # physically motivated ones.
        f1_cap = 700.0 / max(1.0, (R_max - R_min) + 20.0)
        x0 = np.array([s.f_inf, log_amp0, s.f_1, s.R_f, *list(s.coeffs)], dtype=np.float64)
        lo = [-np.inf, log_amp0 - 40.0, 1e-6, -np.inf] + [-np.inf] * n_c
        hi = [np.inf, log_amp0 + 40.0, f1_cap, np.inf] + [np.inf] * n_c

        def decode(x: np.ndarray) -> SmoothR:
            f0 = sign0 * np.exp(x[1])
            return replace(
                s,
                f_inf=float(x[0]),
                f_0=float(f0),
                f_1=float(x[2]),
                R_f=float(x[3]),
                coeffs=tuple(float(c) for c in x[4 : 4 + n_c]),
            )
    else:
        x0 = np.array([s.f_inf, *list(s.coeffs)], dtype=np.float64)
        lo = [-np.inf] * (1 + n_c)
        hi = [np.inf] * (1 + n_c)

        def decode(x: np.ndarray) -> SmoothR:
            return replace(s, f_inf=float(x[0]), coeffs=tuple(float(c) for c in x[1 : 1 + n_c]))

    return x0, lo, hi, decode


def _smooth_grad(s: SmoothR, R: float) -> np.ndarray:
    """`d f(R) / d x` for a `SmoothR` decoded at the current `x`, in exactly
    `_smooth_reparam`'s layout -- `[f_inf, u=log|f_0|, f_1, R_f, *coeffs]`
    when `f_0 != 0`, else `[f_inf, *coeffs]`.

    With `t = f_1 (R - R_f)`, `sig = f_0 / (1 + e^t)` and
    `poly = 1 + sum_i coeffs[i] y_p(R)^(i+1)`:

        df/df_inf = 1
        df/du     = sig * poly          (since d f_0/du = f_0)
        df/df_1   = -sig * poly * (R - R_f) * expit(t)
        df/dR_f   = +sig * poly * f_1 * expit(t)
        df/dc_i   = sig * y_p(R)^(i+1)

    `expit` (`1/(1+e^-t)`) is used for BOTH `e^t/(1+e^t)` and `1/(1+e^t)`
    rather than forming `e^t` directly: `f_1 (R - R_f)` is unbounded during
    an optimizer's excursion, and the naive quotient is `inf/inf -> nan`
    there while the logistic form saturates cleanly at 0 or 1.

    In the degenerate `f_0 == 0` branch `sig == 0`, so the `coeffs` columns
    are legitimately zero -- the model genuinely does not depend on them,
    which is exactly why `_smooth_reparam` holds the sigmoid constants fixed
    there.
    """
    n_c = len(s.coeffs)
    poly = 1.0
    y = 0.0
    if n_c:
        y = float(y_p(R, s.R_e, s.p).real)
        poly = 1.0 + sum(c * y ** (i + 1) for i, c in enumerate(s.coeffs))
    t = s.f_1 * (R - s.R_f)
    sig = s.f_0 * float(expit(-t))
    poly_terms = [sig * y ** (i + 1) for i in range(n_c)]
    if s.f_0 == 0.0:
        return np.array([1.0, *poly_terms], dtype=np.float64)
    e_over = float(expit(t))  # e^t / (1 + e^t)
    return np.array(
        [
            1.0,
            sig * poly,
            -sig * poly * (R - s.R_f) * e_over,
            sig * poly * s.f_1 * e_over,
            *poly_terms,
        ],
        dtype=np.float64,
    )


def _fit_smooth(
    model: FlexibleDiatomicModel,
    prefix: str,
    R: np.ndarray,
    y: np.ndarray,
    weights: np.ndarray | None = None,
) -> tuple[FlexibleDiatomicModel, float]:
    """Fit `SmoothR`'s sigmoid (+ optional polynomial `coeffs`) through node
    data `(R, y)`, returning the fitted model and the relative node residual
    `_smooth_rel_rms` at the optimum.

    The optional `weights` (the alpha pass's
    caller passes `gamma_target(R_j)`, normalized to max 1) down-weight
    unreliable nodes in the OPTIMIZER's objective only (`sqrt(weights) *
    residual`, since `least_squares` minimizes a sum of squares) -- the
    REPORTED/selection metric (`_smooth_rel_rms`, both the multi-start
    winner and the returned `rel_rms`) always uses the RAW, unweighted
    residual, so it keeps answering "how well does the fitted curve
    reproduce the data" regardless of how the optimizer got there.
    `d(chi^2)/d(alpha)` vanishes as `Gamma -> 0` at fixed E_res residual, so
    a barely-open node's tracked `alpha` is intrinsically noisier than a
    fully-open node's and must not carry equal weight in the alpha(R) fit
    (confirmed on F2: R=2.42, Gamma_target=0.027, carried a 3.4% alpha
    error that an unweighted fit could not average out). The `lam` pass
    passes no `weights` and is unaffected.

    The RAW `(f_0, f_1, R_f)` parametrization
    is catastrophically conditioned whenever the seed's amplitude `f_0` is
    astronomically large and its inflection `R_f` sits far outside the
    physical `R` range -- exactly N2's published `lam(R)`: `f_0 ~ -7.4e13`,
    `R_f ~ -28` against a physical `R` in [1.6, 3.0]. There,
    `sig(R) ~ f_0 * exp(-f_1*(R-R_f))`, so a few-percent relative change in
    `f_1` times `|R-R_f| ~ 30` moves the exponent by O(1) -- i.e. the curve by
    a factor of `e` -- even though the tracked pole data feeding this fit
    agrees with the true curve to ~1e-7 (confirmed by direct comparison,
    measured). A plain `least_squares` call over the raw
    constants (the pre-fix code) converged to a WORSE point on that ridge
    (N2: fitted `f_0` doubled to `-1.48e14`, `R_f` drifted to `-29.4`) while
    still leaving a 4-7% relative residual against the very data it was
    fitting -- not overfitting noise, just a badly scaled Jacobian.

    Fixed by fitting in `u = log|f_0|` (sign frozen from the seed, since the
    amplitude never needs to change sign within one tier) with
    `x_scale="jac"` (each parameter's step follows its own residual
    sensitivity instead of raw units) and `f_1` bounded `> 0` via `method=
    "trf"`. When the seed's `f_0` is exactly 0 (an R-INDEPENDENT term -- e.g.
    every one of N2/NO/F2's published `alpha(R)`), `f_1`/`R_f` are degenerate
    (`sig(R) == 0` identically regardless of their value) and are held fixed;
    only `f_inf` and any polynomial `coeffs` are fit. If the primary fit's
    relative residual still exceeds 1e-4, a small multi-start over
    `f_1 in {0.5, 1, 2, 4} x seed` and `R_f in {R_min-2, R_mid, R_max+2}`
    reruns the same reparametrized fit from each combination and keeps
    whichever converges to the lowest relative residual -- the physical curve
    is what's being validated, not any particular point on the degenerate
    ridge of constants that reproduces it.
    """
    s = getattr(model, prefix)
    R_min, R_max = float(np.min(R)), float(np.max(R))
    x0, lo, hi, decode = _smooth_reparam(s, R_min, R_max)
    fit_sigmoid = s.f_0 != 0.0
    sqrt_w = (
        np.ones_like(y)
        if weights is None
        else np.sqrt(np.maximum(np.asarray(weights, dtype=np.float64), 0.0))
    )

    def raw_resid(x: np.ndarray) -> np.ndarray:
        """RAW (unweighted) residual -- used for reporting and multi-start
        selection, never handed straight to `least_squares`."""
        m = replace(model, **{prefix: decode(x)})
        return np.asarray(getattr(m, prefix)(R).real - y, dtype=np.float64)

    def resid(x: np.ndarray) -> np.ndarray:
        return sqrt_w * raw_resid(x)

    if fit_sigmoid:
        # Bound the log-amplitude and f_1 well inside float64's exponent
        # range: an unbounded multi-start candidate can otherwise walk
        # `f1 * (R - R_f)` past ~700 and overflow `exp` (RuntimeWarning,
        # `sig` -> inf/nan) -- see `_smooth_reparam`, these bounds are
        # enormously generous relative to any physical lam(R)/alpha(R) scale
        # so they never bind a genuine fit.
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            sol = least_squares(
                resid,
                x0,
                bounds=(lo, hi),
                method="trf",
                x_scale="jac",
                xtol=1e-14,
                ftol=1e-14,
                max_nfev=5000,
            )
            best_x, best_rms = sol.x, _smooth_rel_rms(raw_resid(sol.x), y)
            if best_rms > 1e-4:
                R_mid = 0.5 * (R_min + R_max)
                for f1_0 in (0.5 * s.f_1, s.f_1, 2.0 * s.f_1, 4.0 * s.f_1):
                    for rf_0 in (R_min - 2.0, R_mid, R_max + 2.0):
                        x0_i = x0.copy()
                        x0_i[2] = float(np.clip(f1_0, lo[2], hi[2]))
                        x0_i[3] = rf_0
                        try:
                            sol = least_squares(
                                resid,
                                x0_i,
                                bounds=(lo, hi),
                                method="trf",
                                x_scale="jac",
                                xtol=1e-14,
                                ftol=1e-14,
                                max_nfev=5000,
                            )
                        except (ValueError, RuntimeError):
                            continue
                        rms = _smooth_rel_rms(raw_resid(sol.x), y)
                        if rms < best_rms:
                            best_x, best_rms = sol.x, rms
    else:
        # f_0 == 0: the sigmoid term is degenerate (contributes nothing
        # regardless of f_1/R_f) -- `_smooth_reparam` already holds it fixed
        # and returns a flat vector of just f_inf + coeffs.
        sol = least_squares(resid, x0, xtol=1e-14, ftol=1e-14, max_nfev=5000)
        best_x, best_rms = sol.x, _smooth_rel_rms(raw_resid(sol.x), y)

    return replace(model, **{prefix: decode(best_x)}), best_rms


def _shell_extra_of_R(model: FlexibleDiatomicModel):
    """The barrier/shell term, as a function of
    `R`, on the SAME scale `track_curve` expects (`extra` added to the bare
    well before the centrifugal term -- see `well_potential`). `None` when
    the model has no shell, so tracking sees exactly the bare-well-only
    potential it always did. `walk_t1`'s re-verify already includes the
    shell (it evaluates `model.surface`, which always adds it); without this,
    tracking and re-verify would see DIFFERENT potentials whenever a shell is
    present.
    """
    if model.shell is None:
        return None

    def extra_of_R(R: float):
        def extra(r):
            rr = np.asarray(r, dtype=np.complex128)
            return np.asarray(
                model.shell_R(np.asarray(R, dtype=np.complex128))
                * np.exp(-model.alpha_b * (rr - model.r_b) ** 2),
                dtype=np.complex128,
            )

        return extra

    return extra_of_R


def _apply_ea_constraint(
    model: FlexibleDiatomicModel, pair: ElectronicPair, ea: float
) -> tuple[FlexibleDiatomicModel, bool]:
    """Returns `(model, applied)`. This must
    NEVER raise out of the fit -- `anion_electronic_states` can raise when
    the well unbinds at a probed `lam.f_inf` (e.g. `0.5 * f0`), which is
    exactly the kind of failure a bracket search is expected to run into.
    Both the endpoint probes AND `brentq`'s own interior evaluations are
    guarded; any failure is treated as "no bracket" and the model is
    returned unchanged.
    """

    def g(f_inf: float) -> float:
        m = with_params(model, {"lam.f_inf": f_inf})
        eps_e, _ = anion_electronic_states(pair.grid_a, m, _R_INF, 1)
        return float(eps_e[0] - m.v0(_R_INF).real + ea)

    def safe_g(f_inf: float) -> float:
        # an UNBOUND anion sits above every negative -ea: count it as "positive",
        # so the unbinding edge is just where the sign of g changes
        try:
            return g(f_inf)
        except (ValueError, ConvergenceError):
            return abs(ea) if ea != 0.0 else 1.0

    f0 = model.lam.f_inf
    lo = 0.5 * f0
    if safe_g(lo) <= 0.0:  # still bound below -ea at half depth: no bracket on this side
        return model, False
    hi = 1.5 * f0
    while safe_g(hi) > 0.0 and hi < 4.0 * f0:
        hi *= 1.5
    if safe_g(hi) > 0.0:
        return model, False
    # bisection (g may be discontinuous at the unbinding edge, so not brentq)
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if safe_g(mid) > 0.0:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-9 * f0:
            break
    return with_params(model, {"lam.f_inf": hi}), True


def _polish_window_about(E: complex, bound: bool) -> Window:
    """A per-node window recentred on `E`, rebuilt here from PUBLIC building
    blocks (`MAX_STEP`, `DEFAULT_BOUND_WINDOW`) since `tracker.py`'s
    `_window_about` is not exported.

    Both branches recentre. The resonant branch is `_window_about`'s formula
    verbatim. The BOUND branch recentres only the REAL interval and keeps
    `DEFAULT_BOUND_WINDOW`'s narrow imaginary band (`im_lo = -1e-6`,
    `im_hi = 1e-6`): a bound state sits on the real axis by definition, so
    opening the imaginary side to `-MAX_STEP` would put the whole resonant
    branch inside a bound node's window for `find_resonance_pole` -- a global
    residual-argmin -- to return instead. `im_hi` is `+1e-6` rather than `0`
    so a pole sitting exactly on the axis stays visible to the window.
    `tracker._window_about` returns the unrecentred `DEFAULT_BOUND_WINDOW`
    here; recentring is strictly narrower, and the wide window is still tried
    as the fallback below.
    """
    if bound:
        re_lo, re_hi, im_lo, im_hi = DEFAULT_BOUND_WINDOW
        return (max(E.real - MAX_STEP, re_lo), min(E.real + MAX_STEP, re_hi), im_lo, im_hi)
    return (E.real - MAX_STEP, E.real + MAX_STEP, min(E.imag - MAX_STEP, -MAX_STEP), 1e-6)


def _joint_polish(
    model: FlexibleDiatomicModel,
    pair: ElectronicPair,
    R_nodes: np.ndarray,
    target_of_R: Callable[[float], complex],
    gamma_t: np.ndarray,
    resonant_mask: np.ndarray,
    seed_poles: list[Pole | None],
    extra_of_R: Callable[[float], object] | None,
    max_nfev: int = 100,
) -> tuple[FlexibleDiatomicModel, str]:
    """A JOINT least-squares polish of `lam`
    AND `alpha` together, minimizing the residual directly on the POLES
    themselves (`E_pole_j(lam(R_j), alpha(R_j)) - E_target_j`) rather than on
    the per-node TRACKED `(lam_j, alpha_j)` sample `_fit_smooth` was fit
    through. `_fit_smooth`'s loss lives on that per-node sample, which can be
    mildly contaminated (a near-threshold node's tracked alpha, a bound-
    branch node's lam biased by an imperfect alpha(R)) enough for a flexible
    4-parameter sigmoid to overfit a wrong-shape curve that nonetheless fits
    the contaminated sample well (measured on F2's lam(R): a
    fitted curve with `rel_rms=0.0048` against the tracked sample scored
    WORSE, 0.0092-0.0124, when evaluated at the TRUE parameters). Putting the
    loss back on the poles removes that indirection: the residual IS the
    same quantity `track_curve`'s own Newton solve drives to zero, on the
    SAME reparametrization (`_smooth_reparam`) `_fit_smooth` uses,
    so both routes describe the model in the same well-conditioned
    coordinates.

    Residual vector: `Re(E_pole_j - E_target_j)` for every node in
    `R_nodes`, followed by `Im(E_pole_j - E_target_j)` for RESONANT nodes
    only (`resonant_mask`), weighted `sqrt(gamma_target_j / max)` exactly as
    `_fit_smooth` weights the alpha fit -- a barely-open node's pole position is
    less reliable, and (unlike `_fit_smooth`'s per-node table) here Im is
    only physically meaningful on the resonant branch at all.

    A node with no gated pole at the trial `(lam(R_j), alpha(R_j))` (or one
    whose `lam`/`alpha` has gone unphysical, i.e. non-positive) contributes
    `_POLISH_MISS_PENALTY` to each of its residual entries -- NOT 0, and not
    dropped (`least_squares` needs a fixed-length vector across
    evaluations). Contributing 0 made "find
    no pole anywhere" a GLOBAL MINIMUM of the sum of squares, and on F2's
    wider node set the optimizer walked straight into it -- one accepted
    step to `lam.f_inf = -0.018`, `alpha.f_inf = -0.090` (a well of negative
    depth and negative width, i.e. no bound or resonant state anywhere),
    reported as `polish_rms=0.00e+00, skipped 7/7, status=1` after nfev=2,
    and left the model 100% wrong at every node. The penalty (1 Ha, ~1000x a
    converged per-node residual) makes such a step ruinously expensive
    instead, so the trust region rejects it and shrinks. It is deliberately
    a CONSTANT: it carries no gradient, so it never steers the fit, it only
    walls off the region where the residual stops being defined.
    `polish_rms` in the returned detail is reported over the nodes that DID
    yield a pole, so the penalty inflates neither the reported accuracy nor
    its comparison with `e_res_rms`; the skip count is reported next to it.

    Each node's pole search recentres on that SAME node's pole from the
    PREVIOUS residual evaluation (`_polish_window_about`, falling back to
    the default bound/resonant window on the first call or whenever nothing
    was found there last time) -- `track_curve`'s own continuity idea,
    applied here across optimizer iterations at a fixed R rather than across
    R at a fixed iterate. Once a node HAS a previous pole, any candidate
    (recentred or fallback) must land within `MAX_STEP` of it to be
    accepted, exactly as `solve_pole_params` and `walk_t1` require -- see the
    guard's own comment in `evaluate`.

    The Jacobian is ANALYTIC: Hellmann-Feynman under the
    c-product gives `dE/dV_i = psi_i^2` (`pole_sensitivity`), hence per node
    `dE/dlam = sum_i s_i (-e^{-alpha r_i^2})` and `dE/dalpha = sum_i s_i
    (lam r_i^2 e^{-alpha r_i^2})` -- the SAME two derivatives
    `solve_pole_params` builds its Newton step from -- chained through
    `_smooth_grad`, the derivative of `lam(R_j)`/`alpha(R_j)` with respect
    to `_smooth_reparam`'s coordinates. It costs ONE extra dense
    eigendecomposition per node per Jacobian, against `(n_params + 1)` full
    pole searches (two eigendecompositions each) for the numeric
    alternative, which is what makes the raised `max_nfev` affordable. A
    node with no pole (see above) has a genuinely zero Jacobian row, since
    its residual entries are the constant penalty.
    """
    ell = model.ell
    R_min, R_max = float(np.min(R_nodes)), float(np.max(R_nodes))
    x0_lam, lo_lam, hi_lam, decode_lam = _smooth_reparam(model.lam, R_min, R_max)
    x0_alpha, lo_alpha, hi_alpha, decode_alpha = _smooth_reparam(model.alpha, R_min, R_max)
    n_lam = x0_lam.size
    x0 = np.concatenate([x0_lam, x0_alpha])
    lo = lo_lam + lo_alpha
    hi = hi_lam + hi_alpha

    targets = np.array([target_of_R(float(Rj)) for Rj in R_nodes], dtype=np.complex128)
    w = np.zeros(R_nodes.size)
    if resonant_mask.any():
        g = np.maximum(gamma_t[resonant_mask], 0.0)
        w[resonant_mask] = g / max(float(np.max(g)), 1e-300)
    sqrt_w_im = np.sqrt(w[resonant_mask])

    last_pole: list[Pole | None] = list(seed_poles)
    # The state of the MOST RECENT residual evaluation, so `jac` can reuse
    # that evaluation's poles instead of searching for them a second time.
    cache: dict[str, object] = {"x": None}

    def evaluate(x: np.ndarray) -> None:
        lam_s = decode_lam(x[:n_lam])
        alpha_s = decode_alpha(x[n_lam:])
        poles: list[Pole | None] = [None] * R_nodes.size
        wells = np.zeros((R_nodes.size, 2))
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            for j, Rj_arr in enumerate(R_nodes):
                Rj = float(Rj_arr)
                lam_j = float(lam_s(Rj).real)
                alpha_j = float(alpha_s(Rj).real)
                wells[j] = (lam_j, alpha_j)
                if not (lam_j > 0.0 and alpha_j > 0.0):
                    continue
                extra = extra_of_R(Rj) if extra_of_R is not None else None
                v_fn = well_potential(ell, lam_j, alpha_j, extra)
                bound_j = not bool(resonant_mask[j])
                windows: list[Window] = []
                seen = last_pole[j]
                if seen is not None:
                    windows.append(_polish_window_about(seen.energy, bound_j))
                windows.append(DEFAULT_BOUND_WINDOW if bound_j else DEFAULT_RES_WINDOW)
                for win in windows:
                    cand = pair.pole(v_fn, win)
                    # CONTINUITY GUARD, the same one `solve_pole_params` and
                    # `walk_t1` carry: `find_resonance_pole` is a GLOBAL
                    # residual-argmin over the window with no preference for
                    # the state being tracked, so a wide fallback window can
                    # return a different angle-stable state. Once this node
                    # has a previous pole, a candidate further than one
                    # recentred-window step from it is discarded as if
                    # nothing had been found (the miss penalty below), rather
                    # than latched into `last_pole[j]` for every later
                    # iteration to recentre on.
                    if cand is not None and (
                        seen is None or abs(cand.energy - seen.energy) <= MAX_STEP
                    ):
                        poles[j] = cand
                        last_pole[j] = cand
                        break
        cache["x"] = np.array(x, dtype=np.float64)
        cache["poles"] = poles
        cache["wells"] = wells
        cache["lam_s"] = lam_s
        cache["alpha_s"] = alpha_s

    def _ensure(x: np.ndarray) -> None:
        cached_x = cache["x"]
        if cached_x is None or not np.array_equal(cached_x, x):
            evaluate(x)

    def residual(x: np.ndarray) -> np.ndarray:
        _ensure(x)
        poles = cache["poles"]
        assert isinstance(poles, list)
        re = np.full(R_nodes.size, _POLISH_MISS_PENALTY)
        im_list: list[float] = []
        for j, pole in enumerate(poles):
            if pole is None:
                if resonant_mask[j]:
                    im_list.append(_POLISH_MISS_PENALTY)
                continue
            f = pole.energy - targets[j]
            re[j] = f.real
            if resonant_mask[j]:
                im_list.append(f.imag)
        # The miss penalty must NOT be scaled by the Gamma weight: a
        # down-weighted node would otherwise be cheap to lose, which is
        # precisely the failure the penalty exists to prevent.
        im = np.asarray(im_list, dtype=np.float64)
        found_im = np.array([poles[j] is not None for j in np.flatnonzero(resonant_mask)])
        if im.size:
            im = np.where(found_im, im * sqrt_w_im, im)
        return np.concatenate([re, im])

    def jac(x: np.ndarray) -> np.ndarray:
        _ensure(x)
        poles = cache["poles"]
        wells = cache["wells"]
        lam_s, alpha_s = cache["lam_s"], cache["alpha_s"]
        assert isinstance(poles, list) and isinstance(wells, np.ndarray)
        assert isinstance(lam_s, SmoothR) and isinstance(alpha_s, SmoothR)
        r_pts = pair.grid_a.points
        n_res = int(resonant_mask.sum())
        out = np.zeros((R_nodes.size + n_res, x.size))
        k = -1
        for j, pole in enumerate(poles):
            if resonant_mask[j]:
                k += 1
            if pole is None:
                continue  # constant penalty -> zero row
            Rj = float(R_nodes[j])
            lam_j, alpha_j = float(wells[j, 0]), float(wells[j, 1])
            extra = extra_of_R(Rj) if extra_of_R is not None else None
            Ha, _ = pair.hamiltonians(well_potential(ell, lam_j, alpha_j, extra))
            try:
                s = pole_sensitivity(Ha, pole.energy)
            except ValueError:
                continue  # c-product self-orthogonal: no usable sensitivity here
            g = np.exp(-alpha_j * r_pts**2)
            dE_dlam = complex(np.sum(s * (-g)))
            dE_dalpha = complex(np.sum(s * (lam_j * r_pts**2 * g)))
            row = np.concatenate(
                [dE_dlam * _smooth_grad(lam_s, Rj), dE_dalpha * _smooth_grad(alpha_s, Rj)]
            )
            out[j] = row.real
            if resonant_mask[j]:
                out[R_nodes.size + k] = sqrt_w_im[k] * row.imag
        return out

    sol = least_squares(
        residual, x0, jac=jac, bounds=(lo, hi), method="trf", x_scale="jac", max_nfev=max_nfev
    )
    r_final = residual(sol.x)  # also refreshes `last_pole`/the cache to the FINAL trial
    final_poles = cache["poles"]
    assert isinstance(final_poles, list)
    n_skip = sum(p is None for p in final_poles)
    kept = np.array([p is not None for p in final_poles])
    kept_res = kept[resonant_mask] if int(resonant_mask.sum()) else np.zeros(0, dtype=bool)
    found = np.concatenate([kept, kept_res])
    r_found = r_final[found]
    polish_rms = float(np.sqrt(np.mean(r_found**2))) if r_found.size else float("nan")
    fitted = replace(model, lam=decode_lam(sol.x[:n_lam]), alpha=decode_alpha(sol.x[n_lam:]))
    detail = (
        f"polish_rms={polish_rms:.2e} Ha (analytic Jacobian); "
        f"polish skipped {n_skip}/{R_nodes.size} nodes; "
        f"least_squares status={sol.status} nfev={sol.nfev}/{max_nfev}"
    )
    return fitted, detail


def _t1_nodes(target: ResonanceTarget, n_nodes: int) -> tuple[np.ndarray, str]:
    """The descending `R` grid `fit_resonance` tracks and re-verifies on, plus
    a short label for `TierResult.detail`.

    When `target.v_ion` came from a TABLE, use the table's own nodes (inside
    `R_range`, evenly subsampled to at most `n_nodes`) rather than an
    independent `linspace`. `Curve.from_table` is a `CubicSpline` THROUGH
    those nodes, so on them the target is exact; between them it is a
    polynomial interpolant, and `E_res(R)` has a BRANCH POINT at the
    resonance threshold `R_c`, which no polynomial interpolates accurately.
    Measured on F2 by evaluating the re-verify residual at the PUBLISHED
    parameters -- the floor no fit can beat -- an off-node grid put 1.39e-3 /
    1.67e-3 Ha of spurious target error into the residual, 140-170% of the
    `e_res_rms = 1e-3` tolerance, BEFORE any fitting happened. `fit_neutral`
    (T0) and `_coupling_eval_grid` (T3) already read their targets on-node
    for the same reason; this is the third.

    Falls back to `linspace` for a target built from callables, which has no
    nodes and no interpolation error to avoid.
    """
    lo, hi = target.R_range
    nodes = target.v_ion.x
    if nodes is None:
        return np.linspace(hi, lo, n_nodes), f"linspace {n_nodes}"
    sel = nodes[(nodes >= lo) & (nodes <= hi)][::-1]  # descending
    if sel.size > n_nodes:
        sel = sel[np.linspace(0, sel.size - 1, n_nodes).astype(int)]
    return np.asarray(sel, dtype=np.float64), f"table {sel.size}/{nodes.size}"


def fit_resonance(
    target: ResonanceTarget,
    model: FlexibleDiatomicModel,
    *,
    pair: ElectronicPair,
    n_nodes: int = 24,
    tol: Tolerances,
    lam_coeffs: int = 0,
    alpha_coeffs: int = 0,
    polish_nfev: int = 100,
) -> tuple[FlexibleDiatomicModel, TierResult]:
    R_desc, node_source = _t1_nodes(target, n_nodes)
    # Every denominator below counts the nodes ACTUALLY used (a table target
    # may hold fewer than `n_nodes`), not the requested `n_nodes`.
    n_used = R_desc.size
    model = replace(
        model,
        lam=_grow_coeffs(model.lam, lam_coeffs),
        alpha=_grow_coeffs(model.alpha, alpha_coeffs),
    )

    def pole_target(R: float) -> complex:
        # The target curve is NaN outside its table (or in a gap);
        # a node with no target is "no target here", not an error -- let
        # `track_curve` flag-and-freeze it via a NaN real part.
        s = float(target.v_ion(R) - model.v0(R).real)
        g = float(target.gamma(R))
        if np.isnan(s) or np.isnan(g):
            return complex(np.nan, np.nan)
        # Classify BOUND iff gamma_target < tol.gamma_floor (2e-3 Ha
        # default), not some far smaller threshold.
        # gamma_floor is already the residual mask's own noise floor below --
        # a node whose true width is under it contributes only its E_res on
        # either branch, so solving it as a real (bound) target rather than
        # asking Newton for a resonant pole with Im ~ 1e-5..1e-3 (numerically
        # negligible, effectively noise) is both cheaper and consistent. F2's
        # first-pass Newton failures land EXACTLY in the old gap between
        # 1e-6 and gamma_floor (gamma_target = 1.05e-5, 1.26e-4, 1.12e-3).
        if g < tol.gamma_floor:
            # A width under the floor with E_res ABOVE the neutral is the crossing
            # slice: neither a bound state (no real eigenvalue exists there) nor a
            # resolvable resonance (its Im is below what the gate can certify) --
            # no target at this node, exactly as a NaN table point.
            return complex(s, 0.0) if s < 0.0 else complex(np.nan, np.nan)
        return complex(s, -0.5 * g)

    # Tracking and re-verify must see the same
    # potential. `walk_t1` re-verifies against `model.surface`, which already
    # includes the shell/barrier term when present; without `extra_of_R` here
    # `track_curve` would tune (lam, alpha) against the bare well alone.
    extra_of_R = _shell_extra_of_R(model)

    seed = WellParams(float(model.lam_R(R_desc[0]).real), float(model.alpha_R(R_desc[0]).real))
    tr = track_curve(pair, model.ell, R_desc, pole_target, seed, extra_of_R=extra_of_R)
    ok = tr.converged
    if int(ok.sum()) < 4:
        return model, TierResult(
            "T1",
            "not met",
            np.inf,
            np.inf,
            f"nodes: {node_source}; only {int(ok.sum())}/{n_used} tracked",
        )

    # `solve_pole_params` (tracker.py) moves ONLY `lam` for a bound target --
    # a single real bound-state energy is one equation in `lam` alone, and
    # `alpha` is degenerate (held at whatever the seed carried). Over the
    # bound branch (R >= the crossing) every tracked `alpha` is therefore
    # just the seed's `alpha`, not physically informative -- fitting the
    # "alpha" SmoothR against ALL converged nodes lets those uninformative
    # points drag the curve away from the true (here R-independent) alpha(R),
    # which then also corrupts the paired `lam` at those same nodes (Newton
    # compensates the wrong `alpha` with the wrong `lam`). Fit "alpha" from
    # the RESONANT nodes only, where the full 2-D Newton solve pins down both
    # parameters uniquely; since alpha(R) is one global smooth functional
    # form (not a per-node table), the fit correctly extrapolates into the
    # bound branch too. Then re-track with that corrected alpha(R) as the new
    # seed at EVERY node (`alpha_of_R`), not just the first, so the
    # bound branch's frozen alpha is the accurate one and the re-tracked
    # `lam` there is no longer contaminated.
    gamma_t = np.array([float(target.gamma(Rj)) for Rj in tr.R])
    # This mask must match `pole_target`'s own bound/resonant
    # boundary (`tol.gamma_floor`, not the old 1e-6) -- a node Newton solved
    # as BOUND (alpha frozen at the seed, not physically informative) must
    # not be classified "resonant" here just because its noise-floor gamma
    # happens to clear a looser threshold.
    resonant = ok & ~np.isnan(gamma_t) & (gamma_t >= tol.gamma_floor)
    # The minimum resonant-node count needed
    # before trusting the resonant-only subset depends on how many DEGREES OF
    # FREEDOM `_fit_smooth` actually has for this SmoothR. When the seed's
    # `alpha.f_0 == 0` (every published molecule here -- alpha(R) is
    # R-independent), `_fit_smooth` now holds f_0/f_1/R_f fixed and fits ONLY
    # f_inf (+ coeffs): a single clean resonant point is informative there,
    # so diluting it with uninformative frozen bound-branch points (see the
    # comment above) is worse, not better -- confirmed on F2 (only 2/10
    # nodes resonant): falling back to all 7 converged nodes pulled the
    # fitted alpha to the MEAN of 5 uninformative points near the seed's
    # (wrong) 3.9 and 2 informative points near the true 3.0, landing at
    # ~3.64 (22% off) and cascading into 6/10 nodes surviving re-tracking.
    # The general sigmoid case (f_0 != 0) keeps the original >=4 threshold
    # (4 shape parameters need >=4 points to be well-posed).
    n_free_alpha = (1 if model.alpha.f_0 == 0.0 else 4) + len(model.alpha.coeffs)
    alpha_nodes = resonant if int(resonant.sum()) >= n_free_alpha else ok
    # Weight the alpha fit by gamma_target(R_j), normalized to
    # max 1 -- `d(chi^2)/d(alpha)` vanishes as Gamma -> 0 at fixed E_res
    # residual, so a barely-open node's tracked alpha is intrinsically
    # noisier than a fully-open node's and must not pull the fit as hard.
    # Floored at 1e-6 (not literally 0) so a node is heavily down-weighted
    # rather than dropped outright -- in the degenerate (single-parameter)
    # branch a zero-weighted point simply vanishes from the objective,
    # which is fine when other points remain, but a hard floor keeps this
    # robust if it's ever the ONLY node available.
    alpha_w = np.maximum(gamma_t[alpha_nodes], 0.0)
    alpha_w = np.maximum(alpha_w / max(float(np.max(alpha_w)), 1e-300), 1e-6)
    model, alpha_rel_rms = _fit_smooth(
        model, "alpha", tr.R[alpha_nodes], tr.alpha[alpha_nodes], weights=alpha_w
    )

    seed2 = WellParams(float(model.lam_R(R_desc[0]).real), float(model.alpha_R(R_desc[0]).real))
    tr2 = track_curve(
        pair,
        model.ell,
        R_desc,
        pole_target,
        seed2,
        extra_of_R=extra_of_R,
        alpha_of_R=lambda R: float(model.alpha_R(R).real),
    )
    ok2 = tr2.converged
    if int(ok2.sum()) < 4:
        return model, TierResult(
            "T1",
            "not met",
            np.inf,
            np.inf,
            f"nodes: {node_source}; only {int(ok2.sum())}/{n_used} nodes re-tracked "
            f"with alpha(R) corrected",
        )
    model, lam_rel_rms = _fit_smooth(model, "lam", tr2.R[ok2], tr2.lam[ok2])

    # The electron-affinity asymptote goes in BEFORE the polish, and then INTO
    # it as one more (bound) pole node at R_inf: applied afterwards it would
    # move `lam.f_inf`, the sigmoid's asymptote, and with it lam(R) at every
    # R -- undoing the polish (measured on O2: a 1 eV drift of V_ion on the
    # bound side). Exactly one of the three statuses describes what happened.
    reaches_asymptote = (
        abs(float(target.v_ion(R_desc[0]) - model.v0(R_desc[0]).real) + target.ea) < 0.05
    )
    ea_status = "skipped (table short of asymptote)"
    ea_node = False
    if reaches_asymptote:
        model, applied = _apply_ea_constraint(model, pair, target.ea)
        ea_status = "applied, held by the polish" if applied else "skipped (no bracket)"
        ea_node = applied

    # Joint polish of lam(R)/alpha(R) TOGETHER
    # against the poles themselves (not the per-node tracked sample -- see
    # `_joint_polish`'s docstring), on exactly the second-pass converged
    # nodes, using the SAME resonant/gamma_target weighting as the alpha fit
    # above. `tr2.poles` is index-aligned with `tr2.R`/`ok2`.
    R_pol = tr2.R[ok2]
    gamma_t2 = np.array([float(target.gamma(Rj)) for Rj in R_pol])
    resonant2 = gamma_t2 >= tol.gamma_floor
    seed_poles = [p for p, keep in zip(tr2.poles, ok2, strict=True) if keep]
    polish_target = pole_target
    if ea_node:
        extra_inf = extra_of_R(_R_INF) if extra_of_R is not None else None
        p_inf = pair.pole(
            well_potential(
                model.ell,
                float(model.lam_R(_R_INF).real),
                float(model.alpha_R(_R_INF).real),
                extra_inf,
            ),
            DEFAULT_BOUND_WINDOW,
        )
        if p_inf is not None:
            R_pol = np.concatenate([[_R_INF], R_pol])
            gamma_t2 = np.concatenate([[0.0], gamma_t2])
            resonant2 = np.concatenate([[False], resonant2])
            seed_poles = [p_inf, *seed_poles]
            ea = float(target.ea)

            def polish_target(R: float, _ea: float = ea) -> complex:
                return complex(-_ea, 0.0) if R == _R_INF else pole_target(R)

    model, polish_detail = _joint_polish(
        model,
        pair,
        R_pol,
        polish_target,
        gamma_t2,
        resonant2,
        seed_poles,
        extra_of_R,
        max_nfev=polish_nfev,
    )

    # Re-verify on the SMOOTHED model with the package's own gated per-node
    # walk: `qscat.core.lcp.resonance_pole_walk` freezes at the
    # crossing on this grid and silently zeroes Gamma over the whole
    # resonant region (measured). `walk_t1` drops (not freezes) a node
    # with no gated pole, exactly like the tracking step above.
    try:
        R_w, shift_w, gamma_w = walk_t1(
            model, pair, R_desc, seed_energy=pole_target(float(R_desc[0]))
        )
    except ValueError as e:
        return model, TierResult(
            "T1",
            "not met",
            np.inf,
            np.inf,
            f"nodes: {node_source}; "
            f"tracked {int(ok.sum())}/{n_used} ({int(resonant.sum())} resonant), "
            f"re-tracked {int(ok2.sum())}/{n_used} nodes; re-walk failed: {e}",
        )

    v_t = target.v_ion(R_w)
    g_t = target.gamma(R_w)
    # Residuals only on nodes where BOTH the target curve is
    # defined (not NaN) and the re-walk found a gated pole (already true of
    # every R_w, since walk_t1 only returns survivors).
    valid = ~np.isnan(v_t) & ~np.isnan(g_t)
    n_valid = int(valid.sum())
    if n_valid < 4:
        return model, TierResult(
            "T1",
            "not met",
            np.inf,
            np.inf,
            f"nodes: {node_source}; "
            f"tracked {int(ok.sum())}/{n_used} ({int(resonant.sum())} resonant), "
            f"re-tracked {int(ok2.sum())}/{n_used} nodes; re-walk verified only "
            f"{n_valid}/{n_used} nodes against the target",
        )

    e_err = (model.v0(R_w[valid]).real + shift_w[valid]) - v_t[valid]
    g_target = g_t[valid]
    g_gate = g_target > tol.gamma_floor
    g_rel = (
        (gamma_w[valid][g_gate] - g_target[g_gate]) / g_target[g_gate]
        if g_gate.any()
        else np.zeros(0)
    )

    e_rms, e_max = float(np.sqrt(np.mean(e_err**2))), float(np.max(np.abs(e_err)))
    g_rms = float(np.sqrt(np.mean(g_rel**2))) if g_rel.size else 0.0
    g_max = float(np.max(np.abs(g_rel))) if g_rel.size else 0.0
    # "met" also requires COVERAGE, not just accuracy where nodes
    # survived -- a curve verified on a shrinking sliver of the requested
    # grid (nodes silently gated out along the way) must not read as "met"
    # just because the sliver that's left agrees well.
    coverage_ok = int(ok2.sum()) >= 0.75 * n_used and n_valid >= 0.75 * n_used
    # "met" is decided by the RE-VERIFY
    # residuals + coverage alone -- the joint polish above puts the loss on
    # the poles themselves, so `smooth_rms_lam` (the per-node-sample
    # residual) no longer needs to be independently gated; it is still
    # reported in `detail` as a diagnostic, alongside `polish_rms`/status.
    met = e_rms <= tol.e_res_rms and g_max <= tol.gamma_rel and coverage_ok
    detail = (
        f"nodes: {node_source}; "
        f"tracked {int(ok.sum())}/{n_used} ({int(resonant.sum())} resonant), "
        f"re-tracked {int(ok2.sum())}/{n_used} nodes; "
        f"re-walk verified {n_valid}/{n_used} nodes; "
        f"E_res rms={e_rms:.2e} max={e_max:.2e} Ha; Gamma rel rms={g_rms:.3f} max={g_max:.3f}; "
        f"smooth_rms_lam={lam_rel_rms:.2e} smooth_rms_alpha={alpha_rel_rms:.2e}; "
        f"{polish_detail}; "
        f"EA constraint {ea_status}"
    )
    return model, TierResult("T1", "met" if met else "not met", e_rms, max(e_max, g_max), detail)


def _coupling_eval_grid(
    target: CouplingTarget, n_eps: int, n_R: int
) -> tuple[np.ndarray, np.ndarray, str]:
    """The `(eps, R)` grid `fit_coupling`'s residual is evaluated on, plus a
    short label saying which of the two grids was used (for `TierResult.detail`).

    When `target` was built `from_table`, its `gamma_tilde` is a log-log
    interpolant that is exact to round-off ON its own build nodes but has
    real curvature error BETWEEN them: `Gamma_tilde(eps, R) ~ eps^a A(R)
    exp(-B(R) eps)`'s `exp(-B(R) eps)` factor is not linear in log(eps), so
    linear-in-log(eps) interpolation between build nodes picks up spurious
    curvature error that would otherwise leak into the fit loss as if it
    were a real model/target mismatch. Evaluating on an even subsample of
    the table's own nodes instead keeps that interpolation error out of the
    loss entirely. Falls back to a geomspace/linspace grid for a target with
    no nodes (`from_alt_houfek`, an analytic form with no interpolation
    error to avoid).
    """
    if target.eps_nodes is None or target.R_nodes is None:
        eps = np.geomspace(target.eps_window[0], target.eps_window[1], n_eps)
        R_asc = np.linspace(target.R_range[0], target.R_range[1], n_R)
        return eps, R_asc, "geomspace/linspace"

    def subsample(nodes: np.ndarray, window: tuple[float, float], k: int) -> np.ndarray:
        sel = nodes[(nodes >= window[0]) & (nodes <= window[1])]
        if sel.size <= k:
            return sel
        idx = np.linspace(0, sel.size - 1, k).astype(int)
        return sel[idx]

    eps = subsample(target.eps_nodes, target.eps_window, n_eps)
    R_asc = subsample(target.R_nodes, target.R_range, n_R)
    return eps, R_asc, "table nodes"


def model_gamma_tilde(
    model: FlexibleDiatomicModel, pair: ElectronicPair, eps: np.ndarray, R_asc: np.ndarray
) -> np.ndarray:
    """The model's own 2 pi |V_dk+(eps, R)|^2 with the R-independent discrete state."""
    phi_d = AsymptoticDiscreteState(pair.grid_a, model, _R_INF)
    out = np.empty((eps.size, R_asc.size))
    for i, e in enumerate(eps):
        out[i] = gamma_from_coupling(v_dk_plus(pair.grid_a, model, phi_d, R_asc, float(e)))
    return out


def _threshold_exponent_mismatch(
    target: CouplingTarget, model: FlexibleDiatomicModel
) -> str | None:
    """Why `target`'s near-threshold power law is out of this ansatz's scope,
    or `None` if it is not.

    The exponent is a PROPERTY of the ansatz, never a fitted parameter: `ell`
    is fixed, so `Gamma~ ~ eps^(l+1/2)` is exact as `eps -> 0` and no shell
    term can change it. A target carrying a different power -- a POLAR
    molecule, where the long-range dipole gives `alpha = sqrt(d + 1/4)`
    instead of `l + 1/2` -- is not describable here at all, so T3 must SAY so
    rather than fit a shell against a law the model cannot produce and report
    whatever residual that leaves.
    """
    want = model.ell + 0.5
    if abs(target.alpha_exponent - want) <= 1e-9:
        return None
    return (
        f"threshold exponent alpha={target.alpha_exponent:.6g} differs from the ansatz's "
        f"l+1/2={want:.6g}: polar/dipole targets are out of scope for this ansatz"
    )


def fit_coupling(
    target: CouplingTarget,
    model: FlexibleDiatomicModel,
    *,
    pair: ElectronicPair,
    n_eps: int = 8,
    n_R: int = 8,
    tol: Tolerances,
    r_b: float = 3.0,
    alpha_b: float = 2.0,
) -> tuple[FlexibleDiatomicModel, TierResult]:
    out_of_scope = _threshold_exponent_mismatch(target, model)
    if out_of_scope is not None:
        return model, TierResult("T3", "not met", np.nan, np.nan, out_of_scope)
    eps, R_asc, grid_source = _coupling_eval_grid(target, n_eps, n_R)
    y = np.log(target.gamma_tilde(eps[:, None], R_asc[None, :]))
    grid_detail = f"log Gamma~ on {grid_source} {eps.size}x{R_asc.size} (eps, R)"
    orig_model = model
    if model.shell is None:
        R_mid = 0.5 * (target.R_range[0] + target.R_range[1])
        shell = SmoothR(f_inf=0.0, f_0=0.0, f_1=1.0, R_f=R_mid, R_e=model.R_e)
        model = model.with_shell(shell, alpha_b, r_b)
    names = ["shell.f_inf", "shell.f_0", "r_b"]

    def resid(x):
        m = unpack(model, names, x)
        try:
            g = model_gamma_tilde(m, pair, eps, R_asc)
        except ValueError:
            # a trial shell that unbinds the anion at R_inf has no discrete
            # state; a large residual sends the optimizer back, no exception
            return np.full(y.size, 10.0)
        return (np.log(np.maximum(g, 1e-300)) - y).ravel()

    x0 = pack(model, names)
    try:
        r0 = resid(x0)
    except ValueError as err:  # no bound anion at R_inf: T1 missed the asymptote
        return model, TierResult(
            "T3", "not met", np.nan, np.nan, f"no discrete state: {str(err)[:100]}"
        )
    rms0 = float(np.sqrt(np.mean(r0**2)))
    if rms0 <= 0.5 * tol.coupling_log_rms:
        # No fit was needed: hand back the CALLER's model, not the
        # zero-strength shell installed above only to probe `rms0` -- an
        # installed-but-unfit shell is not a change the caller asked for.
        return orig_model, TierResult(
            "T3",
            "met",
            rms0,
            float(np.max(np.abs(r0))),
            f"{grid_detail}; shell not needed; seed_rms={rms0:.3e}",
        )
    # The shell is a REPULSIVE barrier by design: amplitudes >= 0 keep it from
    # turning into a second well that binds a state of its own.
    lower, upper = [0.0, 0.0, 0.5], [5.0, 5.0, 12.0]
    x0 = np.clip(x0, np.array(lower) + 1e-9, np.array(upper) - 1e-9)
    sol = least_squares(
        resid, x0, bounds=(lower, upper), diff_step=1e-3, xtol=1e-8, ftol=1e-8, max_nfev=60
    )
    fitted = unpack(model, names, sol.x)
    r = resid(sol.x)
    rms, mx = float(np.sqrt(np.mean(r**2))), float(np.max(np.abs(r)))
    status = "met" if rms <= tol.coupling_log_rms else "not met"
    return fitted, TierResult(
        "T3", status, rms, mx, f"{grid_detail}; shell fitted; seed_rms={rms0:.3e}"
    )


def _crossing(target: Target, model: FlexibleDiatomicModel) -> float | None:
    """The sign change of `v_ion - v0` (the target's ionic curve, relative to
    the fitted neutral curve) over the T1 target's own R_range -- the R at
    which the anion curve crosses the neutral one. `None` if no crossing is
    found in range, or if the target carries no T1 data at all.
    """
    if target.resonance is None:
        return None
    R = np.linspace(target.resonance.R_range[0], target.resonance.R_range[1], 400)
    d = target.resonance.v_ion(R) - model.v0(R).real
    # The target's v_ion is NaN outside its own table (or in a gap);
    # np.sign(nan) is nan, and a nan participating in the product test would
    # itself compare False to 0 either way -- but guard explicitly so a NaN
    # neighbour never masquerades as "no sign change" via a stray true from
    # nan*finite comparisons.
    s = np.sign(d)
    valid = ~np.isnan(s[:-1]) & ~np.isnan(s[1:])
    prod = np.where(valid, s[:-1] * s[1:], 1.0)
    idx = np.flatnonzero(prod < 0)
    if idx.size == 0:
        return None
    i = int(idx[0])
    return float(R[i] - d[i] * (R[i + 1] - R[i]) / (d[i + 1] - d[i]))


def _da_sign(target: Target, model: FlexibleDiatomicModel) -> int | None:
    """`+1` endothermic (N2), `-1` exothermic (F2): the sign of the
    dissociative-attachment threshold energy, `(-ea) - eps_0`, measured from
    the fitted model's own v=0 vibrational level on the package's standard
    nuclear grid."""
    if target.resonance is None:
        return None
    g = nuclear_grid()
    eps, _ = vibrational_states(g, model.mu, 1, model.v0)
    threshold = (-target.resonance.ea) - eps[0]  # DA threshold energy from v=0
    return int(np.sign(threshold)) if threshold != 0.0 else 0


def _t1_recheck(
    target: Target,
    model: FlexibleDiatomicModel,
    pair: ElectronicPair,
    n_nodes: int,
    tol: Tolerances,
) -> tuple[bool, str]:
    """Re-walk the T1 nodes on the FINAL model (after T3's shell) and report
    whether T1's tolerances still hold -- the shell moves the poles, so a T3
    that breaks T1 is not a met tier."""
    assert target.resonance is not None
    R_desc, _ = _t1_nodes(target.resonance, n_nodes)
    s0 = float(target.resonance.v_ion(R_desc[0]) - model.v0(R_desc[0]).real)
    seed = complex(s0, -0.5 * float(target.resonance.gamma(R_desc[0])))
    try:
        R_w, shift_w, gamma_w = walk_t1(model, pair, R_desc, seed_energy=seed)
    except ValueError as err:
        return False, f"post-T3 T1 re-walk failed: {err}"
    v_t, g_t = target.resonance.v_ion(R_w), target.resonance.gamma(R_w)
    valid = ~np.isnan(v_t) & ~np.isnan(g_t)
    e_err = (model.v0(R_w[valid]).real + shift_w[valid]) - v_t[valid]
    mask = g_t[valid] > tol.gamma_floor
    g_rel = (
        np.abs(gamma_w[valid][mask] - g_t[valid][mask]) / g_t[valid][mask]
        if mask.any()
        else np.zeros(0)
    )
    e_rms = float(np.sqrt(np.mean(e_err**2))) if e_err.size else np.inf
    g_max = float(np.max(g_rel)) if g_rel.size else 0.0
    ok = e_rms <= tol.e_res_rms and g_max <= tol.gamma_rel and valid.sum() >= 0.75 * R_desc.size
    return ok, (
        f"post-T3 T1 re-walk {int(valid.sum())}/{R_desc.size} nodes: "
        f"E_res rms={e_rms:.2e} Ha, Gamma rel max={g_max:.3f}"
    )


def fit(
    target: Target,
    seed: FlexibleDiatomicModel,
    *,
    pair: ElectronicPair,
    tol: Tolerances = Tolerances(),  # noqa: B008 (frozen dataclass; safe as a default)
    n_beta: int = 1,
    n_nodes: int = 24,
    lam_coeffs: int = 0,
    alpha_coeffs: int = 0,
    continue_on_miss: bool = False,
    polish_nfev: int = 100,
) -> tuple[FlexibleDiatomicModel, FitReport]:
    """T0 -> T1 -> T3 in order over the tiers present in `target`, stopping
    after the first "not met" (later tiers report "not attempted" without
    running), then filling in the crossing radius and DA-threshold sign from
    the fitted model."""
    model = seed
    tiers: list[TierResult] = []
    halted = False
    stages = [
        (
            "T0",
            target.neutral,
            lambda m: fit_neutral(target.neutral, m, n_beta=n_beta, tol=tol),
        ),
        (
            "T1",
            target.resonance,
            lambda m: fit_resonance(
                target.resonance,
                m,
                pair=pair,
                n_nodes=n_nodes,
                tol=tol,
                lam_coeffs=lam_coeffs,
                alpha_coeffs=alpha_coeffs,
                polish_nfev=polish_nfev,
            ),
        ),
        ("T3", target.coupling, lambda m: fit_coupling(target.coupling, m, pair=pair, tol=tol)),
    ]
    for name, present, stage in stages:
        if present is None:
            tiers.append(TierResult(name, "not attempted", np.nan, np.nan, "no target data"))
            continue
        if halted:
            tiers.append(
                TierResult(name, "not attempted", np.nan, np.nan, "an earlier tier was not met")
            )
            continue
        model, res = stage(model)
        if name == "T3" and target.resonance is not None and res.status != "not attempted":
            ok_t1, note = _t1_recheck(target, model, pair, n_nodes, tol)
            res = TierResult(
                res.name,
                res.status if ok_t1 else "not met",
                res.rms,
                res.max,
                f"{res.detail}; {note}" + ("" if ok_t1 else " -- the shell breaks T1"),
            )
        tiers.append(res)
        halted = res.status == "not met" and not continue_on_miss
    # The nuclear probe tail: a pivot at R = 12 bohr with a straight ECS ray
    # off it. The angle is named rather than inlined so the tail and the
    # angle `ecs_bounded` reports as PROBED cannot drift apart.
    nuclear_deg = 35.0
    R_tail = 12.0 + np.linspace(0.1, 6.0, 8) * np.exp(1j * np.deg2rad(nuclear_deg))
    report = FitReport(
        target_name=target.name,
        parameters=params(model),
        tiers=tiers,
        ecs_bounds_deg=ecs_bounded(model, pair, R_tail, nuclear_deg=nuclear_deg),
        crossing_R=_crossing(target, model),
        da_threshold_sign=_da_sign(target, model),
        provenance={
            k: {"source": v.source, "locator": v.locator} for k, v in target.provenance.items()
        },
    )
    return model, report
