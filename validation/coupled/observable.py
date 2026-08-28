"""Turn the screen's curve difference into a cross-section difference, and
decide the gate.

The coupled and fixed-l curves are fed through `qscat.core.lcp`'s
VIBRATIONAL-EXCITATION route. VE and never DA: NO's dissociative-attachment
channel opens at +0.172 Ha, above the resonance at 0.02-0.05 Ha, so sigma_DA
is a 1e-19 bohr^2 tail on this sweep and the LCP is documented to miss it by
five to seven orders -- a quantity that wrong cannot discriminate a
few-percent change in Gamma(R).

The coupled curve reaches the LCP as a DIFFERENCE applied to the shipped
`local_complex_potential` output, interpolated from the screen's R sample and
zero outside it. That reuses the shipped tail handling, freezing and clamping
untouched, so the comparison cannot pick up a tail artefact -- at the price
of confining the measured effect to R in [1.6, 6.0], which is where the
resonance lives and where the screen sampled.

This module makes no decision of its own beyond the three criteria the spec
declared before the campaign ran. A SHUT gate is a result: it says the
fixed-l reduction is sound for a NO-like model over the full geometric range
of the anisotropy, and it must be reported as prominently as an open one.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import numpy.typing as npt
from qscat.core.grids import electronic_grid
from qscat.core.lcp import lcp_ve_cross_section, local_complex_potential
from qscat.core.vibrational import vibrational_states
from qscat.dvr.grid import FemDvrEcsGrid
from qscat.model import NO

from validation.coupled.screen import ANGLES, NO_ELECTRONIC, RESULTS, CoupledCurve
from validation.diatomic.config import CONFIGS

__all__ = [
    "E_SWEEP",
    "GAMMA_TOL",
    "PERTURBATION_MAX",
    "SIGMA_TOL",
    "gate_decision",
    "lcp_from_curve",
    "main",
]

# 41 energies across NO's resonance (0.02-0.10 Ha), well below the +0.172 Ha
# DA threshold.
E_SWEEP = np.arange(0.020, 0.1001, 0.002)
# Declared in the spec BEFORE the campaign ran. 5 % sits far above the curves'
# 1e-9..1e-7 Ha convergence floor and far below the factor-level departures the
# approximations already in production show.
GAMMA_TOL = 0.05
SIGMA_TOL = 0.05
# The difference-on-shipped construction is a PERTURBATION of the shipped curve
# and is only meaningful while it stays one. MEASURED across the campaign's own
# ladder, max|dGamma| as a fraction of max(Gamma_shipped): 0.01 at s = 0.1,
# 0.62 at s = 0.2, 1.88 at s = 0.3, 2.10 at s = 0.5 -- and the coupled sigma
# collapses with it, 1.6e2 -> 7.4e1 -> 6.6e0 -> 1.5e-5 against a fixed-l sigma
# that stays at 1.2e2. Past this fraction the relative shift pins to 1.0000
# because sigma_full has gone to ZERO, not because the cross section changed by
# 100 %: the construction reports its own failure as a physical effect. The
# threshold admits s = 0.1 and rejects s = 0.2 on this campaign.
PERTURBATION_MAX = 0.25
N_VIB = 4
V_INIT = 0
# Elastic and first inelastic -- the two channels the spec's criterion (c)
# names.
VPRIMES = [0, 1]


def lcp_from_curve(
    full: CoupledCurve,
    fixed: CoupledCurve,
    nuclear_grid: FemDvrEcsGrid,
    vd_shipped: npt.NDArray[np.complex128],
    gamma_shipped: npt.NDArray[np.float64],
) -> dict[str, npt.NDArray[np.float64]]:
    """VE cross sections for the coupled and the fixed-l curve.

    The coupled curve enters as a difference applied to the shipped
    `local_complex_potential` output; see the module docstring.
    """
    R_grid = np.asarray(nuclear_grid.real_points, dtype=np.float64)

    # No-target points are a normal feature of a resonance curve, and `np.interp`
    # would propagate their nan across the whole interpolated range and from
    # there into every cross section. Interpolate the difference over the points
    # where BOTH curves have a pole, and leave it zero elsewhere -- outside that
    # span the shipped LCP curve stands unmodified, which is the same
    # conservative choice already made beyond the sampled R.
    ok = (
        np.isfinite(full.v_d)
        & np.isfinite(fixed.v_d)
        & np.isfinite(full.gamma)
        & np.isfinite(fixed.gamma)
    )
    if int(ok.sum()) < 2:
        raise ValueError(
            f"only {int(ok.sum())} R have a pole in both curves; "
            "cannot interpolate the coupling difference"
        )
    R_ok = full.R[ok]
    inside = (R_grid >= R_ok[0]) & (R_grid <= R_ok[-1])

    d_vd = np.zeros_like(R_grid)
    d_gamma = np.zeros_like(R_grid)
    d_vd[inside] = np.interp(R_grid[inside], R_ok, (full.v_d - fixed.v_d)[ok])
    d_gamma[inside] = np.interp(R_grid[inside], R_ok, (full.gamma - fixed.gamma)[ok])

    eps, chi = vibrational_states(nuclear_grid, NO.mu, N_VIB, NO.v0)
    out: dict[str, npt.NDArray[np.float64]] = {}
    zero = np.zeros_like(d_vd)
    for label, dv, dg in (("fixed", zero, zero), ("full", d_vd, d_gamma)):
        sigma = lcp_ve_cross_section(
            nuclear_grid,
            NO.mu,
            np.asarray(vd_shipped + dv, dtype=np.complex128),
            np.asarray(np.maximum(0.0, gamma_shipped + dg), dtype=np.float64),
            eps,
            chi,
            V_INIT,
            VPRIMES,
            E_SWEEP,
        )
        out[label] = np.asarray(sigma, dtype=np.float64)
    return out


def gate_decision(summary: dict[str, float]) -> dict[str, object]:
    """Apply the three criteria the spec declared. Returns the verdict."""
    reasons: list[str] = []
    if summary["max_n_poles"] > 1:
        reasons.append("a second genuine pole entered the window")
    if summary["max_relative_gamma_shift"] > GAMMA_TOL:
        reasons.append(f"Gamma moved {summary['max_relative_gamma_shift']:.1%} > {GAMMA_TOL:.0%}")
    if summary["median_relative_sigma_shift"] > SIGMA_TOL:
        reasons.append(
            f"sigma_VE moved {summary['median_relative_sigma_shift']:.1%} "
            f"(median) > {SIGMA_TOL:.0%}"
        )
    if reasons:
        return {"open": True, "reason": "; ".join(reasons)}
    return {
        "open": False,
        "reason": (
            "no criterion met -- Phase 2 is deliberately not run. The fixed-l "
            "reduction is sound for a NO-like model over the full geometric "
            "range of the anisotropy."
        ),
    }


def _curve_from_payload(
    R: npt.NDArray[np.float64], payload: dict[str, list[float]]
) -> CoupledCurve:
    """Rebuild a `CoupledCurve` from its JSON payload."""
    v_d = np.asarray(payload["v_d"], dtype=np.float64)
    gamma = np.asarray(payload["gamma"], dtype=np.float64)
    return CoupledCurve(
        R=R,
        E_res=np.asarray(v_d - 0.5j * gamma, dtype=np.complex128),
        residual=np.asarray(payload["residual"], dtype=np.float64),
        n_stable=np.asarray(payload["n_stable"], dtype=np.intp),
        n_poles=np.asarray(payload["n_poles"], dtype=np.intp),
    )


def _perturbation_fraction(
    full: CoupledCurve,
    fixed: CoupledCurve,
    nuclear_grid: FemDvrEcsGrid,
    gamma_shipped: npt.NDArray[np.float64],
) -> float:
    """`max|dGamma| / max(Gamma_shipped)` -- how far the difference is from
    being a perturbation of the curve it is applied to.

    Above 1 the difference exceeds the width it modifies, `Gamma` clamps to
    zero across the doorway, and the resulting cross section is not a smaller
    cross section but no cross section at all.
    """
    R_grid = np.asarray(nuclear_grid.real_points, dtype=np.float64)
    ok = (
        np.isfinite(full.v_d)
        & np.isfinite(fixed.v_d)
        & np.isfinite(full.gamma)
        & np.isfinite(fixed.gamma)
    )
    if int(ok.sum()) < 2:
        return float("inf")
    R_ok = full.R[ok]
    inside = (R_grid >= R_ok[0]) & (R_grid <= R_ok[-1])
    d_gamma = np.zeros_like(R_grid)
    d_gamma[inside] = np.interp(R_grid[inside], R_ok, (full.gamma - fixed.gamma)[ok])
    scale = float(np.max(gamma_shipped))
    return float(np.max(np.abs(d_gamma)) / scale) if scale > 0 else float("inf")


def _summarize(report: dict) -> dict[str, float]:
    """Reduce the campaign JSON to the three numbers the gate consumes.

    The comparison point is `kappa = 0.5` at the largest `s` BOTH branches
    reached -- the biggest anisotropy at which the comparison is still matched.
    The walk stops where Gamma exceeds eps, and the full and fixed-l models need
    not stop together; comparing their two endpoints would silently compare
    different `s`. `full` is the widest channel set the campaign ran, `fixed` is
    the one-channel model, and both come from the SAME campaign on the same
    grids and R sample.
    """
    n_max = max(report["n_channels"])
    full_walk = report["kappa_curves"][str(n_max)]["0.5"]
    fixed_walk = report["kappa_curves"]["1"]["0.5"]
    shared = sorted(set(full_walk) & set(fixed_walk), key=float)
    if not shared:
        raise ValueError("full and fixed-l walks share no s value at kappa = 0.5")
    s_common = shared[-1]
    full_payload = full_walk[s_common]
    fixed_payload = fixed_walk[s_common]

    g_full = np.asarray(full_payload["gamma"], dtype=np.float64)
    g_fixed = np.asarray(fixed_payload["gamma"], dtype=np.float64)
    # `np.max` over a nan returns nan, and `nan > tol` is False -- a curve with
    # one no-target point would silently hold the gate SHUT. Compare only where
    # both curves have a pole, and refuse to report at all if none do, rather
    # than returning a number that means "no evidence" but reads as "no effect".
    both = np.isfinite(g_full) & np.isfinite(g_fixed)
    if not both.any():
        raise ValueError(
            f"no R has a pole in both the full and fixed-l curves at s = {s_common}; "
            "the gate cannot be evaluated"
        )
    rel_gamma = np.abs(g_full[both] - g_fixed[both]) / np.maximum(g_fixed[both], 1e-12)
    max_gamma = float(np.max(rel_gamma))

    # The two groups are nested to DIFFERENT depths and must be walked
    # separately: `s_curves[n_ch]` is {s: payload}, while `kappa_curves[n_ch]`
    # is {kappa: {s: payload}} -- the kappa sweep stores whole walks, because
    # full and fixed-l need not stop at the same s and the comparison has to be
    # made at a matched one. Treating them alike reads a dict where a payload
    # belongs.
    # `n_poles`, NOT `n_stable`. `n_stable` counts every angle-stable state and
    # the spurious near-threshold state is present at every R and every s, so it
    # is 2 everywhere -- a gate reading it would fire on every campaign that ever
    # runs. `n_poles` counts the states that pass the residual cut, which is what
    # "a second pole appeared" has to mean.
    n_poles = 1

    def _note(payload: dict[str, list[float]]) -> None:
        nonlocal n_poles
        n_poles = max(n_poles, int(np.max(payload["n_poles"])))

    # `s_curves[n_ch]` is {s: payload}; `kappa_curves[n_ch]` is
    # {kappa: {s: payload}}; `n_channels_5_check` is {s: payload}. Three
    # structures at two depths -- walk each at its own.
    for per_channel in report["s_curves"].values():
        for payload in per_channel.values():
            _note(payload)
    for per_channel in report["kappa_curves"].values():
        for per_kappa in per_channel.values():
            for payload in per_kappa.values():
                _note(payload)
    for payload in report["n_channels_5_check"].values():
        _note(payload)

    nuclear = CONFIGS["NO"].da_grid().grids[1]
    ga, gb = (electronic_grid(angle_deg=a, **NO_ELECTRONIC) for a in ANGLES)
    vd_shipped, gamma_shipped = local_complex_potential(NO, nuclear, ga, gb)

    R = np.asarray(report["R"], dtype=np.float64)

    # Criterion (c) is evaluated at the largest shared s where the construction
    # is still a perturbation -- which is NOT the largest shared s. Walk down
    # the shared ladder until the difference is small enough to trust.
    sigma_s, sigmas = None, None
    for s in reversed(shared):
        cf = _curve_from_payload(R, full_walk[s])
        cx = _curve_from_payload(R, fixed_walk[s])
        if _perturbation_fraction(cf, cx, nuclear, gamma_shipped) <= PERTURBATION_MAX:
            sigma_s = s
            sigmas = lcp_from_curve(cf, cx, nuclear, vd_shipped, gamma_shipped)
            break
    if sigmas is None:
        raise ValueError(
            "the curve difference exceeds "
            f"{PERTURBATION_MAX:.0%} of the shipped width at every shared s; "
            "the LCP route cannot evaluate criterion (c) for this campaign"
        )

    rel_sigma = np.abs(sigmas["full"] - sigmas["fixed"]) / np.maximum(sigmas["fixed"], 1e-30)
    if not np.isfinite(rel_sigma).any():
        raise ValueError("no finite sigma comparison; the gate cannot be evaluated")
    # `sigma(E)` is a resonance profile of width ~0.006 Ha on a 0.002 Ha energy
    # mesh -- three samples across a peak -- so a small pole shift puts the two
    # profiles locally out of step and a single flank sample can differ by 30x
    # while the curves as a whole differ by a third. A gate reading that
    # pointwise MAXIMUM measures the mesh, not the coupling, and would fire as
    # soon as a peak moved at all, by however little (this repository has paid
    # for that lesson before, on an O2 mesh that read every peak height at
    # 0.69 of its converged value). The MEDIAN is the criterion; the max is
    # kept as a diagnostic so the record still shows the peaks moved.
    median_sigma = float(np.nanmedian(rel_sigma))
    max_sigma = float(np.nanmax(rel_sigma))

    return {
        "max_relative_gamma_shift": max_gamma,
        "median_relative_sigma_shift": median_sigma,
        # Diagnostic only -- mesh-dominated, see the comment above. Recorded so
        # the record shows the peaks moved, not just that the curves differ.
        "max_relative_sigma_shift": max_sigma,
        "max_n_poles": float(n_poles),
        # The two criteria are evaluated at DIFFERENT s and the report must say
        # so: (b) compares two computed curves directly and is sound wherever
        # they exist; (c) rides a construction that is only valid while the
        # difference stays a perturbation.
        "gamma_s": float(s_common),
        "sigma_s": float(sigma_s),
    }


def main(results: Path = RESULTS) -> dict[str, object]:
    """Read the screen report, build the observable, decide, and write it.

    `results` defaults to the committed `RESULTS` directory (reading
    `screen.json` and writing `gate.json` there); a caller that does not
    want to overwrite the tracked gate decision (a test, in particular)
    should pass a scratch directory instead.
    """
    report = json.loads((results / "screen.json").read_text())
    summary = _summarize(report)
    verdict = gate_decision(summary)
    out = {"summary": summary, "verdict": verdict}
    (results / "gate.json").write_text(json.dumps(out, indent=1))
    state = "OPEN" if verdict["open"] else "SHUT"
    print(f"[coupled] gate {state}: {verdict['reason']}")
    return out


if __name__ == "__main__":
    main()
