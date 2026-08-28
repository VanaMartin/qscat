"""The fixed-R electronic screen: E_res(R) and Gamma(R) with the partial
waves coupled, against the same quantities with the fixed-l reduction.

The pole at each R is located by ECS ANGLE STABILITY -- two spectra of the
same coupled Hamiltonian at two electronic ECS angles, matched by
`qscat.ecs.match_angle_stable`. The multi-state matcher rather than
`find_resonance_pole` on purpose: the screen exists to notice when the single
resonance becomes more than one, so the count of stable states in the window
is recorded alongside the pole.

Both branches of the comparison -- "full" (n_channels = N_l) and "fixed-l"
(n_channels = 1) -- run through this one function on the same grids and the
same R sample, so the difference between them is the coupling and nothing
else. In particular the monopole shift from moving the wells apart, which is
NOT coupling, is present in both and cancels.
"""

from __future__ import annotations

import json
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import numpy.typing as npt
from qscat.core.grids import electronic_grid
from qscat.dvr.grid import FemDvrEcsGrid
from qscat.dvr.operators import eigen
from qscat.ecs import match_angle_stable
from qscat.model import NO

from projects.no_coupled_channels.anisotropy import TwoCentreWell
from projects.no_coupled_channels.model import CoupledModel

__all__ = [
    "ANGLES",
    "KAPPA_REFERENCE",
    "KAPPA_VALUES",
    "NO_ELECTRONIC",
    "N_CHANNEL_VALUES",
    "RESULTS",
    "R_SAMPLE",
    "S_VALUES",
    "CoupledCurve",
    "coupled_resonance_curve",
    "main",
    "run_continuation",
]

# NOT the eMoScat NO electronic deck, and not the (35, 44) angle pair the
# published curves use. Both were measured to be inadequate HERE: the
# anisotropy broadens the resonance by an order of magnitude over the
# continuation, and a 6-element tail cannot represent it. At the same physical
# point (s = 0.4, R = 2.0) the two-angle residual is 6e-4 on the published deck
# and 3.0e-9 on this one -- five orders, with no physics changed. The published
# deck loses the pole at s = 0.5; this one follows it to s = 1.
#
# The angle exceeds 45 degrees, which is unusual and is safe HERE for a reason
# worth stating: the constraint is that the Gaussian interaction stay bounded on
# the contour, i.e. Re(z^2) >= 0 along it, and that is a JOINT condition on angle
# and tail extent. This contour keeps min Re(z^2) = 0. Pushing the tail out
# instead (50 degrees, n_complex=12, tail_alpha=0.5, |z| = 1111) overflows the
# potential to 1.75e259.
NO_ELECTRONIC = {"r_max": 16.0, "order": 8, "n_complex": 8}
ANGLES = (44.0, 52.0)


@dataclass(frozen=True)
class CoupledCurve:
    """A resonance curve, plus the evidence needed to trust it."""

    R: npt.NDArray[np.float64]
    E_res: npt.NDArray[np.complex128]
    residual: npt.NDArray[np.float64]
    n_stable: npt.NDArray[np.intp]
    # Angle-stable states that also pass the residual cut. `n_stable` counts
    # every angle-stable state and is DIAGNOSTIC only: measured, the spurious
    # near-threshold state is present at every R and every s, so `n_stable` is
    # 2 everywhere and can never indicate a split. `n_poles` is what a gate can
    # use -- it counts states that could actually be poles.
    n_poles: npt.NDArray[np.intp]

    @property
    def v_d(self) -> npt.NDArray[np.float64]:
        """`Re E_res(R)` -- the LCP's real curve."""
        return np.asarray(self.E_res.real, dtype=np.float64)

    @property
    def gamma(self) -> npt.NDArray[np.float64]:
        """`Gamma(R) = max(0, -2 Im E_res(R))`, clamped as the LCP clamps it."""
        return np.asarray(np.maximum(0.0, -2.0 * self.E_res.imag), dtype=np.float64)


def coupled_resonance_curve(
    model: CoupledModel,
    R_values: npt.ArrayLike,
    grid_a: FemDvrEcsGrid,
    grid_b: FemDvrEcsGrid,
    *,
    seeds: npt.ArrayLike,
    half_width: float = 0.15,
    rel_tol: float = 1e-4,
    atol: float = 1e-3,
    resid_max: float = 1e-5,
) -> CoupledCurve:
    """`E_res(R)` for the coupled electronic problem, by two-angle stability.

    `seeds[i]` centres the search window at `R_values[i]`; the pole returned is
    the angle-stable state nearest that seed. A seed is a WINDOW CENTRE, not a
    result: `test_the_pole_is_insensitive_to_the_seed` locks that.

    `atol = 1e-3` matches the `resid_tol` the shipped
    `qscat.core.lcp.local_complex_potential` walks with, because the s = 0 gate
    demands agreement with that walk and an acceptance tighter than the
    reference's own cannot reproduce it. `match_angle_stable` accepts at
    `max(rel_tol*|E|, atol)`, so at these energies (|E| ~ 0.2 Ha) the absolute
    floor is what binds. CONSEQUENCE: `n_stable` counts states accepted at that
    threshold, which is loose enough to admit a near-degenerate discretised
    continuum state; a count above 1 is a candidate split, not a proven one,
    and must be read together with `residual` before it is believed.

    `resid_max` is what makes the seed a safe DISAMBIGUATOR rather than a
    SELECTOR. A residual IDENTIFIES a state as a pole; proximity to the seed
    only chooses between poles once identity is settled -- so the residual
    cut is applied FIRST, and the nearest-to-seed pick runs only over the
    survivors. A spurious near-threshold state (`eps` ~ +0.001, `Gamma` ~
    0.006, residual ~7e-4) can sit closer to a seed than the genuine pole
    (residual ~1e-9) while still failing the residual cut; picking nearest
    first and guarding afterwards then reports NO TARGET at a point where the
    genuine pole was present in the window and merely outranked by distance,
    rather than recovering it. Filtering first fixes exactly that case:
    measured on the NO screen, two points that a nearest-first order returned
    as `nan` recover the genuine pole at residual 6-8e-9 once the order is
    reversed. If NO candidate in the window survives the residual cut, the
    point genuinely has no target: `E_res` is `nan`, following the
    repository's existing convention for a slice with no resonance, and
    `residual` records the smallest residual seen (not the seed-nearest one)
    so a caller can tell "nothing here" from "something here that failed the
    cut".

    `resid_max = 1e-5` rather than something tighter because a pole crossing
    from bound to resonant is marginal by nature and its residual rises with it:
    measured across the full committed campaign (all 6692 points where a pole
    was actually recorded, over `s_curves`, `kappa_curves` and
    `n_channels_5_check`), the worst GENUINE residual accepted is 9.8e-6,
    right at that crossing -- within about 2% of the cut itself, so the
    accepted population reaches essentially to `resid_max` and a pole only
    marginally broader than the ones accepted here would have been recorded
    as no-target rather than as a pole. A 1e-6 cut punches holes in the curve
    exactly where the physics is most interesting, while 1e-5 still clears the
    artefact's 7e-4 by two orders.

    `half_width = 0.15`, wider than a single-`R` measurement might suggest,
    because `eps_res(R) = E_res(R) - v0(R)` is not close to constant over the
    `R` range this screen walks -- on the NO screen it dips to ~0.035 Ha near
    R=2.6 and rises to ~0.088 Ha at both ends of a typical sample (R=2.0 and
    R=4.0). A seed built from one part of that trend (e.g. `_seeds` in
    `test_screen.py`, which tracks only `v0(R)` and adds a constant offset)
    can miss the pole at another `R` entirely inside a narrower window; 0.15
    keeps it inside the window at every `R` sampled while still excluding the
    bulk of the discretised continuum.
    """
    R_arr = np.asarray(R_values, dtype=np.float64)
    seed_arr = np.asarray(seeds, dtype=np.complex128)
    if seed_arr.shape != R_arr.shape:
        raise ValueError(f"seeds has shape {seed_arr.shape}, expected {R_arr.shape}")

    poles = np.empty(R_arr.size, dtype=np.complex128)
    resids = np.empty(R_arr.size, dtype=np.float64)
    counts = np.empty(R_arr.size, dtype=np.intp)
    survivors = np.empty(R_arr.size, dtype=np.intp)

    for i, (R, seed) in enumerate(zip(R_arr, seed_arr, strict=True)):
        ea, _ = eigen(model.electronic_hamiltonian(grid_a, complex(R)).toarray())
        eb, _ = eigen(model.electronic_hamiltonian(grid_b, complex(R)).toarray())
        window = (
            seed.real - half_width,
            seed.real + half_width,
            seed.imag - half_width,
            seed.imag + half_width,
        )
        energies, residuals, _idx = match_angle_stable(ea, eb, window, rel_tol=rel_tol, atol=atol)
        if energies.size == 0:
            # NOT an error. `match_angle_stable` documents an empty result as a
            # normal outcome -- eigenvalues near the window, none of them
            # angle-stable -- and that is precisely the "no resonance at this R"
            # case, which belongs in the curve as a no-target point exactly like
            # a failed residual cut two branches below. Raising here would take
            # the whole curve down over one bad R, and the campaign that consumes
            # this reads `nan` as its stop signal. (`match_angle_stable` still
            # raises on its own account when the window catches NOTHING in one of
            # the two spectra; that is a seeding failure, not a physics result,
            # and is deliberately left to propagate.)
            poles[i] = np.nan + 1j * np.nan
            resids[i] = np.inf
            counts[i] = 0
            survivors[i] = 0
            continue
        counts[i] = energies.size
        keep = residuals <= resid_max
        survivors[i] = int(keep.sum())
        if not keep.any():
            # No target: nothing in the window is a pole (see `resid_max`).
            # Record the best residual seen, so the campaign can tell "nothing
            # here" from "something here that failed the cut".
            poles[i] = np.nan + 1j * np.nan
            resids[i] = float(np.min(residuals))
            continue
        surviving, surviving_resid = energies[keep], residuals[keep]
        pick = int(np.argmin(np.abs(surviving - seed)))
        poles[i] = surviving[pick]
        resids[i] = surviving_resid[pick]

    return CoupledCurve(R=R_arr, E_res=poles, residual=resids, n_stable=counts, n_poles=survivors)


S_VALUES = tuple(round(0.1 * i, 1) for i in range(11))
KAPPA_VALUES = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5)
KAPPA_REFERENCE = 0.3
N_CHANNEL_VALUES = (1, 2, 3, 4)
# The walk stops when the width exceeds the resonance energy NOWHERE ON THE
# CURVE is it still a resonance -- the MINIMUM of Gamma/eps over the points that
# are actually resonant, not the maximum.
#
# The maximum cannot be the test: measured, the shipped s = 0 model ALREADY has
# Gamma/eps = 1.032 at R = 1.60, where the anion curve is high on the repulsive
# wall and the state is legitimately broad and short-lived. A max rule stops the
# walk at s = 0 on the reference model itself.
#
# The minimum behaves (N_l = 2, kappa = 0.3 walk): 0.180 (s=0, at R=2.26;
# 0.448 at R=2.15 is the SECOND-smallest, not the minimum) -> 0.548 (0.3, at
# R=2.59) -> 0.872 (0.4, at R=2.48) -> 1.177 (0.5, at R=2.37). The narrowest
# point sits in the R = 2.15-2.59 crossing region only from s=0.3 onward --
# where the vibrational wavefunctions live and the cross section is built --
# NOT at s=0.1 (R=5.89 for N_l=1) or s=0.2 (R=2.81-3.14), where the minimum
# briefly wanders before the walk settles into that region. So the walk ends
# at s = 0.5, and it ends because the state has stopped being a resonance
# everywhere that matters rather than because a pole finder gave up -- the
# published deck loses the pole at s = 0.5 for want of tail elements, which
# would have been the deck talking, not the physics.
GAMMA_OVER_EPS_MAX = 1.0
# A width at round-off is not a width. Bound points sit at 1e-13; genuine widths
# are 1e-2. Four orders of clearance either side.
GAMMA_FLOOR = 1e-6
R_SAMPLE = np.linspace(1.6, 6.0, 41)
# NO's resonance sits 0.02-0.05 Ha ABOVE the neutral curve
# (docs/physics/diatomic-ve-cross-sections.md), and the electronic
# Hamiltonian's diagonal contains v0(R) -- so the pole is at v0(R) + this
# offset, and a constant seed would be adrift by up to 0.25 Ha. A window
# centre, not a result.
SEED_OFFSET = 0.03 - 0.01j
RESULTS = Path("validation/coupled/results")


def _electronic_grids() -> tuple[FemDvrEcsGrid, FemDvrEcsGrid]:
    return tuple(electronic_grid(angle_deg=a, **NO_ELECTRONIC) for a in ANGLES)


def run_continuation(
    *,
    n_channels: int,
    kappa: float,
    s_values: Sequence[float] = S_VALUES,
    R_sample: npt.NDArray[np.float64] = R_SAMPLE,
) -> dict[float, CoupledCurve]:
    """The `s` walk at fixed `kappa`, each step seeded from the previous one.

    The chain is anchored at `s = 0`, where the model IS `qscat.model.NO` --
    so no seed anywhere on the trajectory comes from the approximation the
    campaign is measuring.
    """
    ga, gb = _electronic_grids()
    out: dict[float, CoupledCurve] = {}
    v0 = np.asarray(NO.v0(R_sample).real, dtype=np.float64)
    seeds = np.asarray(v0 + SEED_OFFSET, dtype=np.complex128)
    for s in s_values:
        model = CoupledModel(
            well=TwoCentreWell(base=NO, s=float(s), kappa=float(kappa)),
            n_channels=n_channels,
        )
        curve = coupled_resonance_curve(model, R_sample, ga, gb, seeds=seeds)
        out[float(s)] = curve  # record BEFORE testing, so the crossing is kept
        finite = np.isfinite(curve.E_res)
        eps = curve.v_d - v0
        # Only a point that IS a resonance can testify that the resonance has
        # gone. A bound point (eps <= 0) is not broad, it is the opposite -- and
        # there are 34 of them at s = 0, because at R >= ~2.3 the anion starts
        # bound and the anisotropy is what makes it resonant.
        active = finite & (eps > 0.0) & (curve.gamma > GAMMA_FLOOR)
        if active.any() and float(np.min(curve.gamma[active] / eps[active])) >= (
            GAMMA_OVER_EPS_MAX
        ):
            break
        # A no-target point is a NORMAL feature of a resonance curve, not a
        # failure of the walk -- the repository already treats a crossing slice
        # that way. Carry the analytic guess where the last curve had none,
        # rather than propagating a nan seed into the next step.
        seeds = np.asarray(np.where(finite, curve.E_res, v0 + SEED_OFFSET), dtype=np.complex128)
    return out


def main(results: Path = RESULTS) -> dict[str, object]:
    """Run the full continuation and write `results/screen.json`.

    `results` defaults to the committed `RESULTS` directory; a caller that
    does not want to overwrite the tracked campaign data (a test, in
    particular) should pass a scratch directory instead.
    """
    report: dict[str, object] = {
        "n_channels": list(N_CHANNEL_VALUES),
        "s_values": list(S_VALUES),
        "kappa_values": list(KAPPA_VALUES),
        "kappa_reference": KAPPA_REFERENCE,
        "R": R_SAMPLE.tolist(),
        "s_curves": {},
        "kappa_curves": {},
    }
    for n_ch in N_CHANNEL_VALUES:
        s_walk = run_continuation(n_channels=n_ch, kappa=KAPPA_REFERENCE)
        report["s_curves"][str(n_ch)] = {  # type: ignore[index]
            str(s): _curve_payload(c) for s, c in s_walk.items()
        }
        # The whole walk per kappa, not just its endpoint: the walk stops where
        # the pole stops being a resonance, and full and fixed-l need not stop at
        # the same s. The comparison must still be made at a MATCHED s, so the
        # consumer needs both ladders, not two endpoints.
        #
        # KAPPA_REFERENCE is itself in KAPPA_VALUES, so its walk is the one
        # already computed above as `s_walk`. Reuse it rather than repeating
        # it: on a ~36-minute campaign that is one redundant walk per
        # n_channels value.
        per_kappa: dict[str, dict[str, dict[str, list[float]]]] = {}
        for k in KAPPA_VALUES:
            walk = s_walk if k == KAPPA_REFERENCE else run_continuation(n_channels=n_ch, kappa=k)
            per_kappa[str(k)] = {str(s): _curve_payload(c) for s, c in walk.items()}
        report["kappa_curves"][str(n_ch)] = per_kappa  # type: ignore[index]
    # Spec gate 6: N_l convergence. The WHOLE walk on the SAME ladder as every
    # other curve, because a convergence check must be read at a matched s and
    # the walks stop where the physics says, not where a hardcoded index says.
    # An earlier version walked a coarse (0, 0.5, 1) ladder and took [1.0],
    # which both skipped the stop condition and produced a curve at an s no
    # other curve reached -- uncomparable, and it silently looked like a 106%
    # convergence failure when read against N_l = 4 at a different s.
    report["n_channels_5_check"] = {  # type: ignore[assignment]
        str(s): _curve_payload(c) for s, c in run_continuation(n_channels=5, kappa=0.5).items()
    }
    results.mkdir(parents=True, exist_ok=True)
    (results / "screen.json").write_text(json.dumps(report, indent=1))
    print(f"[coupled] wrote {results / 'screen.json'}")
    return report


def _curve_payload(curve: CoupledCurve) -> dict[str, list[float]]:
    return {
        "v_d": curve.v_d.tolist(),
        "gamma": curve.gamma.tolist(),
        "residual": curve.residual.tolist(),
        "n_stable": curve.n_stable.tolist(),
        "n_poles": curve.n_poles.tolist(),
    }


if __name__ == "__main__":
    main()
