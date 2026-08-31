"""Preserving the dissociation limit: how many channels, and how much `lam`.

The two-centre split hands the deeper well only `(1+kappa)/2` of `lam`, so the
`R -> infinity` binding is lost and the anion unbinds beyond `R ~ 0.7/s` for
ANY `s > 0`. Two things must be fixed before the model's asymptote means
anything, and they are not independent:

1. **The channel cutoff.** A single-centre partial-wave expansion about the
   MOLECULAR centre has to represent an electron localised on one nucleus a
   distance `d = sR/2` away, so the count needed grows with `d`. A cutoff
   validated at small `R` says nothing at large `R`.

2. **The `lam` renormalisation.** `f(R)` is solved at every `R` so the
   two-centre resonance POSITION reproduces the shipped model's `E_res(R)`.
   That pins the curve, the crossing and the asymptote by construction and
   leaves `Gamma` free -- which is the physics under study: the anisotropy
   should change how the electron COUPLES, not where the anion state sits.
   With a converged basis the large-`R` limit of `f` must be the analytic
   `2/(1+kappa)`, the inverse of the deeper well's share; that is the check,
   not the input.

Both are measured on `f`, never on a raw eigenvalue. The reason is a trap:
this grid carries a spurious near-threshold state at `+0.00129-0.00283i` that
is present at every `s`, `R` and `lam`, so a search seeded near a vanished
bound state snaps onto it and reports a converged number that is pure grid.
`_responds_to_lam` rejects it on the one property no physical state has --
indifference to the depth of the well binding it.

Run: `python -m validation.coupled.renormalise`
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt
from qscat.core.grids import electronic_grid
from qscat.linalg import ShiftInvertEigs
from qscat.model import NO
from scipy.optimize import brentq

__all__ = [
    "RESULTS",
    "ScaledLam",
    "channel_cutoff",
    "e_res",
    "e_res_pole",
    "lam_factor",
    "lam_factor_pole",
    "main",
]

RESULTS = Path("validation/coupled/results")
KAPPA = 0.5
# The anisotropy screen's electronic deck, not NO's published one. The
# published deck's tail cannot hold the pole once anisotropy broadens the
# resonance -- it loses it outright by s = 0.5, where this one follows it to
# s = 1. A single angle suffices on the bound-state branch: a bound state is
# ECS-invariant, so it needs no angle-stability test to be told from continuum.
GRID = {"angle_deg": 44.0, "r_max": 16.0, "order": 8, "n_complex": 8}
F_TOL = 1e-3  # cutoff criterion on f, which is dimensionless and O(1)
LADDER = (4, 8, 12, 16, 20, 24)


class ScaledLam:
    """`base` with `lam(R)` multiplied by `f`; everything else delegated.

    A wrapper, not a new model: `v0`, `mu`, `alpha_c`, `ell` and the charge
    must stay bit-identical to the shipped model or the `s = 0`, `f = 1`
    embedding stops being an identity.
    """

    def __init__(self, base: Any, f: float) -> None:
        self._b, self._f = base, float(f)

    def __getattr__(self, name: str) -> Any:
        return getattr(self._b, name)

    def lam(self, R: npt.ArrayLike) -> Any:
        return self._f * self._b.lam(R)


def _hamiltonian(s: float, R: float, f: float, n_ch: int, kappa: float):
    from projects.no_coupled_channels.anisotropy import TwoCentreWell
    from projects.no_coupled_channels.model import CoupledModel

    model = CoupledModel(
        well=TwoCentreWell(base=ScaledLam(NO, f), s=s, kappa=kappa), n_channels=n_ch
    )
    return model.electronic_hamiltonian(electronic_grid(**GRID), complex(R))


def e_res(
    s: float, R: float, f: float, n_ch: int, seed: complex, *, kappa: float = KAPPA
) -> complex:
    """The LOWEST bound state near `seed`, measured from the neutral.

    Lowest, not nearest, and the distinction is load-bearing. Deepening the
    well pushes the ground state monotonically DOWN, which is what makes a
    bracketed root-find on `f` valid; the state merely *nearest* the target is
    not monotonic in `f` at all, because excited states sweep down through the
    target as the well deepens. Seeding on nearest gave roots at `f ~ 5-6`,
    an excited state impersonating the anion.

    Returns `nan` when no bound state is present -- the caller needs to see
    that rather than be handed the nearest continuum point.

    Two shifts are pooled, not one: the tracking seed, and a probe just below
    threshold. A single shift is not enough in either direction -- tracking
    alone misses a state that has just appeared at `E -> 0^-`, and the
    threshold probe alone latches onto whichever EXCITED state happens to sit
    near threshold once the well is deep (measured: it returned `f = 2.15`
    for the `s = 0` case whose answer is exactly 1).
    """
    v0 = complex(np.asarray(NO.v0(R)).ravel()[0])
    H = _hamiltonian(s, R, f, n_ch, kappa)
    solver = ShiftInvertEigs(H, k=12)
    pooled: list[complex] = []
    for shift in (seed, complex(-1e-3)):
        out = solver.near(shift + v0)
        vals = np.asarray(out[0] if isinstance(out, tuple) else out).ravel() - v0
        pooled.extend(vals.tolist())
    arr = np.asarray(pooled)
    bound = arr[(np.abs(arr.imag) < 1e-7) & (arr.real < 0.0)]
    return complex(bound[np.argmin(bound.real)]) if bound.size else complex("nan")


def _responds_to_lam(s: float, R: float, f: float, n_ch: int, seed: complex, kappa: float) -> bool:
    """Does this state care how deep the well is? A physical state does.

    The grid's spurious near-threshold state does not move when `lam` changes;
    a bound state or resonance of the well moves a lot (measured: 5 % in `lam`
    is the difference between bound at -0.0598 and not bound at all).
    """
    a = e_res(s, R, f, n_ch, seed, kappa=kappa)
    b = e_res(s, R, f * 1.02, n_ch, seed, kappa=kappa)
    return bool(abs(a - b) > 1e-6)


def lam_factor(
    s: float,
    R: float,
    target: float,
    n_ch: int,
    *,
    kappa: float = KAPPA,
    bracket: tuple[float, float] = (1.0, 3.0),
    step: float = 0.02,
) -> float | None:
    """`f` with `Re E_res = target`, by CONTINUATION in `f`.

    A fixed seed cannot do this job. Deepening the well drags the ground state
    far below the target, out of the set shift-invert returns near a fixed
    shift, so the search silently loses the state it is tracking (measured:
    `nan` at every `f >= 1.5`). Instead `f` is walked in small steps with the
    shift re-seeded on the previous step's eigenvalue, exactly as the
    resonance curve is walked in `R`. The state is then always found near
    where it actually is.

    The walk is seeded AT the target, which is where the state sits when the
    well is already the right depth (the `s = 0`, `f = 1` case, whose answer
    is exactly 1). `e_res` pairs that with its own near-threshold probe, so a
    state that has not appeared yet is still caught the moment it binds.
    """
    seed = complex(target)
    prev_f: float | None = None
    prev_e: float | None = None
    f = bracket[0]
    while f <= bracket[1] + 1e-12:
        e = e_res(s, R, f, n_ch, seed, kappa=kappa)
        if np.isfinite(e.real):
            seed = e
            if e.real <= target:  # crossed from above
                if prev_f is None:
                    return None  # already past it at the first step
                # linear in f between the bracketing pair, then polish
                lo, hi = prev_f, f

                def g(x: float, _seed: complex = seed) -> float:
                    return e_res(s, R, x, n_ch, _seed, kappa=kappa).real - target

                try:
                    root = float(brentq(g, lo, hi, xtol=1e-6, rtol=1e-10))
                except Exception:
                    root = float(np.interp(target, [e.real, prev_e], [hi, lo]))
                return root if _responds_to_lam(s, R, root, n_ch, seed, kappa) else None
            prev_f, prev_e = f, e.real
        f += step
    return None


def channel_cutoff(s: float, R: float, target: float) -> tuple[int | None, list[float | None]]:
    """Smallest `N_l` whose `f*` is within `F_TOL` of the next rung.

    The ladder is returned with the cutoff because a cutoff quoted without its
    ladder cannot be checked.
    """
    vals = [lam_factor(s, R, target, n) for n in LADDER]
    for i in range(len(LADDER) - 1):
        a, b = vals[i], vals[i + 1]
        if a is not None and b is not None and abs(a - b) < F_TOL:
            return LADDER[i], vals
    return None, vals


def e_res_pole(
    s: float,
    R: float,
    f: float,
    n_ch: int,
    seed: complex,
    *,
    kappa: float = KAPPA,
    half_width: float = 0.35,
) -> complex:
    """`E_res` in the RESONANCE region, by two-angle ECS stability.

    Below the crossing (`R < ~2.25`) the anion is not bound, so the bound-state
    filter in `e_res` finds nothing and the renormalisation cannot be closed
    there -- which matters, because that is exactly where the VE cross section
    gets its physics. This delegates to the screen's own pole walk, which
    already applies the residual cut BEFORE the nearest-to-seed pick and so
    rejects the spurious near-threshold state on identity rather than distance.

    Cheap where it is needed: the resonance region is small `R`, hence small
    `d = sR/2`, hence `N_l = 4` by the cutoff rule -- the expensive large-`N_l`
    regime is at large `R`, where `e_res` already has the answer.
    """
    from projects.no_coupled_channels.anisotropy import TwoCentreWell
    from projects.no_coupled_channels.model import CoupledModel
    from validation.coupled.screen import _electronic_grids, coupled_resonance_curve

    ga, gb = _electronic_grids()
    v0 = complex(np.asarray(NO.v0(R)).ravel()[0])
    model = CoupledModel(
        well=TwoCentreWell(base=ScaledLam(NO, f), s=s, kappa=kappa), n_channels=n_ch
    )
    curve = coupled_resonance_curve(model, [R], ga, gb, seeds=[v0 + seed], half_width=half_width)
    return complex(curve.E_res[0] - v0)


def lam_factor_pole(
    s: float,
    R: float,
    target: float,
    n_ch: int,
    *,
    kappa: float = KAPPA,
    bracket: tuple[float, float] = (0.95, 1.25),
) -> float | None:
    """`f` with `Re E_res = target` in the resonance region.

    The bracket is NARROW on purpose, and widening it breaks the search rather
    than strengthening it: the pole is sought in a finite window about the
    seed, so a large `f` drives the state clean out of that window, the walk
    returns `nan`, and the sign change brentq needs never appears. Measured,
    `f` never leaves 1.0-1.1 anywhere in the resonance region.

    A plain bracketed root-find is safe here: `Re E_res` is monotonically
    decreasing in `f` (measured at `R = 2.0`, `s = 0.1`: +0.0567, +0.0195,
    -0.0242 for `f` = 1.00, 1.05, 1.10) and the pole walk returns `nan` rather
    than a wrong state when the window holds no pole.
    """

    def g(x: float) -> float:
        e = e_res_pole(s, R, x, n_ch, complex(target), kappa=kappa)
        return 1.0 if not np.isfinite(e.real) else e.real - target

    try:
        return float(brentq(g, *bracket, xtol=1e-6, rtol=1e-10))
    except Exception:
        return None


def _shipped_target(results: Path) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """The shipped model's own `E_res(R) = v_d - v0`, from the committed screen."""
    d = json.loads((results / "screen.json").read_text())
    R = np.asarray(d["R"], dtype=float)
    vd = np.asarray(d["s_curves"]["1"]["0.0"]["v_d"], dtype=float)
    return R, vd - np.real(np.asarray(NO.v0(R)))


def main(results: Path = RESULTS) -> dict[str, object]:
    R_ship, target_curve = _shipped_target(results)
    analytic = 2.0 / (1.0 + KAPPA)
    report: dict[str, object] = {
        "kappa": KAPPA,
        "f_tol": F_TOL,
        "grid": GRID,
        "analytic_factor": analytic,
    }

    print(f"Analytic large-separation factor 2/(1+kappa) = {analytic:.4f}")
    print(f"Cutoff criterion: |f(N_l) - f(2 N_l)| < {F_TOL}\n")
    print("=== 1. Channel cutoff, against well separation d = sR/2 ===")
    print("   s      R      d      N_l*     " + "   ".join(f"f({n})" for n in LADDER))
    cutoffs = []
    for s in (0.1, 0.2, 0.3, 0.5):
        for R in (3.0, 6.0, 12.0, 20.0):
            tgt = float(np.interp(R, R_ship, target_curve))
            cut, ladder = channel_cutoff(s, R, tgt)
            txt = "  ".join(f"{v:.4f}" if v is not None else "  --  " for v in ladder)
            print(f"  {s:.1f}  {R:5.1f}  {s * R / 2:5.2f}   {cut!s:>4}    {txt}")
            cutoffs.append({"s": s, "R": R, "d": s * R / 2, "n_l": cut, "ladder": ladder})
    report["cutoffs"] = cutoffs

    found = [int(c["n_l"]) for c in cutoffs if c["n_l"] is not None]
    n_ch = max(found) if found else LADDER[-1]
    print(f"\n=== 2. Renormalised lam(R), at the worst-case cutoff N_l = {n_ch} ===")
    curves = {}
    for s in (0.1, 0.2, 0.3):
        print(f"\n  s = {s}")
        print("      R     target E_res     f(R)      Re E_res got      Gamma got")
        rows = []
        for R in (2.5, 3.0, 4.0, 6.0, 8.0, 12.0, 20.0):
            tgt = float(np.interp(R, R_ship, target_curve))
            f = lam_factor(s, R, tgt, n_ch)
            if f is None:
                print(f"    {R:5.1f}   {tgt:+.6f}     no root")
                continue
            got = e_res(s, R, f, n_ch, complex(tgt))
            print(
                f"    {R:5.1f}   {tgt:+.6f}   {f:.5f}    {got.real:+.6f}     {-2 * got.imag:+.6f}"
            )
            rows.append({"R": R, "target": tgt, "f": f, "re": got.real, "gamma": -2 * got.imag})
        curves[str(s)] = rows
    report["renormalised"] = curves
    report["n_channels"] = n_ch

    results.mkdir(parents=True, exist_ok=True)
    (results / "renormalise.json").write_text(json.dumps(report, indent=1))
    print(f"\nwrote {results / 'renormalise.json'}")
    return report


if __name__ == "__main__":
    main()
