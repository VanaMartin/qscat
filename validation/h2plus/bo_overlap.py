"""Is an exact H2+ pole a resonance at all? Ask its overlap with the BO states.

The machinery this module used to carry now lives in the library --
`qscat.core.assignment` (overlap pairing, the verdicts, the closed-channel
admissibility check) and `qscat.core.bo` (the BO product basis). What stays here
is the H2+ campaign: which curves to build, which windows to solve, and the
measured result below.

Energy proximity cannot tell a genuine resonance from a rotated-continuum
eigenvalue that happens to land nearby: both are just numbers on the same axis.
The published work assigns a cross-section feature to a quasi-bound state by
the overlap `|<omega_i^j | Psi>|` instead (M. Vana, doctoral thesis, Charles
University 2017, Table 4.2, p. 73), and the same quantity answers the prior
question of whether there is a state there to assign.

The discriminator is that the two candidates behave differently under overlap,
not just under energy:

- A **genuine resonance** is, to the extent the Born-Oppenheimer picture holds
  at all, a product `phi_Ryj(r; R) chi_v(R)` plus a correction. It therefore has
  substantial overlap with ONE BO state and little with the rest, and the node
  and lobe counts agree -- which is a strong constraint, because neighbouring
  levels differ in Rydberg index or vibrational quantum number and so in node
  count.
- A **discretized-continuum state** rotated onto the real axis by ECS has no
  such structure. It overlaps everything weakly and nothing strongly.

This matters because the campaign found more poles than there are BO levels,
and some of them sit where the published cross section shows no peak at all.
Those split into two classes, and only one is benign:

- **On a BO level, no peak.** The high-`n` Rydberg states (`Ry_6..Ry_11` at
  `v=1,2,3`): E = 0.001563, 0.010781, 0.014779, 0.023513, 0.024206 Ha, each
  within 0.7 resonance widths of its BO level. `Gamma` falls as `n^-3`, so
  these are far narrower than the published sweep's 1e-5 Ha sampling and
  CANNOT produce a resolvable peak. Their absence from the cross section is
  expected.
- **Far from a peak AND far from any BO level.** E = 0.004479 (6.1 widths from
  the nearest peak, 5.3 from `omega_1^8`) and E = 0.021702 / 0.021782 (36 and
  32 widths from a peak, 20 and 24 from a level). Nothing explains these on
  either axis, and they are what this module exists to test.

## Measured, window 0 (r_max = 300)

The two classes separate by three orders of magnitude, so this is not a
marginal judgement:

| E (Ha) | best BO | overlap | second | |
|---|---|---|---|---|
| 0.001563 | `w^6_1` | 0.9870 | `w^4_2` 0.10 | genuine, no peak (narrow) |
| 0.003924 | `w^4_2` | 0.8777 | `w^3_3` 0.33 | genuine, mixed |
| **0.004479** | `w^2_5` | **0.0006** | 0.0001 | **NOT a resonance** |
| 0.004607 | `w^8_1` | 0.9769 | `w^4_2` 0.14 | genuine |
| 0.005661 | `w^3_3` | 0.7831 | `w^2_5` 0.46 | genuine, strongly mixed |
| 0.006316 | `w^2_5` | 0.8700 | `w^3_3` 0.44 | genuine, strongly mixed |
| 0.007171 | `w^12_1` | 0.9901 | `w^5_2` 0.07 | genuine |
| 0.008180 | `w^16_1` | 0.6828 | `w^5_2` 0.65 | genuine, MAXIMALLY mixed |
| 0.008260 | `w^16_1` | 0.7252 | `w^5_2` 0.60 | genuine, MAXIMALLY mixed |

**E = 0.004479 is spurious.** Its best overlap against a locally complete basis
is 6e-4 where real states score 0.97-0.99. It passed the two-angle stability
test that `exact_resonance_states` applies, which is the point: **angle
stability is necessary and not sufficient**, exactly as
`docs/physics/lcp-resonance-levels.md` warns for the 1-D case. Overlap is the
check that catches it.

**The basis has to be deep enough or the test lies.** With curves only to
`Ry_11`, the poles at 0.0072-0.0085 scored 0.02-0.09 and looked spurious too.
They are `Ry_12..Ry_16 v=1` at 0.97-0.99, and their partners simply were not in
the basis. A low overlap means "no partner HERE", never "no partner". The
series accumulates at the threshold, so any finite basis runs out eventually --
0.008395 and 0.008530 still score 0.21 and 0.11 and need `Ry_17+`.

**Two results the overlap gives that energy alone cannot.** The high-`n`
Rydberg states are nearly PURE BO products (0.97-0.99) while the compact low-`n`
ones are strongly MIXED (0.78-0.88, with 0.33-0.46 of a second level) -- the
same regime split the energy shifts show, seen structurally. And at 0.0082-0.0083
the `Ry_16 v=1` series crosses `omega_2^5`: two exact states come out ~0.68/0.65
and ~0.73/0.60 of the two BO levels, i.e. neither is either one. That is the
sharpest Born-Oppenheimer breakdown in the set -- past shifting a level, the
identities themselves stop being well defined.

Run as::

    uv run python -m validation.h2plus.bo_overlap

Each window's poles are cached to a git-ignored npz (delete to recompute).
Building the BO basis is one pass of dense electronic eigensolves over the
nuclear grid -- a few minutes -- shared across all `(curve, vib)` pairs rather
than repeated per pair.
"""

from __future__ import annotations

import pathlib
import time

import numpy as np
import numpy.typing as npt
from qscat.core import (
    BoBasis,
    bo_basis,
    electronic_curves,
    exact_resonance_states,
    pair_by_overlap,
)
from qscat.dvr import TensorGrid
from qscat.model import H2P

__all__ = ["CURVES", "N_VIB", "bo_basis_for", "solve_window", "main"]

# Deep enough to reach the `Ry_16` the top of window 0 needs. The series
# accumulates at the threshold, so this is a CUTOFF rather than a limit: poles
# above it come back `basis-limited`, which is the honest verdict and not a
# claim that they are spurious.
CURVES = tuple(range(2, 17))
N_VIB = 8


def bo_basis_for(tgrid: TensorGrid, *, curves=CURVES, n_vib: int = N_VIB) -> BoBasis:
    """The BO product basis for one window's grid.

    One pass of dense electronic eigensolves over the nuclear grid, shared
    across every `(curve, vib)` pair rather than repeated per pair. The
    machinery is `qscat.core.bo`; H2+-specific here is only how deep to build.
    """
    g_r, g_R = tgrid.grids
    cur = electronic_curves(H2P, g_r, g_R, n_curves=max(curves) + 1, with_states=True)
    return bo_basis(cur, g_R, H2P.mu, n_vib=n_vib, allow_partial=True)


def solve_window(
    window: int, *, r_max: float = 300.0, cache_dir: pathlib.Path | None = None
) -> tuple[npt.NDArray[np.complex128], npt.NDArray[np.complex128], TensorGrid]:
    """Exact poles WITH states for one DR window, cached to a git-ignored npz.

    The cache stores energies and states ONLY, not a whole
    `ExactResonanceStates`. That is deliberate rather than legacy: overlap
    pairing needs exactly those two, and the widths and angle residuals would
    make each of these three files larger without a consumer. A caller that
    wants the full record should use `exact_poles.exact_poles`, which keeps it.
    """
    from validation.h2plus.exact_poles import K_SEARCH, find_seeds, grid_family

    base, moved_el, moved_nu = grid_family(r_max)

    cache = (cache_dir or pathlib.Path(__file__).parent) / f"bo_overlap.w{window}.npz"
    if cache.exists():
        z = np.load(cache)
        return z["energies"], z["states"], base

    seeds = [s for s in find_seeds()[0] if s.window == window]
    lo = min(s.e_tot for s in seeds) - 0.01
    hi = max(s.e_tot for s in seeds) + 0.01
    t0 = time.perf_counter()
    res = exact_resonance_states(
        H2P,
        base,
        moved_el,
        moved_nu,
        shifts=[complex(s.e_tot, -1e-4) for s in seeds],
        window=(lo, hi, -0.01, 0.0),
        k=K_SEARCH,
    )
    print(f"window {window}: {res.energies.size} poles in {time.perf_counter() - t0:.0f}s")
    np.savez(cache, energies=res.energies, states=res.states)
    return res.energies, res.states, base


def main(windows: tuple[int, ...] = (0, 1, 2)) -> None:
    from validation.h2plus.exact_poles import EPS0, THRESHOLDS

    tally: dict[str, int] = {}
    for w in windows:
        energies, states, base = solve_window(w)
        print(f"\n=== window {w}: {energies.size} poles ===", flush=True)
        print(f"building BO basis, curves {CURVES[0]}..{CURVES[-1]} ...", flush=True)
        basis = bo_basis_for(base)

        print(
            f"{'E (Ha)':>10} {'level':>9} {'overlap':>8} {'2nd':>9} "
            f"{'2nd val':>8} {'shift(meV)':>11}  verdict"
        )
        for i in np.argsort(energies.real):
            p = pair_by_overlap(energies[i], states[:, i], basis, THRESHOLDS)
            tally[p.verdict] = tally.get(p.verdict, 0) + 1
            lvl = "-" if p.level is None else f"w^{p.level[0]}_{p.level[1]}"
            second = (
                "-" if p.second_level is None else f"w^{p.second_level[0]}_{p.second_level[1]}"
            )
            print(
                f"{p.pole_energy - EPS0:>10.6f} {lvl:>9} {p.overlap:>8.4f} "
                f"{second:>9} {p.second_overlap:>8.4f} {p.shift_mev:>11.3f}  {p.verdict}"
            )
    print(f"\nverdict tally across {len(windows)} window(s): {tally}")


if __name__ == "__main__":
    main()
