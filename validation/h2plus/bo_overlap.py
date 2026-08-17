"""Is an exact pole a resonance at all? Ask its overlap with the BO states.

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

Reuses the cached window-0 states from `resonance_state_figures` (delete that
cache to recompute). Building the BO basis is one pass of dense electronic
eigensolves over the nuclear grid -- a few minutes -- shared across all
`(curve, vib)` pairs rather than repeated per pair.
"""

from __future__ import annotations

import pathlib
import time
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from qscat.core import exact_resonance_states, vibrational_states
from qscat.dvr import TensorGrid, eigen, kinetic
from qscat.linalg import c_product
from qscat.model import H2P
from qscat.units import HARTREE_TO_EV

__all__ = ["BoState", "OverlapPair", "bo_basis", "overlap", "pair_by_overlap", "main"]

# An overlap below this is "no partner in this basis". It does NOT distinguish a
# spurious state from a real one whose partner was never built -- see the module
# docstring; both look identical here and only a deeper basis separates them.
NO_PARTNER = 0.10

# Between NO_PARTNER and this, the best match exists but is not an
# identification: a genuine assignment scores 0.87-0.99 and even the strongly
# mixed states reach 0.63-0.78, so anything under half is a state this basis
# does not really describe.
WEAK = 0.50

# `second/best` above this means the state is a blend, not an assignment: at the
# Ry_16 v=1 / omega_2^5 crossing the ratio reaches 0.95 and neither label fits.
MIXED_RATIO = 0.70

# A partner further than this in energy is reported even when the overlap is
# high. Overlap alone would happily pair a state with a level several meV away,
# and past a few meV that is a statement about mixing rather than a shift worth
# quoting. Calibrated on the measured set: every clean pairing sits well inside
# it, and the rows that exceed it are the known crossings.
MAX_SHIFT_MEV = 20.0


@dataclass(frozen=True)
class BoState:
    """One Born-Oppenheimer product state and the energy its curve gives it."""

    psi: npt.NDArray[np.complex128]
    energy: float


@dataclass(frozen=True)
class OverlapPair:
    """One exact pole paired to a BO level BY OVERLAP, with the checks."""

    pole_energy: float  # absolute, Ha
    level: tuple[int, int] | None  # (curve, vib), None when nothing matched
    overlap: float
    second_level: tuple[int, int] | None
    second_overlap: float
    shift_mev: float
    # "ok" | "spurious" | "basis-limited" | "weak" | "mixed" | "distant"
    verdict: str

    @property
    def label(self) -> str:
        return "-" if self.level is None else f"w^{self.level[0]}_{self.level[1]}"


def bo_basis(tgrid: TensorGrid, curves: list[int], vibs: int) -> dict[tuple[int, int], BoState]:
    """`{(curve, vib): phi_Rycurve(r; R) chi_vib(R)}` on `tgrid`, c-normalized.

    One pass over the nuclear grid builds every requested curve's electronic
    eigenvector, so the cost is a single sweep of dense eigensolves rather than
    one per `(curve, vib)` pair.

    Each curve's eigenvector is phase-aligned across `R` by continuity: the
    phase is arbitrary at every nuclear point independently, and without
    alignment the product flips sign at random `R`, which destroys the overlap
    this module measures (the sign flips integrate to near zero against a smooth
    partner).
    """
    g_r, g_R = tgrid.grids
    top = max(curves) + 1
    phi = {j: np.empty((g_r.n, g_R.n), dtype=np.complex128) for j in curves}
    vals_all = np.empty((top, g_R.n), dtype=np.complex128)
    prev: dict[int, npt.NDArray[np.complex128] | None] = {j: None for j in curves}
    for k, R in enumerate(g_R.points):
        H_el = kinetic(g_r, 1.0) + np.diag(H2P.surface(g_r.points, complex(R)))
        vals, vecs = eigen(H_el)
        vals_all[:, k] = vals[:top]
        for j in curves:
            vec = vecs[:, j]
            p = prev[j]
            if p is not None:
                ov = complex(np.vdot(p, vec))
                if ov != 0:
                    vec = vec * (abs(ov) / ov)
            prev[j] = vec
            phi[j][:, k] = vec

    out: dict[tuple[int, int], BoState] = {}
    for j in curves:

        def v_n(
            _R: npt.ArrayLike, _c: npt.NDArray[np.complex128] = vals_all[j]
        ) -> npt.NDArray[np.complex128]:
            return np.asarray(_c, dtype=np.complex128)

        basis = vibrational_states(g_R, H2P.mu, vibs, v_n)
        for v in range(vibs):
            psi = (phi[j] * basis.chi[v][None, :]).ravel()
            nrm = complex(c_product(psi, psi))
            out[(j, v)] = BoState(
                psi=np.asarray(psi / np.sqrt(nrm), dtype=np.complex128),
                energy=float(basis.eps[v]),
            )
    return out


def overlap(a: npt.NDArray[np.complex128], b: npt.NDArray[np.complex128]) -> float:
    """`|<a|b>|` under the c-product, with both sides c-normalized.

    The c-product (bilinear, not conjugated) is the ECS-correct pairing -- the
    same convention `qscat.core`'s cross sections use. A conjugated dot would
    weight the exponentially growing ECS tail instead of cancelling it.
    """
    na = complex(c_product(a, a))
    nb = complex(c_product(b, b))
    if na == 0 or nb == 0:
        return 0.0
    return float(abs(complex(c_product(a, b)) / np.sqrt(na * nb)))


def admissible_levels(e_tot: float, thresholds: npt.NDArray[np.float64]) -> list[tuple[int, int]]:
    """The `(curve, vib)` levels that can exist AT this energy, from energy alone.

    A Rydberg series is attached to a CLOSED channel: above `eps[v]` that
    vibrational channel is open (the VE channel), and states attached to it are
    continuum rather than bound. So only thresholds above `e_tot` contribute,
    and each contributes one index:

        binding = eps[v] - e_tot,  n_eff = 1/sqrt(2*binding),  curve ~ n_eff - 1

    (the last from the measured series, where `Ry_j` has `n_eff ~ j+1`).

    The consequence is a strong constraint and the reason this function exists:
    at fixed energy a HIGHER vibrational level needs a LARGER binding and so a
    LOWER Rydberg index. The admissible set is therefore finite and small, and
    a basis can be checked for covering it -- which turns "is this state
    spurious, or is its partner merely missing?" from a judgement call into a
    computation. See `basis_covers`.

    The set is finite only away from an accumulation region: as `e_tot ->
    eps[v]` the binding tends to zero and the admissible index diverges. That
    happens in the last ~1 mHa below each threshold, which the DR windows are
    trimmed to exclude and the `n_eff <= 12` cutoff excludes again.
    """
    out: list[tuple[int, int]] = []
    for v, thr in enumerate(thresholds):
        if thr <= e_tot:
            continue
        n_eff = 1.0 / np.sqrt(2.0 * (thr - e_tot))
        j = int(round(n_eff - 1.0))
        if j >= 0:
            out.append((j, v))
    return out


def basis_covers(
    e_tot: float,
    thresholds: npt.NDArray[np.float64],
    basis: dict[tuple[int, int], BoState],
) -> bool:
    """Does `basis` contain every level energetically admissible at `e_tot`?

    When it does, a low overlap means the state has no BO partner at all --
    spurious. When it does not, a low overlap is uninformative. Without this
    distinction the two are indistinguishable, and conflating them once nearly
    cost eight genuine `Ry_12..Ry_16` states (see the module docstring).

    Curves are checked with a +/-1 tolerance because the `n_eff ~ j+1` mapping
    is the asymptotic hydrogenic one and the low curves depart from it.
    """
    for j, v in admissible_levels(e_tot, thresholds):
        if not any((jj, v) in basis for jj in (j - 1, j, j + 1)):
            return False
    return True


def pair_by_overlap(
    pole_energy: complex,
    pole_state: npt.NDArray[np.complex128],
    basis: dict[tuple[int, int], BoState],
    thresholds: npt.NDArray[np.float64] | None = None,
) -> OverlapPair:
    """Pair one exact pole to the BO level it most resembles, with checks.

    Overlap is the PRIMARY criterion -- it tests physical identity, where energy
    proximity only tests that two numbers are close and cannot tell a resonance
    from a rotated-continuum state. Energy then enters as a CHECK rather than as
    the assignment, because a large overlap with an energetically distant level
    is real information: it means the state is mixing across a gap, not that it
    is that level shifted.

    Three verdicts other than "ok", in priority order:

    - `no-partner` -- nothing in the basis reaches `NO_PARTNER`. Either the
      state is spurious or its partner was never built; **this function cannot
      tell those apart** and deliberately does not pretend to.
    - `mixed` -- the second-best overlap is within `MIXED_RATIO` of the best, so
      the state is a blend and neither label describes it.
    - `distant` -- the identification is clear but the partner lies more than
      `MAX_SHIFT_MEV` away, so the "shift" is a mixing statement.
    """
    ranked = sorted(
        ((overlap(b.psi, pole_state), key) for key, b in basis.items()), key=lambda t: -t[0]
    )
    best_v, best_k = ranked[0]
    second_v, second_k = ranked[1] if len(ranked) > 1 else (0.0, None)
    shift = (float(pole_energy.real) - basis[best_k].energy) * 1000.0 * HARTREE_TO_EV

    # Energy decides which reading a poor overlap supports. If every level that
    # could exist at this energy is in the basis, a poor overlap is a fact about
    # the STATE; otherwise it is a fact about the BASIS. This is checked for any
    # weak identification, not only for a vanishing one: a state whose true
    # partner is absent still scores moderately against a wrong partner --
    # E = 0.008395 reaches 0.21 against `omega_2^5` while the `Ry_18 v=1` it
    # actually is was never built.
    covered = (
        basis_covers(float(pole_energy.real), thresholds, basis)
        if thresholds is not None
        else False
    )
    if best_v < WEAK and not covered:
        verdict, level = "basis-limited", None
    elif best_v < NO_PARTNER:
        verdict, level = "spurious", None
    elif best_v < WEAK:
        verdict, level = "weak", best_k
    elif second_v / best_v > MIXED_RATIO:
        verdict, level = "mixed", best_k
    elif abs(shift) > MAX_SHIFT_MEV:
        verdict, level = "distant", best_k
    else:
        verdict, level = "ok", best_k

    return OverlapPair(
        pole_energy=float(pole_energy.real),
        level=level,
        overlap=best_v,
        second_level=second_k,
        second_overlap=second_v,
        shift_mev=shift,
        verdict=verdict,
    )


def solve_window(
    window: int, *, r_max: float = 300.0, cache_dir: pathlib.Path | None = None
) -> tuple[npt.NDArray[np.complex128], npt.NDArray[np.complex128], TensorGrid]:
    """Exact poles WITH states for one DR window, cached to a git-ignored npz."""
    from validation.h2plus.exact_poles import (
        _EL_BASE_DEG,
        _EL_MOVED_DEG,
        _NUC_BASE_DEG,
        _NUC_MOVED_DEG,
        K_SEARCH,
        _electronic_box,
        _nuclear_box,
        find_seeds,
    )

    el_a, el_b = _electronic_box(r_max, _EL_BASE_DEG), _electronic_box(r_max, _EL_MOVED_DEG)
    nu_a, nu_b = _nuclear_box(_NUC_BASE_DEG), _nuclear_box(_NUC_MOVED_DEG)
    base = TensorGrid([el_a, nu_a])

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
        TensorGrid([el_b, nu_a]),
        TensorGrid([el_a, nu_b]),
        shifts=[complex(s.e_tot, -1e-4) for s in seeds],
        window=(lo, hi, -0.01, 0.0),
        k=K_SEARCH,
    )
    print(f"window {window}: {res.energies.size} poles in {time.perf_counter() - t0:.0f}s")
    np.savez(cache, energies=res.energies, states=res.states)
    return res.energies, res.states, base


def main(windows: tuple[int, ...] = (0, 1, 2)) -> None:
    from validation.h2plus.exact_poles import EPS0, THRESHOLDS

    # Deep enough to reach the Ry_16 the top of window 0 needs; the series
    # accumulates, so this is a cutoff rather than a limit (see the docstring).
    curves = list(range(2, 17))
    tally: dict[str, int] = {}
    for w in windows:
        energies, states, base = solve_window(w)
        print(f"\n=== window {w}: {energies.size} poles ===", flush=True)
        print(f"building BO basis, curves {curves[0]}..{curves[-1]} ...", flush=True)
        basis = bo_basis(base, curves, 8)

        print(
            f"{'E (Ha)':>10} {'level':>9} {'overlap':>8} {'2nd':>9} "
            f"{'2nd val':>8} {'shift(meV)':>11}  verdict"
        )
        for i in np.argsort(energies.real):
            p = pair_by_overlap(energies[i], states[:, i], basis, THRESHOLDS)
            tally[p.verdict] = tally.get(p.verdict, 0) + 1
            second = "-" if p.second_level is None else f"w^{p.second_level[0]}_{p.second_level[1]}"
            print(
                f"{p.pole_energy - EPS0:>10.6f} {p.label:>9} {p.overlap:>8.4f} "
                f"{second:>9} {p.second_overlap:>8.4f} {p.shift_mev:>11.3f}  {p.verdict}"
            )
    print(f"\nverdict tally across {len(windows)} window(s): {tally}")


if __name__ == "__main__":
    main()
