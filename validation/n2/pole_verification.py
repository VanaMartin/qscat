"""Are N2's exact 2-D poles resonances? The overlap check they never had.

`docs/physics/exact-2d-resonances.md` reports five exact poles for N2 and pairs
them to the BO/LCP levels BY SORTED INDEX. Both halves of that were untested:
nothing confirmed the poles are resonances rather than rotated-continuum
eigenvalues that survived the angle test, and nothing confirmed the pairing.

On H2+ the same procedure produced four non-resonances out of 57, so the
question is not hypothetical. This module answers it for N2 with the machinery
promoted out of that campaign (`qscat.core.assignment`, `qscat.core.bo`).

## The neutral basis

H2+'s reference states are Rydberg products `phi_Ryj(r; R) chi_v(R)`: the
electron is BOUND in the cation's field. N2's are not -- the anion state is a
resonance, so the electronic factor comes from `resonance_curve` (the two-angle
pole walk, keeping the eigenvector `local_complex_potential` discards) and the
nuclear factor from `resonance_levels` (which solves the quasi-bound problem in
the complex curve). `bo_basis_from_levels` combines them. That is the same
comparator serving a different builder, which is what the promotion was for.

`basis_complete=True` is asserted here and the assertion is defensible: an anion
resonance curve carries ONE electronic state, so the basis holds every state the
BO picture admits once its vibrational ladder is built. There is no Rydberg
series to run out of, which is exactly why `admissible_levels` has nothing to
say about a neutral.

Run as::

    uv run python -m validation.n2.pole_verification
"""

from __future__ import annotations

import time

import numpy as np
from qscat.core import (
    anion_electronic_states,
    bo_basis_from_levels,
    ecs_angle_family,
    exact_resonance_states,
    pair_by_overlap,
    pair_one_to_one,
    real_weight,
    resonance_curve,
    resonance_levels,
)
from qscat.core.grids import electronic_grid, nuclear_grid
from qscat.model import N2
from qscat.units import HARTREE_TO_EV

# The same deck and seeds as `exact_resonance_figures.py`, so this verifies the
# poles that figure and docs/physics/exact-2d-resonances.md actually report.
EL = dict(r_max=24.0, order=8, n_complex=6)
NU = dict(r_max=24.0, n_complex=6, quadrature=12)
EL_ANGLES = (35.0, 44.0)
NU_ANGLES = (35.0, 30.0)
SEEDS = [
    -0.673960 - 0.0025j,
    -0.664231 - 0.0027j,
    -0.654640 - 0.0030j,
    -0.645187 - 0.0032j,
    -0.635863 - 0.0034j,
]
WINDOW = (-0.85, -0.45, -0.08, 0.0)

__all__ = ["main", "verify"]


def verify(n_levels: int = 8):
    """Pair every N2 exact pole to a BO level by overlap. Returns the pairs."""
    el_a, el_b = (
        electronic_grid(angle_deg=EL_ANGLES[0], **EL),
        electronic_grid(angle_deg=EL_ANGLES[1], **EL),
    )
    nu_a, nu_b = (
        nuclear_grid(angle_deg=NU_ANGLES[0], **NU),
        nuclear_grid(angle_deg=NU_ANGLES[1], **NU),
    )
    base, moved_el, moved_nu = ecs_angle_family(
        lambda a: electronic_grid(angle_deg=a, **EL),
        lambda a: nuclear_grid(angle_deg=a, **NU),
        electronic_angles=EL_ANGLES,
        nuclear_angles=NU_ANGLES,
    )

    t0 = time.perf_counter()
    bo = resonance_levels(N2, nu_a, nu_b, el_a, el_b, n_levels=n_levels)
    print(f"BO/LCP levels: {bo.energies.size} in {time.perf_counter() - t0:.0f}s", flush=True)

    # The electronic factor of the SAME curve those levels live in, seeded from
    # the bound anion state at the dissociation limit -- the seed
    # `local_complex_potential` uses, so both walks follow one pole.
    eps_e, _ = anion_electronic_states(el_a, N2, nu_a.R0, 1)
    seed = (eps_e[0] - 0.05, eps_e[0] + 0.05, -0.05, 0.05)
    t0 = time.perf_counter()
    cur = resonance_curve(N2, el_a, el_b, nu_a, seed, with_states=True)
    print(f"resonance curve + states in {time.perf_counter() - t0:.0f}s", flush=True)

    basis = bo_basis_from_levels(cur, bo.energies, bo.states)

    t0 = time.perf_counter()
    res = exact_resonance_states(N2, base, moved_el, moved_nu, shifts=SEEDS, k=8, window=WINDOW)
    print(
        f"exact poles: {res.energies.size} in {time.perf_counter() - t0:.0f}s",
        flush=True,
    )

    order = np.argsort(res.energies.real)
    # `localization` is supplied because the overlap cannot see a state that has
    # left the box -- the c-product cancels the rotated tail by construction. On
    # H2+ that blindness hid 18 poles whose orbitals are larger than their grid.
    pairs = [
        pair_by_overlap(
            res.energies[i],
            res.states[i],
            basis,
            basis_complete=True,
            localization=real_weight(res.states[i], base),
        )
        for i in order
    ]
    return res, bo, basis, pairs, order


def main() -> None:
    res, _bo, basis, pairs, order = verify()

    base, _, _ = ecs_angle_family(
        lambda a: electronic_grid(angle_deg=a, **EL),
        lambda a: nuclear_grid(angle_deg=a, **NU),
        electronic_angles=EL_ANGLES,
        nuclear_angles=NU_ANGLES,
    )
    print(
        f"\n{'E_r (Ha)':>12} {'Gamma (Ha)':>11} {'level':>8} {'overlap':>8} "
        f"{'2nd':>8} {'2nd val':>8} {'shift(meV)':>11} {'real_wt':>8}  verdict"
    )
    tally: dict[str, int] = {}
    for p, i in zip(pairs, order, strict=True):
        tally[p.verdict] = tally.get(p.verdict, 0) + 1
        lvl = "-" if p.level is None else f"w_{p.level[1]}"
        second = "-" if p.second_level is None else f"w_{p.second_level[1]}"
        print(
            f"{p.pole_energy:>12.6f} {res.widths[i]:>11.3e} {lvl:>8} {p.overlap:>8.4f} "
            f"{second:>8} {p.second_overlap:>8.4f} {p.shift_mev:>11.3f} "
            f"{real_weight(res.states[i], base):>8.4f}  {p.verdict}"
        )
    print(f"\nverdict tally: {tally}")

    # The claim under test: does overlap agree with the sorted-index pairing
    # `exact_resonance_figures.py` and docs/physics/exact-2d-resonances.md use?
    level_e, keys = basis.flat()
    n = min(res.energies.size, level_e.size)
    print("\nindex pairing vs overlap pairing (the doc's table rests on the first):")
    disagreements = 0
    for k in range(n):
        by_index = keys[k]
        by_overlap = pairs[k].level
        mark = "" if by_index == by_overlap else "   <-- DISAGREE"
        if mark:
            disagreements += 1
        print(f"  pole {k}: index -> {by_index}, overlap -> {by_overlap}{mark}")
    print(f"disagreements: {disagreements}/{n}")

    by_energy = pair_one_to_one(res.energies.real, level_e.real, max_distance=0.05)
    print(f"\nHungarian (energy) pairing: {by_energy}")

    print("\nexact vs BO, on the OVERLAP pairing:")
    for p, i in zip(pairs, order, strict=True):
        if p.level is None:
            continue
        e_bo = basis[p.level].energy
        d_e = (res.energies[i].real - e_bo.real) * HARTREE_TO_EV * 1000.0
        d_g = (res.widths[i] - max(0.0, -2.0 * e_bo.imag)) * HARTREE_TO_EV * 1000.0
        print(f"  w_{p.level[1]}: dE = {d_e:+8.3f} meV   dGamma = {d_g:+8.3f} meV")


if __name__ == "__main__":
    main()
