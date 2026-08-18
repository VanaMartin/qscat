"""Does the H2+ pole search converge on a small electronic box?

The cross-section deck needs 1300 bohr for the incident Coulomb wave. A
Rydberg resonance is closed-channel and localized, so it may not. Measured
here rather than assumed, because the answer sets the cost of the whole pole
campaign.

Result (measured 2026-08-17, laptop SuperLU, no MUMPS): **converges for the COMPACT
states probed here; NOT established for the states the pole campaign actually needs.**

Both seeds, and both converged poles, sit well below the ion's `v=0`
threshold `eps[0] = -0.097604` Ha (`vibrational_states(full_grid().grids[1],
H2P.mu, 5, H2P.v0).eps[0]`). Reading that binding as a hydrogenic Rydberg
series (`binding = 1/(2 n_eff^2)`, `<r> ~ 1.5 n_eff^2` bohr) gives:

| pole | E (Ha) | binding below eps[0] | n_eff | <r> ~ (bohr) |
|---|---|---|---|---|
| 1 | -0.113959946 | 0.016356 Ha | 5.53 | 46 |
| 2 | -0.112423647 | 0.014820 Ha | 5.81 | 51 |

A 150-bohr box is ~3x either state's Rydberg extent, so their being box-
independent across 150/300/600 bohr is EXPECTED, not informative: it shows
the probe correctly recovers a compact state on a box already several times
its size, not that a box this small suffices for the whole campaign.

The pole campaign's actual targets sit AT OR ABOVE `eps[0]`, attached to the `v >= 1`
cation thresholds -- but "binding" for THOSE states is measured against the
NEXT ion vibrational threshold above the state's `e_tot`, not against
`eps[0]`, and the threshold spacing is NOT the ~2.8 mHa an earlier version
of this docstring assumed -- that number is wrong by about 3x; the measured
spacings follow. The six lowest ion-core
levels (`vibrational_states(full_grid().grids[1], H2P.mu, 6, H2P.v0).eps`):

    -0.097604  -0.087802  -0.078519  -0.069754  -0.061507  -0.053780  Ha

spaced 9.80, 9.28, 8.77, 8.25, 7.73 mHa apart -- decreasing (anharmonic),
averaging ~8.8 mHa, about 3x the earlier assumed value. Locating actual
levels inside the published DR windows via `rydberg_levels(H2P,
proxy_grid().grids[0], full_grid().grids[1], n_curves=4, n_vib=5)` (electron
energy `E = e_tot - eps[0]`):

| window | curve | vib | E (Ha) | next threshold above `e_tot` | binding | n_eff | <r> ~ (bohr) |
|---|---|---|---|---|---|---|---|
| 0 | 3 | 3 | 0.005529 | `eps[1]=-0.087802` | 4.27 mHa | 10.82 | 176 |
| 1 | 3 | 4 | 0.013454 | `eps[2]=-0.078519` | 5.63 mHa | 9.42 | 133 |

So DR-window states are `n_eff` ~ 9-11, `<r>` ~ 133-176 bohr -- not the
13-25 / 270-940 bohr this docstring previously claimed. (Only two levels
appear because `n_vib=5` is capped by the bound spectrum on the REDUCED
electronic deck used for this cheap lookup; the count is a floor, but the
`n_eff` scale is representative and it is the scale that sets the box.)

**Governing scope decision (repository owner):** states close to a
vibrational threshold are expected to be
poorly described and are DELIBERATELY OMITTED from the pole campaign via a
Rydberg-index cutoff -- not because they are intractable, but because they
are not isolated resonances in the first place: level spacing falls as
`n^-3`, so peaks merge before widths do, and a "resonance" identity stops
being meaningful before the numerics would. Concretely: **cut at
`n_eff <= 12`**, equivalently exclude any level bound by less than ~3.5 mHa
to the next threshold. At that cutoff the extreme case is `<r> ~ 216` bohr
(`1.5 * 12^2`), so a converged box wants roughly 650-900 bohr (~3x the
extreme `<r>`, the same margin `r_max=150` gave the compact states above).

**This module does not establish that `r_max=150` (or 600) suffices for
the pole campaign.** What would: seeding INSIDE the DR window BELOW the `n_eff <= 12`
cutoff (not the deeper curve-1/2 levels used here) and sweeping
300/600/1200 bohr -- which brackets the ~650-900-bohr requirement, so 600
is EXPECTED to suffice, though that expectation is not itself measured --
WITH the production nuclear grid (not this probe's reduced one; a
near-threshold state's width is nuclear-continuum-coupling-sensitive, see
the width caveat below). That is almost certainly a MUMPS/compute-host job: the
production nuclear grid alone (n=818 vs. this probe's n=181) already pushes
the 2-D size well past what ran here, on top of a larger electronic box.

Residuals for the two converged poles, across all three boxes -- the numbers
the width discussion below hinges on:

| box (bohr) | E (Ha) | Gamma (Ha) | res_electronic | res_nuclear |
|---|---|---|---|---|
| 150 | -0.113959946 | 2.936e-06 | 1.8e-17 | 5.9e-08 |
| 300 | -0.113959946 | 2.936e-06 | 4.1e-20 | 5.9e-08 |
| 600 | -0.113959946 | 2.936e-06 | 2.3e-19 | 5.9e-08 |
| 150 | -0.112423647 | 5.380e-10 | 2.6e-17 | 8.6e-12 |
| 300 | -0.112423647 | 5.380e-10 | 1.4e-18 | 8.6e-12 |
| 600 | -0.112423647 | 5.380e-10 | 6.1e-19 | 8.6e-12 |

**Positions (Re E) are well established**: matched to the 9 decimals printed
across a 4x box range, residuals 8-11 orders of magnitude below the
acceptance tolerance. **Widths are not.** `exact_resonance_states`'s default
acceptance is `|E_a - E_b| < max(rel_tol*|E_a|, atol)`, `rel_tol=1e-4`,
`atol=1e-8`; at `|E| ~ 0.112-0.114` that ceiling is ~1.1e-5 Ha. Pole 1's
`Gamma = 2.936e-06` Ha is only ~4x below that ceiling; pole 2's
`Gamma = 5.380e-10` Ha is ~2e4x below it. A state whose TRUE width sat
anywhere under that ~1.1e-5 Ha ceiling would pass the identical acceptance
test, so passing alone does not certify a width at the 1e-6 or 1e-10 level --
only the actual residuals (table above) argue for it, and those come from a
nuclear grid this probe explicitly does NOT claim converged (`n_complex=3`
vs. `full_grid()`'s production 25 -- see "Cost" below): width is a nuclear-
continuum-coupling-sensitive quantity, so a claim resting on an unconverged
nuclear tail is not a settled physical width. This repo's own convention
elsewhere (`docs/physics/lcp-resonance-levels.md`) is that widths below
~1e-6 Ha are noise, which brackets pole 2 squarely and leaves pole 1 only
marginally above the floor. (`Gamma` is also only matched to the 3
significant figures printed above, not to 9 decimals like `E` -- that alone
is not evidence of width convergence at the residual's own precision.)
**Read this module as establishing POSITIONS for the compact states probed,
with widths reported at or below this probe's resolution limit -- not as
confirmed physical widths, and not as evidence for the pole-campaign boxes.**

Two corrections to this probe's original design, both load-bearing:

1. `fem_grid_exp_tail`'s real signature takes `real_segments` (a list of
   `(n_elements, endpoint)` pairs), not `r_max`/`angle_deg` alone -- see
   `validation/h2plus/config.py`. `_electronic_box` below reuses
   `full_grid()`'s inner segments (0.1/0.3/1.0/4.0 bohr resolution out to
   100 bohr, where `v_int` lives) and only varies where the outer segment
   /ECS pivot ends, so the box sweep changes ONLY `r_max`, nothing else.

2. **The originally specified `WINDOW = (0.0, 0.05, -0.01, 0.0)` is the wrong energy
   frame and would have searched pure electron continuum.** `model.hamiltonian`
   returns the FULL system's absolute energy (electron + nuclear, in the
   same units `v0`/`surface` are built from); its bound/quasi-bound spectrum
   is negative throughout, because the Morse ion-core curve `v0(R)` <= 0
   everywhere and the electron-core interaction (`v_int` + the `charge/r`
   Coulomb tail) is attractive. Measured directly: `rydberg_levels`'s curve 1
   (n=2 Rydberg) vibrational levels run -0.161..-0.135 Ha and curve 2 (n=3
   Rydberg) -0.132..-0.098 Ha (`v = 0..4`, on `proxy_grid()`'s electronic grid
   x `full_grid()`'s nuclear grid); the ion-core-alone levels (`v0` diagonalized
   by itself) are -0.098..-0.062 Ha. The config module's `E_LO`/`E_HI` =
   [0.001, 0.050] Ha is the electron's KINETIC energy above the ion's `v=0`
   threshold (`E = e_tot - eps[v_init]`, team note), so in the SAME absolute
   frame as `model.hamiltonian` that is `e_tot` in about [-0.097, -0.048] Ha --
   negative, not positive. `WINDOW`/`_seeds()` below are set from the
   measured curve energies instead of that placeholder.

Cost (measured, this reduced deck -- see `_nuclear_box`'s `n_complex=3`,
`quadrature=8` below, NOT `full_grid()`'s production nuclear grid): one
`near()` call at `r_max=150` (n=108781) took 30.3s; at `r_max=600` (n=165796),
41.9s -- far below the naive projection from N2's own factorization-time
scaling (a first, now-abandoned timing probe at a WRONG seed inside the
electron continuum, where ARPACK has to fight much higher local spectral
density, took 447s at a comparable size -- the fix was the seed, not the
grid). The reduced nuclear grid (n=181, vs. `full_grid()`'s production
nuclear grid at n=818) is deliberate for this probe: it targets the
ELECTRONIC box only, so the nuclear grid is held fixed and cheap throughout,
not claimed converged in its own right. The full 3-box sweep (18
factorizations total) ran in 671s (~11.2 minutes): 188.9s (150 bohr, n2d=
108781) + 219.3s (300 bohr, n2d=127786) + 262.9s (600 bohr, n2d=165796).

Seed 1 (curve 1 v=4, -0.135194 Ha) returned only near-real eigenvalues
(`|Im E|` ~ 1e-16 to 1e-18) -- true BOUND states of the full 2-D Hamiltonian:
this part of curve 1's ladder sits well below any open dissociation or
ionization channel, so there is nothing for those states to decay into.
Seed 2 (curve 2 v=2, -0.113954 Ha) also turned up the two states tabulated
above, whose imaginary parts are nonzero and angle-stable to high relative
precision -- worth further investigation on a converged nuclear grid, but
not (per the width caveat above) a certified physical width from this probe
alone.

The near-degenerate cluster around both seeds SHRANK as the box grew (8
states at 150 bohr, 6 at 300, 3 at 600), even though positions that do
recur across boxes match closely. Read the direction correctly: as `r_max`
grows, `n_r` grows, and MORE states crowd into the fixed `k=8` window
around each seed -- the true near-seed states get pushed out of a budget
that does not grow with them. That is a WARNING that the search is
under-resourced, not a harmless wobble, and it warns exactly where the pole campaign
will feel it most: Rydberg state density diverges approaching the
ionization/dissociation threshold, so a fixed small `k` will saturate
FASTER there than it did in this probe's deeper, less-crowded region.
The pole campaign should budget `k` generously (or scan it) rather than reuse `k=8`
unexamined.

Reproducibility: the seeds below are NOT hardcoded numbers -- `_seeds()`
calls `rydberg_levels` directly, so
re-running this module regenerates them (~1 extra minute on top of the
18-factorization sweep) rather than trusting values nobody can check without
re-deriving them by hand. Re-running `_seeds()` alone (cheap, no
`exact_resonance_states` factorizations) reproduces `-0.13519369779747126`
and `-0.11395404995874556` bit-for-bit.

Three readings this module once carried are WITHDRAWN, recorded so they are
not re-derived from the same output:

1. **"Converges already at the smallest box tried", i.e. `r_max=150` suffices.**
   It does not follow: the measurement covers COMPACT states well below the
   ion's `v=0` threshold, and says nothing about the diffuse near-threshold
   states the DR windows are made of. The sections above replace it.
2. **The two nonzero-width poles as "actual poles of the full 2-D S-matrix".**
   Their widths sit at or below this probe's resolution limit, so the positions
   are established and the widths are not.
3. **`n_eff` = 13.4/18.3/25.0 and `<r>` ~ 270/500/940 bohr for the DR-window
   states.** This assumed a ~2.8 mHa ion vibrational threshold spacing that was
   never measured and is wrong by about 3x; measured directly, the spacings
   average ~8.8 mHa, which puts the same states at `n_eff` ~ 9-11 and
   `<r>` ~ 133-176 bohr. The estimate above uses the measured thresholds.

All three were corrected from the SAME 671 s of `exact_resonance_states`
output -- no sweep was re-run, only the interpretation and the residuals that
had been omitted. The corrected `n_eff`/`<r>` numbers were re-derived from the
cheap `vibrational_states`/`rydberg_levels` calls before being written down.
"""

from __future__ import annotations

import time

import numpy as np
from qscat.core import ecs_angle_family, exact_resonance_states
from qscat.core.grids import fem_grid_exp_tail, nuclear_grid
from qscat.dvr import FemDvrEcsGrid
from qscat.model import H2P

from validation.h2plus.config import full_grid, proxy_grid
from validation.h2plus.rydberg_levels import rydberg_levels

# WINDOW brackets the measured curve-1/curve-2 vibrational ladder (-0.161 to
# -0.098 Ha) plus margin either side, in the model's absolute (electron +
# nuclear) energy frame -- see the module docstring for why this is negative,
# not the originally specified placeholder [0.0, 0.05].
WINDOW = (-0.16, -0.03, -0.01, 0.0)


def _seeds() -> list[complex]:
    """Curve 1 v=4 and curve 2 v=2, from `rydberg_levels` directly (`rydberg_levels`) --
    NOT hardcoded, so re-running this module regenerates them. Spread across
    two curves rather than clustered in one, so the k=8 shift-invert windows
    sample different parts of the ladder instead of mostly overlapping. See
    the module docstring's "Reproducibility" note for the exact values this
    reproduces.
    """
    lv = rydberg_levels(H2P, proxy_grid().grids[0], full_grid().grids[1], n_curves=3, n_vib=5)
    return [complex(lv.energies[1, 4], -1e-4), complex(lv.energies[2, 2], -1e-4)]


# Electronic ECS angle pair (base/moved): 30/40 deg, both comfortably under
# the model's own 45 deg bound (`exp(-r^2/3)` in `v_int`, `2*theta < pi/2`).
_EL_BASE_DEG = 30.0
_EL_MOVED_DEG = 40.0

# Nuclear ECS angle pair (base/moved): 18/12 deg, both under H2P's
# `max_nuclear_ecs_angle_deg = 22.5` (the `a3*R**4` term in `v_int` diverges
# at or above that under rotation). Held fixed across the r_max sweep below
# -- this probe is about the ELECTRONIC box only.
_NUC_BASE_DEG = 18.0
_NUC_MOVED_DEG = 12.0


def _electronic_box(r_max: float, angle_deg: float) -> FemDvrEcsGrid:
    """H2+ electronic grid at box `r_max`, `full_grid()`'s inner-segment
    resolution (0.1/0.3/1.0/4.0 bohr out to 100 bohr) with only the outer
    segment (100 -> r_max) and ECS pivot varying. `r_max` must be a multiple
    of 10 bohr above 100 so the outer segment divides exactly, matching
    `full_grid()`'s own 10-bohr outer element length.
    """
    n_outer = round((r_max - 100.0) / 10.0)
    return fem_grid_exp_tail(
        [(10, 1.0), (10, 4.0), (16, 20.0), (20, 100.0), (n_outer, r_max)],
        angle_deg=angle_deg,
        quadrature=8,
        tail_n=25,
    )


def _nuclear_box(angle_deg: float) -> FemDvrEcsGrid:
    """Reduced nuclear grid, fixed across the r_max sweep -- see the module
    docstring for why this is deliberately cheaper than `full_grid()`'s
    production nuclear deck (n=818): this probe targets the electronic box
    only, so the nuclear side just needs to be present and angle-movable, not
    itself converged.
    """
    return nuclear_grid(angle_deg=angle_deg, r_max=14.0, n_complex=3, quadrature=8)


def main() -> None:
    seeds = _seeds()
    print(f"seeds: {[f'{s:.6f}' for s in seeds]}", flush=True)

    results: dict[float, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = {}

    for r_max in (150.0, 300.0, 600.0):
        # `ecs_angle_family` enforces the invariant this probe rests on: each
        # partner grid moves exactly one ECS angle and shares every real node.
        grid_base, grid_electronic, grid_nuclear = ecs_angle_family(
            lambda a, _r=r_max: _electronic_box(_r, a),
            _nuclear_box,
            electronic_angles=(_EL_BASE_DEG, _EL_MOVED_DEG),
            nuclear_angles=(_NUC_BASE_DEG, _NUC_MOVED_DEG),
        )

        t0 = time.perf_counter()
        res = exact_resonance_states(
            H2P,
            grid_base,
            grid_electronic,
            grid_nuclear,
            shifts=seeds,
            window=WINDOW,
            k=8,
        )
        dt = time.perf_counter() - t0

        el_base, nu_base = grid_base.grids
        n2d = el_base.n * nu_base.n
        print(
            f"r_max={r_max:6.0f}  n_r={el_base.n:5d}  n_R={nu_base.n:4d}  "
            f"n2d={n2d:8d}  found={res.energies.size}  {dt:.1f}s",
            flush=True,
        )
        for e, g, re_, rn in zip(
            res.energies,
            res.widths,
            res.residual_electronic,
            res.residual_nuclear,
            strict=True,
        ):
            print(
                f"    E={e.real:+.9f}  G={g:.3e}  res_el={re_:.1e}  res_nuc={rn:.1e}",
                flush=True,
            )
        results[r_max] = (
            res.energies,
            res.widths,
            res.residual_electronic,
            res.residual_nuclear,
        )

    # Boxes need not find the SAME NUMBER of states (a fixed k=8/2-seed
    # shift-invert window can pick up a different subset of a
    # near-degenerate, essentially-zero-width bound cluster as the box grows
    # -- see the module docstring's saturation warning). Match by nearest
    # energy instead of requiring equal set sizes, chaining 150 -> 300 -> 600
    # so all three boxes (not just the endpoints) inform the convergence
    # read -- the genuinely recurring states, in particular any narrow,
    # nonzero-width resonance, still get a convergence number either way.
    e150, e300, e600 = results[150.0][0], results[300.0][0], results[600.0][0]
    if e150.size and e300.size and e600.size:
        for e in e150:
            j300 = int(np.argmin(np.abs(e300 - e)))
            e_300 = e300[j300]
            j600 = int(np.argmin(np.abs(e600 - e_300)))
            e_600 = e600[j600]
            print(
                f"E(150)={e.real:+.9f}  nearest E(300)={e_300.real:+.9f} "
                f"(|diff|={abs(e - e_300):.3e})  nearest E(600)={e_600.real:+.9f} "
                f"(|diff|={abs(e_300 - e_600):.3e})",
                flush=True,
            )


if __name__ == "__main__":
    main()
