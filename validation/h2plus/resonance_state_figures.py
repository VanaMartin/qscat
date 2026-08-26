"""Exact 2-D resonance states of the H2+ model, and one against its BO twin.

The shift table in `exact_poles.py` says the Born-Oppenheimer levels sit up to
a few meV away from the exact poles. These figures show WHAT that shift is made
of: the exact state is a genuine 2-D object, and the Born-Oppenheimer
approximation asserts it factorizes as `phi_Ryj(r; R) * chi_v(R)`. Where the two
pictures disagree in position, they disagree in shape too.

Three states from window 0 (electron energy `E = e_tot - EPS0`):

- **The near-degenerate pair at E ~ 0.0055 Ha.** Two BO levels sit 20 uHa
  (0.5 meV) apart there -- `omega_1^9` (`Ry_9 v=1`) at 0.005509 and
  `omega_3^3` (`Ry_3 v=3`) at 0.005529 -- and the exact solver returns two
  poles 154 uHa (4.2 meV) apart, at 0.005507 and 0.005661. **The exact
  treatment splits the BO near-degeneracy about eightfold, and does it
  asymmetrically**: `omega_1^9` barely moves (-0.04 meV) while `omega_3^3` is
  pushed up 3.60 meV.

  The two panels show why. `omega_1^9` is a DIFFUSE Rydberg orbital -- ~9
  radial lobes reaching past 250 bohr -- with a single node in `R`. `omega_3^3`
  is COMPACT in `r` (under ~70 bohr, ~4 lobes) with three nodes in `R`,
  spreading instead along the nuclear coordinate to R ~ 4 bohr. A distant
  Rydberg electron follows the nuclei adiabatically and is nearly BO-exact; a
  compact low-`n` state overlapping the dissociative channel is not. Here the
  two regimes sit at the same energy and repel.

- **The state at E ~ 0.0039 Ha** -- pole at `e_tot` = -0.093680172, paired to
  the BO level `Ry_4 v=2` at -0.093564, shift **-3.154 meV**. This one gets a
  side-by-side with its BO product state, and the difference is directly
  visible: both carry the same core (5 radial lobes, two nodes in `R`), but the
  exact state has additional amplitude near `r ~ 130-160` bohr that the product
  state has nowhere -- admixture of a diffuse high-`n` Rydberg component. That
  extra lobe is what the -3.154 meV is made of.

Note the pairing here comes from the FULL BO level set (12 Rydberg curves on a
300-bohr electronic box), not the 5-curve enumeration `exact_poles.py`'s
campaign ran with, which could not hold an `n_eff >= 6` orbital and so did not
contain `omega_1^9` at all.

**The BO comparison is a construction, not a solver output.** `Psi_BO` is built
here as `phi_Ry4(r; R) * chi_2(R)`: the 5th electronic eigenvector at each
nuclear point, phase-aligned across R by continuity (the eigenvector's phase is
arbitrary per point, so a naive loop produces a state that flips sign at random
R and is meaningless to look at), times the vibrational function that curve
supports. That product is exactly what the Born-Oppenheimer approximation
claims the true state is.

Reduced mass: these are the repository's own results, so they use the shipped
`H2P.mu = 918.076` -- NOT the `918.25` the published sweep used. The
figures in `dr_levels_figure.py` are the ones that must match published curves;
these do not.

Run as::

    uv run python -m validation.h2plus.resonance_state_figures

Writes `docs/physics/figures/h2p-exact-2d-resonance-state-*.png`. The 2-D solve
is ~15 minutes at `r_max=300`; it is cached to a git-ignored `.npz` beside this
file, so restyling a figure does not re-solve. Delete that file to recompute.
"""

from __future__ import annotations

import pathlib
import time

import numpy as np
import numpy.typing as npt
from qscat.core import (
    ExactResonanceStates,
    bo_basis,
    electronic_curves,
    exact_resonance_states,
)
from qscat.dvr import TensorGrid
from qscat.model import H2P
from qscat.viz import EquidistantProjector, plot_wavefunction_2d

from validation.h2plus.exact_poles import EPS0, K_SEARCH, find_seeds, grid_family

FIGURES = pathlib.Path(__file__).resolve().parents[2] / "docs" / "physics" / "figures"
CACHE = pathlib.Path(__file__).with_suffix(".cache.npz")

R_MAX = 300.0

# The states are diffuse Rydberg orbitals (`n_eff` ~ 9-11, `<r>` ~ 133-176
# bohr), so the electronic window has to be wide enough to contain one; a
# 10-bohr view like N2's would show an empty core and nothing else. The nuclear
# window is the ion-core well around its R_e = 2.0 bohr equilibrium.
STATE_EXTENT = ((0.0, 260.0), (1.0, 5.0))

# Targets, as absolute energies (Ha). Chosen for what they show, not roundness.
PAIR = (-0.092096961, -0.091942982)  # E ~ 0.0055 Ha, 154 uHa apart
BO_PARTNERED = -0.093680172  # E ~ 0.0039 Ha, pairs with Ry_4 v=2 (-3.161 meV)
BO_CURVE, BO_VIB = 4, 2


def _solve() -> ExactResonanceStates:
    """Window-0 poles WITH their states, cached."""
    if CACHE.exists():
        print(f"loaded cached states from {CACHE}", flush=True)
        return ExactResonanceStates.load(CACHE)

    seeds = [s for s in find_seeds()[0] if s.window == 0]
    base, moved_el, moved_nu = grid_family(R_MAX)
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
    print(f"{res.energies.size} poles in {time.perf_counter() - t0:.0f}s", flush=True)
    res.save(CACHE)
    return res


def bo_product_state(tgrid: TensorGrid, curve: int, vib: int) -> npt.NDArray[np.complex128]:
    """`phi_Ry<curve>(r; R) * chi_<vib>(R)` -- what the BO picture asserts the
    exact state is.

    A thin call into `qscat.core.bo`, which owns the phase alignment this
    picture depends on: the electronic eigenvector's phase is arbitrary at each
    nuclear point, and without alignment the product flips sign at random `R`
    and the figure is noise rather than a wavefunction.
    """
    g_r, g_R = tgrid.grids
    cur = electronic_curves(H2P, g_r, g_R, n_curves=curve + 1, with_states=True)
    return bo_basis(cur, g_R, H2P.mu, n_vib=vib + 1, allow_partial=True)[(curve, vib)].psi


def main() -> None:
    res = _solve()
    base, _, _ = grid_family(R_MAX)
    projector = EquidistantProjector(base, samples=(520, 320), extent=STATE_EXTENT)
    FIGURES.mkdir(parents=True, exist_ok=True)

    def draw(psi: npt.NDArray[np.complex128], name: str, title: str, level: float) -> None:
        plot_wavefunction_2d(
            projector,
            psi,
            mag=float(np.abs(projector.project(psi)).max()),
            path=FIGURES / name,
            title=title,
            xlabel=r"$R$ (bohr)",
            ylabel=r"$r$ (bohr)",
            inverse=True,  # print style: white ground
            contours=5,
            contour_field="magnitude",
            potential=lambda r, R_: np.asarray(H2P.surface(r, R_)).real,
            # The state's OWN turning surface -- where the 2-D potential equals
            # its energy, i.e. the boundary of the region it is classically
            # allowed to occupy. Both `potential` and `potential_levels` are
            # required for the overlay to draw at all; passing the potential
            # alone silently does nothing.
            potential_levels=[level],
        )
        print(f"wrote {name}", flush=True)

    for tag, e_tot, label, bo_e, shift in (
        ("a", PAIR[0], r"$\omega_1^9$ (Ry$_9$, $v=1$)", 0.005509, -0.040),
        ("b", PAIR[1], r"$\omega_3^3$ (Ry$_3$, $v=3$)", 0.005529, +3.596),
    ):
        i = int(np.argmin(np.abs(res.energies.real - e_tot)))
        e, g = res.energies[i], res.widths[i]
        draw(
            res.states[i],
            f"h2p-exact-2d-resonance-state-pair-{tag}.png",
            rf"{label}: exact pole at $E={e.real - EPS0:.6f}$ Ha, "
            rf"BO at ${bo_e:.6f}$ — shift ${shift:+.2f}$ meV"
            "\n"
            rf"($\Gamma={g:.2e}$ Ha; the two BO levels are 0.5 meV apart, "
            rf"the two exact poles 4.2 meV)",
            float(e.real),
        )

    i = int(np.argmin(np.abs(res.energies.real - BO_PARTNERED)))
    e, g = res.energies[i], res.widths[i]
    draw(
        res.states[i],
        "h2p-exact-2d-resonance-state-vs-bo-exact.png",
        rf"EXACT 2-D pole, $E_{{\rm tot}}={e.real:.6f}$ Ha, $\Gamma={g:.2e}$ Ha"
        "\n"
        r"(note the extra amplitude near $r\sim130$–$160$ bohr — absent below)",
        float(e.real),
    )
    draw(
        bo_product_state(base, BO_CURVE, BO_VIB),
        "h2p-exact-2d-resonance-state-vs-bo-product.png",
        rf"BO product $\varphi_{{\rm Ry{BO_CURVE}}}(r;R)\,\chi_{{v={BO_VIB}}}(R)$"
        "\n"
        r"($E_{\rm BO}=-0.093564$ Ha — same core, no large-$r$ lobe: "
        r"that difference is the $-3.15$ meV)",
        -0.093564,
    )


if __name__ == "__main__":
    main()
