"""N2 figures for the exact (non-BO) 2-D resonance states.

The N2-specific driver behind two committed figures, mirroring `ti_curve.py`'s
role for the cross-section figure: the generic plotting lives in the library
(`qscat.core.plot_resonance_levels`, `qscat.viz.plot_wavefunction_2d`) and only
the molecule, the grids and the seeds are chosen here.

Both figures follow the conventions of the published versions of this model
(M. Vana, doctoral thesis, Charles University 2017 — the level diagram of
Fig. 3.2 and the wave-function panels of Fig. 3.3), so they can be read against
them directly:

- levels: `V0(R)` and `E_res(R)` with the `Gamma(R)/2` envelope shaded, and the
  levels drawn as horizontal lines in the `E_res` well. The published figure
  shows the LCP levels alone; here the exact 2-D levels are overlaid on them,
  which is the comparison this capability exists to make.
- states: `R` horizontal, `r` vertical increasing downward, complex phase as
  hue and magnitude as brightness, with potential contours in light grey.

Run as::

    uv run python -m validation.n2.exact_resonance_figures

Costs a few minutes: three 2-D solves per seed on a ~47k-unknown grid.
"""

from __future__ import annotations

import pathlib
import time

import numpy as np
from qscat.core import (
    ExactResonanceStates,
    exact_resonance_states,
    plot_resonance_levels,
    resonance_levels,
)
from qscat.core.grids import electronic_grid, nuclear_grid
from qscat.dvr import TensorGrid
from qscat.model import N2
from qscat.viz import EquidistantProjector, plot_wavefunction_2d

FIGURES = pathlib.Path(__file__).resolve().parents[2] / "docs" / "physics" / "figures"

# Converged settings from docs/physics/exact-2d-resonances.md: the polynomial
# order converges at 8 and the nuclear grid at quadrature 12 / r_max 24. The
# electronic r_max = 24 used here keeps the figure affordable, and is NOT the
# converged box: sweeping it to 72 bohr leaves only v = 0's difference from the
# BO levels converged. The figure says so in its title.
EL = dict(r_max=24.0, order=8, n_complex=6)
NU = dict(r_max=24.0, n_complex=6, quadrature=12)
SEEDS = [
    -0.673960 - 0.0025j,
    -0.664231 - 0.0027j,
    -0.654640 - 0.0030j,
    -0.645187 - 0.0032j,
    -0.635863 - 0.0034j,
]
WINDOW = (-0.85, -0.45, -0.08, 0.0)

R_VIEW = (1.6, 3.0)  # the well, as published
E_VIEW = (-0.76, -0.55)
# The interaction zone: electron close in, R around equilibrium. The full
# electronic box would be mostly empty space plus the ECS tail.
STATE_EXTENT = ((0.0, 10.0), (1.7, 2.8))


def main(cache: pathlib.Path | None = None) -> None:
    el_a = electronic_grid(angle_deg=35.0, **EL)
    el_b = electronic_grid(angle_deg=44.0, **EL)
    nu_a = nuclear_grid(angle_deg=35.0, **NU)
    nu_b = nuclear_grid(angle_deg=30.0, **NU)
    base = TensorGrid([el_a, nu_a])

    t0 = time.perf_counter()
    bo, v_d, gamma = resonance_levels(
        N2, nu_a, nu_b, el_a, el_b, n_levels=len(SEEDS), return_curve=True
    )
    # The 2-D search is minutes of linear algebra and the figures get iterated
    # on; caching it keeps a restyle from costing a resolve. Delete the file to
    # force a recompute.
    if cache is not None and cache.exists():
        res = ExactResonanceStates.load(cache)
        print(f"loaded cached states from {cache}", flush=True)
    else:
        res = exact_resonance_states(
            N2,
            base,
            TensorGrid([el_b, nu_a]),
            TensorGrid([el_a, nu_b]),
            shifts=SEEDS,
            k=8,
            window=WINDOW,
        )
        if cache is not None:
            res.save(cache)
    print(
        f"n2d={el_a.n * nu_a.n}  exact={res.energies.size} levels  "
        f"BO={bo.energies.size}  {time.perf_counter() - t0:.0f}s",
        flush=True,
    )
    for e, g in zip(res.energies, res.widths, strict=True):
        print(f"  exact  E_r={e.real:+.9f}  Gamma={g:.9f}", flush=True)

    # The curves live on the nuclear grid's REAL region; the ECS tail is not
    # part of the picture (and its points are complex).
    pts = nu_a.points
    real = np.flatnonzero(np.abs(pts.imag) < 1e-12)
    R = pts[real].real
    keep = (R >= R_VIEW[0]) & (R <= R_VIEW[1])
    R = R[keep]
    e_res = v_d[real][keep].real
    half = 0.5 * gamma[real][keep]

    FIGURES.mkdir(parents=True, exist_ok=True)
    n = min(res.energies.size, bo.energies.size)
    plot_resonance_levels(
        {
            "exact 2-D": res.energies[:n],
            "BO / LCP": bo.energies[:n].real - 0.5j * bo.widths[:n],
        },
        curves={
            r"$V_0(R)$": (R, np.asarray(N2.v0(R), dtype=np.complex128).real),
            r"$E_{\rm res}(R)$": (R, e_res),
        },
        band=(R, e_res, half),
        xlim=R_VIEW,
        ylim=E_VIEW,
        baseline="BO / LCP",
        title=(
            "N$_2^-$ resonance levels $\\omega_i$: exact 2-D vs BO/LCP\n"
            "(electronic $r_{\\max}=24$; only $v=0$ is box-converged — see text)"
        ),
        path=FIGURES / "n2-exact-2d-resonance-levels.png",
    )

    projector = EquidistantProjector(base, samples=(420, 300), extent=STATE_EXTENT)
    for i in range(min(3, res.energies.size)):
        psi = res.states[i]
        e, g = res.energies[i], res.widths[i]
        plot_wavefunction_2d(
            projector,
            psi,
            # One global scale. Per-region scaling (`region_magnitudes`) does
            # reveal the outgoing tail, but at the price of visible seams at the
            # region boundaries and a brightness that no longer means the same
            # thing across the panel -- too much cost for a static figure whose
            # job is the shape of the trapped state. The tail's behaviour is
            # reported numerically in the note instead.
            mag=float(np.abs(projector.project(psi)).max()),
            path=FIGURES / f"n2-exact-2d-resonance-state-v{i}.png",
            title=rf"$v={i}$:  $E_r={e.real:.5f}$ Ha,  $\Gamma={g:.5f}$ Ha",
            xlabel=r"$R$ (bohr)",
            ylabel=r"$r$ (bohr)",
            inverse=True,  # print style: white ground, as published
            contours=6,
            potential=lambda r, R_: np.asarray(N2.surface(r, R_)).real,
            contour_field="magnitude",
            # The turning surface of each level: where the 2-D potential equals
            # that state's own energy, which is the boundary of the region it is
            # classically allowed to occupy.
            potential_levels=[float(w.real) for w in res.energies[:n]],
        )
    print(f"figures written to {FIGURES}", flush=True)


if __name__ == "__main__":
    main(cache=pathlib.Path(__file__).with_suffix(".cache.npz"))
