"""Task 3 (sub-project #9): the dense exact-2-D sigma_{0->v'}(E) curve for N2,
displayed against Houfek's golden `CSVE.V00.J00` via the generic
`projects.n2_2d_cross_section.cross_section_plot.plot_cross_sections`.

Sub-projects #3-#8 validated single ANCHOR points; this module is the first
to sweep the exact 2-D solver's OWN symbolic-reuse machinery (Tasks 1-2 of
this sub-project: `ve_cross_section_2d` factors once and refactors per
energy) across a dense energy grid and look at the resulting curve, not
just isolated coordinates.

`compute_ti_curve` mirrors `validation.n2.exact2d._build_system()` exactly:
`working_tgrid()` (the converged N~26857 grid) + `vibrational_states` with
the same `N_VIB` the anchor tables use, so results here are directly
comparable to `exact2d.compute_exact2d_results()` at any shared
`(energy, channel)` coordinate.

The layering rule is one-directional: `validation/` may import `projects/`
(several modules here do, `exact2d.py` among them), while `projects/` must
not import `validation/`. That is why the Houfek data stays on this side --
this module reads the golden file (`houfek_reference`) and hands the
resulting plain arrays to the generic, Houfek-agnostic
`plot_cross_sections`, which knows nothing about the reference.

Run: `uv run python -m validation.n2.ti_curve` (Task 3, Step 3) -- full
density over the working grid costs MINUTES (one sparse LU per energy on
the scipy/SuperLU backend; `SparseLU.refactor` re-runs `splu` from scratch
on this backend, see its docstring, so there is no per-energy speedup here
beyond skipping the Python-level setup) and is NOT run by the test suite.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from qscat.dvr import TensorGrid

from projects.n2_2d_cross_section.convergence import working_tgrid
from projects.n2_2d_cross_section.cross_section_2d import ve_cross_section_2d
from projects.n2_2d_cross_section.cross_section_plot import plot_cross_sections
from projects.n2_2d_cross_section.hamiltonian2d import MU
from projects.n2_ti_cross_section.vibrational import vibrational_states
from validation.figures import FIGURE_DIR
from validation.n2 import loader
from validation.n2.cross_section import N_VIB

__all__ = ["FIGURE_PATH", "NPZ_PATH", "compute_ti_curve", "houfek_reference", "main"]

FIGURE_PATH = FIGURE_DIR / "n2-2d-ti-cross-section.png"
NPZ_PATH = FIGURE_PATH.with_suffix(".npz")

# The channels the 6 `reference.anchors()` exercise (0=elastic, 1..3=VE):
# plotting all four lets `main()`'s figure show every anchor coordinate.
DEFAULT_VPRIMES = [0, 1, 2, 3]


def compute_ti_curve(
    E_grid: npt.NDArray[np.float64],
    vprimes: list[int],
    *,
    tgrid: TensorGrid | None = None,
) -> npt.NDArray[np.float64]:
    """sigma_{0->v'}(E) over `E_grid`, exact 2-D solver, v_init=0 (neutral
    N2's vibrational ground state).

    Builds (or accepts, via `tgrid`) the working 2-D tensor grid and its
    matching neutral-N2 vibrational states exactly like
    `exact2d._build_system()`, then makes ONE call to `ve_cross_section_2d`
    with the whole `E_grid` array so Tasks 1-2's analyze-once/refactor-per-
    energy reuse applies across the entire sweep, not per-point.

    Returns `sigma[len(E_grid), len(vprimes)]`, bohr^2, real, >=0.
    """
    if tgrid is None:
        tgrid = working_tgrid()
    eps, chi = vibrational_states(tgrid.grids[1], MU, N_VIB)
    sigma = ve_cross_section_2d(tgrid, eps, chi, 0, vprimes, E_grid)
    return np.asarray(sigma, dtype=np.float64)


def houfek_reference(
    vprimes: list[int],
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Houfek's `CSVE.V00.J00`, sliced to the requested v' columns.

    Returns `(E_ref, sigma_ref)`: `E_ref` is `(N,)` Hartree,
    `sigma_ref` is `(N, len(vprimes))` bohr^2, same column order as
    `vprimes` -- directly usable as `plot_cross_sections`'s `reference=`.
    """
    d = loader.load()
    return d.energy, d.sigma[:, vprimes]


def main() -> None:
    """Generate the committed figure + `.npz` (Task 3, Step 3).

    Dense `E_grid` over (0, 0.2] Ha at the full working-grid density -- one
    sparse LU factorization per energy point, ~3-4s each on this machine's
    scipy backend, so this is a several-minute run. `E` starts just above 0
    (not at exactly 0) since `ve_cross_section_2d` short-circuits `E<=0` to
    zero (no driven-equation solve below threshold, see its docstring) --
    there is no scattering solution to plot at E=0 itself.
    """
    vprimes = DEFAULT_VPRIMES
    E_grid = np.linspace(0.005, 0.2, 60)
    sigma = compute_ti_curve(E_grid, vprimes)
    e_ref, sigma_ref = houfek_reference(vprimes)

    FIGURE_PATH.parent.mkdir(parents=True, exist_ok=True)
    np.savez(NPZ_PATH, E_grid=E_grid, sigma=sigma, channels=np.array(vprimes))

    plot_cross_sections(
        E_grid,
        sigma,
        channels=vprimes,
        reference=(e_ref, sigma_ref),
        title="N2 exact 2-D VE cross section sigma_{0->v'}(E) vs Houfek CSVE.V00.J00",
        path=FIGURE_PATH,
    )
    print(f"wrote {FIGURE_PATH}")
    print(f"wrote {NPZ_PATH}")


if __name__ == "__main__":
    main()
