"""Nuclear-coordinate density: exact 2-D driven solution vs the 1-D LCP
driven solution, at the same collision energy (sub-project #6, Task 6 --
EXPLORATORY, no pass/fail on the physics).

An integrated cross section averages away much of what distinguishes an
exact treatment from an approximate one. This module compares the
*nuclear-coordinate density* of the driven solution instead:

  rho_exact(R)  = integral |Psi^(+)(r,R)|^2 dr     (2-D exact solver,
                  sub-project #6's `ve_cross_section_2d`, projected over the
                  electronic coordinate)
  rho_lcp(R)    = |xi(R)|^2                        (1-D Local Complex
                  Potential driven solution, sub-project #3's
                  `ve_cross_section`, exposed here as `xi` itself)

Both are masked to the UNSCALED region in every ECS coordinate before use:
under exterior complex scaling the complex tail carries outgoing flux, not
probability density, so summing/plotting it as if it were density is
physically meaningless (same masking discipline `TensorGrid.real_mask()`
and `channel_vector` already enforce for the driving/channel terms).

`rho_exact` and `rho_lcp` are NOT the same object dimensionally: one is a
2-D projected density (units of 1/length after the r-integral), the other a
raw 1-D driven-solution intensity. `compare_to_lcp` therefore normalizes
both to unit area over R before comparing shapes (centroid, RMS width) --
only the SHAPE comparison is meaningful, never an absolute-scale one.

DVR-coefficient convention (see `cross_section_2d.py`'s and
`vibrational.py`'s module docstrings): basis functions are pre-normalized by
`1/sqrt(w)`, so the coefficient magnitude squared `|Psi^(+)_ij|^2` already
IS the physical density weight at grid point `(r_i, R_j)` -- no extra
quadrature-weight factor is applied when summing over `i`. This was
sanity-checked here by requiring `nuclear_density`'s output to integrate
(via `np.trapezoid` over the real `R` points) to a finite, strictly positive
number for a genuine driven solution -- see `test_nuclear_density.py`.
`np.trapz` was deprecated in numpy 2.0 and is unavailable here (numpy
2.5.1 pinned); `np.trapezoid` is used throughout instead.

The 1-D LCP driven solution `xi(R)` faithfully reproduces sub-project #3's
`ve_cross_section` (`projects/n2_ti_cross_section/cross_section.py`): that
module's `_sigma_at_one_energy` solves `M @ xi = doorway[v_init]` internally
but returns only `sigma`, never `xi`. No helper anywhere already exposes
`xi`, so `lcp_driven_solution` below inlines the SAME minimal driven solve
using the SAME building blocks (`vres_on_grid`, `kinetic`, the doorway
`sqrt(Gamma/2pi) * chi_v`) documented in that module's own docstring --
`cross_section.py` itself is not modified.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from qscat.dvr import FemDvrEcsGrid, TensorGrid, kinetic

from projects.n2_2d_cross_section.convergence import working_tgrid
from projects.n2_2d_cross_section.cross_section_2d import ve_cross_section_2d
from projects.n2_2d_cross_section.hamiltonian2d import MU
from projects.n2_ti_cross_section.vibrational import vibrational_states
from projects.n2_ti_cross_section.vres import vres_on_grid

__all__ = ["compare_to_lcp", "lcp_driven_solution", "nuclear_density"]

# Enough bound vibrational states to cover any reasonable v_init, matching
# the N_VIB used elsewhere in this sub-project (`convergence.STUDY_VP`'s
# sibling constant); grown if a caller asks for a higher v_init.
_N_VIB_DEFAULT = 4


def nuclear_density(
    tgrid: TensorGrid, psi: npt.ArrayLike
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """`rho(R_j) = sum_i |Psi_ij|^2` over electronic index `i`, restricted to
    the UNSCALED region in BOTH coordinates.

    `psi` is the flat driven solution on `tgrid` (C order; axis 0 =
    electronic r, axis 1 = nuclear R -- `ve_cross_section_2d`'s
    `return_wavefunction=True` convention). Returns `(R_real, density)`:
    `R_real` the unscaled nuclear grid points (ascending), `density` the
    real, non-negative projected density there. The electronic index is
    restricted to unscaled `r` BEFORE summing (the complex tail there
    carries outgoing flux, not density); the nuclear index is restricted to
    unscaled `R` AFTER summing (same reason, the other axis).
    """
    g_r, g_R = tgrid.grids
    psi2d = np.asarray(psi, dtype=np.complex128).reshape(tgrid.shape)

    r_mask = g_r.real_points <= g_r.R0
    R_mask = g_R.real_points <= g_R.R0

    density_all_R = np.sum(np.abs(psi2d[r_mask, :]) ** 2, axis=0)
    density = np.asarray(density_all_R[R_mask], dtype=np.float64)
    R_real = np.asarray(g_R.real_points[R_mask], dtype=np.float64)
    return R_real, density


def lcp_driven_solution(
    grid: FemDvrEcsGrid,
    mu: float,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    E: float,
) -> npt.NDArray[np.complex128]:
    """The 1-D LCP driven solution `xi(R)` at collision energy `E`, channel
    `v_init` -- the same `xi` sub-project #3's `ve_cross_section` solves
    internally (`M @ xi = doorway[v_init]`, `M = E_tot*I - H_res`) but never
    returns. Inlined here from the SAME building blocks
    (`vres_on_grid`, `kinetic`, the doorway `sqrt(Gamma/2pi) * chi_v`)
    documented in `cross_section.py`'s module docstring, so this is a
    faithful re-exposure, not a reimplementation with different physics.
    """
    Vd, Gamma = vres_on_grid(grid)
    doorway_v_init = np.sqrt(Gamma / (2.0 * np.pi)) * chi[v_init]

    T = kinetic(grid, mu)
    H_res = T + np.diag(Vd - 1j * Gamma / 2.0)

    e_tot = E + eps[v_init]
    M = e_tot * np.eye(grid.n, dtype=np.complex128) - H_res
    xi = np.linalg.solve(M, doorway_v_init)
    return np.asarray(xi, dtype=np.complex128)


def _normalize(R: npt.NDArray[np.float64], rho: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """`rho / integral(rho dR)`, unit area over the (ascending) real grid `R`."""
    area = np.trapezoid(rho, R)
    return np.asarray(rho / area, dtype=np.float64)


def _centroid_and_width(
    R: npt.NDArray[np.float64], rho: npt.NDArray[np.float64]
) -> tuple[float, float]:
    """`<R>` and RMS width `sqrt(<(R - <R>)^2>)`, `rho` assumed unit area."""
    centroid = float(np.trapezoid(R * rho, R))
    variance = float(np.trapezoid((R - centroid) ** 2 * rho, R))
    return centroid, float(np.sqrt(variance))


def compare_to_lcp(E: float, v_init: int = 0) -> dict[str, object]:
    """Compare the exact 2-D nuclear density to the 1-D LCP `|xi(R)|^2` at
    the same collision energy `E` (Hartree) and initial channel `v_init`.

    Both densities are normalized to unit area over `R` first -- they are
    not the same object dimensionally (a 2-D projected density vs. a raw 1-D
    driven-solution intensity), so only their SHAPES are comparable. This is
    EXPLORATORY: the returned centroid/width numbers are reported as data,
    with no pass/fail threshold on their difference.

    Returns a dict with both `R` grids (identical arrays here, since both
    solutions are evaluated on the SAME nuclear grid, `working_tgrid()`'s),
    both unit-area-normalized densities, and each one's centroid `<R>` and
    RMS width in `R`.
    """
    tgrid = working_tgrid()
    n_vib = max(_N_VIB_DEFAULT, v_init + 1)
    eps, chi = vibrational_states(tgrid.grids[1], MU, n_vib)

    _, psi_plus = ve_cross_section_2d(
        tgrid, eps, chi, v_init, [v_init], E, return_wavefunction=True
    )
    if not isinstance(psi_plus, np.ndarray):
        # `E` is a scalar here by construction, so `psi_plus` is either the
        # single wavefunction array or `None` (below threshold) -- never the
        # `list[...]` shape `ve_cross_section_2d` uses for array `E`. This
        # also narrows the type for the call below.
        raise ValueError(
            f"compare_to_lcp: E={E} Ha is at or below threshold; there is no "
            "driven-equation solution to project a density from"
        )

    R_exact, rho_exact = nuclear_density(tgrid, psi_plus)
    rho_exact_n = _normalize(R_exact, rho_exact)
    centroid_exact, width_exact = _centroid_and_width(R_exact, rho_exact_n)

    g_R = tgrid.grids[1]
    xi = lcp_driven_solution(g_R, MU, eps, chi, v_init, E)
    R_mask = g_R.real_points <= g_R.R0
    R_lcp = np.asarray(g_R.real_points[R_mask], dtype=np.float64)
    rho_lcp = np.asarray(np.abs(xi[R_mask]) ** 2, dtype=np.float64)
    rho_lcp_n = _normalize(R_lcp, rho_lcp)
    centroid_lcp, width_lcp = _centroid_and_width(R_lcp, rho_lcp_n)

    return {
        "E": E,
        "v_init": v_init,
        "R_exact": R_exact,
        "density_exact": rho_exact_n,
        "R_lcp": R_lcp,
        "density_lcp": rho_lcp_n,
        "centroid_exact": centroid_exact,
        "width_exact": width_exact,
        "centroid_lcp": centroid_lcp,
        "width_lcp": width_lcp,
    }
