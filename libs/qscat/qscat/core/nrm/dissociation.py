"""Dissociative attachment in the nonlocal resonance model.

The nuclear wave equation, Eq. (52), in DVR coefficient space:

    [E I - T_R - diag(V_d(R))] psi - F(E) psi = V_dk+(R) chi_vi(R)

and the cross section, Eq. (54):

    sigma_DA(E) = (2 pi^2 / k_i^2) (K_DA / mu) |Psi_d(R -> inf)|^2

`Psi_d(R -> inf)` is the wavefunction VALUE at the outermost real node,
`psi[b]/sqrt(w_b)`, not the DVR coefficient. That extraction is algebraically
identical to the one `qscat.core.lcp.lcp_da_cross_section` performs
(`sigma = 4 pi^3 |S_DA|^2 / 2E` with `S_DA = sqrt(K/2 pi mu) psi(X)` expands to
`pi^2 K |psi|^2 / (mu E)`, and so does Eq. 54 with `k_i^2 = 2E`), so the two
methods are compared on exactly the same footing.

`solve_nuclear` and `da_sigma_from_psi` are public precisely so the local limit
can be driven directly: `F -> diag(-(i/2) Gamma)` with the LCP's own doorway
must reproduce `lcp_da_cross_section`, which is this package's differential
oracle against already-validated code.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

from qscat.dvr import FemDvrEcsGrid, kinetic

from .coupling import v_dk_plus
from .discrete_state import DiscreteState
from .ingredients import NrmIngredients, nrm_ingredients
from .nonlocal_potential import continue_to_tail, nonlocal_operator

if TYPE_CHECKING:
    from qscat.model import ResonanceModel

__all__ = ["da_sigma_from_psi", "nrm_da_cross_section", "solve_nuclear"]


def _boundary_node(grid: FemDvrEcsGrid) -> int:
    """Index of the outermost REAL node -- the flux surface `X`."""
    real_idx = np.flatnonzero(grid.points.imag == 0.0)
    return int(real_idx[np.argmax(grid.points[real_idx].real)])


def solve_nuclear(
    nuclear_grid: FemDvrEcsGrid,
    mu: float,
    v_d_full: npt.NDArray[np.complex128],
    f_matrix: npt.NDArray[np.complex128],
    rhs: npt.NDArray[np.complex128],
    e_total: float,
) -> npt.NDArray[np.complex128]:
    """Solve Eq. (52) for `Psi_d+` in nuclear DVR coefficient space.

    Parameters
    ----------
    nuclear_grid : FemDvrEcsGrid
        The nuclear radial grid.
    mu : float
        Nuclear reduced mass.
    v_d_full : ndarray
        `V_d(R)` (Eq. 20) on the full nuclear grid, complex.
    f_matrix : ndarray
        `F(E)` from `nonlocal_operator`, or any operator standing in for it
        (a diagonal `-(i/2) Gamma` reproduces the LCP).
    rhs : ndarray
        The right-hand side `V_dk+(R) chi_vi(R)`.
    e_total : float
        Total energy `E` (hartree).

    Returns
    -------
    ndarray
        `Psi_d+` as DVR coefficients.
    """
    a = (
        e_total * np.eye(nuclear_grid.n, dtype=np.complex128)
        - kinetic(nuclear_grid, mu)
        - np.diag(v_d_full)
        - f_matrix
    )
    out: npt.NDArray[np.complex128] = np.linalg.solve(a, rhs)
    return out


def da_sigma_from_psi(
    nuclear_grid: FemDvrEcsGrid,
    mu: float,
    psi: npt.NDArray[np.complex128],
    e_total: float,
    eps_e: float,
    e_kin: float,
) -> float:
    """`sigma_DA` (bohr^2) from `Psi_d+` -- Eq. (54).

    Returns `0.0` for a closed channel (`e_kin <= 0` or `e_total <= eps_e`),
    matching `lcp_da_cross_section`'s convention.
    """
    if e_kin <= 0.0:
        return 0.0
    e_dr = e_total - eps_e
    if e_dr <= 0.0:
        return 0.0
    b = _boundary_node(nuclear_grid)
    value = psi[b] / np.sqrt(complex(nuclear_grid.weights[b]))
    k_da = float(np.sqrt(2.0 * mu * e_dr))
    k_i2 = 2.0 * e_kin
    return float(2.0 * np.pi**2 / k_i2 * (k_da / mu) * abs(value) ** 2)


def nrm_da_cross_section(
    nuclear_grid: FemDvrEcsGrid,
    elec_grid: FemDvrEcsGrid,
    model: ResonanceModel,
    phi_d: DiscreteState,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    E: float | npt.ArrayLike,
    *,
    ingredients: NrmIngredients | None = None,
    n_states: int | None = None,
) -> npt.NDArray[np.float64]:
    """`sigma_DA(E)` in the nonlocal resonance model (bohr^2).

    Parameters
    ----------
    nuclear_grid, elec_grid : FemDvrEcsGrid
        The nuclear and electronic radial grids.
    model : ResonanceModel
        The molecule.
    phi_d : DiscreteState
        The discrete-state choice under test.
    eps, chi : ndarray
        Neutral vibrational energies and states (`qscat.core.vibrational`).
    v_init : int
        Initial vibrational level.
    E : float or array
        Incident electron kinetic energy or energies (hartree).
    ingredients : NrmIngredients, optional
        Precomputed ingredients; built here if omitted. Pass them in when
        sweeping energies or comparing discrete-state choices -- they are
        energy-independent and dominate the cost.
    n_states : int, optional
        Truncate the sum over projected electronic states. `None` uses all.

    Returns
    -------
    ndarray
        `sigma_DA` per energy; scalar-shaped for a scalar `E`.
    """
    real = nuclear_grid.points.imag == 0.0
    R_desc = np.sort(nuclear_grid.points[real].real)[::-1]
    ing = ingredients or nrm_ingredients(elec_grid, model, phi_d, R_desc)

    v_d_full = continue_to_tail(ing.v_d_discrete, ing.R, nuclear_grid)
    eps_e = float(v_d_full[_boundary_node(nuclear_grid)].real)

    e_arr = np.atleast_1d(np.asarray(E, dtype=np.float64))
    out = np.zeros(e_arr.size, dtype=np.float64)
    for ie, e_kin in enumerate(e_arr):
        if float(e_kin) <= 0.0:
            continue
        e_total = float(e_kin) + float(eps[v_init])
        v_dk = v_dk_plus(elec_grid, model, phi_d, ing.R, float(e_kin))
        rhs = continue_to_tail(v_dk, ing.R, nuclear_grid) * chi[v_init]
        f = nonlocal_operator(ing, nuclear_grid, model, e_total, n_states=n_states)
        psi = solve_nuclear(nuclear_grid, model.mu, v_d_full, f, rhs, e_total)
        out[ie] = da_sigma_from_psi(nuclear_grid, model.mu, psi, e_total, eps_e, float(e_kin))

    scalar = np.isscalar(E) or (isinstance(E, np.ndarray) and np.ndim(E) == 0)
    return np.asarray(out[0] if scalar else out, dtype=np.float64)
