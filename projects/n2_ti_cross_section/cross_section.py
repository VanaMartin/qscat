"""Vibrationally-elastic/inelastic (VE) cross section via the
resolvent/driven-equation (sub-project #3, Task 3 -- THE CRUX).

`.superpowers/sdd/ti-cross-section-extraction.md` sections 3-4:

- Doorway: `d_v(R) = sqrt(Gamma(R)/(2*pi)) * chi_v(R)`.
- Driven equation, for collision energy `E` and initial vibrational channel
  `v_init`: `E_tot = E + eps_{v_init}`;
  `H_res = T_nuc(mu) + diag(V_d(R) - i*Gamma(R)/2)`;
  `M = E_tot*I - H_res`; solve `M @ xi = d_{v_init}` (`np.linalg.solve`).
- S-matrix: `S_{v'<-v_init} = <d_{v'} | xi> = sum_j d_{v'}[j] * xi[j]` --
  the **c-product** (no conjugate), matching `chi`/`eigen()`'s convention
  (see `vibrational.py`'s docstring): the DVR basis is already
  1/sqrt(weight)-normalized, so this inner product is a plain coefficient
  dot product with no extra quadrature weights, and no conjugation because
  `xi` is a genuinely complex ECS-driven solution (not an eigenvector pair
  needing a Hermitian norm). This was verified empirically here: the
  c-product convention gives real, non-negative `sigma`; the Hermitian
  (conjugated) convention does not (see `.superpowers/sdd/task-3-report.md`).
- Cross section: `sigma_{v_init->v'}(E) = 4*pi**3*|S|**2/(2*E)`, set to 0
  if the final channel is energetically closed (`E_tot - eps_{v'} <= 0`).

Efficiency: `xi` depends only on `(E, v_init)`, not on `v'`, so it is solved
for once per `E` and then reused (a single dot product per `v'`) for every
channel in `vprimes` -- the "reuse the LU across vprimes at fixed E" of the
task brief.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from qscat.dvr import FemDvrEcsGrid, kinetic

__all__ = ["ve_cross_section"]


def _sigma_at_one_energy(
    grid: FemDvrEcsGrid,
    mu: float,
    Vd: npt.NDArray[np.complex128],
    Gamma: npt.NDArray[np.float64],
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    vprimes: list[int],
    E: float,
) -> npt.NDArray[np.float64]:
    n = grid.n
    doorway = np.sqrt(Gamma / (2.0 * np.pi))[None, :] * chi  # (n_vib, n)

    T = kinetic(grid, mu)
    H_res = T + np.diag(Vd - 1j * Gamma / 2.0)

    E_tot = E + eps[v_init]
    M = E_tot * np.eye(n, dtype=np.complex128) - H_res
    xi = np.linalg.solve(M, doorway[v_init])

    sigma = np.zeros(len(vprimes), dtype=np.float64)
    if E <= 0.0:
        return sigma
    for k, vp in enumerate(vprimes):
        if E_tot - eps[vp] <= 0.0:
            continue  # closed channel
        S = np.dot(doorway[vp], xi)  # c-product: no conjugate
        sigma[k] = 4.0 * np.pi**3 * np.abs(S) ** 2 / (2.0 * E)
    return sigma


def ve_cross_section(
    grid: FemDvrEcsGrid,
    mu: float,
    Vd: npt.NDArray[np.complex128],
    Gamma: npt.NDArray[np.float64],
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    vprimes: list[int],
    E: float | npt.ArrayLike,
) -> npt.NDArray[np.float64]:
    """sigma_{v_init->v'}(E) (bohr^2) for each `v'` in `vprimes`.

    `E` (collision energy, Hartree) may be a scalar or an array; scalar `E`
    returns shape `(len(vprimes),)`, array `E` returns shape
    `(len(E), len(vprimes))`.
    """
    E_arr = np.atleast_1d(np.asarray(E, dtype=np.float64))
    out = np.stack(
        [
            _sigma_at_one_energy(grid, mu, Vd, Gamma, eps, chi, v_init, vprimes, float(e))
            for e in E_arr
        ]
    )
    if np.isscalar(E) or (isinstance(E, np.ndarray) and E.ndim == 0):
        return np.asarray(out[0], dtype=np.float64)
    return out
