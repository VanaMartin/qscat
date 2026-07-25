"""Time-domain propagation engine: propagate Psi(0) under H_2D and sample.

Two cadences (the numeric-output design):
  * `c_{v'}(t_n)` -- recorded at EVERY step; the raw material of the
    Tannor-Weeks energy transform and the literal "formation from the
    correlation functions".
  * density/norm snapshots -- recorded on a COARSE schedule (`sample_period`
    steps, or explicit `snapshot_times`), so the wavefunction is observed at
    static points without storing every step.

`H_2D` is time-independent, so the sparse Crank-Nicolson factorization is built
once and reused; under the ECS contour `||Psi||` decays as outgoing flux is
absorbed (the resonance depletes). Norm here is the Hermitian L2 norm
`np.linalg.norm(psi)` -- the physical remaining-probability diagnostic, real
and provably monotone non-increasing under CN with an absorbing (ECS) H. The
c-product is a different object, reserved for the correlations `c_{v'}(t)`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from qscat.dvr import TensorGrid
from qscat.evolution import make_sparse_cn_stepper
from qscat.linalg import c_product

from projects.n2_2d_cross_section.hamiltonian2d import build_h2d

__all__ = ["Snapshot", "PropagationResult", "propagate"]


@dataclass(frozen=True)
class Snapshot:
    time: float
    rho_R: npt.NDArray[np.float64]  # nuclear density, sum_r |Psi|^2 (unscaled)
    rho_r: npt.NDArray[np.float64]  # electronic density, sum_R |Psi|^2 (unscaled)
    psi: npt.NDArray[np.complex128] | None  # full state, only if requested


@dataclass(frozen=True)
class PropagationResult:
    t: npt.NDArray[np.float64]  # (n_t,)  sample times n*dt
    c: npt.NDArray[np.complex128]  # (n_t, n_channels)  c_{v'}(t_n)
    norm: npt.NDArray[np.float64]  # (n_t,)  np.linalg.norm(psi) -- Hermitian L2
    snapshots: list[Snapshot]


def _densities(
    tgrid: TensorGrid, psi: npt.NDArray[np.complex128]
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    block = psi.reshape(tgrid.shape)
    dens = np.abs(block) ** 2
    r_real = tgrid.grids[0].real_points <= tgrid.grids[0].R0
    R_real = tgrid.grids[1].real_points <= tgrid.grids[1].R0
    rho_R = dens[r_real, :].sum(axis=0)
    rho_r = dens[:, R_real].sum(axis=1)
    return rho_R.astype(np.float64), rho_r.astype(np.float64)


def propagate(
    tgrid: TensorGrid,
    psi0: npt.NDArray[np.complex128],
    out_channels: list[npt.NDArray[np.complex128]],
    *,
    dt: float,
    n_steps: int,
    sample_period: int = 0,
    snapshot_times: list[float] | None = None,
    keep_psi_at: list[float] | None = None,
) -> PropagationResult:
    """Propagate and sample. See module docstring for the two cadences."""
    H = build_h2d(tgrid)
    step = make_sparse_cn_stepper(H, dt)

    n_t = n_steps + 1
    t = np.arange(n_t, dtype=np.float64) * dt
    n_ch = len(out_channels)
    c = np.empty((n_t, n_ch), dtype=np.complex128)
    norm = np.empty(n_t, dtype=np.float64)

    if snapshot_times is not None:
        snap_set = {round(x / dt) for x in snapshot_times}
    elif sample_period > 0:
        snap_set = set(range(0, n_t, sample_period)) | {n_t - 1}
    else:
        snap_set = {0, n_t - 1}
    keep_set = {round(x / dt) for x in (keep_psi_at or [])}
    snap_set |= keep_set  # a requested keep_psi_at time always gets a snapshot

    psi = psi0.astype(np.complex128).copy()
    snapshots: list[Snapshot] = []
    for n in range(n_t):
        for k in range(n_ch):
            c[n, k] = c_product(out_channels[k], psi)  # correlation: c-product
        norm[n] = float(np.linalg.norm(psi))  # Hermitian L2: physical, monotone
        if n in snap_set:
            rho_R, rho_r = _densities(tgrid, psi)
            snapshots.append(
                Snapshot(float(t[n]), rho_R, rho_r, psi.copy() if n in keep_set else None)
            )
        if n < n_steps:
            psi = step(psi)

    return PropagationResult(t=t, c=c, norm=norm, snapshots=snapshots)
