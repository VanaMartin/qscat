"""The coupled-channel driven solve: sigma_VE(E) with the partial waves mixed.

Mirrors `qscat.core.driven.ve_cross_section`. Three things differ, and only
three:

- the interaction is a MATRIX (`CoupledModel.interaction_matrix`), so the
  Lippmann-Schwinger step is `psi_i + lu.solve(V @ psi_i)` rather than an
  elementwise product;
- the entrance is a single-channel vector embedded in one block;
- the exit is summed over blocks, because the coupling lets the electron leave
  in a partial wave it did not enter on.

The post-form T-matrix, the non-conjugated c-product and the
`4 pi^3 |T|^2 / 2E` normalisation are unchanged.

This duplicates roughly forty lines of sweep boilerplate from `driven.py`. That
is deliberate at toy stage -- generalising a shipped solver before the coupled
shape has been used twice is what the lifecycle exists to prevent -- and
`test_scattering.py` gates the duplicate against the original at s = 0.

See docs/physics/coupled-partial-waves.md.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp
from qscat.core.channels import channel_vector
from qscat.dvr import TensorGrid
from qscat.linalg import Ordering, SparseLU, c_product

from projects.no_coupled_channels.model import CoupledModel

__all__ = ["CoupledSigma", "channel_block", "coupled_channel_vector", "coupled_ve_cross_section"]


def channel_block(
    vec: npt.NDArray[np.complex128], channel: int, n: int
) -> npt.NDArray[np.complex128]:
    """The `channel`-th length-`n` slice of a channel-outermost vector."""
    return np.asarray(vec[channel * n : (channel + 1) * n], dtype=np.complex128)


def coupled_channel_vector(
    tgrid: TensorGrid,
    k: float,
    chi_v: npt.NDArray[np.complex128],
    ells: tuple[int, ...],
    channel: int,
    *,
    charge: int = 0,
) -> npt.NDArray[np.complex128]:
    """`F_{E,l}(r) chi_v(R)` for one partial wave, embedded in its block.

    `k` is shared across channels: it depends on the vibrational level alone,
    since the partial wave changes the Bessel ORDER in `F_{E,l}` rather than
    the momentum. `ells[channel]` is what selects that order.
    """
    n = tgrid.size
    out = np.zeros(len(ells) * n, dtype=np.complex128)
    out[channel * n : (channel + 1) * n] = channel_vector(
        tgrid, k, chi_v, ells[channel], charge=charge
    )
    return out


@dataclass(frozen=True)
class CoupledSigma:
    """`sigma_{v_init->v'}(E)` in bohr^2, two ways.

    `total` sums over EXIT partial waves -- what an angle-integrated
    measurement sees, and the like-for-like partner of a fixed-l model's single
    exit. `restricted` keeps only the exit channel equal to the entrance, which
    isolates how the coupling changes the entrance amplitude through virtual
    excursions into other waves. Their difference is the flux the coupling
    redistributes, and it costs nothing extra: both come from one solve.
    """

    E: npt.NDArray[np.float64]
    total: npt.NDArray[np.float64]
    restricted: npt.NDArray[np.float64]


def coupled_ve_cross_section(
    tgrid: TensorGrid,
    model: CoupledModel,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    vprimes: list[int],
    E: npt.ArrayLike,
    *,
    entrance: int = 0,
    ordering: Ordering = "COLAMD",
) -> CoupledSigma:
    """Exact 2-D coupled-channel `sigma_{v_init->v'}(E)`.

    The analysis is done ONCE and reused: `SparseLU` is built at the first
    energy that needs a solve and `refactor`ed per subsequent energy, since
    `E_tot * I - H` keeps one sparsity pattern across the sweep. On the
    SuperLU backend that reuse does not happen (it re-runs `splu`), which is
    why the production sweep needs MUMPS.

    Energies at or below threshold return zeros without any factorisation.
    """
    E_arr = np.atleast_1d(np.asarray(E, dtype=np.float64))
    ells = model.channel_ells()
    n = tgrid.size
    charge = model.charge

    H = model.hamiltonian(tgrid)
    V = model.interaction_matrix(tgrid)
    ident = sp.identity(H.shape[0], dtype=np.complex128, format="csr")

    total = np.zeros((E_arr.size, len(vprimes)), dtype=np.float64)
    restricted = np.zeros((E_arr.size, len(vprimes)), dtype=np.float64)
    lu: SparseLU | None = None

    for i, e in enumerate(E_arr):
        if e <= 0.0:
            continue
        e_tot = float(e) + eps[v_init]
        a = sp.csc_matrix(e_tot * ident - H)
        if lu is None:
            lu = SparseLU(a, ordering=ordering)
        else:
            lu.refactor(a)

        k = float(np.sqrt(2.0 * e))
        psi_i = coupled_channel_vector(tgrid, k, chi[v_init], ells, entrance, charge=charge)
        psi_plus = psi_i + lu.solve(V @ psi_i)
        v_psi = V @ psi_plus

        for j, vp in enumerate(vprimes):
            excess = e_tot - eps[vp]
            if excess <= 0.0:
                continue  # closed channel
            kp = float(np.sqrt(2.0 * excess))
            for c in range(len(ells)):
                # Only block `c` of the projection is non-zero, so take that
                # block rather than building a full-length vector of zeros.
                phi = channel_vector(tgrid, kp, chi[vp], ells[c], charge=charge)
                t = c_product(phi, channel_block(v_psi, c, n))
                s = 4.0 * np.pi**3 * abs(t) ** 2 / (2.0 * float(e))
                total[i, j] += s
                if c == entrance:
                    restricted[i, j] = s

    return CoupledSigma(E=E_arr, total=total, restricted=restricted)
