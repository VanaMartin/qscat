"""Thin N2-binding shim over the model-agnostic, promoted implementation in
`qscat.core.time_dependent` (sub-project #A, Task 5) -- `_propagate`,
`_s_vector_one_energy`, `_sigma_one_energy`, `sigma_from_correlations`, and
`td_ve_cross_section_2d` here just fix `model=N2` and delegate (the
`qscat.core.time_dependent` names they wrap are `propagate_wavepacket`,
`s_vector_one_energy`, `sigma_one_energy` -- promoted to public, `eps`
dropped from `propagate_wavepacket`, lib-m7/lib-M16). See
`qscat.core.time_dependent`'s module docstring for the full physics/
convention writeup (the Tannor-Weeks energy transform, the elastic
free-reference fix) -- unchanged by this promotion.

Kept as a module (not deleted) so existing callers/imports in this project
(and its tests) are unaffected by the move; no new physics or numerics live
here.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from qscat.core.time_dependent import PropagationResult
from qscat.core.time_dependent import propagate_wavepacket as _core_propagate
from qscat.core.time_dependent import s_vector_one_energy as _core_s_vector_one_energy
from qscat.core.time_dependent import sigma_from_correlations as _core_sigma_from_correlations
from qscat.core.time_dependent import sigma_one_energy as _core_sigma_one_energy
from qscat.core.time_dependent import td_ve_cross_section as _core_td_ve_cross_section
from qscat.dvr import TensorGrid
from qscat.model import N2

__all__ = ["sigma_from_correlations", "td_ve_cross_section_2d"]

# Wavepacket parameter dict keys `initial_state`/`outgoing_channel` accept
# (r0/p0/sigma for the incident packet; r0_out/p0_out/sigma_out for the
# outgoing test function).
_WpIn = dict[str, float]
_WpOut = dict[str, float]


def _propagate(
    tgrid: TensorGrid,
    chi: npt.NDArray[np.complex128],
    v_init: int,
    vprimes: list[int],
    *,
    dt: float,
    n_steps: int,
    wp_in: _WpIn,
    wp_out: _WpOut,
    free: bool = False,
    order: int = 3,
) -> PropagationResult:
    """N2-binding shim over `qscat.core.time_dependent.propagate_wavepacket`;
    see there for the full docstring."""
    return _core_propagate(
        tgrid,
        N2,
        chi,
        v_init,
        vprimes,
        dt=dt,
        n_steps=n_steps,
        wp_in=wp_in,
        wp_out=wp_out,
        free=free,
        order=order,
    )


def _s_vector_one_energy(
    tgrid: TensorGrid,
    result: PropagationResult,
    eps: npt.NDArray[np.float64],
    v_init: int,
    vprimes: list[int],
    E: float,
    dt: float,
    wp_in: _WpIn,
    wp_out: _WpOut,
) -> npt.NDArray[np.complex128]:
    """N2-binding shim over `qscat.core.time_dependent.s_vector_one_energy`;
    see there for the full docstring."""
    return _core_s_vector_one_energy(tgrid, N2, result, eps, v_init, vprimes, E, dt, wp_in, wp_out)


def _sigma_one_energy(
    tgrid: TensorGrid,
    result: PropagationResult,
    eps: npt.NDArray[np.float64],
    v_init: int,
    vprimes: list[int],
    E: float,
    dt: float,
    wp_in: _WpIn,
    wp_out: _WpOut,
    free_result: PropagationResult | None = None,
) -> npt.NDArray[np.float64]:
    """N2-binding shim over `qscat.core.time_dependent.sigma_one_energy`;
    see there for the full docstring."""
    return _core_sigma_one_energy(
        tgrid, N2, result, eps, v_init, vprimes, E, dt, wp_in, wp_out, free_result
    )


def sigma_from_correlations(
    tgrid: TensorGrid,
    result: PropagationResult,
    eps: npt.NDArray[np.float64],
    v_init: int,
    vprimes: list[int],
    E: float | npt.ArrayLike,
    *,
    dt: float,
    wp_in: _WpIn,
    wp_out: _WpOut,
    free_result: PropagationResult | None = None,
) -> npt.NDArray[np.float64]:
    """N2-binding shim over `qscat.core.time_dependent.sigma_from_correlations`;
    see there for the full docstring."""
    return _core_sigma_from_correlations(
        tgrid,
        N2,
        result,
        eps,
        v_init,
        vprimes,
        E,
        dt=dt,
        wp_in=wp_in,
        wp_out=wp_out,
        free_result=free_result,
    )


def td_ve_cross_section_2d(
    tgrid: TensorGrid,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    vprimes: list[int],
    E: float | npt.ArrayLike,
    *,
    dt: float,
    n_steps: int,
    wp_in: _WpIn,
    wp_out: _WpOut,
    subtract_free_reference: bool = True,
    order: int = 3,
) -> npt.NDArray[np.float64]:
    """N2-binding shim over `qscat.core.time_dependent.td_ve_cross_section`;
    see there for the full docstring (scalar/array `E` contract,
    `subtract_free_reference` semantics)."""
    return _core_td_ve_cross_section(
        tgrid,
        N2,
        eps,
        chi,
        v_init,
        vprimes,
        E,
        dt=dt,
        n_steps=n_steps,
        wp_in=wp_in,
        wp_out=wp_out,
        order=order,
        subtract_free_reference=subtract_free_reference,
    )
