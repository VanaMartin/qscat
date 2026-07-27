"""Exact 2-D VE cross section by the driven Lippmann-Schwinger equation.

Thin N2-binding shim over the model-agnostic, promoted implementation in
`qscat.core.channels`/`qscat.core.driven` (sub-project #A, Task 4) --
`channel_vector` and `ve_cross_section_2d` here just fix `l=ELL`/`model=N2`
and delegate. See `qscat.core.driven`'s module docstring for the full
physics/convention writeup (driven-equation formula, c-product masking,
elastic/inelastic T-matrix unification) -- unchanged by this promotion.

Kept as a module (not deleted) so existing callers/imports in this project
(and its tests) are unaffected by the move; no new physics or numerics live
here.
"""

from __future__ import annotations

from typing import Literal, overload

import numpy as np
import numpy.typing as npt
from qscat.core.channels import channel_vector as _channel_vector
from qscat.core.driven import ve_cross_section as _ve_cross_section
from qscat.dvr import TensorGrid
from qscat.model import N2

from projects.n2_2d_cross_section.hamiltonian2d import ELL

__all__ = ["channel_vector", "ve_cross_section_2d"]


def channel_vector(
    tgrid: TensorGrid,
    k: float,
    chi_v: npt.NDArray[np.complex128],
    *,
    l: int = ELL,
) -> npt.NDArray[np.complex128]:
    """DVR coefficients of `F_{E,l}(r) chi_v(R)`, masked to the unscaled region.

    N2-binding shim over `qscat.core.channels.channel_vector`; see there for
    the full docstring.
    """
    return _channel_vector(tgrid, k, chi_v, l)


# `return_wavefunction=False`'s return type (the common case).
_Sigma = npt.NDArray[np.float64]
_Psi = npt.NDArray[np.complex128] | None
_PsiOut = _Psi | list[_Psi]


@overload
def ve_cross_section_2d(
    tgrid: TensorGrid,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    vprimes: list[int],
    E: float | npt.ArrayLike,
    *,
    ordering: str = ...,
    lam_scale: float = ...,
    return_wavefunction: Literal[False] = ...,
) -> _Sigma: ...


@overload
def ve_cross_section_2d(
    tgrid: TensorGrid,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    vprimes: list[int],
    E: float | npt.ArrayLike,
    *,
    ordering: str = ...,
    lam_scale: float = ...,
    return_wavefunction: Literal[True],
) -> tuple[_Sigma, _PsiOut]: ...


def ve_cross_section_2d(
    tgrid: TensorGrid,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    vprimes: list[int],
    E: float | npt.ArrayLike,
    *,
    ordering: str = "COLAMD",
    lam_scale: float = 1.0,
    return_wavefunction: bool = False,
) -> _Sigma | tuple[_Sigma, _PsiOut]:
    """sigma_{v_init->v'}(E) in bohr^2, exact 2-D driven-equation solution.

    N2-binding shim over `qscat.core.driven.ve_cross_section`; see there for
    the full docstring (scalar/array `E` contract, `lam_scale` test lever,
    `return_wavefunction` semantics).
    """
    if return_wavefunction:
        return _ve_cross_section(
            tgrid,
            N2,
            eps,
            chi,
            v_init,
            vprimes,
            E,
            ordering=ordering,
            lam_scale=lam_scale,
            return_wavefunction=True,
        )
    return _ve_cross_section(
        tgrid,
        N2,
        eps,
        chi,
        v_init,
        vprimes,
        E,
        ordering=ordering,
        lam_scale=lam_scale,
        return_wavefunction=False,
    )
