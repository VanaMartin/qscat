"""Asymptotic channel functions for the model-independent VE-scattering engine.

The entrance/exit channel of a VE transition is a free electron of momentum
`k` in partial wave `l`, times a neutral vibrational state. The electronic
factor is the ENERGY-NORMALIZED regular free radial solution
`riccati_bessel_en` (`qscat.special`).

Promoted from `projects/n2_2d_cross_section/cross_section_2d.py`'s
`channel_vector` (sub-project #A, Task 4). ONE change from the original: `l`
is now a REQUIRED parameter (the original defaulted it to N2's fixed
partial wave `ELL` -- `qscat.core` must never import anything N2-specific).
Callers now pass `model.ell` explicitly.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from qscat.dvr import TensorGrid
from qscat.exceptions import GridError
from qscat.linalg import c_product
from qscat.special import coulomb_f_en, riccati_bessel_en

__all__ = ["channel_vector"]

# `channel_vector` divides by `sqrt(c_product(chi, chi))`; guard against a
# (near-)null vibrational vector producing a divide-by-(near-)zero rather
# than a clear error. In practice `c_product(chi, chi)` is within ~7e-15 of
# 1.0 for every vibrational state this repo uses (see `vibrational_states`'s
# docstring), so this threshold is cheap insurance, not a normal code path.
_MIN_NORM2 = 1e-12


def channel_vector(
    tgrid: TensorGrid,
    k: float,
    chi_v: npt.NDArray[np.complex128],
    l: int,
    *,
    charge: int = 0,
) -> npt.NDArray[np.complex128]:
    """DVR coefficients of `F_{E,l}(r) chi_v(R)`, masked to the unscaled region.

    `chi_v` is already a coefficient vector; `F` is a function and picks up
    `sqrt(w_r)`. `charge=0` (default) is the neutral-target case: the
    electronic factor is the fast `riccati_bessel_en` free radial function,
    bit-for-bit unchanged from before `charge` existed (every VE/DA caller).
    `charge != 0` is the ionic-target case (e.g. H2+): the electronic factor
    is the energy-normalized Coulomb function `coulomb_f_en` instead (mass-1
    electron).
    """
    g_r = tgrid.grids[0]
    f_vals = (
        riccati_bessel_en(g_r.real_points, k, l)
        if charge == 0
        else coulomb_f_en(g_r.real_points, k, float(charge), 1.0, l)
    )
    # sqrt_weights() is per-axis and broadcast-shaped ((n_r, 1) at D=2); ravel
    # it to pair elementwise with the 1-D electronic function values.
    sqrt_w_r = tgrid.sqrt_weights()[0].ravel()
    f_coeff = f_vals * sqrt_w_r

    chi = np.asarray(chi_v, dtype=np.complex128)
    norm2 = c_product(chi, chi)
    if abs(norm2) < _MIN_NORM2:
        raise GridError(
            f"channel_vector: c-product norm^2 of chi_v is ~0 ({norm2!r}); "
            "cannot normalize a (near-)null vibrational vector"
        )
    chi = chi / np.sqrt(norm2)  # c-product normalization, not Hermitian

    psi = tgrid.outer([f_coeff, chi])
    psi[~tgrid.real_mask()] = 0.0
    return psi
