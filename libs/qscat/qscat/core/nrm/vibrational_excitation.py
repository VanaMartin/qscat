"""Vibrational excitation in the nonlocal resonance model.

PRA 77 Sec. III B splits the VE T-matrix into a background and a resonant
part (Eq. 28):

    T^VE_{vi->vf} = T^bg + T^res
    T^res = <chi_vf| V^{-*}_{dk_f} |Psi_d^+>                     Eq. (31)
    T^bg  = <chi_vf phi^-_{k_f}| V_int |chi_vi J_{k_i}>
            - <chi_vf| V^{-*}_{dk_f} J_{dk_i} |chi_vi>           Eq. (37)
    J_{dk}(R) = Int dr phi_d*(r;R) J_k(r)                        Eq. (38)
    sigma = 4 pi^3 |T^VE|^2 / k_i^2

`Psi_d^+` is the SAME Eq. (52) solution `dissociation.nrm_da_cross_section`
computes, so the resonant part is a contraction over already-gated machinery;
the background part is what is new here, and it is what PRA 77 shows a bare
LCP curve is missing (largest for the broadest resonance).

CONJUGATION -- three citations, three separate jobs; do not let one stand in
for another:

- **Eq. (34)** is the CONDITION selecting the implementable branch: for the
  radial case with a real discrete state, `phi_k^- = (phi_k^+)^*`.
- **Eq. (35)** is why that condition drops the conjugation on the COUPLING
  FACTOR itself, `V^{-*}_{dk} -> V^+_{dk}` (non-conjugated) -- it follows
  from `H_el` Hermiticity in the exact, unscaled theory, combined with
  Eq. (34). The paper's own emphasis: "we can use the matrix element
  `Vdk^+` but *without complex conjugation*."
- **p. 012710-6** is the separate reason the SCALAR PRODUCT itself (the
  c-product pairing of any two DVR coefficient vectors,
  `qscat.linalg.c_product`) carries no conjugation -- under exterior
  complex scaling `P H_el P` is complex symmetric, not Hermitian.

Eq. (38) is printed with `phi_d*(r;R)`, conjugating the discrete state. That
is a no-op here: Eq. (34)/(35) already restrict this model to a real
`phi_d`, so `phi_d* = phi_d` and `j_dk` below conjugates nothing.

WEIGHTS. DVR vectors are coefficients (`c_j = f(x_j) sqrt(w_j)`), so a
c-product of two coefficient vectors already IS the quadrature integral and a
function multiplying them is elementwise. No explicit weights appear here.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from qscat.dvr import FemDvrEcsGrid
from qscat.linalg import c_product

from .discrete_state import DiscreteState
from .scattering import incident_coefficients

__all__ = ["j_dk", "t_resonant"]


def j_dk(
    elec_grid: FemDvrEcsGrid,
    phi_d: DiscreteState,
    R_values: npt.NDArray[np.float64],
    energy: float,
    ell: int,
) -> npt.NDArray[np.complex128]:
    """`J_dk(R)` of Eq. (38): the overlap of `phi_d(.;R)` with the incident wave.

    Parameters
    ----------
    elec_grid : FemDvrEcsGrid
        The electronic radial grid.
    phi_d : DiscreteState
        The discrete-state choice; consumed only through `phi_d.phi_d(R)`.
    R_values : ndarray
        Nuclear coordinates at which to evaluate the overlap.
    energy : float
        The real electron energy `E = k^2/2` (hartree), positive.
    ell : int
        Partial wave.

    Returns
    -------
    ndarray
        Complex `J_dk(R)`, shape `(R_values.size,)`.
    """
    if energy <= 0.0:
        raise ValueError(f"energy must be positive, got {energy}")
    k = float(np.sqrt(2.0 * energy))
    inc = incident_coefficients(elec_grid, k, ell)
    R = np.asarray(R_values, dtype=np.float64)
    out = np.empty(R.size, dtype=np.complex128)
    for j in range(R.size):
        out[j] = c_product(phi_d.phi_d(float(R[j])), inc)
    return out


def t_resonant(
    chi_f: npt.NDArray[np.complex128],
    v_dk_f: npt.NDArray[np.complex128],
    psi_d: npt.NDArray[np.complex128],
) -> complex:
    """`T^res` of Eq. (31): `<chi_vf| V^{-*}_{dk_f} |Psi_d^+>`.

    All three arguments live on the full nuclear grid. `chi_f` and `psi_d`
    are DVR coefficient vectors; `v_dk_f` is the coupling as a function of
    `R` (`nonlocal_potential.continue_to_tail`'s output), so it multiplies
    elementwise -- no explicit weights, per the module docstring.

    The paper's `V^{-*}_{dk_f}` is not conjugated here: Eq. (34) is the
    CONDITION (radial case, real discrete state, so `phi_k^- = (phi_k^+)^*`);
    Eq. (35) is why that condition drops the conjugation off the COUPLING
    FACTOR, collapsing it to `V^+_{dk_f}` by `H_el` Hermiticity ("we can use
    the matrix element `Vdk^+` but *without complex conjugation*", p.
    012710-4); p. 012710-6 is the separate reason the SCALAR PRODUCT itself
    (`qscat.linalg.c_product`) carries no conjugation, since `P H_el P` is
    complex symmetric, not Hermitian, under exterior complex scaling. The
    paper warns this distinction "becomes important when the background
    terms ... are added to the resonant T matrix, since `V^±_dk` are in
    general complex even when the discrete state is real" (p. 012710-4) --
    Eq. (37)'s background term must use the same non-conjugated convention.

    Parameters
    ----------
    chi_f : ndarray
        Final vibrational state, DVR coefficients on the nuclear grid.
    v_dk_f : ndarray
        `V^+_{dk_f}(R)`, the coupling continued to the tail, on the same
        nuclear grid.
    psi_d : ndarray
        `Psi_d^+(R)` of Eq. (52) (the same solution the DA path computes),
        DVR coefficients on the same nuclear grid.

    Returns
    -------
    complex
        `T^res_{vi->vf}`.
    """
    if not (chi_f.size == v_dk_f.size == psi_d.size):
        raise ValueError(
            "chi_f, v_dk_f and psi_d must have the same length, got "
            f"{chi_f.size}, {v_dk_f.size}, {psi_d.size}"
        )
    return complex(c_product(chi_f, v_dk_f * psi_d))
