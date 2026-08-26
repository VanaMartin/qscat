"""Exact 2-D VE cross section by the driven Lippmann-Schwinger equation.

    Psi_i    = F_{E,l}(r) chi_v(R)                     [masked to unscaled region]
    Psi_sc   = (E_tot - H_2D)^{-1} V_int Psi_i         [one sparse LU per energy]
    Psi^(+)  = Psi_i + Psi_sc
    T_{v->v'} = <chi_v' F_{E',l} | V_int | Psi^(+)>    [c-product, masked]
    sigma     = 4 pi^3 |T|^2 / k^2                     [bohr^2]

Conventions that must not be gotten wrong (each has bitten this repo or the
reference implementation):

- `chi` from `vibrational_states` is ALREADY a DVR coefficient vector; `F` is
  a FUNCTION and must be converted with `c_j = F(r_j) sqrt(w_j)` using the
  bridge-summed complex weight (`TensorGrid.sqrt_weights()`). Mixing the two
  rescales every cross section.
- Everything is paired with the C-PRODUCT (no conjugate): under ECS `H = H^T`,
  not `H^dagger`. With both sides in coefficient form the c-product IS the
  quadrature integral.
- `Psi_i` and every `Phi_f` are masked to the unscaled region. eMoScat uses a
  Hermitian dot here and is saved only by doing the same masking.

Elastic and inelastic share one formula: with `S = 1 - 2 pi i T`,
`|S - 1|^2 = 4 pi^2 |T|^2`, so `pi |S-1|^2 / k^2` and our `4 pi^3 |T|^2 / k^2`
are the same expression. The normalization is Vana & Houfek, Phys. Rev. A 95,
022714 (2017), Eqs. (38)-(39) -- NOT Houfek et al., Phys. Rev. A 73, 032721
(2006), whose Eqs. (25)-(26) use a differently normalized T. Unlike the 1-D LCP
model, this elastic T-matrix DOES contain the non-resonant background
scattering.

Promoted from `projects/n2_2d_cross_section/cross_section_2d.py`'s
`ve_cross_section_2d`. The only change is that the
Hamiltonian and interaction diagonal now come from a caller-supplied
`model: qscat.model.ResonanceModel` (`model.hamiltonian(tgrid)` /
`model.interaction_diag(tgrid)` / `model.ell`) instead of a hardcoded N2
`build_h2d`/`interaction_diag`/`ELL` import -- this is what makes the solver
model-agnostic. `qscat.core` never imports `qscat.model` at runtime (only
under `TYPE_CHECKING`, for the annotation); the solver depends on the
*protocol*, not the concrete class.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal, overload

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp

from qscat.dvr import TensorGrid
from qscat.linalg import Ordering, SparseLU, c_product

from .channels import channel_vector

if TYPE_CHECKING:
    from qscat.model import ResonanceModel

__all__ = ["ve_cross_section"]


def _sigma_at_one_energy(
    tgrid: TensorGrid,
    lu: SparseLU | None,
    v_diag: npt.NDArray[np.complex128],
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    vprimes: list[int],
    E: float,
    l: int,
    *,
    charge: int = 0,
    want_psi: bool,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.complex128] | None]:
    sigma = np.zeros(len(vprimes), dtype=np.float64)
    if E <= 0.0:
        # No driven-equation solve happens below threshold (there is no
        # scattering wavefunction to compute): `psi` is `None` here
        # regardless of `want_psi`, by design -- not a `want_psi`-dependent
        # omission.
        return sigma, None

    e_tot = E + eps[v_init]
    k = float(np.sqrt(2.0 * E))

    # `lu` is never `None` here: the caller builds (or refactors) the solver
    # for every `E > 0` before this call, and the `E <= 0` case returned
    # above. The assert makes that cross-function invariant explicit (and
    # narrows the `SparseLU | None` type for the `.solve` below).
    assert lu is not None
    psi_i = channel_vector(tgrid, k, chi[v_init], l, charge=charge)
    psi_plus = psi_i + lu.solve(v_diag * psi_i)
    v_psi = v_diag * psi_plus

    for j, vp in enumerate(vprimes):
        excess = e_tot - eps[vp]
        if excess <= 0.0:
            continue  # closed channel
        kp = float(np.sqrt(2.0 * excess))
        phi_f = channel_vector(tgrid, kp, chi[vp], l, charge=charge)
        t = c_product(phi_f, v_psi)
        sigma[j] = 4.0 * np.pi**3 * abs(t) ** 2 / (2.0 * E)

    return sigma, (psi_plus if want_psi else None)


# `return_wavefunction=False`'s return type (the common case).
_Sigma = npt.NDArray[np.float64]
# `return_wavefunction=True`'s wavefunction slot: `None` below threshold,
# else the `(psi_plus, ...)` for a scalar `E`, or one such entry per energy
# for an array `E`.
_Psi = npt.NDArray[np.complex128] | None
_PsiOut = _Psi | list[_Psi]


@overload
def ve_cross_section(
    tgrid: TensorGrid,
    model: ResonanceModel,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    vprimes: list[int],
    E: float | npt.ArrayLike,
    *,
    ordering: Ordering = ...,
    lam_scale: float = ...,
    return_wavefunction: Literal[False] = ...,
) -> _Sigma: ...


@overload
def ve_cross_section(
    tgrid: TensorGrid,
    model: ResonanceModel,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    vprimes: list[int],
    E: float | npt.ArrayLike,
    *,
    ordering: Ordering = ...,
    lam_scale: float = ...,
    return_wavefunction: Literal[True],
) -> tuple[_Sigma, _PsiOut]: ...


# bool catch-all (open()-style): callers holding a runtime flag forward it
# directly; the union return is narrowed by the Literal overloads above when
# the flag is literal.
@overload
def ve_cross_section(
    tgrid: TensorGrid,
    model: ResonanceModel,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    vprimes: list[int],
    E: float | npt.ArrayLike,
    *,
    ordering: Ordering = ...,
    lam_scale: float = ...,
    return_wavefunction: bool = ...,
) -> _Sigma | tuple[_Sigma, _PsiOut]: ...


def ve_cross_section(
    tgrid: TensorGrid,
    model: ResonanceModel,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    vprimes: list[int],
    E: float | npt.ArrayLike,
    *,
    ordering: Ordering = "COLAMD",
    lam_scale: float = 1.0,
    return_wavefunction: bool = False,
) -> _Sigma | tuple[_Sigma, _PsiOut]:
    """sigma_{v_init->v'}(E) in bohr^2, exact 2-D driven-equation solution.

    *Provisional API* (docs/adr/0004-public-api-stability-policy.md): this wide
    functional signature is the layer the context-object refactor targets and
    may change in a minor release; `ScatteringProblem.ve_cross_section` is the stable route.

    `model` supplies the Hamiltonian (`model.hamiltonian(tgrid)`), the
    interaction diagonal (`model.interaction_diag(tgrid)`) and the fixed
    partial wave (`model.ell`) -- the entire molecule-specific input.

    `model.charge` is forwarded to `channel_vector`, so an IONIC target's
    entrance/exit channels are the energy-normalized Coulomb functions
    (`coulomb_f_en`) rather than the free `riccati_bessel_en`; `charge == 0`
    (every neutral model) takes the identical pre-existing free-function
    branch. This is what lets `dissociation.dr_cross_section` reuse this
    solver's driven sweep for H2+ instead of re-inlining it.

    `E` may be scalar or an array; scalar returns shape `(len(vprimes),)`,
    array returns `(len(E), len(vprimes))`. One sparse LU per energy is
    reused across all `vprimes`.

    `lam_scale` scales `V_int` ONLY, for the free-particle and first-Born
    validation limits. It is a test lever, never a physics knob.

    If `return_wavefunction`, also returns `psi_plus` (or `None` when
    `E <= 0`, since no driven-equation solve happens below threshold): one
    array for scalar `E`, one list entry per energy for array `E`.
    """
    e_arr = np.atleast_1d(np.asarray(E, dtype=np.float64))
    H = model.hamiltonian(tgrid)
    v_diag = lam_scale * model.interaction_diag(tgrid)
    ident = sp.identity(tgrid.size, format="csc", dtype=np.complex128)

    # `A(E) = e_tot*ident - H` has an E-INDEPENDENT sparsity pattern (H fixed;
    # the identity only shifts the already-present diagonal), so the symbolic
    # analysis is done ONCE and reused: build `SparseLU` at the first energy
    # that needs a solve, then `refactor` it per subsequent energy. Energies
    # with `E <= 0` return zeros without any factorization (`_sigma_at_one_energy`
    # short-circuits before touching `lu`), so the solver is built lazily at the
    # first `E > 0` -- never for a below-threshold energy.
    out = []
    psis = []
    lu: SparseLU | None = None
    for e in e_arr:
        if float(e) > 0.0:
            e_tot = float(e) + eps[v_init]
            a = (e_tot * ident - H).tocsc()
            if lu is None:
                lu = SparseLU(a, ordering=ordering)
            else:
                lu.refactor(a)
        s, psi = _sigma_at_one_energy(
            tgrid,
            lu,
            v_diag,
            eps,
            chi,
            v_init,
            vprimes,
            float(e),
            model.ell,
            charge=model.charge,
            want_psi=return_wavefunction,
        )
        out.append(s)
        psis.append(psi)

    sigma = np.stack(out)
    scalar = np.isscalar(E) or (isinstance(E, np.ndarray) and np.ndim(E) == 0)
    result = np.asarray(sigma[0], dtype=np.float64) if scalar else sigma
    if return_wavefunction:
        return result, (psis[0] if scalar else psis)
    return result
