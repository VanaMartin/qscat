"""1-D TI resolvent solvers on the LCP curve: `lcp_da_cross_section`, `lcp_ve_cross_section`."""

from __future__ import annotations

from typing import Literal, overload

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp

from qscat.dvr import FemDvrEcsGrid, kinetic_sparse
from qscat.linalg import Ordering, SparseLU

# `return_wavefunction` output types (same convention as driven/dissociation):
# the 1-D nuclear resolvent `psi_sc(R)` per energy (`None` when the DA channel
# is closed), one array for scalar `E`, one list entry per energy for array `E`.
_Sigma = npt.NDArray[np.float64]
_Psi = npt.NDArray[np.complex128] | None
_PsiOut = _Psi | list[_Psi]


@overload
def lcp_da_cross_section(
    nuclear_grid: FemDvrEcsGrid,
    mu: float,
    Vd: npt.NDArray[np.complex128],
    Gamma: npt.NDArray[np.float64],
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    E: float | npt.ArrayLike,
    *,
    ordering: Ordering = ...,
    return_wavefunction: Literal[False] = ...,
) -> _Sigma: ...


@overload
def lcp_da_cross_section(
    nuclear_grid: FemDvrEcsGrid,
    mu: float,
    Vd: npt.NDArray[np.complex128],
    Gamma: npt.NDArray[np.float64],
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    E: float | npt.ArrayLike,
    *,
    ordering: Ordering = ...,
    return_wavefunction: Literal[True],
) -> tuple[_Sigma, _PsiOut]: ...


# bool catch-all (open()-style): callers holding a runtime flag forward it
# directly; the union return is narrowed by the Literal overloads above when
# the flag is literal.
@overload
def lcp_da_cross_section(
    nuclear_grid: FemDvrEcsGrid,
    mu: float,
    Vd: npt.NDArray[np.complex128],
    Gamma: npt.NDArray[np.float64],
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    E: float | npt.ArrayLike,
    *,
    ordering: Ordering = ...,
    return_wavefunction: bool = ...,
) -> _Sigma | tuple[_Sigma, _PsiOut]: ...


def lcp_da_cross_section(
    nuclear_grid: FemDvrEcsGrid,
    mu: float,
    Vd: npt.NDArray[np.complex128],
    Gamma: npt.NDArray[np.float64],
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    E: float | npt.ArrayLike,
    *,
    ordering: Ordering = "COLAMD",
    return_wavefunction: bool = False,
) -> _Sigma | tuple[_Sigma, _PsiOut]:
    """LCP dissociative-attachment sigma_DA(E) (bohr^2), TI resolvent form.

    *Provisional API* (docs/adr/0004-public-api-stability-policy.md): this wide
    functional signature is the layer the context-object refactor targets and
    may change in a minor release; `ScatteringProblem.lcp_da_cross_section` is the stable route.

    Solve `(E_tot I - H_res) psi_sc = d`, `H_res = T_nuc + diag(V_d - i Gamma/2)`,
    doorway `d = sqrt(Gamma/2pi) chi_{v_init}`; the DA amplitude is the outgoing
    dissociation flux at the boundary `X` (outermost real point):
    `S_DA = sqrt(K/2pi mu) psi_sc(X)`, `psi_sc(X) = psi_sc[b]/sqrt(w_b)` (the
    wavefunction VALUE, not the DVR coefficient), `sigma = 4 pi^3 |S_DA|^2/2E`.
    The DA threshold `eps_e = V_d(R_inf) = Vd[b].real` (open iff `E_tot > eps_e`).

    Requires the FINE per-molecule nuclear grid (the K~58 outgoing wave is
    unresolved on a coarse grid). The T->infty limit of eMoScat's TD
    `ModelLCP/SMatrix.cpp`. The approximation under test vs the exact-2D
    `da_cross_section` oracle -- validated at sigma_DA(F2,0.03)=1.47 vs ~1.66.

    Argument-order note (docs/adr/0007): this solver deliberately takes
    `(nuclear_grid, mu, Vd, Gamma, ...)` rather than a `model` -- the LCP
    equation contains no model; its physics input IS the curve, which may come
    from `resonance_levels(return_curve=True)`, a fit, or a file.
    `ScatteringProblem.lcp_da_cross_section` supplies `mu`/`eps`/`chi`/`v_init`
    from its bundle.

    If `return_wavefunction`, also returns the 1-D nuclear resolvent
    `psi_sc(R) = (E_tot I - H_res)^-1 d` per energy (`None` when the DA channel
    is closed -- `E <= 0` or `E_tot <= eps_e`): one array for scalar `E`, one
    list entry per energy for array `E`, same convention as
    `driven`/`dissociation`. `psi_sc` is the DVR-coefficient vector on the full
    nuclear grid (length `nuclear_grid.n`).
    """
    pts = nuclear_grid.points
    real_idx = np.flatnonzero(pts.imag == 0.0)
    b = int(real_idx[np.argmax(pts[real_idx].real)])
    eps_e = float(Vd[b].real)
    sqrt_wb = np.sqrt(complex(nuclear_grid.weights[b]))

    doorway = np.sqrt(Gamma / (2.0 * np.pi)).astype(np.complex128) * chi[v_init]
    H_res = (kinetic_sparse(nuclear_grid, mu) + sp.diags(Vd - 0.5j * Gamma)).tocsc()
    ident = sp.identity(nuclear_grid.n, format="csc", dtype=np.complex128)

    e_arr = np.atleast_1d(np.asarray(E, dtype=np.float64))
    out = np.zeros(e_arr.size, dtype=np.float64)
    psi_list: list[_Psi] = [None] * e_arr.size
    lu: SparseLU | None = None
    for ie, e in enumerate(e_arr):
        if float(e) <= 0.0:
            continue
        e_tot = float(e) + eps[v_init]
        e_dr = e_tot - eps_e
        if e_dr <= 0.0:
            continue
        a = (e_tot * ident - H_res).tocsc()
        if lu is None:
            lu = SparseLU(a, ordering=ordering)
        else:
            lu.refactor(a)
        psi_sc = lu.solve(doorway)
        psi_list[ie] = np.asarray(psi_sc, dtype=np.complex128)
        k_r = float(np.sqrt(2.0 * mu * e_dr))
        val = psi_sc[b] / sqrt_wb
        s_da = np.sqrt(k_r / (2.0 * np.pi * mu)) * val
        out[ie] = 4.0 * np.pi**3 * abs(s_da) ** 2 / (2.0 * float(e))

    scalar = np.isscalar(E) or (isinstance(E, np.ndarray) and np.ndim(E) == 0)
    sigma = np.asarray(out[0] if scalar else out, dtype=np.float64)
    if return_wavefunction:
        return sigma, (psi_list[0] if scalar else psi_list)
    return sigma


@overload
def lcp_ve_cross_section(
    nuclear_grid: FemDvrEcsGrid,
    mu: float,
    Vd: npt.NDArray[np.complex128],
    Gamma: npt.NDArray[np.float64],
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    vprimes: list[int],
    E: float | npt.ArrayLike,
    *,
    ordering: Ordering = ...,
    return_wavefunction: Literal[False] = ...,
) -> _Sigma: ...


@overload
def lcp_ve_cross_section(
    nuclear_grid: FemDvrEcsGrid,
    mu: float,
    Vd: npt.NDArray[np.complex128],
    Gamma: npt.NDArray[np.float64],
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    vprimes: list[int],
    E: float | npt.ArrayLike,
    *,
    ordering: Ordering = ...,
    return_wavefunction: Literal[True],
) -> tuple[_Sigma, _PsiOut]: ...


# bool catch-all (open()-style): callers holding a runtime flag forward it
# directly; the union return is narrowed by the Literal overloads above when
# the flag is literal.
@overload
def lcp_ve_cross_section(
    nuclear_grid: FemDvrEcsGrid,
    mu: float,
    Vd: npt.NDArray[np.complex128],
    Gamma: npt.NDArray[np.float64],
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    vprimes: list[int],
    E: float | npt.ArrayLike,
    *,
    ordering: Ordering = ...,
    return_wavefunction: bool = ...,
) -> _Sigma | tuple[_Sigma, _PsiOut]: ...


def lcp_ve_cross_section(
    nuclear_grid: FemDvrEcsGrid,
    mu: float,
    Vd: npt.NDArray[np.complex128],
    Gamma: npt.NDArray[np.float64],
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    vprimes: list[int],
    E: float | npt.ArrayLike,
    *,
    ordering: Ordering = "COLAMD",
    return_wavefunction: bool = False,
) -> _Sigma | tuple[_Sigma, _PsiOut]:
    """LCP vibrational-excitation sigma_{v_init->v'}(E) (bohr^2), TI resolvent form.

    *Provisional API* (docs/adr/0004-public-api-stability-policy.md): this wide
    functional signature is the layer the context-object refactor targets and
    may change in a minor release.

    Solve `(E_tot I - H_res) xi = d_{v_init}`, `H_res = T_nuc(mu) + diag(V_d
    - i Gamma/2)`, doorway `d_v = sqrt(Gamma/2pi) chi_v`; S-matrix element
    `S_{v'<-v_init} = <d_{v'}|xi>` by the DVR c-product (no conjugate);
    `sigma = 4 pi^3 |S|^2 / 2E`, exactly zero for `E <= 0` and for a closed
    final channel (`E_tot - eps[v'] <= 0`).

    Graduated from `projects/n2_ti_cross_section/cross_section.py`'s
    `ve_cross_section` (the deliberately dense 1-D toy model). This version
    is SPARSE and sweep-reusing: `A(E) = E_tot I - H_res` has an
    E-independent sparsity pattern, so the symbolic analysis is done once
    and `SparseLU.refactor` re-runs only the numeric factor per energy --
    the same structure as `lcp_da_cross_section` and `driven.ve_cross_section`.
    `xi` depends only on `(E, v_init)`, so one solve per energy serves every
    channel in `vprimes`.

    If `return_wavefunction`, also returns `xi(R)` per energy (`None` when
    `E <= 0`): one array for scalar `E`, one list entry per energy for array
    `E` -- the driven solution `nuclear_density.lcp_driven_solution` consumes.
    """
    doorway = np.sqrt(Gamma / (2.0 * np.pi)).astype(np.complex128)[None, :] * chi
    H_res = (kinetic_sparse(nuclear_grid, mu) + sp.diags(Vd - 0.5j * Gamma)).tocsc()
    ident = sp.identity(nuclear_grid.n, format="csc", dtype=np.complex128)

    e_arr = np.atleast_1d(np.asarray(E, dtype=np.float64))
    out = np.zeros((e_arr.size, len(vprimes)), dtype=np.float64)
    psi_list: list[_Psi] = [None] * e_arr.size
    lu: SparseLU | None = None
    for ie, e in enumerate(e_arr):
        if float(e) <= 0.0:
            continue
        e_tot = float(e) + eps[v_init]
        a = (e_tot * ident - H_res).tocsc()
        if lu is None:
            lu = SparseLU(a, ordering=ordering)
        else:
            lu.refactor(a)
        xi = lu.solve(doorway[v_init])
        psi_list[ie] = np.asarray(xi, dtype=np.complex128)
        for k, vp in enumerate(vprimes):
            if e_tot - eps[vp] <= 0.0:
                continue  # closed channel
            s_el = np.dot(doorway[vp], xi)  # c-product: no conjugate
            out[ie, k] = 4.0 * np.pi**3 * np.abs(s_el) ** 2 / (2.0 * float(e))

    scalar = np.isscalar(E) or (isinstance(E, np.ndarray) and np.ndim(E) == 0)
    sigma = np.asarray(out[0] if scalar else out, dtype=np.float64)
    if return_wavefunction:
        return sigma, (psi_list[0] if scalar else psi_list)
    return sigma
