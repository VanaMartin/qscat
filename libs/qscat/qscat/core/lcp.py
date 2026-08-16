"""Model-independent local complex potential `(V_d(R), Gamma(R))`.

Sub-project B: the local-complex-potential (LCP) approximation of
dissociative attachment (and vibrational excitation) reduces the full 2-D
(electronic r x nuclear R) resonance problem to a 1-D nuclear problem by
replacing the fixed-R electronic resonance with a single complex number at
each R -- eMoScat's `ModelLCP` (`v0(R) + E_res(R)` real part, width
`-2*Im(E_res(R))`). This is the RESEARCH-PROGRAM "approximation under test":
the exact 2-D solver (`qscat.core.dissociation`/`driven`) is the oracle, and
`local_complex_potential` is the reduction whose accuracy against that oracle
is what sub-project B is built to measure -- not a description of the "real"
physics.

`V_d(R) = Re(E_pole(R))`, `Gamma(R) = max(0, -2*Im(E_pole(R)))`, where
`E_pole(R)` is the two-ECS-angle-matched resonance pole (`qscat.ecs.
find_resonance_pole`) of the fixed-R electronic Hamiltonian
`H_el(R) = -1/2 d^2/dr^2 + model.surface(r, R)`. Because `model.surface`
ALREADY INCLUDES `v0(R)` (the neutral-molecule channel potential), `V_d =
Re(E_pole)` directly -- adding `v0(R)` again would double-count it. (Contrast
`projects/n2_ti_cross_section/vres.py`, whose `v_eff_el` EXCLUDES `v0`, so
that code adds `v0(R) + E_res(R)` separately; the two bookkeeping schemes
must and do agree on the observable `V_d(R)`, checked in
`test_lcp.py::test_matches_n2_vres_oracle`.)

The continuation is seeded from the bound anion state at `R_inf =
nuclear_grid.R0` (`qscat.core.dissociation.anion_electronic_states`, the
same electronic bound state the exact-2D DA solver uses for its exit
channel) and walked INWARD (decreasing R), matching
`projects/n2_resonance/pole.py`'s `resonance_curve` /
`projects/n2_ti_cross_section/vres.py`'s `vres_on_grid` continuation, made
model-independent. The first rejected step (pole finder raises, or residual
>= `resid_tol`) freezes the electronic shift `s = V_d - v0(R)` at its last
accepted value for all remaining smaller-R points (`vres.py`'s documented
small-R breakdown -- physically irrelevant there since `v0(R)` is already
Ha above threshold and the vibrational wavefunction has negligible density).
The complex ECS tail (R.imag != 0) clamps `Gamma = 0` and freezes the shift
at its outermost-real-R (asymptotic) value.

See `docs/superpowers/specs/2026-07-27-da-cross-sections-design.md`.

`qscat.core` never imports `qscat.model`/`projects.*` at runtime: `model` is
typed against the `ResonanceModel` protocol under `TYPE_CHECKING` only,
exactly like `driven.py`/`dissociation.py`.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal, overload

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp

from qscat.dvr import FemDvrEcsGrid, eigen, kinetic, kinetic_sparse
from qscat.ecs import find_resonance_pole
from qscat.linalg import SparseLU, c_product

from .dissociation import anion_electronic_states

if TYPE_CHECKING:
    from qscat.model import ResonanceModel

__all__ = [
    "local_complex_potential",
    "lcp_da_cross_section",
    "resonance_pole_walk",
    "resonance_eigenstate",
    "resonance_eigenstate_at_peak_width",
]

# `return_wavefunction` output types (same convention as driven/dissociation):
# the 1-D nuclear resolvent `psi_sc(R)` per energy (`None` when the DA channel
# is closed), one array for scalar `E`, one list entry per energy for array `E`.
_Sigma = npt.NDArray[np.float64]
_Psi = npt.NDArray[np.complex128] | None
_PsiOut = _Psi | list[_Psi]

# `resonance_eigenstate_at_peak_width`: max |Re(E_pole_fresh) - V_d_walk| for a
# real R to count as GENUINE (walk-accepted) rather than a frozen-tail point. At
# an accepted R the fresh find and the walk are the same computation (agree to
# round-off); on the frozen plateau V_d is stale and the fresh pole drifts O(1e-2).
_FROZEN_TOL = 1e-3


def _h_el(model: ResonanceModel, R: complex, g: FemDvrEcsGrid) -> npt.NDArray[np.complex128]:
    return kinetic(g, 1.0) + np.diag(model.surface(g.points, R))


def resonance_pole_walk(
    model: ResonanceModel,
    R_descending: npt.NDArray[np.float64],
    elec_grid_a: FemDvrEcsGrid,
    elec_grid_b: FemDvrEcsGrid,
    seed_window: tuple[float, float, float, float],
    *,
    re_half_width: float = 0.05,
    im_half_width: float = 0.05,
    resid_tol: float = 1e-3,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Continuation walk of the resonance pole over descending real `R`.

    Returns `(shift, gamma)` aligned with `R_descending`, where
    `shift[j] = Re(E_pole(R_j)) - v0(R_j)` and `gamma[j] = max(0, -2
    Im(E_pole(R_j)))`. The pole finder (`qscat.ecs.find_resonance_pole`) is
    seeded from `seed_window` at the FIRST (outermost) `R` and recentered on
    each accepted pole; on breakdown (residual `>= resid_tol` or a solver
    error) the LAST accepted `(shift, gamma)` is FROZEN for all remaining
    (inner) `R`. The freeze holds the electronic SHIFT `s = V_d - v0(R)`
    constant, not the absolute pole. Raises `RuntimeError` if the finder
    fails already at the seed edge (no accepted pole to freeze).
    """
    window = seed_window
    shift = np.empty(R_descending.size, dtype=np.float64)  # s = V_d - v0(R)
    gamma_w = np.empty(R_descending.size, dtype=np.float64)
    last_s: float | None = None
    last_g = 0.0
    broken = False
    for j in range(R_descending.size):
        R = float(R_descending[j])
        if not broken:
            try:
                E_pole, resid = find_resonance_pole(
                    eigen(_h_el(model, R, elec_grid_a))[0],
                    eigen(_h_el(model, R, elec_grid_b))[0],
                    window,
                )
            except (ValueError, np.linalg.LinAlgError):
                resid = np.inf
            else:
                if resid < resid_tol:
                    v0R = float(np.real(model.v0(np.asarray(R))))
                    last_s = E_pole.real - v0R
                    last_g = max(0.0, -2.0 * E_pole.imag)
                    window = (
                        E_pole.real - re_half_width,
                        E_pole.real + re_half_width,
                        E_pole.imag - im_half_width,
                        E_pole.imag + im_half_width,
                    )
                    shift[j], gamma_w[j] = last_s, last_g
                    continue
            broken = True
        if last_s is None:
            raise RuntimeError("resonance_pole_walk: pole finder failed at the seed edge")
        shift[j], gamma_w[j] = last_s, last_g
    return shift, gamma_w


def resonance_eigenstate(
    model: ResonanceModel,
    elec_grid_a: FemDvrEcsGrid,
    elec_grid_b: FemDvrEcsGrid,
    R: float,
    window: tuple[float, float, float, float],
) -> tuple[complex, npt.NDArray[np.complex128]]:
    """The resonance pole energy + its electronic eigenfunction at nuclear `R`.

    Diagonalizes the fixed-`R` electronic Hamiltonian
    `H_el(R) = -1/2 d^2/dr^2 + model.surface(r, R)` at the two ECS angles,
    matches the angle-stable pole (`qscat.ecs.find_resonance_pole`, restricted to
    `window = (re_lo, re_hi, im_lo, im_hi)`), and returns `(E_pole, phi_res)`:
    `E_pole = V_d - i*Gamma/2` (complex; `Gamma = -2*Im`), and `phi_res` the
    angle-`a` eigenvector nearest `E_pole`, c-product-normalized over the
    electronic real region. The eigenstate counterpart of `local_complex_potential`
    (which keeps only the pole energy). Raises whatever `find_resonance_pole`
    raises if no pole lies in `window`.
    """
    E_a, V_a = eigen(_h_el(model, R, elec_grid_a))
    E_b, _ = eigen(_h_el(model, R, elec_grid_b))
    E_pole, _resid = find_resonance_pole(E_a, E_b, window)
    idx = int(np.argmin(np.abs(E_a - E_pole)))
    phi = V_a[:, idx].astype(np.complex128)
    real = elec_grid_a.real_points <= elec_grid_a.R0
    p = phi.copy()
    p[~real] = 0.0
    phi = phi / np.sqrt(c_product(p, p))
    return complex(E_pole), np.asarray(phi, dtype=np.complex128)


def resonance_eigenstate_at_peak_width(
    model: ResonanceModel,
    nuclear_grid: FemDvrEcsGrid,
    elec_grid_a: FemDvrEcsGrid,
    elec_grid_b: FemDvrEcsGrid,
    *,
    re_half_width: float = 0.05,
    im_half_width: float = 0.05,
) -> tuple[float, complex, npt.NDArray[np.complex128]]:
    """`(R_star, E_pole, phi_res)` at the nuclear geometry of maximum resonance width.

    Runs `local_complex_potential` to find `Gamma(R)`, then re-solves the resonance
    eigenstate (`resonance_eigenstate`) at the real-`R` points in DESCENDING order
    of `Gamma`, returning the first GENUINE (non-frozen) resolve. A point is
    genuine iff the fresh single-`R` pole reproduces `local_complex_potential`'s
    `V_d(R)` there (to `_FROZEN_TOL`): at a walk-accepted `R` the two are the same
    computation and agree to round-off, whereas on the frozen small-`R`
    continuation tail (where the walk holds the shift constant) `V_d` is stale and
    the fresh pole drifts away -- so this both excludes the unphysical frozen
    plateau from the width search AND guarantees the returned `E_pole` is
    consistent with the reported `V_d(R_star)`. Lands on the genuinely
    most-resonant geometry -- the natural single representative resonance state.

    Raises `RuntimeError` if no real-`R` point has a resolvable, genuine width
    (`Gamma` ~ 0 everywhere, or every wide point is frozen).
    """
    Vd, gamma = local_complex_potential(
        model,
        nuclear_grid,
        elec_grid_a,
        elec_grid_b,
        re_half_width=re_half_width,
        im_half_width=im_half_width,
    )
    pts = nuclear_grid.points
    real = np.flatnonzero(pts.imag == 0.0)
    order = real[np.argsort(gamma[real])[::-1]]  # widest first
    for j in order:
        if gamma[j] < 1e-4:
            break
        R = float(pts[j].real)
        e_re, g = float(Vd[j].real), float(gamma[j])
        window = (
            e_re - re_half_width,
            e_re + re_half_width,
            -0.5 * g - im_half_width,
            -0.5 * g + im_half_width,
        )
        try:
            E_pole, phi = resonance_eigenstate(model, elec_grid_a, elec_grid_b, R, window)
        except (ValueError, np.linalg.LinAlgError):
            continue  # unresolvable at this R -- try the next-widest
        if abs(E_pole.real - e_re) > _FROZEN_TOL:
            continue  # frozen point: fresh pole disagrees with the stale V_d -- skip
        return R, E_pole, phi
    raise RuntimeError(
        "resonance_eigenstate_at_peak_width: no real-R point has a resolvable, "
        "genuine (non-frozen) resonance width"
    )


def _assemble_lcp(
    model: ResonanceModel,
    grid: FemDvrEcsGrid,
    shift: npt.NDArray[np.float64],
    gamma_w: npt.NDArray[np.float64],
) -> tuple[npt.NDArray[np.complex128], npt.NDArray[np.float64]]:
    """Place one `resonance_pole_walk` result onto `grid` as `(V_d, Gamma)`.

    `shift`/`gamma_w` are indexed by DESCENDING real `R` (the walk order).
    Real nodes get `V_d = v0(R) + shift`; the complex ECS tail gets the
    analytic continuation `v0(z) + shift[0]` (the shift at the largest real `R`,
    i.e. the asymptotic electronic shift) with `Gamma = 0`.

    Factored out of `local_complex_potential` so `resonance_levels` can run the
    expensive electronic walk ONCE and lay the same curve onto two nuclear grids
    that differ only in their ECS tail angle.
    """
    pts = grid.points
    real_idx = np.flatnonzero(pts.imag == 0.0)
    order = np.argsort(pts[real_idx].real)[::-1]  # descending R: outer -> inner
    walk = real_idx[order]

    Vd = np.empty(grid.n, dtype=np.complex128)
    Gamma = np.zeros(grid.n, dtype=np.float64)
    Vd[walk] = model.v0(pts[walk].real) + shift
    Gamma[walk] = gamma_w

    tail = np.flatnonzero(pts.imag != 0.0)
    if tail.size:
        Vd[tail] = model.v0(pts[tail]) + shift[0]
    return Vd, Gamma


def local_complex_potential(
    model: ResonanceModel,
    nuclear_grid: FemDvrEcsGrid,
    elec_grid_a: FemDvrEcsGrid,
    elec_grid_b: FemDvrEcsGrid,
    *,
    re_half_width: float = 0.05,
    im_half_width: float = 0.05,
    resid_tol: float = 1e-3,
) -> tuple[npt.NDArray[np.complex128], npt.NDArray[np.float64]]:
    """Local complex potential `(V_d(R), Gamma(R))` for the LCP DA/VE model.

    `V_d(R) = Re(E_pole(R))`, `Gamma(R) = max(0, -2 Im(E_pole(R)))`, `E_pole`
    the two-angle-matched resonance pole of `-1/2 d^2/dr^2 + model.surface(r,R)`
    (surface includes v0(R), so V_d = Re(E_pole) directly). Seeded from the
    bound anion at `R_inf = nuclear_grid.R0` (`anion_electronic_states`) and
    continued INWARD; small-R breakdown freezes the electronic shift; the
    complex tail clamps Gamma=0. See module docstring.
    """
    R_inf = nuclear_grid.R0
    eps_e, _ = anion_electronic_states(elec_grid_a, model, R_inf, 1)
    seed_window = (
        eps_e[0] - re_half_width,
        eps_e[0] + re_half_width,
        -im_half_width,
        im_half_width,
    )

    pts = nuclear_grid.points
    real_idx = np.flatnonzero(pts.imag == 0.0)
    order = np.argsort(pts[real_idx].real)[::-1]  # descending R: outer -> inner
    walk = real_idx[order]
    R_real = pts[walk].real

    shift, gamma_w = resonance_pole_walk(
        model,
        R_real,
        elec_grid_a,
        elec_grid_b,
        seed_window,
        re_half_width=re_half_width,
        im_half_width=im_half_width,
        resid_tol=resid_tol,
    )

    return _assemble_lcp(model, nuclear_grid, shift, gamma_w)


# Mirrors `driven.py`'s (private) `_Ordering` -- scipy's `splu`'s
# `permc_spec`. Not imported directly: that name is an internal detail of
# `SparseLU`, not part of its public API.
_Ordering = Literal["NATURAL", "MMD_ATA", "MMD_AT_PLUS_A", "COLAMD"]


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
    ordering: _Ordering = ...,
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
    ordering: _Ordering = ...,
    return_wavefunction: Literal[True],
) -> tuple[_Sigma, _PsiOut]: ...


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
    ordering: _Ordering = "COLAMD",
    return_wavefunction: bool = False,
) -> _Sigma | tuple[_Sigma, _PsiOut]:
    """LCP dissociative-attachment sigma_DA(E) (bohr^2), TI resolvent form.

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
