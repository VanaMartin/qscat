"""Model-independent local complex potential `(V_d(R), Gamma(R))`.

The local-complex-potential (LCP) approximation of
dissociative attachment (and vibrational excitation) reduces the full 2-D
(electronic r x nuclear R) resonance problem to a 1-D nuclear problem by
replacing the fixed-R electronic resonance with a single complex number at
each R -- eMoScat's `ModelLCP` (`v0(R) + E_res(R)` real part, width
`-2*Im(E_res(R))`). This is the RESEARCH-PROGRAM "approximation under test":
the exact 2-D solver (`qscat.core.dissociation`/`driven`) is the oracle, and
`local_complex_potential` is the reduction whose accuracy against that oracle
is the thing being measured -- not a description of the "real" physics.

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

**Naming caution.** This module's `V_d` is the resonance pole's real part --
the quantity Vana & Houfek, Phys. Rev. A 95, 022714 (2017) Eq. (41) call
`E_res(R)`. It is NOT the `V_d(R)` of Houfek, Rescigno & McCurdy, Phys. Rev. A
77, 012710 (2008) Eq. (20), which is the DISCRETE-STATE potential
`<phi_d|H_el|phi_d>` of the nonlocal theory. The two only *almost* coincide,
and only for that paper's "physical" choice of discrete state; this code never
constructs a `phi_d` at all. See
`reference/literature/houfek-2008-pra77-012710.md` for the terminology map.

`qscat.core` never imports `qscat.model`/`projects.*` at runtime: `model` is
typed against the `ResonanceModel` protocol under `TYPE_CHECKING` only,
exactly like `driven.py`/`dissociation.py`.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal, overload

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp

from qscat.dvr import FemDvrEcsGrid, eigen, kinetic, kinetic_sparse
from qscat.ecs import find_resonance_pole, match_angle_stable
from qscat.exceptions import ConvergenceError
from qscat.linalg import SparseLU, c_product

from .dissociation import anion_electronic_states
from .grids import assert_shared_real_nodes

if TYPE_CHECKING:
    from qscat.model import ResonanceModel

__all__ = [
    "ResonanceLevels",
    "lcp_da_cross_section",
    "lcp_resonance_levels",
    "lcp_ve_cross_section",
    "local_complex_potential",
    "resonance_eigenstate",
    "resonance_eigenstate_at_peak_width",
    "resonance_levels",
    "resonance_pole_walk",
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

# `_assemble_lcp`: largest `Gamma` at the outermost REAL node that still counts as
# "zero" at the real/ECS-tail junction, where the tail's clamped `Gamma = 0` meets
# the walk's value. Above it, `W = V_d - i Gamma/2` STEPS at the junction and the
# rotated tail reflects instead of absorbing.
_JUNCTION_GAMMA_TOL = 1e-8

# `_levels_from`: smallest |sum_i c_i^2| a state may have and still be normalized
# by the bilinear c-product. `eigen` returns Euclidean-unit vectors, so this
# quantity is <= 1 and equals 1 for a real vector; a value near 0 means the state
# is (numerically) self-orthogonal and `c / sqrt(sum c^2)` blows it up into noise.
_C_NORM_TOL = 1e-8

# `resonance_eigenstate_at_peak_width`: smallest `Gamma(R)` still worth
# re-solving a pole at. The search walks real `R` widest-first, so once it
# reaches a width this small every remaining point is narrower still, and a
# pole that narrow is not separable from the discretized continuum by the
# two-angle match -- the loop stops rather than grinding through them.
_MIN_RESOLVABLE_GAMMA = 1e-4

# `resonance_levels`: largest `Gamma(R)` tolerated, without warning, in the
# region where the anion curve lies BELOW the neutral one and autodetachment
# is energetically closed (Vana 2017, Sec. 1.5) -- there Gamma must vanish, so
# anything above this is pole-finder noise leaking into a closed region.
_CLOSED_REGION_GAMMA_TOL = 1e-6


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
    constant, not the absolute pole. Raises `ConvergenceError` if the finder
    fails already at the seed edge (no accepted pole to freeze).

    Warns
    -----
    UserWarning
        When the walk freezes, with the breakdown `R` and how many nodes are
        held constant. The freeze is DELIBERATE and usually harmless -- at
        small `R` the neutral curve is already hartrees above threshold and
        the vibrational wavefunction has negligible density there -- but it is
        an EXTRAPOLATION, and how far in it starts depends on the electronic
        grid. Measured on F2's production nuclear deck (2026-08-24): a
        55-point electronic grid breaks down at `R = 2.5033` and holds
        `Gamma = 0.00949256` over the inner 198 of 819 real nodes, while a
        132-point one runs on to `R = 1.8657` and gives `Gamma` = 0.0104 /
        0.140 / 0.539 at R = 2.49 / 2.20 / 1.51 -- 57x the frozen value at the
        innermost of those. Anything that reads
        `Gamma(R)` inside the frozen region (a doorway, a packet's turning
        point) is reading a constant, not the molecule, and silence about
        that is how it goes unnoticed.
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
                raise ConvergenceError("resonance_pole_walk: pole finder failed at the seed edge")
            warnings.warn(
                f"resonance_pole_walk: the pole finder broke down at R = "
                f"{R:.4f}; the electronic shift and Gamma are FROZEN at their "
                f"last accepted values (shift = {last_s:.6g}, Gamma = "
                f"{last_g:.6g}) for the remaining "
                f"{R_descending.size - j} of {R_descending.size} nodes. That "
                "region is an extrapolation, not a computed curve, and how far "
                "in it starts depends on the electronic grid -- refine it if "
                "anything you care about (a doorway, a turning point) lives "
                "inside the frozen range.",
                stacklevel=2,
            )
        # Repeated deliberately: the first check above already raised, so this
        # branch is unreachable at runtime -- it is here so mypy can narrow
        # `last_s` for the return below. Removing it as dead code breaks the
        # type check.
        if last_s is None:
            raise ConvergenceError("resonance_pole_walk: pole finder failed at the seed edge")
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

    Raises `ConvergenceError` if no real-`R` point has a resolvable, genuine width
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
        if gamma[j] < _MIN_RESOLVABLE_GAMMA:
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
    raise ConvergenceError(
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

    Warns if `Gamma` is not ~0 at the outermost real node (`gamma_w[0]`, the
    real/tail junction): the tail is force-zeroed here, so a nonzero value there
    makes `W = V_d - i Gamma/2` STEP discontinuously at the junction, and the ECS
    tail reflects the outgoing dissociative wave instead of absorbing it.
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
        if gamma_w.size and gamma_w[0] > _JUNCTION_GAMMA_TOL:
            warnings.warn(
                f"_assemble_lcp: Gamma = {float(gamma_w[0]):.3e} at the outermost real "
                f"node R = {float(pts[walk[0]].real):.4f}, but the ECS tail clamps "
                "Gamma = 0 -- the local complex potential STEPS at the real/tail "
                "junction and the tail will reflect rather than absorb. Extend the "
                "real region outward until the autodetachment width has died off.",
                UserWarning,
                stacklevel=3,
            )
    return Vd, Gamma


def _walk_from_anion_seed(
    model: ResonanceModel,
    nuclear_grid: FemDvrEcsGrid,
    elec_grid_a: FemDvrEcsGrid,
    elec_grid_b: FemDvrEcsGrid,
    *,
    re_half_width: float,
    im_half_width: float,
    resid_tol: float,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """`resonance_pole_walk` over `nuclear_grid`'s real nodes, anion-seeded.

    The seed window is centered on the bound anion electronic state at
    `R_inf = nuclear_grid.R0` (`anion_electronic_states`) and the walk runs
    INWARD over the real nodes in descending `R`. Returns `(shift, gamma_w)`
    in that walk order, ready for `_assemble_lcp`.

    The one implementation shared by `local_complex_potential` and
    `resonance_levels`; they differ only in how many nuclear grids they lay
    the resulting curve onto.
    """
    eps_e, _ = anion_electronic_states(elec_grid_a, model, nuclear_grid.R0, 1)
    seed_window = (
        eps_e[0] - re_half_width,
        eps_e[0] + re_half_width,
        -im_half_width,
        im_half_width,
    )
    pts = nuclear_grid.points
    real_idx = np.flatnonzero(pts.imag == 0.0)
    walk = real_idx[np.argsort(pts[real_idx].real)[::-1]]  # descending R: outer -> inner
    return resonance_pole_walk(
        model,
        pts[walk].real,
        elec_grid_a,
        elec_grid_b,
        seed_window,
        re_half_width=re_half_width,
        im_half_width=im_half_width,
        resid_tol=resid_tol,
    )


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
    shift, gamma_w = _walk_from_anion_seed(
        model,
        nuclear_grid,
        elec_grid_a,
        elec_grid_b,
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
    ordering: _Ordering = ...,
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
    ordering: _Ordering = "COLAMD",
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
    ordering: _Ordering = ...,
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
    ordering: _Ordering = ...,
    return_wavefunction: Literal[True],
) -> tuple[_Sigma, _PsiOut]: ...


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
    ordering: _Ordering = "COLAMD",
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


@dataclass(frozen=True)
class ResonanceLevels:
    """The quasi-bound vibrational levels of the anion in the LCP curve.

    The thesis's `omega_j` (Vana 2017, Sec. 1.5/3.4), promoted from real levels
    of `Re V_res` to genuine complex eigenvalues. These are complex-scaled (ECS)
    resonance eigenstates -- NOT Siegert pseudostates, which carry an
    outgoing-wave condition at a finite radius and a surface-corrected
    orthogonality relation (Hvizdos et al., Phys. Rev. A 97, 022704 (2018),
    App. A). ECS rotates rather than truncates, so the plain bilinear c-product
    is the complete inner product here.

    - `energies`: `E_v - i Gamma_v/2` (Hartree), ascending in `Re E`.
    - `widths`: `Gamma_v = max(0, -2 Im E_v)` (Hartree). A level below the anion
      dissociation limit carries only the ELECTRONIC autodetachment width; one
      above it also carries a NUCLEAR (dissociative) width. Both come out of the
      one diagonalization.
    - `states`: shape `(n_levels, grid.n)` DVR COEFFICIENTS `c_i`
      (`psi(R_i) = c_i / sqrt(w_i)`), c-product-normalized: `sum_i c_i^2 = 1`.
    - `residuals`: the two-angle ECS-TAIL stability residual per level,
      `|E_a - E_b|` between the matched eigenvalues of the two rotation
      angles (`qscat.ecs.match_angle_stable`). A large residual means the
      level is contaminated by the rotated continuum, not a genuine pole.
      It is NOT a real-region convergence diagnostic: `nuclear_grid_a` and
      `nuclear_grid_b` are required to share every real node and
      quadrature, so real-region discretization error is common to both
      spectra and cancels out of the difference -- `residuals` stays near
      machine precision even on a badly under-resolved real grid. Judge
      real-region convergence separately, by refining the shared real
      nodes and checking that `energies` itself does not move.
    - `real_weight`: fraction of `|c|^2` inside the real region -- a diagnostic,
      not a normalization. Near 1 for a well-localized level.
    - `golden_rule`: `E_v^(0) - i Re<chi_v|Gamma|chi_v>/2`, the perturbative
      comparator (the `Gamma = 0` levels plus the first-order width). The
      expectation is taken with the bilinear c-product, so on an ECS grid it is
      complex in general and only its REAL part is a width; the discarded
      imaginary part is a tail-amplitude residue, negligible for a level
      localized in the real region and a sign that the comparator is
      inapplicable for one that is not. This is
      what eMoScat and the thesis actually computed. Agreement with `energies`
      means the level is perturbative; divergence means it is genuinely broad
      and the non-perturbative treatment is load-bearing. `nan` where no
      comparator level could be paired, and all-`nan` when `golden_rule=False`.
      The distance guard also produces `nan` when the Gamma-induced real
      shift exceeds half the local level spacing (the strongly
      non-perturbative regime, where `nan` is the honest answer) and when
      output levels are near-degenerate.
    """

    energies: npt.NDArray[np.complex128]
    widths: npt.NDArray[np.float64]
    states: npt.NDArray[np.complex128]
    residuals: npt.NDArray[np.float64]
    real_weight: npt.NDArray[np.float64]
    golden_rule: npt.NDArray[np.complex128]


def _check_shared_real_nodes(grid_a: FemDvrEcsGrid, grid_b: FemDvrEcsGrid) -> None:
    """Reject two nuclear grids that do not share every real node.

    The two-angle stability test compares eigenvalues of two discretizations
    that must differ ONLY in their ECS tail angle; a different real-region mesh
    makes the residuals meaningless. Called before anything is laid onto either
    grid, so a mismatch surfaces as this message rather than as a downstream
    numpy broadcast error.
    """
    assert_shared_real_nodes(grid_a, grid_b, what="nuclear_grid_a and nuclear_grid_b")


def _default_window(
    Vd: npt.NDArray[np.complex128],
    Gamma: npt.NDArray[np.float64],
    grid: FemDvrEcsGrid,
    atol: float,
) -> tuple[float, float, float, float]:
    """`Re` spanning the anion curve over the REAL nodes, `Im` down to `-max Gamma`.

    The real span is taken from the real nodes only (the ECS tail's continued
    `v0(z)` is complex and says nothing about where levels lie), and covers the
    whole curve so neither the well bottom nor levels above the neutral
    dissociation limit `v0(inf) = 0` are cut.

    **The `Im` band is sized for AUTODETACHMENT widths only.** `Gamma_v =
    <chi_v|Gamma|chi_v> <= max_R Gamma(R)`, so `-max Gamma` is a correct floor
    for a level below the anion dissociation limit. It is NOT a bound on the
    NUCLEAR (dissociative) width of a level ABOVE that limit: that width is
    generated by the ECS rotation of the tail and bears no relation to
    `Gamma(R)` at all -- it can be orders of magnitude larger (a barrierless
    dissociative width is ~1e-3 Ha). Such levels fall outside this window and
    are silently absent from the result. Pass an explicit `window` with a low
    enough `im_lo` to look for them. When `Gamma` is ~0 over the whole grid the
    band degenerates to `+-atol` and NO dissociative level whatsoever can be
    represented; that case warns.
    """
    real = grid.points.imag == 0.0
    v = Vd[real].real
    gmax = float(Gamma.max()) if Gamma.size else 0.0
    if gmax <= atol:
        warnings.warn(
            f"lcp_resonance_levels: the default window's Im band is [-{atol:.1e}, "
            f"{atol:.1e}] because max Gamma(R) = {gmax:.3e} <= atol. That band can "
            "represent only real (bound) levels: any DISSOCIATIVE level -- one above "
            "the anion dissociation limit, whose width comes from the ECS tail and "
            "is unrelated to Gamma(R) -- is excluded from the result. Pass an "
            "explicit `window` with a lower `im_lo` if you are looking for those.",
            UserWarning,
            stacklevel=3,
        )
    return (float(v.min()), float(v.max()), -float(max(gmax, atol)), atol)


def _levels_from(
    grid_a: FemDvrEcsGrid,
    grid_b: FemDvrEcsGrid,
    mu: float,
    W_a: npt.NDArray[np.complex128],
    W_b: npt.NDArray[np.complex128],
    window: tuple[float, float, float, float],
    rel_tol: float,
    atol: float,
) -> tuple[npt.NDArray[np.complex128], npt.NDArray[np.float64], npt.NDArray[np.complex128]]:
    """Diagonalize `T(mu) + diag(W)` on both grids, keep the angle-stable states.

    Returns `(energies, residuals, states)` with `states` shape
    `(n_levels, grid_a.n)`, c-product-normalized, taken from grid `a`.

    Warns (and leaves the state Euclidean-normalized instead) if a state is
    numerically SELF-ORTHOGONAL under the bilinear c-product,
    `|sum_i c_i^2| < _C_NORM_TOL`: dividing by that square root would amplify
    round-off into a meaningless vector while silently reporting a normalized
    one. `energies`/`residuals` are unaffected -- normalization is a property
    of the eigenvector, not of the eigenvalue.
    """
    E_a, V_a = eigen(kinetic(grid_a, mu) + np.diag(W_a))
    E_b, _ = eigen(kinetic(grid_b, mu) + np.diag(W_b))
    energies, residuals, idx = match_angle_stable(E_a, E_b, window, rel_tol=rel_tol, atol=atol)
    states = np.empty((idx.size, grid_a.n), dtype=np.complex128)
    for k, j in enumerate(idx):
        c = V_a[:, j].astype(np.complex128)
        norm2 = complex(c_product(c, c))
        if abs(norm2) < _C_NORM_TOL:
            warnings.warn(
                f"lcp_resonance_levels: level {k} (Re E = {energies[k].real:.6g}) is "
                f"numerically self-orthogonal under the c-product (|sum c^2| = "
                f"{abs(norm2):.3e} < {_C_NORM_TOL:.0e}); leaving it Euclidean-"
                "normalized instead. Its `states` row (and any overlap computed "
                "from it, including `golden_rule`) is not trustworthy.",
                UserWarning,
                stacklevel=3,
            )
            states[k] = c
            continue
        states[k] = c / np.sqrt(norm2)
    return energies, residuals, states


def lcp_resonance_levels(
    nuclear_grid_a: FemDvrEcsGrid,
    nuclear_grid_b: FemDvrEcsGrid,
    mu: float,
    Vd_a: npt.NDArray[np.complex128],
    Vd_b: npt.NDArray[np.complex128],
    Gamma: npt.NDArray[np.float64],
    *,
    window: tuple[float, float, float, float] | None = None,
    n_levels: int | None = None,
    rel_tol: float = 1e-4,
    atol: float = 1e-8,
    golden_rule: bool = True,
) -> ResonanceLevels:
    """Quasi-bound levels of `H_N = T(mu) + V_d(R) - i Gamma(R)/2`.

    The Born-Oppenheimer approximation to the 2-D model's resonance energies:
    step 1 (the fixed-`R` electronic pole, `local_complex_potential`) supplies
    the complex curve; this is step 2, the nuclear eigenvalue problem in it. The
    thesis's `H_LCP` (Vana 2017 Eq. 1.65).

    `nuclear_grid_a`/`nuclear_grid_b` must share every real node and differ only
    in their ECS tail angle -- that is what makes the two spectra comparable.
    `Vd_a`/`Vd_b` are the curve laid onto each grid (identical on the real
    nodes, differing in the continued tail); `Gamma` is real and tail-zero, so
    the same array serves both.

    Physical levels are selected by two-angle stability (`match_angle_stable`);
    the rotated dissociative continuum fails that test and drops out. Levels
    with `Im E > atol` are unphysical and are dropped with a warning.

    See `docs/physics/lcp-resonance-levels.md`.
    """
    if mu <= 0.0:
        raise ValueError(f"mu must be positive, got {mu}")
    for name, arr, grid in (
        ("Vd_a", Vd_a, nuclear_grid_a),
        ("Vd_b", Vd_b, nuclear_grid_b),
        ("Gamma", Gamma, nuclear_grid_a),
    ):
        if arr.shape != (grid.n,):
            raise ValueError(f"{name} has shape {arr.shape}, expected ({grid.n},)")
    _check_shared_real_nodes(nuclear_grid_a, nuclear_grid_b)

    if window is None:
        window = _default_window(Vd_a, Gamma, nuclear_grid_a, atol)

    half_i_gamma = 0.5j * Gamma
    energies, residuals, states = _levels_from(
        nuclear_grid_a,
        nuclear_grid_b,
        mu,
        Vd_a - half_i_gamma,
        Vd_b - half_i_gamma,
        window,
        rel_tol,
        atol,
    )

    physical = energies.imag <= atol
    if not physical.all():
        warnings.warn(
            f"lcp_resonance_levels: dropped {int((~physical).sum())} level(s) with "
            f"Im E > atol = {atol:.1e} (unphysical: a growing state; the tolerance "
            "admits round-off-level positive Im E). Usually an over-wide window "
            "or an under-resolved grid.",
            UserWarning,
            stacklevel=2,
        )
    energies, residuals, states = energies[physical], residuals[physical], states[physical]

    if n_levels is not None:
        energies, residuals, states = (energies[:n_levels], residuals[:n_levels], states[:n_levels])

    widths = np.maximum(0.0, -2.0 * energies.imag)
    real_mask = nuclear_grid_a.points.imag == 0.0
    dens = np.abs(states) ** 2
    total = dens.sum(axis=1)
    real_weight = np.divide(
        dens[:, real_mask].sum(axis=1),
        total,
        out=np.zeros_like(total),
        where=total > 0.0,
    )

    gr = np.full(energies.size, np.nan + 1j * np.nan, dtype=np.complex128)
    if golden_rule and energies.size:
        try:
            E0, _resid0, chi0 = _levels_from(
                nuclear_grid_a,
                nuclear_grid_b,
                mu,
                Vd_a,
                Vd_b,
                (window[0], window[1], -atol, atol),
                rel_tol,
                atol,
            )
        except ValueError:
            # The Gamma=0 comparator problem can genuinely have no angle-
            # stable state in this window: a level near/above the
            # dissociation limit already carries a nonzero Im E from V_d's
            # own complex ECS-tail continuation (no Gamma needed), which the
            # tight [-atol, atol] comparator band excludes even though the
            # primary (wider-window) solve above correctly kept it as
            # physical. This is a failure of the DIAGNOSTIC comparator, not
            # of the primary result -- leave `gr` all-nan and move on.
            E0 = np.empty(0, dtype=np.complex128)
            chi0 = np.empty((0, nuclear_grid_a.n), dtype=np.complex128)
        if E0.size:
            g1 = np.array([c_product(c, Gamma * c).real for c in chi0])
            # Pair each complex level to the nearest comparator level in Re
            # E, but only accept a pairing within a physically plausible
            # distance -- otherwise a level whose true comparator is simply
            # missing (dropped by the window above, or never existed) gets
            # silently glued to an unrelated one. The natural distance scale
            # is half the local Re-E spacing between NEIGHBORING levels in
            # this same output spectrum (the vibrational quantum): a
            # comparator farther than that is closer to some other level's
            # true partner than to this one. With fewer than two levels
            # there is no such spacing to measure, so fall back to half the
            # window's Re-span (the whole region a comparator could
            # plausibly belong to).
            if energies.size >= 2:
                gaps = np.diff(energies.real)  # energies is ascending in Re
                local_spacing = np.empty(energies.size, dtype=np.float64)
                local_spacing[0] = gaps[0]
                local_spacing[-1] = gaps[-1]
                if energies.size > 2:
                    local_spacing[1:-1] = np.minimum(gaps[:-1], gaps[1:])
            else:
                local_spacing = np.full(energies.size, window[1] - window[0])
            max_dist = 0.5 * local_spacing

            dist = np.abs(energies.real[:, None] - E0.real[None, :])
            near = np.argmin(dist, axis=1)
            nearest_dist = dist[np.arange(energies.size), near]
            paired = nearest_dist <= max_dist
            gr[paired] = E0[near[paired]].real - 0.5j * g1[near[paired]]

    return ResonanceLevels(
        energies=np.asarray(energies, dtype=np.complex128),
        widths=np.asarray(widths, dtype=np.float64),
        states=np.asarray(states, dtype=np.complex128),
        residuals=np.asarray(residuals, dtype=np.float64),
        real_weight=np.asarray(real_weight, dtype=np.float64),
        golden_rule=np.asarray(gr, dtype=np.complex128),
    )


def _check_angle_bound(model: ResonanceModel, *grids: FemDvrEcsGrid) -> None:
    """Reject nuclear grids whose ECS tail angle reaches or exceeds the model's bound.

    Strict rejection at the boundary itself: the derivation (Hvizdos et al.
    2018, Sec. II) requires `4*theta < pi/2`, so `theta == max_nuclear_ecs_
    angle_deg` (`4*theta == pi/2`) is already the marginal, non-decaying
    case, not a safe edge.
    """
    bound = getattr(model, "max_nuclear_ecs_angle_deg", None)
    if bound is None:
        return
    for g in grids:
        worst = max((el.angle_deg for el in g.spec.elements), default=0.0)
        if worst >= bound:
            raise ValueError(
                f"nuclear grid ECS angle {worst} deg reaches or exceeds this "
                f"model's max_nuclear_ecs_angle_deg = {bound} deg; at or "
                "beyond it the interaction potential diverges under the "
                "rotation (Hvizdos et al., Phys. Rev. A 97, 022704 (2018), "
                "Sec. II)"
            )


@overload
def resonance_levels(
    model: ResonanceModel,
    nuclear_grid_a: FemDvrEcsGrid,
    nuclear_grid_b: FemDvrEcsGrid,
    elec_grid_a: FemDvrEcsGrid,
    elec_grid_b: FemDvrEcsGrid,
    *,
    re_half_width: float = ...,
    im_half_width: float = ...,
    resid_tol: float = ...,
    window: tuple[float, float, float, float] | None = ...,
    n_levels: int | None = ...,
    rel_tol: float = ...,
    atol: float = ...,
    golden_rule: bool = ...,
    return_curve: Literal[False] = ...,
) -> ResonanceLevels: ...


@overload
def resonance_levels(
    model: ResonanceModel,
    nuclear_grid_a: FemDvrEcsGrid,
    nuclear_grid_b: FemDvrEcsGrid,
    elec_grid_a: FemDvrEcsGrid,
    elec_grid_b: FemDvrEcsGrid,
    *,
    re_half_width: float = ...,
    im_half_width: float = ...,
    resid_tol: float = ...,
    window: tuple[float, float, float, float] | None = ...,
    n_levels: int | None = ...,
    rel_tol: float = ...,
    atol: float = ...,
    golden_rule: bool = ...,
    return_curve: Literal[True],
) -> tuple[ResonanceLevels, npt.NDArray[np.complex128], npt.NDArray[np.float64]]: ...


# bool catch-all (open()-style): callers holding a runtime flag forward it
# directly; the union return is narrowed by the Literal overloads above when
# the flag is literal.
@overload
def resonance_levels(
    model: ResonanceModel,
    nuclear_grid_a: FemDvrEcsGrid,
    nuclear_grid_b: FemDvrEcsGrid,
    elec_grid_a: FemDvrEcsGrid,
    elec_grid_b: FemDvrEcsGrid,
    *,
    re_half_width: float = ...,
    im_half_width: float = ...,
    resid_tol: float = ...,
    window: tuple[float, float, float, float] | None = ...,
    n_levels: int | None = ...,
    rel_tol: float = ...,
    atol: float = ...,
    golden_rule: bool = ...,
    return_curve: bool = ...,
) -> (
    ResonanceLevels | tuple[ResonanceLevels, npt.NDArray[np.complex128], npt.NDArray[np.float64]]
): ...


def resonance_levels(
    model: ResonanceModel,
    nuclear_grid_a: FemDvrEcsGrid,
    nuclear_grid_b: FemDvrEcsGrid,
    elec_grid_a: FemDvrEcsGrid,
    elec_grid_b: FemDvrEcsGrid,
    *,
    re_half_width: float = 0.05,
    im_half_width: float = 0.05,
    resid_tol: float = 1e-3,
    window: tuple[float, float, float, float] | None = None,
    n_levels: int | None = None,
    rel_tol: float = 1e-4,
    atol: float = 1e-8,
    golden_rule: bool = True,
    return_curve: bool = False,
) -> ResonanceLevels | tuple[ResonanceLevels, npt.NDArray[np.complex128], npt.NDArray[np.float64]]:
    """Quasi-bound levels of `model`'s anion, straight from the model.

    Runs the electronic pole walk ONCE (`resonance_pole_walk`, seeded from the
    asymptotic anion bound state exactly as `local_complex_potential` does),
    lays the resulting curve onto BOTH nuclear grids with `_assemble_lcp`, and
    diagonalizes (`lcp_resonance_levels`). `E_res(R)` at real `R` does not
    depend on the nuclear tail angle, so the second grid costs one extra nuclear
    diagonalization and nothing else.

    `nuclear_grid_b` must share `nuclear_grid_a`'s real segments and quadrature
    and differ only in its ECS tail angle -- conventionally a SMALLER angle,
    which is always safe against the model's divergence bound.

    If `return_curve`, also returns `(Vd_a, Gamma)` -- the very curve the levels
    were computed in, on `nuclear_grid_a`. A caller that wants both the levels
    and the LCP curve (to solve `lcp_da_cross_section` in it, or to plot it
    under the levels) MUST take this route rather than calling
    `local_complex_potential` separately: that would repeat the expensive
    electronic walk AND, if any setting differed, report a curve that is not
    the one the levels came from.
    """
    _check_angle_bound(model, nuclear_grid_a, nuclear_grid_b)
    _check_shared_real_nodes(nuclear_grid_a, nuclear_grid_b)

    shift, gamma_w = _walk_from_anion_seed(
        model,
        nuclear_grid_a,
        elec_grid_a,
        elec_grid_b,
        re_half_width=re_half_width,
        im_half_width=im_half_width,
        resid_tol=resid_tol,
    )

    Vd_a, Gamma = _assemble_lcp(model, nuclear_grid_a, shift, gamma_w)
    Vd_b, _ = _assemble_lcp(model, nuclear_grid_b, shift, gamma_w)

    pts = nuclear_grid_a.points
    real = pts.imag == 0.0
    bound_region = Vd_a[real].real < np.asarray(model.v0(pts[real].real)).real
    if np.any(Gamma[real][bound_region] > _CLOSED_REGION_GAMMA_TOL):
        warnings.warn(
            "resonance_levels: Gamma(R) is nonzero where the anion curve lies "
            "BELOW the neutral (v0 > E_res), where autodetachment is closed "
            "(Vana 2017, Sec. 1.5). The widths downstream are suspect.",
            UserWarning,
            stacklevel=2,
        )

    levels = lcp_resonance_levels(
        nuclear_grid_a,
        nuclear_grid_b,
        model.mu,
        Vd_a,
        Vd_b,
        Gamma,
        window=window,
        n_levels=n_levels,
        rel_tol=rel_tol,
        atol=atol,
        golden_rule=golden_rule,
    )
    if return_curve:
        return levels, Vd_a, Gamma
    return levels
