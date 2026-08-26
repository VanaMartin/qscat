"""Fixed-`R` electronic resonance-pole machinery for the LCP approximation.

`V_d(R) = Re(E_pole(R))`, `Gamma(R) = max(0, -2*Im(E_pole(R)))`, where
`E_pole(R)` is the two-ECS-angle-matched resonance pole (`qscat.ecs.
find_resonance_pole`) of the fixed-R electronic Hamiltonian
`H_el(R) = -1/2 d^2/dr^2 + model.surface(r, R)`. Because `model.surface`
ALREADY INCLUDES `v0(R)` (the neutral-molecule channel potential), `V_d =
Re(E_pole)` directly -- adding `v0(R)` again would double-count it. (Contrast
`projects/n2_ti_cross_section/vres.py`, whose `v_eff_el` EXCLUDES `v0`, so
that code adds `v0(R) + E_res(R)` separately; the two bookkeeping schemes
must and do agree on the observable `V_d(R)`, checked in
`validation/n2/test_lcp_vres_parity.py::test_matches_n2_vres_oracle`.)

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
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

from qscat.dvr import FemDvrEcsGrid, eigen, kinetic
from qscat.ecs import find_resonance_pole
from qscat.exceptions import ConvergenceError
from qscat.linalg import c_product

from ..dissociation import anion_electronic_states

if TYPE_CHECKING:
    from qscat.model import ResonanceModel

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

# `resonance_eigenstate_at_peak_width`: smallest `Gamma(R)` still worth
# re-solving a pole at. The search walks real `R` widest-first, so once it
# reaches a width this small every remaining point is narrower still, and a
# pole that narrow is not separable from the discretized continuum by the
# two-angle match -- the loop stops rather than grinding through them.
_MIN_RESOLVABLE_GAMMA = 1e-4


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
