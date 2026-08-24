"""Cross sections by propagation instead of by a per-energy solve.

`dissociation.nrm_da_cross_section` and
`vibrational_excitation.nrm_ve_cross_section` build `F(E)` and solve PRA 77
Eq. (52) once per energy. This module reaches the SAME `Psi_d^+(R;E)` by the
time-dependent route of PRA 47 Sec. II: launch `initial_packet`'s low-rank
factorization of Eq. (2.5), propagate it under `extended_hamiltonian`'s arrow
block Hamiltonian (`propagation.propagate_nrm`), and half-Fourier-transform
back. The extractions are then the shipped ones --
`dissociation.da_sigma_from_psi` (Eq. 54) for attachment,
`vibrational_excitation.t_resonant` + `t_background` (Eq. 28/31/37) for
vibrational excitation, all UNCHANGED -- so the two routes differ in how
`Psi_d^+` is obtained and in nothing else.

That is the point: `td_nrm_da_cross_section` and `td_nrm_ve_cross_section`
take their time-independent counterparts' arguments in the same order, so
swapping one identifier for the other is the whole comparison. Anything else
(a different `eps_e`, a different flux surface, a different `k_i`, a
conjugation introduced on one side only) would compare two different physical
quantities and still look like a passing test, which is why `eps_e` here is
read from the same `continue_to_tail(ing.v_d_discrete, ...)` at the same
`_boundary_node` rather than recomputed, and why the VE route calls
`t_resonant`/`t_background` rather than re-deriving Eq. (31)/(37).

`markovian=True` runs the LOCAL sibling of all this: PRA 47 Eq. (2.11) collapses
the memory kernel to `i[Delta_L - (i/2)Gamma_L] delta(R-R') delta(t)`, so the LCP
IS the Markovian approximation to the nonlocal model, and Eq. (2.15) is what
remains. `extended.lcp_limit_hamiltonian` + `extended.lcp_initial_packet` replace
`extended_hamiltonian` + `initial_packet`; everything downstream is unchanged. It
takes `qscat.core.lcp.local_complex_potential`'s `(Vd, Gamma)` -- NOT
`ing.v_d_discrete`, which is missing the level shift `Delta_L` and measurably does
not reproduce `lcp_da_cross_section` (`lcp_limit_hamiltonian`'s docstring). With
no arms the propagated matrix is `N_R` square rather than `(1 + n_states) * N_R`,
which on F2's deck is 974 against 53570 and 4 s against ~30 min. For VE the
substitution has a second half: the local doorway `sqrt(Gamma_L/2pi) chi_v`
replaces `V^+_dk chi_v` at the EXIT as well as at the launch, which is what
makes the route reproduce the LCP's own `S_{v'<-v}` rather than a hybrid of
the two models -- and `include_background=True` is refused there, since
Eq. (37)'s background is built from a `phi_d` the local model does not have.

WHY DISSOCIATIVE ATTACHMENT IS THE HARD CHANNEL FOR THIS ROUTE. In the
vibrational-excitation case the packet decays IN PLACE -- autodetachment
empties the discrete state and `S(t)` falls to round-off. In DA the packet
must physically LEAVE: the surviving amplitude travels out along `V_d(R)` to
the ECS absorber, and the transform only converges once it has. Two
consequences, both measured on F2 (2026-08-19, the fixture of
`libs/qscat/tests/test_nrm_td_cross_section.py`):

- THE NUCLEAR GRID MUST RESOLVE THE EXIT WAVE. On F2 at `E = 0.02-0.05 Ha`
  the dissociating wave carries `K_R = 55.6-64.2`, i.e. a wavelength of
  0.098-0.113 bohr. A coarse N2-style nuclear deck cannot represent it: the
  packet's centroid creeps 2.66 -> 2.77 bohr and then oscillates, `S(t)`
  plateaus at 0.67, and the transform is dominated by trapped flux that
  LOOKS like a bound component. On the FINE per-molecule production deck
  (`validation/diatomic/config.py`'s F2 nuclear segments, 974 points, 65
  points/bohr over `R = 3-10.7`, largest node spacing 0.023 bohr) the
  centroid climbs monotonically and `<P>_t` rises to `K_R` -- the signature
  of a packet that is actually leaving. `centroid` growing roughly linearly
  is the diagnostic to look at; oscillation about a fixed `R` means the
  answer is grid, not physics.
- `S(t)` DOES NOT REACH ZERO ON F2. Its `V_d(R)` has a 0.0223 Ha well
  (minimum -0.149264 Ha at R = 3.363 against the F+F- asymptote -0.126931)
  supporting >= 24 near-real modes (`|Im E| = 1.5e-7 ... 7.7e-6`) that the
  launch populates with 5.08e-3 of its real-region norm. Those modes need
  `T >~ 1e7` to decay and no affordable propagation removes them, so `S(t)`
  PLATEAUS rather than decays -- measured, it falls 0.94 -> 0.007 as the
  dissociating wave crosses into the absorber and then flattens at
  0.006-0.009 from `T ~ 12000` on. `unabsorbed` is reported and warned on,
  not gated: a convergence criterion keyed to an absolute `S(T)/S(0)` floor
  would warn forever on this molecule. What converges instead is the
  observable -- `sigma_DA` reads the wavefunction VALUE at the outermost
  real node, where those well-localized modes have almost no amplitude, and
  it goes stationary (ratio inside [0.977, 1.024], no trend) over
  `T = 12000-15000` while `S(t)` is doing nothing at all. CONVERGENCE HERE
  MEANS `sigma_DA` STATIONARY IN `T`, NOT `S(T)` SMALL.
"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

from qscat.dvr import FemDvrEcsGrid

from .coupling import v_dk_plus
from .discrete_state import DiscreteState
from .dissociation import _boundary_node, da_sigma_from_psi
from .extended import (
    LaunchBasis,
    extended_hamiltonian,
    initial_packet,
    lcp_initial_packet,
    lcp_limit_hamiltonian,
)
from .ingredients import NrmIngredients, nrm_ingredients
from .nonlocal_potential import continue_to_tail
from .propagation import TdNrmResult, propagate_nrm
from .vibrational_excitation import t_background, t_resonant

if TYPE_CHECKING:
    import scipy.sparse as sp

    from qscat.model import ResonanceModel

__all__ = ["td_nrm_da_cross_section", "td_nrm_ve_cross_section"]

#: Default `unabsorbed_tol`: warn above 1% of the initial real-region norm.
#: Not arbitrary -- it is where the alarm discriminates on the one molecule
#: measured. On F2 (2026-08-19, `test_nrm_td_cross_section.py`'s fixture) the
#: worst-energy `sigma_TD/sigma_TI` runs 1.13 at `S(T)/S(0) = 0.045`, 1.065 at
#: 0.018, and 1.014 at 0.009, so 1e-2 fires on the runs whose cross section is
#: still visibly truncated and stays quiet on the converged one. It sits just
#: ABOVE F2's own spectral floor (`S(t)` plateaus at 0.006-0.009, the bound
#: modes of the module docstring), which is the tightest a threshold can be on
#: this molecule without warning forever.
_UNABSORBED_TOL = 1e-2


def _check_markovian_arguments(
    markovian: bool,
    Vd: npt.NDArray[np.complex128] | None,
    Gamma: npt.NDArray[np.float64] | None,
    n_states: int | None,
) -> None:
    """Reject argument sets whose failure would otherwise be a silent wrong answer.

    Passing `Vd`/`Gamma` WITHOUT `markovian` is the dangerous direction: the
    nonlocal branch would ignore them and return a perfectly plausible
    nonlocal cross section under a call that reads as a local one. The
    reverse (`markovian` without them) cannot be defaulted -- the local
    complex potential is an electronic-structure computation
    (`qscat.core.lcp.local_complex_potential`), not something this module
    can reconstruct from `ing`. `n_states` is rejected rather than ignored
    because the local Hamiltonian has no arms at all, so any value for it
    describes something this call does not do.
    """
    if markovian:
        if Vd is None or Gamma is None:
            raise ValueError(
                "markovian=True needs both Vd and Gamma (qscat.core.lcp."
                "local_complex_potential's two returns) -- the local complex "
                "potential cannot be derived from the nonlocal ingredients"
            )
        if n_states is not None:
            raise ValueError(
                f"n_states={n_states} is meaningless with markovian=True: the "
                "local limit has no projected-state arms to truncate"
            )
    elif Vd is not None or Gamma is not None:
        raise ValueError(
            "Vd/Gamma are the markovian=True arguments; passing them to the "
            "nonlocal route would silently ignore them"
        )


def _setup_propagation(
    nuclear_grid: FemDvrEcsGrid,
    elec_grid: FemDvrEcsGrid,
    model: ResonanceModel,
    phi_d: DiscreteState,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    energies: npt.NDArray[np.float64],
    *,
    ingredients: NrmIngredients | None,
    n_states: int | None,
    markovian: bool,
    Vd: npt.NDArray[np.complex128] | None,
    Gamma: npt.NDArray[np.float64] | None,
    rank_tol: float,
) -> tuple[LaunchBasis, sp.csr_matrix, NrmIngredients | None, npt.NDArray[np.complex128]]:
    """The launch state, the Hamiltonian, and the nuclear curve behind them.

    Shared by both entry points, so the nonlocal/Markovian branch is written
    ONCE: `sigma_DA` and `sigma_VE` differ in what they contract `Psi_d`
    with, not in how `Psi_d` is obtained.

    Returns `(launch, h_ext, ing, v_d_full)`. In the Markovian branch `ing`
    is `None` -- no ingredients are built, which is the whole cost argument
    for that route -- and `v_d_full` is the LOCAL complex potential
    `Vd - (i/2)Gamma`; in the nonlocal branch it is Eq. (20)'s
    `v_d_discrete` continued to the tail. Either way it is the curve whose
    outermost real node is the DA channel threshold `eps_e`.
    """
    if markovian:
        assert Vd is not None and Gamma is not None  # narrowed by _check_markovian_arguments
        gamma = np.asarray(Gamma, dtype=np.float64)
        v_res = np.asarray(Vd, dtype=np.complex128) - 0.5j * gamma
        launch = lcp_initial_packet(nuclear_grid, gamma, eps, chi, v_init, energies)
        return launch, lcp_limit_hamiltonian(nuclear_grid, model, v_res), None, v_res

    real = nuclear_grid.points.imag == 0.0
    R_desc = np.sort(nuclear_grid.points[real].real)[::-1]
    ing = ingredients or nrm_ingredients(elec_grid, model, phi_d, R_desc)
    v_d_full = continue_to_tail(ing.v_d_discrete, ing.R, nuclear_grid)
    launch = initial_packet(
        nuclear_grid,
        elec_grid,
        model,
        phi_d,
        ing,
        eps,
        chi,
        v_init,
        energies,
        n_states=n_states,
        rank_tol=rank_tol,
    )
    h_ext = extended_hamiltonian(ing, nuclear_grid, model, n_states=n_states)
    return launch, h_ext, ing, v_d_full


def _warn_if_truncated(
    res: TdNrmResult, caller: str, dt: float, n_steps: int, unabsorbed_tol: float
) -> None:
    """Warn when a packet has not finished leaving the real region.

    The failure mode is not a blow-up: a truncated half-Fourier transform
    returns a finite, positive, plausible cross section that is simply
    wrong, and nothing else in the returned array says so.

    THE TEST IS A DA TEST, and on VE it is CONSERVATIVE rather than wrong.
    `sigma_DA` reads the wavefunction VALUE at the outermost real node, so
    "has the packet left?" is very nearly the convergence question itself;
    `sigma_VE` contracts against `chi_f` in the interaction region, and
    converges once the amplitude THERE has decayed, whether or not anything
    has crossed the box. Measured on F2 (2026-08-24,
    `test_nrm_td_cross_section.py`'s fixture, `dt=2, T=2000`): this warns at
    `unabsorbed/S(0) = 0.938` while `sigma_VE` is already within 6e-5 of the
    time-independent route. It is kept at one threshold anyway -- it reports
    a fact ("the packet is still in the box") that is true in both cases,
    and a second, looser VE threshold would only make the DA calibration
    ambiguous.
    """
    left = res.unabsorbed / np.where(res.survival[0] > 0.0, res.survival[0], np.inf)
    if bool(np.any(left > unabsorbed_tol)):
        warnings.warn(
            f"{caller}: the packet still holds "
            f"{float(np.max(left)):.3g} of its initial real-region norm at "
            f"T={dt * n_steps:g} (tolerance {unabsorbed_tol:g}), so the "
            "half-Fourier transform is truncated and the cross section may "
            "be confidently wrong rather than merely noisy. Propagate "
            "longer, or -- if S(t) has plateaued rather than still falling "
            "-- establish convergence by showing the cross section is "
            "stationary in T.",
            stacklevel=3,
        )


def td_nrm_da_cross_section(
    nuclear_grid: FemDvrEcsGrid,
    elec_grid: FemDvrEcsGrid,
    model: ResonanceModel,
    phi_d: DiscreteState,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    E: float | npt.ArrayLike,
    *,
    ingredients: NrmIngredients | None = None,
    n_states: int | None = None,
    markovian: bool = False,
    Vd: npt.NDArray[np.complex128] | None = None,
    Gamma: npt.NDArray[np.float64] | None = None,
    dt: float,
    n_steps: int,
    order: int = 3,
    rank_tol: float = 1e-6,
    unabsorbed_tol: float = _UNABSORBED_TOL,
) -> npt.NDArray[np.float64]:
    """`sigma_DA(E)` (bohr^2) by propagation -- the TD route to Eq. (54).

    The arguments up to `n_states` are `nrm_da_cross_section`'s, in its
    order; the rest are the propagation knobs.

    ONE propagation serves the whole energy batch: `H_ext` is
    energy-independent, so `initial_packet`'s `r` singular-vector columns are
    stepped once and every energy is reconstructed from them
    (`propagate_nrm`'s reconstruction), rather than one propagation per
    energy.

    Parameters
    ----------
    nuclear_grid, elec_grid : FemDvrEcsGrid
        The nuclear and electronic radial grids. The nuclear grid must
        resolve the exit wave (module docstring) -- a coarse deck returns a
        confidently wrong number.
    model : ResonanceModel
        The molecule.
    phi_d : DiscreteState
        The discrete-state choice under test.
    eps, chi : ndarray
        Neutral vibrational energies and states (`qscat.core.vibrational`).
    v_init : int
        Initial vibrational level.
    E : float or array
        Incident electron kinetic energy or energies (hartree).
        Non-positive entries give `0.0`, as they do on the TI route.
    ingredients : NrmIngredients, optional
        Precomputed ingredients; built here if omitted.
    n_states : int, optional
        Number of projected electronic states. `None` (the default, and the
        only value that should be used in production) keeps the COMPLETE
        sum. This is a CORRECTNESS knob, not a performance one: `V_dn` and
        `E_n` are complex, so a truncated arm set can leave `H_ext` with a
        growing eigenmode, the transform's premise fails, and `psi_d` comes
        back exponentially wrong rather than under-converged.
        `propagate_nrm` warns at runtime when that happens. Rejected with
        `markovian=True`, which has no arms to truncate.
    markovian : bool, optional
        Propagate the LOCAL (LCP) limit instead -- PRA 47 Eq. (2.11)/(2.15),
        `lcp_limit_hamiltonian` plus `lcp_initial_packet`'s doorway. Requires
        `Vd`/`Gamma` and ignores `elec_grid`, `phi_d`, `ingredients`,
        `n_states` and `rank_tol` (none of them enter the local model; no
        ingredients are built). Default `False`.
    Vd, Gamma : ndarray, optional
        The local complex potential, on the full nuclear grid -- exactly
        `qscat.core.lcp.local_complex_potential`'s two returns, under the
        names `lcp_da_cross_section` gives them. Required when `markovian`,
        rejected otherwise. `Vd` is the RESONANCE POSITION `E_res(R)` (plus
        `V_0`), NOT `NrmIngredients.v_d_discrete` -- see
        `lcp_limit_hamiltonian` for why Eq. (2.14) forces that reading and
        what the other one measures.
    dt, n_steps, order : float, int, int
        Propagation step, number of steps, and diagonal-Pade order
        (`qscat.evolution.make_pade_stepper`). The total propagation error
        is `T`- and `dt`-separable and the two parts add IN QUADRATURE:
        `truncation(T) ~ 0.40*sqrt(S(T)/S(0))` and, at order 3,
        `propagation(dt=1) = 1.43e-4` falling as `dt^6`
        (`test_nrm_propagation.py`'s measurements).
    rank_tol : float, optional
        Launch-basis truncation (`initial_packet`); `0.0` keeps every mode.
    unabsorbed_tol : float, optional
        Warn when any energy column still holds more than this fraction of
        its initial real-region norm at `T`. Default `1e-2`.

    Returns
    -------
    ndarray
        `sigma_DA` per energy; scalar-shaped for a scalar `E`.

    Raises
    ------
    ValueError
        If `markovian` is set without both `Vd` and `Gamma`, if either is
        passed without `markovian`, or if `markovian` is combined with an
        explicit `n_states`.

    Warns
    -----
    UserWarning
        If any column's `unabsorbed / survival[0]` exceeds `unabsorbed_tol`.
        A truncated transform is the failure mode that produces plausible
        wrong numbers -- it does not blow up, it quietly reports the
        cross section of a packet that has not finished leaving. On a
        molecule whose `S(t)` plateaus above the tolerance for spectral
        reasons (F2, module docstring), this warning is expected and the
        convergence evidence has to come from `sigma_DA` itself being
        stationary in `T`, not from `unabsorbed`.
    """
    _check_markovian_arguments(markovian, Vd, Gamma, n_states)

    e_arr = np.atleast_1d(np.asarray(E, dtype=np.float64))
    out = np.zeros(e_arr.size, dtype=np.float64)
    scalar = np.isscalar(E) or (isinstance(E, np.ndarray) and np.ndim(E) == 0)
    open_idx = np.flatnonzero(e_arr > 0.0)
    if open_idx.size == 0:
        return np.asarray(out[0] if scalar else out, dtype=np.float64)

    launch, h_ext, _, v_d_full = _setup_propagation(
        nuclear_grid,
        elec_grid,
        model,
        phi_d,
        eps,
        chi,
        v_init,
        e_arr[open_idx],
        ingredients=ingredients,
        n_states=n_states,
        markovian=markovian,
        Vd=Vd,
        Gamma=Gamma,
        rank_tol=rank_tol,
    )
    # In the Markovian branch this reads the LOCAL curve, and is identical to
    # lcp_da_cross_section's own `eps_e = Vd[b].real`: the -(i/2)Gamma term is
    # purely imaginary and cannot move it.
    eps_e = float(v_d_full[_boundary_node(nuclear_grid)].real)

    res = propagate_nrm(h_ext, launch, nuclear_grid, dt=dt, n_steps=n_steps, order=order)
    _warn_if_truncated(res, "td_nrm_da_cross_section", dt, n_steps, unabsorbed_tol)

    for col, ie in enumerate(open_idx):
        out[ie] = da_sigma_from_psi(
            nuclear_grid,
            model.mu,
            res.psi_d[:, col],
            float(launch.e_total[col]),
            eps_e,
            float(e_arr[ie]),
        )

    return np.asarray(out[0] if scalar else out, dtype=np.float64)


def td_nrm_ve_cross_section(
    nuclear_grid: FemDvrEcsGrid,
    elec_grid: FemDvrEcsGrid,
    model: ResonanceModel,
    phi_d: DiscreteState,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    vprimes: list[int],
    E: float | npt.ArrayLike,
    *,
    ingredients: NrmIngredients | None = None,
    n_states: int | None = None,
    include_background: bool = True,
    markovian: bool = False,
    Vd: npt.NDArray[np.complex128] | None = None,
    Gamma: npt.NDArray[np.float64] | None = None,
    dt: float,
    n_steps: int,
    order: int = 3,
    rank_tol: float = 1e-6,
    unabsorbed_tol: float = _UNABSORBED_TOL,
) -> npt.NDArray[np.float64]:
    """`sigma_{v_init->v'}(E)` (bohr^2) by propagation -- the TD route to Eq. (28).

    The arguments up to `include_background` are `nrm_ve_cross_section`'s, in
    its order; the rest are the propagation knobs. Swapping one identifier
    for the other IS the comparison.

    ONE propagation serves every energy AND every final channel: `H_ext` is
    energy-independent (so the whole batch rides on `initial_packet`'s
    singular-vector columns) and `Psi_d^+` does not depend on `v'` at all --
    the final channel enters only through the contraction below, exactly as
    it does on the time-independent route.

    `T^res` is `vibrational_excitation.t_resonant(chi_f, V^+_dk_f, Psi_d)`
    and `T^bg` is `t_background`, both called UNCHANGED: the half-Fourier
    transform returns the same `Psi_d^+(R;E)` the Eq. (52) solve returns, so
    every downstream convention is inherited rather than re-derived. In
    particular PRA 77 Eq. (34)/(35)'s collapse of `V^{-*}_{dk}` to a
    NON-conjugated `V^+_{dk}`, and Eq. (37)'s bra carrying `phi^+` at the
    FINAL channel energy, are properties of those two functions and are not
    restated here. `T^bg` contains no `Psi_d` at all -- it is energy-domain
    and static -- so it is bit-identical between the two routes and any
    disagreement that appears ONLY with `include_background=True` is a
    combination error in the resonant term, not a background error.

    WHY VE IS THE EASY CHANNEL FOR THIS ROUTE, where
    `td_nrm_da_cross_section` documents DA as the hard one. `T^res` is a
    contraction of `Psi_d` against `chi_f`, which lives in the INTERACTION
    region: the packet decays IN PLACE by autodetachment and the transform
    converges as `S(t)` falls, with no requirement that anything traverse
    the grid. DA instead reads the wavefunction VALUE at the outermost real
    node, so its packet must physically cross the box before the transform
    means anything. On N2 that is the difference between a converged run at
    `T = 4000` and one at `T = 12000` on a nuclear deck five times finer.
    A molecule whose DA channel is OPEN (F2 above threshold) pays the DA
    cost even when the observable is VE, because the same `Psi_d` must
    converge either way.

    Parameters
    ----------
    nuclear_grid, elec_grid : FemDvrEcsGrid
        The nuclear and electronic radial grids.
    model : ResonanceModel
        The molecule.
    phi_d : DiscreteState
        The discrete-state choice under test.
    eps, chi : ndarray
        Neutral vibrational energies and states (`qscat.core.vibrational`).
    v_init : int
        Initial vibrational level.
    vprimes : list of int
        Final vibrational levels.
    E : float or array
        Incident electron kinetic energy or energies (hartree).
        Non-positive entries give `0.0`, as they do on the TI route.
    ingredients : NrmIngredients, optional
        Precomputed ingredients; built here if omitted.
    n_states : int, optional
        Number of projected electronic states. `None` (the default, and the
        only value that should be used in production) keeps the COMPLETE
        sum -- a CORRECTNESS requirement of the propagation, not a
        performance knob; see `td_nrm_da_cross_section`. Rejected with
        `markovian=True`, which has no arms to truncate.
    include_background : bool, default True
        Add `T^bg` (Eq. 37) to `T^res` before squaring -- PRA 77's
        "nonlocal + background" curve against its bare "nonlocal" one.
        Must be `False` when `markovian`: see below.
    markovian : bool, optional
        Propagate the LOCAL (LCP) limit instead -- PRA 47 Eq. (2.11)/(2.15).
        The doorway `sqrt(Gamma_L/2pi) chi_v` then replaces `V^+_dk chi_v`
        on BOTH sides: it is the launch state (`lcp_initial_packet`) and, at
        the exit, the coupling `t_resonant` contracts against. Using the
        local kernel on one side and the nonlocal `V^+_dk_f` on the other
        would be a third model, neither the shipped LCP nor the nonlocal
        one. With that substitution this route reproduces the LCP's own
        `S_{v'<-v} = <sqrt(Gamma/2pi) chi_v' | Psi_d>` exactly.
        Requires `Vd`/`Gamma`; ignores `elec_grid`, `phi_d`, `ingredients`,
        `n_states` and `rank_tol` (no ingredients are built). Default
        `False`.
    Vd, Gamma : ndarray, optional
        The local complex potential on the full nuclear grid --
        `qscat.core.lcp.local_complex_potential`'s two returns. Required
        when `markovian`, rejected otherwise. `Vd` is the RESONANCE POSITION
        `E_res(R)` (plus `V_0`), NOT `NrmIngredients.v_d_discrete`; see
        `extended.lcp_limit_hamiltonian`.
    dt, n_steps, order : float, int, int
        Propagation step, number of steps, and diagonal-Pade order
        (`qscat.evolution.make_pade_stepper`).
    rank_tol : float, optional
        Launch-basis truncation (`initial_packet`); `0.0` keeps every mode.
    unabsorbed_tol : float, optional
        Warn when any energy column still holds more than this fraction of
        its initial real-region norm at `T`. Default `1e-2`.

    Returns
    -------
    ndarray
        `sigma_{v_init->v'}` per energy; scalar `E` returns shape
        `(len(vprimes),)`, array `E` returns `(len(E), len(vprimes))` --
        `nrm_ve_cross_section`'s (and `driven.ve_cross_section`'s)
        convention. A closed channel (`E_tot - eps_vf <= 0`, or `E <= 0`)
        contributes `0.0`.

    Raises
    ------
    ValueError
        If `markovian` is set without both `Vd` and `Gamma`, if either is
        passed without `markovian`, if `markovian` is combined with an
        explicit `n_states`, or if `markovian` is combined with
        `include_background=True`.

    Warns
    -----
    UserWarning
        If any column's `unabsorbed / survival[0]` exceeds
        `unabsorbed_tol` -- see `td_nrm_da_cross_section`.
    """
    _check_markovian_arguments(markovian, Vd, Gamma, n_states)
    if markovian and include_background:
        # PRA 77 Eq. (37)'s background is built from `phi_d` and the P-space
        # scattering states -- objects the local model does not have, and
        # which `markovian=True` documents itself as ignoring. The LCP curve
        # this repository ships (and every LCP curve PRA 77 plots) is the
        # bare resonant term; "local resonant + nonlocal background" is a
        # third model, so it is refused rather than silently assembled.
        raise ValueError(
            "markovian=True needs include_background=False: the local model "
            "carries no phi_d, so Eq. (37)'s background is not part of it -- "
            "adding it would be a hybrid of the local and nonlocal models"
        )

    e_arr = np.atleast_1d(np.asarray(E, dtype=np.float64))
    out = np.zeros((e_arr.size, len(vprimes)), dtype=np.float64)
    scalar = np.isscalar(E) or (isinstance(E, np.ndarray) and np.ndim(E) == 0)
    open_idx = np.flatnonzero(e_arr > 0.0)
    if open_idx.size == 0:
        return np.asarray(out[0] if scalar else out, dtype=np.float64)

    launch, h_ext, ing, _ = _setup_propagation(
        nuclear_grid,
        elec_grid,
        model,
        phi_d,
        eps,
        chi,
        v_init,
        e_arr[open_idx],
        ingredients=ingredients,
        n_states=n_states,
        markovian=markovian,
        Vd=Vd,
        Gamma=Gamma,
        rank_tol=rank_tol,
    )

    res = propagate_nrm(h_ext, launch, nuclear_grid, dt=dt, n_steps=n_steps, order=order)
    _warn_if_truncated(res, "td_nrm_ve_cross_section", dt, n_steps, unabsorbed_tol)

    # The LOCAL exit coupling, energy-independent by construction (`Gamma_L(R)
    # = Gamma(E_res(R), R)`), so it is built once rather than per channel.
    local_coupling = (
        np.sqrt(np.asarray(Gamma, dtype=np.float64) / (2.0 * np.pi)).astype(np.complex128)
        if markovian
        else None
    )

    for col, ie in enumerate(open_idx):
        e_kin = float(e_arr[ie])
        e_total = float(launch.e_total[col])
        psi_d = res.psi_d[:, col]
        for jv, vp in enumerate(vprimes):
            excess = e_total - float(eps[vp])
            if excess <= 0.0:
                continue  # closed channel
            if local_coupling is not None:
                v_dk_f = local_coupling
            else:
                assert ing is not None  # the nonlocal branch always builds them
                v_dk_f = continue_to_tail(
                    v_dk_plus(elec_grid, model, phi_d, ing.R, excess), ing.R, nuclear_grid
                )
            t = t_resonant(chi[vp], v_dk_f, psi_d)
            if include_background:
                assert ing is not None  # include_background implies the nonlocal branch
                t += t_background(
                    elec_grid,
                    nuclear_grid,
                    model,
                    phi_d,
                    ing.R,
                    chi[v_init],
                    chi[vp],
                    v_dk_f,
                    e_kin,
                    excess,
                )
            out[ie, jv] = 4.0 * np.pi**3 * abs(t) ** 2 / (2.0 * e_kin)

    return np.asarray(out[0] if scalar else out, dtype=np.float64)
