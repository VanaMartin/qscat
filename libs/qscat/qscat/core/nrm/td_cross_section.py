"""`sigma_DA(E)` by propagation instead of by a per-energy solve.

`dissociation.nrm_da_cross_section` builds `F(E)` and solves PRA 77 Eq. (52)
once per energy. This module reaches the SAME `Psi_d^+(R;E)` by the
time-dependent route of PRA 47 Sec. II: launch `initial_packet`'s low-rank
factorization of Eq. (2.5), propagate it under `extended_hamiltonian`'s arrow
block Hamiltonian (`propagation.propagate_nrm`), and half-Fourier-transform
back. The extraction is then `dissociation.da_sigma_from_psi` -- Eq. (54),
UNCHANGED -- so the two routes differ in how `Psi_d^+` is obtained and in
nothing else.

That is the point: `td_nrm_da_cross_section` takes `nrm_da_cross_section`'s
arguments in the same order, so swapping one identifier for the other is the
whole comparison. Anything else (a different `eps_e`, a different flux
surface, a different `k_i`) would compare two different physical quantities
and still look like a passing test, which is why `eps_e` here is read from
the same `continue_to_tail(ing.v_d_discrete, ...)` at the same
`_boundary_node` rather than recomputed.

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

from .discrete_state import DiscreteState
from .dissociation import _boundary_node, da_sigma_from_psi
from .extended import extended_hamiltonian, initial_packet
from .ingredients import NrmIngredients, nrm_ingredients
from .nonlocal_potential import continue_to_tail
from .propagation import propagate_nrm

if TYPE_CHECKING:
    from qscat.model import ResonanceModel

__all__ = ["td_nrm_da_cross_section"]

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
        `propagate_nrm` warns at runtime when that happens.
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
    real = nuclear_grid.points.imag == 0.0
    R_desc = np.sort(nuclear_grid.points[real].real)[::-1]
    ing = ingredients or nrm_ingredients(elec_grid, model, phi_d, R_desc)

    v_d_full = continue_to_tail(ing.v_d_discrete, ing.R, nuclear_grid)
    eps_e = float(v_d_full[_boundary_node(nuclear_grid)].real)

    e_arr = np.atleast_1d(np.asarray(E, dtype=np.float64))
    out = np.zeros(e_arr.size, dtype=np.float64)
    open_idx = np.flatnonzero(e_arr > 0.0)
    if open_idx.size:
        launch = initial_packet(
            nuclear_grid,
            elec_grid,
            model,
            phi_d,
            ing,
            eps,
            chi,
            v_init,
            e_arr[open_idx],
            n_states=n_states,
            rank_tol=rank_tol,
        )
        h_ext = extended_hamiltonian(ing, nuclear_grid, model, n_states=n_states)
        res = propagate_nrm(h_ext, launch, nuclear_grid, dt=dt, n_steps=n_steps, order=order)

        left = res.unabsorbed / np.where(res.survival[0] > 0.0, res.survival[0], np.inf)
        if bool(np.any(left > unabsorbed_tol)):
            warnings.warn(
                f"td_nrm_da_cross_section: the packet still holds "
                f"{float(np.max(left)):.3g} of its initial real-region norm at "
                f"T={dt * n_steps:g} (tolerance {unabsorbed_tol:g}), so the "
                "half-Fourier transform is truncated and sigma_DA may be "
                "confidently wrong rather than merely noisy. Propagate longer, "
                "or -- if S(t) has plateaued rather than still falling -- "
                "establish convergence by showing sigma_DA is stationary in T.",
                stacklevel=2,
            )

        for col, ie in enumerate(open_idx):
            out[ie] = da_sigma_from_psi(
                nuclear_grid,
                model.mu,
                res.psi_d[:, col],
                float(launch.e_total[col]),
                eps_e,
                float(e_arr[ie]),
            )

    scalar = np.isscalar(E) or (isinstance(E, np.ndarray) and np.ndim(E) == 0)
    return np.asarray(out[0] if scalar else out, dtype=np.float64)
