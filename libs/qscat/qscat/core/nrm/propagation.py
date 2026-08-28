"""Propagate the extended-space packet and half-Fourier-transform it back.

`extended_hamiltonian` (`extended.py`) turns the frequency-domain resolvent
sum PRA 77 Eq. (52) solves into a time-LOCAL propagation problem: launch
`Psi_ext(0) = [Psi_d(R,0); 0; ...; 0]` (the discrete-state block carries the
launch state of `initial_packet`; every arm block starts at zero) and
propagate under

    i d/dt Psi_ext = H_ext Psi_ext.

This module runs that propagation and inverts it -- recovers the SAME
`Psi_d^TI(R;E)` Eq. (52) solves for directly, for a whole batch of energies
at once, from one propagation of `initial_packet`'s low-rank launch basis.
`extended.py`'s module docstring is the RESUMMATION argument
(eliminating the arms from `(E - H_ext)^-1` reproduces `F(E)` exactly); this
module is the other half -- the time-domain route to that same resolvent.

DERIVATION of the `-i` prefactor (module-local; not itself a numbered
equation in either paper -- it is the standard half-Fourier / one-sided
Laplace inversion of `i d/dt Psi = H Psi`, specialized to this extended
Hamiltonian). Define the half-Fourier transform at a real energy `E`,
regularized as `E + i*0+`:

    Psi_hat(E) = Int_0^infty dt e^{i(E + i0)t} Psi(t).

Differentiate `e^{iEt} Psi(t)` and integrate by parts, using
`d/dt Psi = -i H_ext Psi`:

    Int_0^infty dt d/dt(e^{iEt} Psi) = [e^{iEt} Psi(t)]_0^infty
        = iE Psi_hat(E) + Int_0^infty dt e^{iEt} (-i H_ext Psi)
        = iE Psi_hat(E) - i H_ext Psi_hat(E).

The `t -> infty` boundary term vanishes PROVIDED every eigenmode of `H_ext`
decays. That is a condition on the ingredients, not a free gift of ECS: it
holds for the COMPLETE projected-state sum and can fail for a truncated one
(`_warn_if_diverging` carries the measurements and checks it at runtime).
Granting it, the `t = 0` term is `-Psi(0)`, hence

    -Psi(0) = i(E - H_ext) Psi_hat(E)   =>   (E - H_ext) Psi_hat(E) = i Psi(0).

Solving `(E - H_ext) Psi_TI = Psi(0)` directly (the frequency-domain form
this module never actually assembles or factors) would give
`Psi_TI = (E - H_ext)^-1 Psi(0) = -i Psi_hat(E)`, i.e.

    Psi_d^TI(R;E) = -i * Int_0^infty dt e^{iEt} Psi_d(R,t),

restricted to the `d`-block (block 0) since that is the only block
`Psi_d^TI` names. `test_transform_of_a_single_decaying_mode_is_the_resolvent`
is the smallest possible check of this identity -- a 1x1 "Hamiltonian" whose
propagated state and transformed resolvent can both be written in closed
form.

`r` singular-vector columns are stepped per energy WINDOW
(`extended.LaunchBasis`); the per-energy packet
`Psi_d(R,t;E_j) = sum_m coeffs[m,j]
U_m(R,t)` is reconstructed at every step from the shared propagation (legal
because `H_ext` is energy-independent, so the superposition commutes with
the propagator) and fed to both the transform accumulator and the
diagnostics -- no per-step snapshot history is stored.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp

from qscat.core.time_dependent import quadrature_weights
from qscat.dvr import FemDvrEcsGrid, dvr_first_derivative_at_node
from qscat.evolution import make_pade_stepper

from .extended import LaunchBasis

__all__ = ["TdNrmResult", "propagate_nrm"]


@dataclass(frozen=True)
class TdNrmResult:
    """The propagated, transformed packet plus PRA 47 Eq. (4.4)-(4.6)
    diagnostics.

    All arrays are batched over the `n_E` energies of the `LaunchBasis` that
    was propagated. `centroid`/`momentum`/`survival` are computed on the
    REAL nuclear region only (see `_record`'s conjugation note); with
    `nuclear_grid=None` (the diagnostics-off path) `centroid`/`momentum` are
    empty and `survival` falls back to the full `|Psi|^2` sum.
    """

    psi_d: npt.NDArray[np.complex128]  # (N_R, n_E) -- Psi_d(R;E), the transform
    time: npt.NDArray[np.float64]  # (n_steps+1,)
    survival: npt.NDArray[np.float64]  # (n_steps+1, n_E) -- Eq. (4.4) S(t)
    centroid: npt.NDArray[np.float64]  # (n_steps+1, n_E) -- Eq. (4.5) <R>_t
    momentum: npt.NDArray[np.float64]  # (n_steps+1, n_E) -- Eq. (4.6) <P>_t
    unabsorbed: npt.NDArray[np.float64]  # (n_E,) -- survival[-1], packet norm
    # still in the real region at t_max (what ECS absorption has NOT eaten)
    rank: int  # number of columns actually stepped (LaunchBasis.rank)


def _real_derivative_matrix(
    grid: FemDvrEcsGrid, real_idx: npt.NDArray[np.intp]
) -> npt.NDArray[np.complex128]:
    """`(len(real_idx), n)` row-stack of `dvr_first_derivative_at_node` over
    `real_idx` (the grid's real, unscaled nodes). Only those rows are ever
    needed -- `_record` reads `deriv @ d` straight, no further row masking --
    so this returns the thin real-rows-only matrix rather than a dense
    `(n, n)` one with unused ECS-tail rows. Built once per `propagate_nrm`
    call -- `O(n_real * n)`, cheap enough that caching it (e.g. an
    `lru_cache` keyed on `id(grid)`) is not worth the complexity.
    """
    n = grid.n
    mat = np.zeros((real_idx.size, n), dtype=np.complex128)
    for row, i in enumerate(real_idx):
        mat[row, :] = dvr_first_derivative_at_node(grid, int(i))
    return mat


def _diagnostics_setup(
    grid: FemDvrEcsGrid | None,
) -> tuple[
    npt.NDArray[np.bool_] | None,
    npt.NDArray[np.complex128] | None,
    npt.NDArray[np.float64] | None,
    npt.NDArray[np.float64] | None,
]:
    """Everything `_record` needs, computed ONCE (not per step): the
    real-region mask, the real-rows-only derivative matrix, the per-real-node
    `sqrt(w)` DVR weight (see the UNITS NOTE in `_record`'s docstring for why
    this factor matters), and the real node positions. `grid=None` returns
    all-`None` (the diagnostics-off path).
    """
    if grid is None:
        return None, None, None, None
    mask = grid.real_points <= grid.R0
    real_idx = np.flatnonzero(mask)
    deriv = _real_derivative_matrix(grid, real_idx)
    sqrt_w_real = np.sqrt(grid.weights[mask].real)
    r_real = grid.real_points[mask]
    return mask, deriv, sqrt_w_real, r_real


def _record(
    m: int,
    d: npt.NDArray[np.complex128],
    mask: npt.NDArray[np.bool_] | None,
    deriv: npt.NDArray[np.complex128] | None,
    sqrt_w_real: npt.NDArray[np.float64] | None,
    r_real: npt.NDArray[np.float64] | None,
    survival: npt.NDArray[np.float64],
    centroid: npt.NDArray[np.float64],
    momentum: npt.NDArray[np.float64],
) -> None:
    """PRA 47 Eq. (4.4)-(4.6): survival, centroid, and momentum of `Psi_d(R,t)`.

    CONJUGATION NOTE: these are probability densities over the REAL nuclear
    region, not c-products -- `|psi|^2` is the physically meaningful
    quantity there, and the complex ECS tail is an absorber excluded by
    `mask`. This is the one place in this module where `np.conjugate` is
    correct rather than a c-product violation.

    UNITS NOTE (Eq. 4.6): `deriv @ d` (`dvr_first_derivative_at_node`'s own
    contract) returns `d/dR Psi` as a wavefunction VALUE, while `d[mask, :]`
    itself is a DVR COEFFICIENT (`sqrt(w)` already absorbed, module docstring
    convention). `<P>_t = int dR Psi* (-i d/dR) Psi / S(t)` is a
    coefficient-space sum `sum_i conj(coeff_i) * sqrt(w_i) * (-i value'_i)`
    -- mixing a coefficient against a value with no `sqrt(w_i)` between them
    silently drops the DVR quadrature weight from one factor of the
    integrand. Measured impact of the missing factor (2026-08-19, a Gaussian
    with analytic p=1.7 on `qscat.core.grids.nuclear_grid()`): 5.93 without
    `sqrt_w_real`, 1.70 with it -- not a constant rescaling, since GLL
    weights vary within and across elements, so the bug would have looked
    like grid-dependent physics rather than a fixed offset.

    ZERO-SURVIVAL NOTE: when `s` (survival) is exactly zero -- the packet has
    fully decayed out of the real region, the expected end state of a
    converged run (not an error) -- `centroid`/`momentum` come back NaN by
    plain IEEE `0/0` semantics (the numerators are built from the same `d`
    that made `s` zero, so they vanish identically too). This is deliberate,
    not "hold the last valid value": there is no packet left to have a
    position or momentum. `np.errstate` only silences the resulting
    `RuntimeWarning`; it does not change the NaN.
    """
    if mask is None:
        survival[m] = np.sum(np.abs(d) ** 2, axis=0)
        return
    assert deriv is not None
    assert sqrt_w_real is not None
    assert r_real is not None
    d_real = d[mask, :]
    dens = np.abs(d_real) ** 2
    s = dens.sum(axis=0)
    survival[m] = s
    with np.errstate(invalid="ignore", divide="ignore"):
        centroid[m] = (r_real[:, None] * dens).sum(axis=0) / s
        momentum[m] = np.real(
            (np.conjugate(d_real) * sqrt_w_real[:, None] * (-1j * (deriv @ d))).sum(axis=0) / s
        )


def _warn_if_diverging(survival: npt.NDArray[np.float64]) -> None:
    """Warn if the packet ENDS heavier than it started -- `H_ext` had a
    growing eigenmode and `psi_d` is exponentially wrong.

    The half-Fourier transform this module inverts assumes every eigenmode of
    `H_ext` decays. That holds for the COMPLETE projected-state sum
    (`n_states=None`) and can fail for a truncated one: `V_dn` and `E_n` are
    complex (electronic ECS), so the arrow matrix's anti-Hermitian part is
    indefinite and eigenvalues can leave the lower half-plane even though
    every diagonal block on its own is dissipative. Measured 2026-08-19 by
    dense `eig` (see the transform-identity tests): F2's production
    ingredients give `max Im(E)` = +5.29e-05 / +2.48e-04 / +5.84e-04 /
    +1.59e-03 at `n_states` = 1 / 2 / 3 / 6 against +1.07e-14 at
    `n_states=0`, and the complete sum comes back marginally stable
    (+2.4e-12, i.e. zero to solver tolerance). The failure is not monotone in
    `n_states` -- a small fixture measured here loses its growing modes again
    by `n_states=12` -- so a truncated set is not reliably safe OR reliably
    unsafe, which is exactly why this is a runtime check and not a
    precondition on `n_states`.

    Without this, a truncated run returns a plausible-looking `psi_d` whose
    error is O(1) or worse: on F2 at `n_states=3` the relative vector error
    against the time-independent solve walks 1.23 at T=2000 to 1.4e3 at
    T=64000, with nothing in the returned object flagging it.

    TWO CONDITIONS, because neither alone is enough. `survival` legitimately
    RISES early (the arms feed amplitude back into the discrete state), so
    the first test is against the FIRST sample: a converged run ends orders
    below where it started. But a growing mode can still leave the final
    sample below the initial one if the physical packet absorbed first --
    measured on F2 at `n_states=3`, `S(T)/S(0)` reads 1.12 / 0.48 / 0.59 /
    2.92 at T = 2e3 / 4e3 / 8e3 / 16e3, i.e. the T=8000 run is diverging
    (`rel = 2.3`) while ending at 0.59 of its start. So the second test is
    against the MINIMUM over the run: growth after the minimum is the actual
    signature, and it reads 1.47 (T=8000) / 7.22 (T=16000) there against
    exactly 1.0000 for every converging run measured (F2 `n_states=0` at
    T=8000 and T=16000, and the N2 deck of the transform-identity test at
    T=1000 and T=4000 -- all four decay monotonically to their last sample).
    The 1.2 factor is slack against a round-off-noisy floor, sized on that
    gap; a false alarm costs a warning, a missed divergence costs a silently
    wrong answer.

    A run can still slip through -- F2 at T=4000 ends exactly at its own
    minimum, having neither absorbed nor yet visibly grown -- so a silent
    return is not a certificate. It is a cheap alarm on the failure that
    would otherwise be invisible, not a proof of stability.
    """
    first, last = survival[0], survival[-1]
    floor = survival.min(axis=0)
    grew = (last > first) | (last > 1.2 * floor)
    if bool(np.any(grew)):
        safe_first = np.where(first > 0.0, first, np.inf)
        safe_floor = np.where(floor > 0.0, floor, np.inf)
        warnings.warn(
            f"propagate_nrm: the packet grew rather than decayed in at least "
            f"one energy column (final/initial survival "
            f"{float(np.max(last / safe_first)):.3g}, final/minimum "
            f"{float(np.max(last / safe_floor)):.3g}), so H_ext has a growing "
            "eigenmode and psi_d is exponentially wrong rather than merely "
            "under-converged. The usual cause is a TRUNCATED arm set: "
            "n_states must be None (the complete projected-state sum) for "
            "the time-dependent route -- any truncation needs its own dense "
            "eig check that max Im(eigenvalue) <= 0.",
            stacklevel=2,
        )


def propagate_nrm(
    h_ext: sp.spmatrix,
    launch: LaunchBasis,
    nuclear_grid: FemDvrEcsGrid | None,
    *,
    dt: float,
    n_steps: int,
    order: int = 3,
) -> TdNrmResult:
    """Propagate `launch.vectors` under `h_ext` and transform back to energy.

    Parameters
    ----------
    h_ext : sparse matrix
        The extended block Hamiltonian (`extended_hamiltonian`), or (in
        tests) any square matrix the launch vectors' first block matches.
    launch : LaunchBasis
        The low-rank launch state (`initial_packet`): `vectors` are the `r`
        columns actually propagated, `coeffs` reconstructs each of the
        `n_E` energies from them, `e_total` is the transform frequency for
        each energy.
    nuclear_grid : FemDvrEcsGrid or None
        The nuclear grid the `d`-block lives on. `None` disables the
        real-region diagnostics (`centroid`/`momentum` come back empty,
        `survival` falls back to the full `|Psi|^2` sum, and `psi_d` is the
        whole propagated state rather than just the `d`-block) -- the path
        the synthetic unit tests exercise; real callers always pass a grid.
    dt, n_steps, order
        Passed to `qscat.evolution.make_pade_stepper`; `n_steps + 1` samples
        (`t = 0, dt, ..., n_steps*dt`) are taken.

    Returns
    -------
    TdNrmResult

    Raises
    ------
    ValueError
        If `h_ext` is not square, if `launch.vectors` does not have one row
        per `h_ext` column, or if `nuclear_grid.n` does not evenly divide
        `h_ext`'s size -- a cheap guard against propagating a `LaunchBasis`
        against a mismatched `h_ext`/`nuclear_grid` pair (e.g. a refined
        grid against a stale cached `h_ext`), rather than silently slicing
        the wrong block out of `psi`.
    """
    e = launch.e_total
    psi = np.asarray(launch.vectors, dtype=np.complex128)  # (N_ext, r)
    coeffs = launch.coeffs  # (r, n_E)
    if h_ext.shape[0] != h_ext.shape[1]:
        raise ValueError(f"h_ext must be square, got shape {h_ext.shape}")
    if psi.shape[0] != h_ext.shape[0]:
        raise ValueError(
            f"launch.vectors has {psi.shape[0]} rows but h_ext is {h_ext.shape[0]}x{h_ext.shape[1]}"
        )
    n_r = h_ext.shape[0] if nuclear_grid is None else nuclear_grid.n
    if nuclear_grid is not None and h_ext.shape[0] % n_r != 0:
        raise ValueError(
            f"h_ext size {h_ext.shape[0]} is not a multiple of "
            f"nuclear_grid.n={n_r} -- nuclear_grid does not match the grid "
            "h_ext was built from"
        )

    step = make_pade_stepper(h_ext, dt, order)
    w = quadrature_weights(n_steps + 1)
    t = dt * np.arange(n_steps + 1, dtype=np.float64)

    acc = np.zeros((n_r, e.size), dtype=np.complex128)
    survival = np.empty((n_steps + 1, e.size))
    centroid = np.empty_like(survival) if nuclear_grid is not None else np.empty((0, 0))
    momentum = np.empty_like(centroid)
    mask, deriv, sqrt_w_real, r_real = _diagnostics_setup(nuclear_grid)

    for m in range(n_steps + 1):
        # Reconstruct the per-energy packet from the r propagated columns --
        # legal because H_ext is energy-independent, so the superposition
        # commutes with the propagator.
        d = psi[:n_r, :] @ coeffs  # (N_R, n_E)
        acc += (w[m] * dt) * np.exp(1j * e * t[m])[None, :] * d
        _record(m, d, mask, deriv, sqrt_w_real, r_real, survival, centroid, momentum)
        if m < n_steps:
            psi = step(psi)

    _warn_if_diverging(survival)
    return TdNrmResult(
        psi_d=-1j * acc,
        time=t,
        survival=survival,
        centroid=centroid,
        momentum=momentum,
        unabsorbed=survival[-1].copy(),
        rank=launch.rank,
    )
