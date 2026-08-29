"""Memory observables of the time-dependent nonlocal resonance model.

The nonlocal model's memory is not a stored history: `extended.py` resums
Gertitschke & Domcke, Phys. Rev. A 47, 1031 (1993) Eq. (2.1)'s convolution into
one auxiliary nuclear packet `phi_n` per projected electronic state, so the
memory is STATE and can be looked at while it is being made. This module holds
what is looked at.

`local_width` is the Markovian reference the nonlocal observables are read
against: the local limit of the SAME `F(E)` the propagation is built from,
rather than `qscat.core.lcp`'s pole-walk `Gamma`. Both sides of every
comparison then come from one object, and no pole walk is involved -- which
also removes two ways of being quietly wrong. The walk freezes silently at
small `R`; and `local_complex_potential` takes a PAIR of electronic grids at
two DIFFERENT ECS angles, the two-angle match being how the physical pole is
selected, so passing the same grid twice -- `local_complex_potential(model,
nuc, elec, elec)` -- discriminates nothing and reports a freeze radius that is
an artefact of the degenerate pair. That is where the R = 1.5783 quoted in this
sub-project's first design draft came from; it is struck rather than corrected,
because nothing here needs a walk at all. A genuine 35/40-degree pair on the
same deck freezes at R = 2.055 instead.

WHAT `‖phi_n‖^2` MAY BE CALLED, MEASURED. It is NOT a population, so the
partition observable built on it -- `arm_norm` and `arm_norm_by_channel` -- must
be reported as a RELATIVE channel decomposition. Under ECS `H_ext` is complex
SYMMETRIC, so the conjugating norm is conserved by nothing; the one exact
statement left is that what the coupling removes from `Psi_d` it adds to the
arms, and that holds only to the extent `V_dn` is real. The two rates are
`2 Im<Psi_d|sum_n V_dn phi_n>` and `2 Im<phi_n|V_dn Psi_d>`, and their sum is
`4 sum_n Re[conj(Psi_d) phi_n] Im[V_dn]`. Measured 2026-08-27 on
`test_nrm_memory.py`'s N2 gate deck (`n_states=None`, E_kin = 0.10 Ha, 73 arms,
`H_ext` = 13246): that residual is **median 0.82 and max 1.06 of the LARGER of
the two rates** over the shipped 200-step window (median 0.90, same max, over
the full 4000-step run), i.e. the same size as the exchange itself -- the arms
gain conjugating norm several times faster than `Psi_d` loses it. It sits
entirely in the real nuclear region (the ECS tail carries 6e-71 of it over 200
steps and 2e-55 over 4000, `V_dn` being ~1e-13 there by Eq. 67) and within that
in the autodetachment region, not the absorber: R in [2.0, 2.2] 68.6%,
[1.8, 2.0] 17.2%, [2.2, 2.5] 14.2% over the 200-step window, moving outward
with the packet to 57% / 11% / 32% by 4000 steps, with under 0.5% outside
[1.8, 2.5] at either. The electronic rotation leaks into nuclear-space
bookkeeping exactly where the physics is.

The exchange RATE is untouched by this. It is a rate, not an amplitude, and it
never needed the population reading -- a positive `2 Im<Psi_d|sum_n V_dn
phi_n>` is amplitude returning to the discrete state whatever the arms' norm
means.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import cast

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp

from qscat.dvr import FemDvrEcsGrid

__all__ = ["MemoryRecorder", "MemorySpec", "local_width"]


def local_width(
    f_matrix: npt.NDArray[np.complex128], nuclear_grid: FemDvrEcsGrid
) -> npt.NDArray[np.float64]:
    """Local limit of `F(E)`: `Gamma_loc(R_i) = -2 Im[(F sqrt(w))_i / sqrt(w)_i]`.

    `F` (`nonlocal_potential.nonlocal_operator`, PRA 77 Eq. 60-61) is a
    nuclear-DVR COEFFICIENT-space matrix. A local potential in that
    representation is `diag(V(R_i))`, so the local limit of a nonlocal kernel is
    what it does to a function that is CONSTANT over the kernel's span: a
    constant of value 1 has coefficients `sqrt(w_i)`, `(F sqrt(w))_i` is the
    resulting coefficient, and dividing by `sqrt(w_i)` turns it back into a
    value. Those two `sqrt(w)` factors are the coefficient<->value conversion
    and nothing else -- Eq. (60)'s own `sqrt(W)` factors are already absorbed
    (see `nonlocal_potential`'s module docstring) and are NOT reapplied here.

    `diag F` is NOT the local limit. The kernel spans ~10 nodes, so most of the
    row sits off the diagonal: measured on the N2 deck of
    `test_nrm_memory.py`, `diag F` is 0.16x this row sum (0.14x on the
    time-independent decks of `docs/physics/nonlocal-resonance-model.md` Sec.
    9).

    WHAT THIS IS A WIDTH AT, and why it is not `qscat.core.lcp`'s `Gamma`.
    `F(E)` carries ONE total energy, so its local limit is the width at the
    LOCAL electron energy `eps_loc(R) = E - V_0(R)`. That is the quantity
    `nonlocal-resonance-model.md` Sec. 9 records as reproducing Eq. (68)'s
    `Gamma(eps_loc, R) = 2 pi |V_dk+(R)|^2` to median 0.977 on NO and 1.011 on
    F2, and which `test_nrm_memory.py` re-measures at 0.996 on N2.
    `local_complex_potential`'s `Gamma` is a DIFFERENT width -- the one at the
    resonance position `eps_res(R) = E_res(R) - V_0(R)` -- and on N2 the two
    energies are nowhere near each other (`eps_loc` crosses zero at R = 2.4284
    where `eps_res` is still ~0.07 Ha), so the ratio to it runs 0.12-8.9 across
    the region where `Gamma` is nonzero. That is a difference of energy
    argument, not an error in either. Eq. (68) against the pole width AT
    `eps_res` is gated separately, in `test_nrm_coupling.py`.

    WHERE IT IS SUPPORTED, since a consumer may be tempted to divide by it.
    `V_0(R)` is a well, so `eps_loc` is positive only BETWEEN two crossings --
    R = 1.7426 and R = 2.4284 on `test_nrm_memory.py`'s N2 deck at
    E_kin = 0.10 Ha, 31 of 153 real nodes. Outside that window the electron
    channel is CLOSED and `Gamma_loc` decays smoothly rather than vanishing: it
    is 1.4e-4 one node outside the outer crossing and ~1e-12 by R = 4, and it
    is EXACTLY zero at none of the 122 closed-channel nodes (largest 8.6e-4,
    3.4% of the peak, at the node just inside the INNER crossing; median
    2.5e-10). The ECS tail is 3e-22 of the peak. So `-<Psi_d|Gamma_loc|Psi_d>`
    is well defined everywhere and needs no masking, but any RATIO taken
    against `Gamma_loc` divides by ~1e-10 wherever the packet sits outside the
    open window.

    IT IS NOT CLAMPED, and 32 of 179 nodes come back negative -- all of them
    round-off rather than physics: the most negative is -3.8e-9 (1.5e-7 of the
    peak), the median negative magnitude is 7e-29, and every one sits where the
    true width is zero (R <= 0.79 and the ECS tail). Measured consequence for
    the Markovian rate: over a 200-step N2 run those nodes contribute a
    would-be GAIN of 4.1e-48, i.e. 5e-44 of `|exchange_local|`, which stays
    strictly negative (max -3.5e-6). "The Markovian limit can only lose"
    therefore survives numerically without a clamp -- and a clamp would hide
    the one thing a negative entry would be worth knowing about.

    Parameters
    ----------
    f_matrix : ndarray
        `(N_R, N_R)` complex-symmetric `F(E)` in nuclear DVR coefficient space
        (`nonlocal_operator`).
    nuclear_grid : FemDvrEcsGrid
        The nuclear grid `f_matrix` was built on; supplies the DVR weights.

    Returns
    -------
    ndarray
        `(N_R,)` real `Gamma_loc(R)`, one entry per nuclear node. Not clamped
        at zero: a negative entry is a diagnostic (`F`'s anti-Hermitian part is
        negative semidefinite, so a genuinely negative local width means the
        ingredients are wrong), and clamping would hide it.

    Raises
    ------
    ValueError
        If `f_matrix` is not `(nuclear_grid.n, nuclear_grid.n)`.
    """
    f = np.asarray(f_matrix, dtype=np.complex128)
    n = nuclear_grid.n
    if f.shape != (n, n):
        raise ValueError(
            f"f_matrix has shape {f.shape}, expected ({n}, {n}) -- one row and "
            "column per nuclear DVR node"
        )
    sqrt_w = np.sqrt(np.asarray(nuclear_grid.weights, dtype=np.complex128))
    out: npt.NDArray[np.float64] = -2.0 * np.imag((f @ sqrt_w) / sqrt_w)
    return out


# `eq=False` because `gamma_local` is an ndarray: the generated `__eq__` would
# compare it with `==`, get an array back, and raise on the `bool()`; and the
# frozen dataclass's generated `__hash__` would raise too. Neither is used, and
# identity semantics are the right answer for a spec that carries a grid-sized
# array -- but a caller memoizing on a `MemorySpec` would otherwise meet the
# exception rather than the design decision.
@dataclass(frozen=True, eq=False)
class MemorySpec:
    """What `propagate_nrm` needs before it will record the memory observables.

    Recording is OPT-IN and this is the switch: `memory=None` (the default)
    leaves the propagation loop exactly as it was.

    Attributes
    ----------
    gamma_local : ndarray
        `(N_R,)` real `Gamma_loc(R)` on the FULL nuclear grid -- real nodes and
        ECS tail -- normally `local_width(f_matrix, nuclear_grid)`. It is NOT
        clamped and must not be: see `local_width`'s docstring for what a
        negative entry means and why hiding one would be the wrong trade.
    n_channels : int or None
        How many per-arm series `arm_norm_by_channel` keeps. `None` (the
        default) keeps EVERY channel, which is correct wherever it fits and is
        `n_arm x (n_steps+1) x n_E x 8` bytes -- 2.3 MB for the N2 gate deck,
        ~113 MB for a production F2 one. An integer `k` keeps the FIRST `k` arm
        blocks of `h_ext`, i.e. `ing`'s own adiabatic-tracking order.

        THAT IS POSITIONAL, NOT "the `k` largest by norm", and must not be
        called that: which arms are largest is knowable only after the run, and
        knowing it would need exactly the history the truncation exists to
        avoid. Selecting by `|V_dn|` at `t = 0` would be one-pass but assumes
        coupling strength predicts which channel ends up carrying amplitude --
        a finding the campaign might test, not a rule to bake in here.
        `arm_peak` is what makes a truncated run honest: it is a running
        maximum over EVERY channel, so a caller can see whether the first `k`
        were the right ones instead of assuming it, and re-run pointed at the
        interesting ones if the time series is wanted. `arm_norm` (the
        aggregate) always covers every arm regardless.

        MEASURED on the N2 gate deck (73 arms, 4000 steps): ranked by
        `arm_peak`, the top blocks are 0, 1, 2, 3, 4, 5, 6, 8 -- so the first
        four ARE the four largest here, and the first eight miss only block 7.
        At the arm-norm peak the first four carry 93.2% of the total and the
        first eight 99.3%. That is a property of THIS deck, not a rule; it is
        what `arm_peak` exists to let a caller check per molecule.
    """

    gamma_local: npt.NDArray[np.float64]
    n_channels: int | None = None

    def __post_init__(self) -> None:
        """Validate shape, realness and `n_channels` -- everything checkable
        without a grid; the length against `nuclear_grid.n` is `MemoryRecorder`'s
        to check, since only it has the grid."""
        g = np.asarray(self.gamma_local)
        if g.ndim != 1:
            raise ValueError(f"gamma_local must be 1-D, got shape {g.shape}")
        if np.iscomplexobj(g):
            raise ValueError(
                "gamma_local must be real -- `local_width` returns a real width, "
                "and the Markovian rate -<Psi_d|Gamma_loc|Psi_d> is real by "
                "construction"
            )
        if self.n_channels is not None and self.n_channels < 0:
            raise ValueError(f"n_channels must be >= 0 or None, got {self.n_channels}")


class MemoryRecorder:
    """Per-step memory observables of the extended-space propagation.

    Built once per `propagate_nrm` call and fed the reconstructed per-energy
    blocks at every step. It stores no snapshots -- four `(n_steps+1, n_E)`
    series, one `(n_steps+1, k, n_E)` and one `(n_arm, n_E)` running maximum,
    and nothing of the state itself.

    WHAT IT MEASURES, AND WHICH OPERATOR. The coupling is read out of `h_ext`
    ITSELF -- rows `0:N_R`, columns `N_R:` for the discrete side and the
    transpose block for the arm side -- not rebuilt from `NrmIngredients`. Two
    consequences, both wanted. The observables describe the operator actually
    being propagated, so a change to `extended_hamiltonian`'s assembly shows up
    here; and nothing is assumed about the coupling block's structure, so a
    dense or otherwise non-`sp.diags` block is handled rather than silently
    mis-read. `test_nrm_memory.py` exercises it on both shapes.

    WHAT IT COSTS, MEASURED DIRECTLY RATHER THAN BY WALL CLOCK. `record` takes
    0.097-0.100 ms per call on the N2 gate deck (13246 unknowns, 73 arms, one
    energy), against a 34.2 ms order-3 Pade step -- **+0.29%** -- plus 0.017 ms
    (+0.05%) for the arm-block reconstruction `propagate_nrm` does only when
    recording. It is the SAME at `n_channels=None` (73 per-channel series) as
    at `n_channels=4`: the cost is the two sparse mat-vecs and the real-region
    copy, not the per-channel write.

    That number is timed around the call, because an A/B of two whole
    propagations cannot resolve it: repeated `memory=None` runs of the same
    deck spread by 1.8-3.0% among themselves, and paired on/off runs measured
    +2.8%, -0.6%, +1.5% and +0.1% on the same machine. Any of those quoted
    alone would be noise reported as a measurement.

    CONJUGATING PRODUCT, RESTRICTED TO THE REAL NUCLEAR REGION -- deliberately
    NOT the c-product the model's own algebra uses. These are rates and
    magnitudes of `|psi|^2`, which is the physically meaningful density on the
    real region; the ECS tail is an absorber and is excluded.

    THE RESTRICTION IS ON THE OUTER INDEX ONLY, which "restricted to the real
    region" can be read as denying. In `exchange` the BRA `Psi_d` is restricted
    and the coupling's COLUMNS are not, so `(C phi)_r` sums over every arm node
    including the ECS tail; likewise the arm rate restricts `phi` and not `d`.
    That is the correct quantity -- it is exactly the coupling term of
    `d/dt||d_real||^2`, where the operator acts on the whole state and only the
    projection is real-region -- and a both-sided restriction would be a
    different, wrong one. The two coincide for the real model, whose `V_dn`
    blocks are diagonal in R (a real row has only its own real column), and
    differ for a general coupling block, which is why
    `test_nrm_memory.py`'s `_synthetic_deck` uses a dense one. The module
    docstring records what Task 1 measured about how far that reading can be
    pushed, and the short version is that `arm_norm` is a RELATIVE channel
    decomposition and not a population.

    Attributes
    ----------
    arm_norm : ndarray
        `(n_steps+1, n_E)` `sum_n ||phi_n||^2` over the real nuclear region,
        summed over ALL arms. NOT A POPULATION -- read it against itself over
        time and against the other channels, never as a probability, and never
        as something that sums with `survival` to a conserved total. Measured
        2026-08-27 (module docstring): the norm the coupling adds to the arms
        exceeds what it removes from `Psi_d` by O(1), so there is no transfer
        to book-keep.
    arm_norm_by_channel : ndarray
        `(n_steps+1, k, n_E)`, `k = n_arm` unless `spec.n_channels` truncates
        it -- the same quantity resolved over the FIRST `k` arm blocks of
        `h_ext` (see `MemorySpec.n_channels`). Same caveat: relative, not a
        population.
    arm_peak : ndarray
        `(n_arm, n_E)` running maximum of each channel's `||phi_n||^2` over the
        whole run -- EVERY channel, whatever `n_channels` does to the time
        series. `O(n_arm x n_E)` and one pass, so it costs nothing and is what
        a truncated run should be read against: it answers "which channels
        received the flux" without claiming that the ones kept were those.
    exchange : ndarray
        `(n_steps+1, n_E)` `2 Im<Psi_d|sum_n V_dn phi_n>`, the rate at which
        the coupling feeds the discrete state. This is a RATE and is unaffected
        by the caveat above. Its SIGN is the observable: positive is amplitude
        returning from the continuum, which the local complex potential cannot
        represent at all.

        UNNORMALIZED, deliberately -- which normalization to read it against is
        the consumer's choice and the two in use differ by orders. Measured
        2026-08-27 on `test_nrm_memory.py`'s N2 gate deck (`n_states=None`,
        E_kin = 0.10 Ha, 73 arms, dt = 1, 4000 steps): `exchange` is positive at
        85 of the 4001 steps, first at t = 132, with raw extremes +8.78e-7
        (t = 155) and -2.59e-4 (t = 1). Divided by `S_d(0)` = 3.63e-3 those
        become +2.42e-4 and -7.15e-2; divided by `S_d(t)`, +1.64e-3 and
        -7.51e-2. Quoting a figure without saying which is how +8.78e-7 and
        +2.42e-4 came to be the same measurement. Restricting to the real
        region rather than the full grid changes none of it to four figures --
        `V_dn` in the ECS tail is ~1e-13 by Eq. (67).
    exchange_local : ndarray
        `(n_steps+1, n_E)` `-<Psi_d|Gamma_loc|Psi_d>`, the Markovian limit of
        the SAME `F(E)`. Non-positive wherever `Gamma_loc` is, by construction.
        Compare it to `exchange` as a DIFFERENCE, or normalize both by
        `survival` -- never divide one by the other, and never divide anything
        by `Gamma_loc`, which is ~1e-10 rather than 0 over most of the grid.
    imbalance : ndarray
        `(n_steps+1, n_E)` `exchange` plus the arms' own coupling rate
        `2 Im<phi|V_dn Psi_d>` -- the residual `4 sum_n Re[conj(Psi_d) phi_n]
        Im[V_dn]` of the module docstring, zero only where `V_dn` is real.
        Under ECS it is not, and this is measured at O(1); it is a first-class
        output for that reason, not a footnote.
    """

    def __init__(
        self,
        h_ext: sp.spmatrix,
        nuclear_grid: FemDvrEcsGrid,
        spec: MemorySpec,
        *,
        n_steps: int,
        n_energies: int,
    ) -> None:
        """Slice the coupling blocks and allocate the series.

        Raises
        ------
        ValueError
            If `spec.gamma_local` does not have one entry per nuclear DVR node,
            or if `h_ext`'s size is not a whole number of nuclear blocks.
        """
        n_r = nuclear_grid.n
        gamma = np.asarray(spec.gamma_local, dtype=np.float64)
        if gamma.size != n_r:
            raise ValueError(
                f"gamma_local has {gamma.size} entries, expected {n_r} -- one per "
                "nuclear DVR node of the grid h_ext was built on"
            )
        n_ext = int(h_ext.shape[0])
        if n_ext % n_r != 0:
            raise ValueError(
                f"h_ext size {n_ext} is not a multiple of nuclear_grid.n={n_r} -- "
                "the memory observables need h_ext's block structure"
            )
        n_arm = n_ext // n_r - 1

        # The real-region mask, identical to `_record`'s: `real_points <= R0`
        # is the same set of nodes as `points.imag == 0`, since `ecs_map` is
        # the identity below the pivot.
        real_idx = np.flatnonzero(nuclear_grid.real_points <= nuclear_grid.R0)
        n_real = int(real_idx.size)
        # Channel-major, so a reshape to (n_arm, n_real, n_E) is the per-arm
        # split with no further bookkeeping.
        arm_real_idx = (n_r * np.arange(n_arm)[:, None] + real_idx[None, :]).ravel()

        # Concrete complex128 CSR, for the row slicing below. Via `csc_matrix`
        # because scipy-stubs accepts a generic `spmatrix` there and not in
        # `csr_matrix` -- the same stub-shaped detour `make_pade_stepper` takes
        # on this matrix. One O(nnz) conversion per `propagate_nrm` call.
        h: sp.csr_matrix[np.complex128] = sp.csc_matrix(h_ext, dtype=np.complex128).tocsr()
        # Rows restricted BEFORE the matvec: only the real-region rows are ever
        # read, so the ECS-tail rows are never multiplied.
        # scipy-stubs types a sliced `csr_matrix` as `spmatrix[Any]`; the
        # `cast`s say what it is at runtime (same interaction as
        # `extended.py`'s `bmat`/`block_diag` casts).
        self._c_up: sp.csr_matrix | None = (
            cast("sp.csr_matrix", h[:n_r, n_r:])[real_idx] if n_arm else None
        )
        self._c_lo: sp.csr_matrix | None = (
            cast("sp.csr_matrix", h[n_r:, :n_r])[arm_real_idx] if n_arm else None
        )

        self._n_r = n_r
        self._n_arm = n_arm
        self._n_real = n_real
        self._real_idx = real_idx
        self._arm_real_idx = arm_real_idx
        self._gamma_real = gamma[real_idx]
        self._n_keep = n_arm if spec.n_channels is None else min(int(spec.n_channels), n_arm)

        shape = (n_steps + 1, n_energies)
        self.arm_norm: npt.NDArray[np.float64] = np.empty(shape)
        self.arm_norm_by_channel: npt.NDArray[np.float64] = np.empty(
            (n_steps + 1, self._n_keep, n_energies)
        )
        # Zeros, not empty: this is a running maximum over a nonnegative
        # quantity, so 0.0 is the correct starting value and not a placeholder.
        self.arm_peak: npt.NDArray[np.float64] = np.zeros((n_arm, n_energies))
        self.exchange: npt.NDArray[np.float64] = np.empty(shape)
        self.exchange_local: npt.NDArray[np.float64] = np.empty(shape)
        self.imbalance: npt.NDArray[np.float64] = np.empty(shape)

    def record(
        self,
        m: int,
        d: npt.NDArray[np.complex128],
        arms: npt.NDArray[np.complex128],
    ) -> None:
        """Fill step `m` from the reconstructed per-energy blocks.

        Parameters
        ----------
        m : int
            Step index, `0 <= m <= n_steps`.
        d : ndarray
            `(N_R, n_E)` discrete-state block, DVR coefficients.
        arms : ndarray
            `(n_arm * N_R, n_E)` stacked arm blocks in `h_ext`'s own block
            order, DVR coefficients.

        Raises
        ------
        ValueError
            If either block's row count disagrees with the `h_ext` this was
            built from.
        """
        if d.shape[0] != self._n_r or arms.shape[0] != self._n_arm * self._n_r:
            raise ValueError(
                f"blocks have {d.shape[0]} / {arms.shape[0]} rows, expected "
                f"{self._n_r} / {self._n_arm * self._n_r}"
            )
        n_e = d.shape[1]
        dr = d[self._real_idx]

        # CONJUGATING product, restricted to the real nuclear region: a
        # probability density, not a c-product. See the class docstring.
        self.exchange_local[m] = -(self._gamma_real[:, None] * np.abs(dr) ** 2).sum(axis=0)

        if self._n_arm == 0:
            self.arm_norm[m] = 0.0
            self.exchange[m] = 0.0
            self.imbalance[m] = 0.0
            return

        assert self._c_up is not None
        assert self._c_lo is not None
        ar = arms[self._arm_real_idx]  # (n_arm * n_real, n_E)
        dens = np.abs(ar.reshape(self._n_arm, self._n_real, n_e)) ** 2
        per_channel = dens.sum(axis=1)  # (n_arm, n_E)
        self.arm_norm[m] = per_channel.sum(axis=0)
        self.arm_norm_by_channel[m] = per_channel[: self._n_keep]
        np.maximum(self.arm_peak, per_channel, out=self.arm_peak)

        # `2 Im<Psi_d|sum_n V_dn phi_n>` and its partner `2 Im<phi|V_dn Psi_d>`,
        # both taken from h_ext's own coupling blocks (class docstring).
        self.exchange[m] = 2.0 * np.einsum("ij,ij->j", dr.conj(), self._c_up @ arms).imag
        arm_rate = 2.0 * np.einsum("ij,ij->j", ar.conj(), self._c_lo @ d).imag
        self.imbalance[m] = self.exchange[m] + arm_rate
