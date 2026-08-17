"""The NRM discrete state `phi_d(r;R)` -- PRA 77's central choice.

The nonlocal approximation is completely determined by the discrete state, and
the paper's headline finding is that the results depend on which one is taken
(Sec. VI). Two of its three choices are implemented, behind one protocol:

- `PhysicalDiscreteState` (choice A, Sec. VI A) -- the fixed-nuclei electron
  scattering function at the real resonance energy `Re E_res(R)`, smoothly
  truncated by Eq. (69). This is the "intuitive" choice, and the one whose
  Born-Oppenheimer breakdown the paper documents for DA.
- `AsymptoticDiscreteState` (choice B, Sec. VI B) -- the R-independent bound
  state `phi_b`, the `R -> infinity` limit of choice A. Near-exact in the
  paper's tests.

Both satisfy Eq. (67) (`phi_d -> phi_b` as `R -> infinity`) and both return
DVR COEFFICIENTS, c-product-normalized to 1 and localized in the real
electronic region. Because a DVR coefficient is `d_j = phi(r_j) sqrt(w_j)`,
Eq. (58)'s projector is `outer(d, d)` with no further weights, and the
c-product normalization already supplies the paper's `exp[-i delta(R)]`
phase-fixing (p. 012710-8).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Protocol

import numpy as np
import numpy.typing as npt

from qscat.core.dissociation import anion_electronic_states
from qscat.core.lcp import resonance_pole_walk
from qscat.dvr import FemDvrEcsGrid, eigen, kinetic
from qscat.exceptions import ConvergenceError
from qscat.linalg import c_product

from .scattering import scattering_state

if TYPE_CHECKING:
    from qscat.model import ResonanceModel

__all__ = [
    "AsymptoticDiscreteState",
    "DiscreteState",
    "PhysicalDiscreteState",
    "electronic_hamiltonian",
    "truncate",
]

# Below this c-product norm a vector is self-orthogonal to numerical
# precision and cannot be c-normalized; a real discrete state never is.
_MIN_NORM2 = 1e-12

# `PhysicalDiscreteState`'s bound branch: an eigenvalue of `H_el(R)` counts as
# "genuinely bound" if `Re(E) < 0` and `|Im(E)|` is below this bound. Looser
# than `dissociation.anion_electronic_states`' `_IM_TOL_HA = 1e-6` (calibrated
# for a state resolved at an asymptotic, well-separated R) because this branch
# is evaluated at ARBITRARY R along a continuation walk, including right at
# the bound/resonant crossing, where a state can carry a small but genuine
# width even while `Re(E) < 0` -- a true near-threshold Feshbach state, not
# numerical noise. F2's crossing sampled on the coarse test fixture (R ~ 2.6
# bohr) has such a state with `|Im(E)| ~ 5.6e-6`, three orders of magnitude
# below a genuine resonance's width (~1e-3) there. The gate exists to catch the
# OTHER failure -- no negative-Re eigenvalue at all, i.e. the walk's bound/
# resonant sign decision was simply wrong for this R -- not to demand
# asymptotic-grade purity.
#
# SIZED ON THE PRODUCTION DECKS, not on that fixture -- WHY 3e-5 AND NOT A
# NUMBER SOMEONE LIKED. eMoScat's F2 deck has an ELEMENT BOUNDARY at
# R = 2.596908 (`validation.diatomic.config`'s `nuc_real`, hand-placed on the
# crossing), so it puts a DVR node 0.0005 bohr from it -- far closer than the
# coarse fixture above. NO's deck has no such boundary but its uniform
# 37-element segment still lands a node ~0.005 bohr from its own crossing.
# The near-threshold state at those nodes carries a width just over 1e-5.
# Measured spectra of `electronic_hamiltonian`, four lowest by Re(E):
#
#   F2, R = 2.5974:  -1.432236e-4 - 1.059e-5j   <- the state, |Im|/|Re| = 0.074
#                    +4.527175e-3 - 4.435e-3j
#                    +1.338226e-2 - 1.310e-2j
#                    +2.666222e-2 - 2.609e-2j
#   NO, R = 2.2884:  -6.288891e-4 - 1.272e-5j   <- the state, |Im|/|Re| = 0.020
#                    +4.559187e-3 - 4.551e-3j
#
# The neighbouring F2 nodes bracket the crossing and behave: R = 2.60 gives
# -9.300079e-4 - 5.585e-6j (bound, passes even at 1e-5) and R = 2.59 gives
# +2.069295e-3 - 5.458e-5j (positive shift, scattering branch, never reaches
# this guard). Only the node ON the crossing trips it.
#
# In both molecules the classification is unambiguous: the next eigenvalue up
# sits at Re(E) ~ +4.5e-3 with |Im(E)| ~ 4.4e-3, i.e. 7-32x further out in
# |Re(E)| with a width 358-419x larger. So at 1e-5 the guard was rejecting the
# very state it exists to accept, and `PhysicalDiscreteState` could not be
# built on either production deck at all. 3e-5 keeps a 2.4x margin over both
# measurements while staying two orders below a genuine resonance width there.
#
# KNOWN DESIGN LIMITATION: an ABSOLUTE |Im(E)| bound is intrinsically arbitrary
# near a crossing, where the state passes continuously from bound to resonant
# and both parts shrink together -- which is why it needed re-sizing per deck
# at all. A scale-free criterion (`|Im(E)| < c |Re(E)|`, satisfied by both
# points above at c = 0.074 and 0.020) would not. Not attempted here; see
# docs/physics/nonlocal-resonance-model.md.
_BOUND_IM_TOL = 3e-5


class DiscreteState(Protocol):
    """Anything that supplies `phi_d(r;R)` as c-normalized DVR coefficients."""

    def phi_d(self, R: float) -> npt.NDArray[np.complex128]:
        """The discrete state at nuclear coordinate `R`."""
        ...


def electronic_hamiltonian(
    grid: FemDvrEcsGrid, model: ResonanceModel, R: float
) -> npt.NDArray[np.complex128]:
    """`H_el(r;R) = T_r + ell(ell+1)/2r^2 + V_int(R,r)` (Eq. 17), mass 1.

    Note this EXCLUDES `V_0(R)`: `model.surface` includes it, and the paper's
    `H_el` does not (`V_0` re-enters in Eq. 20's `V_d` and Eq. 61's `M(n)`).
    """
    v = model.surface(grid.points, R) - model.v0(R)
    out: npt.NDArray[np.complex128] = kinetic(grid, 1.0) + np.diag(v)
    return out


def truncate(
    coeffs: npt.NDArray[np.complex128], grid: FemDvrEcsGrid, r_d: float = 10.0
) -> npt.NDArray[np.complex128]:
    """Apply Eq. (69)'s smooth cutoff and zero the ECS tail.

    `f(r) = 1 - 1/(1 + e^{-(r-r_d)}) = 1/(1 + e^{(r-r_d)})`, the logistic
    written in the numerically stable form. `r_d = 10 a_0` is the paper's own
    value (p. 012710-8). The complex tail is zeroed outright: the paper's
    discrete state is "real and localized in the inner region where the
    electronic coordinate is not complex scaled" (p. 012710-6).
    """
    real = grid.points.imag == 0.0
    f = np.zeros(grid.n, dtype=np.complex128)
    r = grid.points[real].real
    # expit(-(r - r_d)) == 1/(1 + exp(r - r_d)), overflow-free.
    f[real] = 1.0 / (1.0 + np.exp(np.clip(r - r_d, -700.0, 700.0)))
    out: npt.NDArray[np.complex128] = coeffs * f
    return out


def _c_normalize(d: npt.NDArray[np.complex128]) -> npt.NDArray[np.complex128]:
    """Scale to `c_product(d, d) == 1`, which also fixes the phase (real, +)."""
    norm2 = c_product(d, d)
    if abs(norm2) < _MIN_NORM2:
        raise ConvergenceError(
            f"discrete state is c-product self-orthogonal (norm^2={norm2!r}); "
            "cannot normalize -- the state is not a usable phi_d"
        )
    out: npt.NDArray[np.complex128] = d / np.sqrt(norm2)
    return out


class AsymptoticDiscreteState:
    """Choice B (Sec. VI B): the R-independent bound state `phi_b`.

    One electronic eigenproblem for the whole calculation. Eq. (67) holds
    trivially. The paper finds this choice essentially exact for DA and for VE
    provided the background T-matrix terms are included (p. 012710-9-10).
    """

    def __init__(self, grid: FemDvrEcsGrid, model: ResonanceModel, R_inf: float) -> None:
        _eps, phi = anion_electronic_states(grid, model, R_inf, n_states=1)
        # Eq. (69)'s cutoff applies only to the SCATTERING-derived state
        # (p. 012710-8): a bound state is already square-integrable, and
        # truncating it breaks the exact eigenrelation `H_el phi_b = E_b
        # phi_b` that Eq. (67)'s decoupling (`V_dk -> 0` by P-space
        # orthogonality) relies on. Use the raw eigenvector.
        d = np.asarray(phi[0], dtype=np.complex128)
        self._d = _c_normalize(d)

    def phi_d(self, R: float) -> npt.NDArray[np.complex128]:
        """The same state at every `R`, by construction."""
        return self._d


class PhysicalDiscreteState:
    """Choice A (Sec. VI A): the scattering function at `Re E_res(R)`.

    At each `R` the resonance energy is taken from the two-angle ECS pole walk
    (`qscat.core.lcp.resonance_pole_walk`, which returns the ELECTRONIC shift
    `s(R) = Re E_pole(R) - v0(R)`). Where `s(R) > 0` the state is the
    fixed-nuclei scattering function at that energy (Eq. 17's `H_el`); where
    `s(R) <= 0` the electron is bound and the state is the lowest GENUINELY
    bound eigenvector of `H_el(R)` instead (p. 012710-8): `Re(E) < 0` and
    `|Im(E)| < _BOUND_IM_TOL`. If no such eigenvalue exists, the walk's sign
    decision does not match the fresh spectrum at this `R` (most likely a
    frozen/stale continuation point from a `R_descending` too coarse for
    `resonance_pole_walk` to track -- see `re_half_width`/`im_half_width`
    below) and `ConvergenceError` is raised rather than silently returning
    whatever eigenvector happens to have the smallest real part -- that
    eigenvector is not necessarily real or bound. Only the SCATTERING branch
    is then truncated by Eq. (69) (the paper's cutoff applies to the
    scattering-derived state only, p. 012710-8); the bound branch is already
    square-integrable and is used as-is, so its exact eigenrelation
    `H_el phi_b = E_b phi_b` survives -- truncating it would spuriously
    prevent `V_dk -> 0` at large `R` (Eq. 67). Both branches are then
    c-normalized, which supplies the paper's `exp[-i delta(R)]`.

    `R_descending` must be descending -- the pole walk is seeded at large `R`
    and continued inward. States are precomputed on those nodes; `phi_d(R)`
    returns the nearest one, which is how the ingredient builder consumes it
    (one call per nuclear DVR point, on the same node set). The per-node
    electronic shift and which branch was taken are recorded as `shift`,
    `gamma`, and `used_scattering` (all aligned with `R_descending`), so a
    caller can inspect what happened rather than infer it from the output
    state alone.

    The pole walk's seed window is centered on the bound ANION electronic
    state at the outermost (largest) `R` in `R_descending`
    (`qscat.core.dissociation.anion_electronic_states`), mirroring
    `qscat.core.lcp._walk_from_anion_seed`: at that `R` the discrete state is
    expected to already be essentially the asymptotic bound state (Eq. 67), so
    its energy is a reliable seed for the pole finder. `seed_window` overrides
    this and is passed straight to `resonance_pole_walk` if given.

    `re_half_width`/`im_half_width` default wider (`0.08` Ha) than
    `qscat.core.lcp`'s own default (`0.05` Ha, tuned for a FINE nuclear grid's
    closely-spaced real nodes): `R_descending` here is an arbitrary,
    caller-supplied array and may be sampled far more coarsely (e.g. eMoScat's
    hand-picked nuclear decks step by 0.4 bohr near F2's bound/resonant
    crossing around R ~ 2.6-3.0 bohr, where the electronic shift moves by
    ~0.09 Ha over that step) -- too coarse for the 0.05 Ha window to
    recenter onto without losing the pole, which freezes the walk exactly
    where the physics is most interesting (the crossing itself).
    """

    def __init__(
        self,
        grid: FemDvrEcsGrid,
        model: ResonanceModel,
        R_descending: npt.NDArray[np.float64],
        elec_grid_b: FemDvrEcsGrid,
        *,
        r_d: float = 10.0,
        re_half_width: float = 0.08,
        im_half_width: float = 0.08,
        resid_tol: float = 1e-3,
        seed_window: tuple[float, float, float, float] | None = None,
    ) -> None:
        R = np.asarray(R_descending, dtype=np.float64)
        if R.size > 1 and np.any(np.diff(R) >= 0.0):
            raise ValueError("R_descending must be strictly descending")
        if seed_window is None:
            eps_e, _phi = anion_electronic_states(grid, model, float(R[0]), n_states=1)
            seed_window = (
                eps_e[0] - re_half_width,
                eps_e[0] + re_half_width,
                -im_half_width,
                im_half_width,
            )
        shift, gamma = resonance_pole_walk(
            model,
            R,
            grid,
            elec_grid_b,
            seed_window,
            re_half_width=re_half_width,
            im_half_width=im_half_width,
            resid_tol=resid_tol,
        )
        self._R = R
        self.shift = shift
        self.gamma = gamma
        self.used_scattering = np.zeros(R.size, dtype=np.bool_)
        self._states = np.empty((R.size, grid.n), dtype=np.complex128)
        for j in range(R.size):
            h_el = electronic_hamiltonian(grid, model, float(R[j]))
            e_res = float(shift[j])
            if e_res > 0.0:
                raw = scattering_state(h_el, grid, e_res, model.ell)
                self.used_scattering[j] = True
                raw = truncate(raw, grid, r_d)
            else:
                evals, evecs = eigen(h_el)  # ascending Re(E)
                bound = np.flatnonzero((np.abs(evals.imag) < _BOUND_IM_TOL) & (evals.real < 0.0))
                if bound.size == 0:
                    raise ConvergenceError(
                        f"PhysicalDiscreteState: resonance_pole_walk reports R="
                        f"{float(R[j]):.4f} bohr as bound (electronic shift "
                        f"{e_res:.6g} Ha <= 0), but H_el(R) has no eigenvalue with "
                        f"Re(E) < 0 and |Im(E)| < {_BOUND_IM_TOL:.0e} Ha -- the "
                        "walk's sign decision does not match the fresh spectrum "
                        "here. Likely R_descending is too coarse for "
                        "resonance_pole_walk to track the pole across a step "
                        "(try a denser sampling or a wider re_half_width/"
                        "im_half_width)."
                    )
                raw = evecs[:, int(bound[0])].astype(np.complex128)
                # Bound branch: no Eq. (69) truncation (p. 012710-8) -- see
                # AsymptoticDiscreteState's docstring note above.
            self._states[j] = _c_normalize(raw)

    def phi_d(self, R: float) -> npt.NDArray[np.complex128]:
        """The discrete state at the precomputed node nearest `R`."""
        j = int(np.argmin(np.abs(self._R - R)))
        out: npt.NDArray[np.complex128] = self._states[j]
        return out
