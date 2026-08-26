"""Born-Oppenheimer reference states: the picture an exact 2-D pole departs from.

`qscat.core.exact_resonance_states` returns approximation-free resonance poles.
On its own a pole is a number in the complex plane -- it carries no label, no
quantum numbers, and no evidence that it is a resonance at all. This module
builds the reference the pole is compared against: the Born-Oppenheimer product
states

    Psi_BO(r, R) = phi_j(r; R) chi_v(R)

where `phi_j` is an electronic eigenstate at frozen `R` and `chi_v` a
vibrational level of that curve. `qscat.core.assignment` then pairs poles to
these states by overlap, which is what supplies both the label and the
evidence.

## Two kinds of electronic curve

The `phi_j` differ by target, and this module builds both:

- **Bound curves** (`electronic_curves`). For an ION the electron is bound in
  the residual cation's field, giving the Rydberg series
  `E_Ry0(R) < E_Ry1(R) < ...`; a dissociative-recombination peak is
  conventionally assigned to a vibrational level of one of these curves. For a
  neutral the same function gives the bound electronic curves where they exist.
- **A resonance curve** (`resonance_curve`). For a NEUTRAL target the state of
  interest is not bound at all -- it is the anion resonance `V_d(R) - i
  Gamma(R)/2` that `qscat.core.lcp` reduces to a local complex potential. Its
  electronic eigenfunction comes from the same two-angle pole walk, and the
  product states built on it are what the N2/NO/F2 exact poles must be checked
  against.

Both return the same `ElectronicCurves`, so `bo_basis` does not care which one
it was handed. That is the seam: one comparator, two builders.

## Phase alignment is not cosmetic

An electronic eigenvector's phase is arbitrary at every nuclear point
independently. Without alignment the product `phi_j(r; R) chi_v(R)` flips sign
at random `R`, and the overlap this basis exists to measure integrates those
flips to near zero against any smooth partner -- a genuine state then scores
like an artefact. Every curve here is phase-aligned across `R` by continuity
(each column rotated to make its overlap with the previous column real and
positive) and callers should not re-normalize in a way that destroys it.

## Levels are not a rectangle

Curves do not share one vibrational capacity: a deep curve may support five
clean bound levels where a shallow Rydberg one supports twelve. `n_vib` is
therefore a per-curve REQUEST. `allow_partial=True` asks each curve for as many
as it can supply and pads the rest of the row with `NaN`; the default keeps the
strict behaviour so a caller assuming a full table finds out when it is wrong.

## The basis has to be deep enough or the overlap test lies

A low overlap means "no partner in THIS basis", never "no partner". Measured on
H2+: with curves only to `Ry_11`, eight genuine `Ry_12..Ry_16` states scored
0.02-0.09 and looked spurious. `n_eff`, `admissible_levels` and `basis_covers`
exist to make that distinction computable rather than a judgement call -- see
`admissible_levels` for the closed-channel argument they rest on.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

from qscat.dvr import FemDvrEcsGrid, eigen, kinetic
from qscat.ecs import find_resonance_pole
from qscat.exceptions import ConvergenceError, GridError
from qscat.linalg import c_product

from .vibrational import vibrational_states

if TYPE_CHECKING:  # pragma: no cover - typing only
    from qscat.model import ResonanceModel

__all__ = [
    "BoBasis",
    "BoState",
    "ElectronicCurves",
    "admissible_levels",
    "basis_covers",
    "bo_basis",
    "bo_basis_from_levels",
    "electronic_curves",
    "n_eff",
    "resonance_curve",
]


# Above this, a curve's imaginary part in the REAL nuclear region is a width
# rather than round-off, and its levels are quasi-bound rather than bound.
_REAL_CURVE_TOL = 1e-10


def _align_phase(
    vec: npt.NDArray[np.complex128], prev: npt.NDArray[np.complex128] | None
) -> npt.NDArray[np.complex128]:
    """Rotate `vec` so its overlap with `prev` is real and positive."""
    if prev is None:
        return vec
    ov = complex(np.vdot(prev, vec))
    if ov == 0:
        return vec
    return np.asarray(vec * (abs(ov) / ov), dtype=np.complex128)


@dataclass(frozen=True)
class ElectronicCurves:
    """Electronic eigen-curves tabulated on a nuclear grid.

    - `energies`: `(n_curves, n_R)` complex. `energies[j, k]` is curve `j` at
      nuclear point `g_R.points[k]`. Complex on the nuclear ECS tail even for a
      genuinely bound curve -- feeding a complex `R` into `model.surface` makes
      the whole electronic problem complex, which is a curve-parametrization
      artifact of the nuclear rotation and NOT a resonance width.
    - `states`: `(n_curves, n_r, n_R)` complex, phase-aligned across `R`, or an
      empty `(0, 0, 0)` array when built with `with_states=False`. Each column
      carries `qscat.dvr.eigen`'s Hermitian normalization; `bo_basis`
      c-normalizes the finished product.
    """

    energies: npt.NDArray[np.complex128]
    states: npt.NDArray[np.complex128]

    @property
    def n_curves(self) -> int:
        """Number of electronic curves (rows of `energies`)."""
        return int(self.energies.shape[0])

    @property
    def has_states(self) -> bool:
        """True when the curves were built `with_states=True`."""
        return bool(self.states.size)


def electronic_curves(
    model: ResonanceModel,
    g_r: FemDvrEcsGrid,
    g_R: FemDvrEcsGrid,
    *,
    n_curves: int,
    with_states: bool = False,
) -> ElectronicCurves:
    """The `n_curves` lowest electronic eigen-curves over `g_R`.

    At every nuclear grid point `R` (including the complex ECS-tail points, not
    just the real region) this diagonalizes the frozen-nucleus electronic
    problem `-1/2 d^2/dr^2 + model.surface(r, R)` (electron mass 1) on `g_r` and
    keeps the `n_curves` lowest-`Re E` eigenvalues. For an ion those are the
    Rydberg series; for a neutral, the bound electronic curves.

    This reproduces `qscat.core.anion_electronic_states` wherever that
    function's "genuinely bound" gate succeeds, and deliberately does NOT reuse
    that gate, because the gate fails on the ECS tail: a complex `R_inf` makes
    every eigenvalue there pick up an `O(Im v0(R))` imaginary shift, already
    past the library's `1e-6` bound-state tolerance a hair's width into the
    tail. Gating on it would raise on perfectly good curve points.

    Parameters
    ----------
    model : ResonanceModel
        Supplies `surface(r, R)`.
    g_r, g_R : FemDvrEcsGrid
        Electronic and nuclear radial grids.
    n_curves : int
        How many lowest-`Re E` curves to keep.
    with_states : bool, optional
        Also return the eigenvectors, phase-aligned across `R`. Costs
        `n_curves * g_r.n * g_R.n` complex128 -- ~100 MB on a production H2+
        deck -- so it is opt-in. Required for `bo_basis` to build product
        states; the energies alone are enough for a level table.

    Returns
    -------
    ElectronicCurves

    Raises
    ------
    GridError
        If `n_curves` exceeds the electronic grid's dimension.
    """
    if n_curves < 1:
        raise GridError(f"n_curves must be >= 1, got {n_curves}")
    if n_curves > g_r.n:
        raise GridError(f"n_curves={n_curves} exceeds the electronic grid dimension {g_r.n}")

    pts = g_R.points
    energies = np.empty((n_curves, pts.size), dtype=np.complex128)
    states = (
        np.empty((n_curves, g_r.n, pts.size), dtype=np.complex128)
        if with_states
        else np.empty((0, 0, 0), dtype=np.complex128)
    )
    prev: list[npt.NDArray[np.complex128] | None] = [None] * n_curves

    T = kinetic(g_r, 1.0)
    for k, R in enumerate(pts):
        vals, vecs = eigen(T + np.diag(model.surface(g_r.points, complex(R))))
        energies[:, k] = vals[:n_curves]
        if not with_states:
            continue
        for j in range(n_curves):
            vec = _align_phase(np.asarray(vecs[:, j], dtype=np.complex128), prev[j])
            prev[j] = vec
            states[j, :, k] = vec
    return ElectronicCurves(energies=energies, states=states)


def resonance_curve(
    model: ResonanceModel,
    g_r_a: FemDvrEcsGrid,
    g_r_b: FemDvrEcsGrid,
    g_R: FemDvrEcsGrid,
    seed_window: tuple[float, float, float, float],
    *,
    re_half_width: float = 0.05,
    im_half_width: float = 0.05,
    resid_tol: float = 1e-3,
    with_states: bool = True,
) -> ElectronicCurves:
    """The anion RESONANCE electronic curve over `g_R`, as a one-curve family.

    The neutral-target counterpart of `electronic_curves`: the state a
    dissociative-attachment resonance is built on is not bound, so it cannot be
    picked out by "lowest `Re E`". It is found instead as the angle-stable pole
    of the fixed-`R` electronic problem at two ECS angles
    (`qscat.ecs.find_resonance_pole`), continued inward over descending real `R`
    exactly as `qscat.core.lcp.resonance_pole_walk` does -- and this function
    keeps the eigenvector that walk discards.

    The grid-sizing sibling `qscat.tuning.resonance.resonance_curve_arrays`
    (renamed from `resonance_curve` in the 2026-08-25 API surface pass) runs
    the same underlying pole walk but returns plain `(R, V_d, Gamma)` arrays
    with no states -- use that one to size a grid, this one to build BO
    basis states.

    `energies[0]` is `V_d(R) - i Gamma(R)/2`, so
    `qscat.core.lcp.local_complex_potential`'s curve and this one are the same
    object viewed with and without its eigenfunctions.

    **The ECS tail is frozen, not walked.** `find_resonance_pole` needs a real
    `R`. Tail points therefore take the outermost real node's eigenvector and
    that node's electronic shift added to `v0(z)` -- the same analytic
    continuation `qscat.core.lcp._assemble_lcp` applies to the potential, for
    the same reason. A product state's tail behaviour is carried by `chi_v`,
    not by `phi_res`, so this is where the approximation belongs.

    **On breakdown the last accepted state is frozen inward**, matching
    `resonance_pole_walk`'s behaviour rather than raising: a pole finder that
    loses the track at small `R` has usually run into a region where the
    resonance has merged with the continuum, and the frozen value is the honest
    continuation. `ConvergenceError` is raised only when the finder fails at the
    seed edge, where there is nothing to freeze.

    Parameters
    ----------
    model : ResonanceModel
    g_r_a, g_r_b : FemDvrEcsGrid
        Electronic grids differing ONLY in ECS tail angle. `g_r_a` is the one
        the returned eigenvectors live on.
    g_R : FemDvrEcsGrid
        Nuclear grid the curve is tabulated on.
    seed_window : tuple of float
        `(re_lo, re_hi, im_lo, im_hi)` for the pole at the OUTERMOST real `R`.
        `qscat.core.anion_electronic_states` at `g_R.R0` is the natural source.
    re_half_width, im_half_width : float, optional
        Half-widths of the window recentred on each accepted pole.
    resid_tol : float, optional
        Angle-stability residual above which the walk is considered broken.
    with_states : bool, optional
        Keep the eigenvectors (default `True` -- unlike `electronic_curves`,
        since a one-curve family is cheap and the states are the point).

    Returns
    -------
    ElectronicCurves
        `n_curves == 1`.

    Raises
    ------
    ConvergenceError
        If no pole is accepted at the seed edge.
    """
    pts = g_R.points
    real_idx = np.flatnonzero(pts.imag == 0.0)
    walk = real_idx[np.argsort(pts[real_idx].real)[::-1]]  # descending R: outer -> inner
    tail = np.flatnonzero(pts.imag != 0.0)

    energies = np.empty((1, pts.size), dtype=np.complex128)
    states = (
        np.empty((1, g_r_a.n, pts.size), dtype=np.complex128)
        if with_states
        else np.empty((0, 0, 0), dtype=np.complex128)
    )

    T = kinetic(g_r_a, 1.0)
    T_b = kinetic(g_r_b, 1.0)
    real_mask = g_r_a.real_points <= g_r_a.R0

    window = seed_window
    last: tuple[complex, float, npt.NDArray[np.complex128]] | None = None  # (E, shift, phi)
    prev: npt.NDArray[np.complex128] | None = None
    broken = False

    for idx in walk:
        R = float(pts[idx].real)
        v0R = float(np.real(model.v0(np.asarray(R))))
        if not broken:
            try:
                E_a, V_a = eigen(T + np.diag(model.surface(g_r_a.points, R)))
                E_b, _ = eigen(T_b + np.diag(model.surface(g_r_b.points, R)))
                E_pole, resid = find_resonance_pole(E_a, E_b, window)
            except (ValueError, np.linalg.LinAlgError):
                resid = np.inf
            else:
                if resid < resid_tol:
                    col = int(np.argmin(np.abs(E_a - E_pole)))
                    phi = np.asarray(V_a[:, col], dtype=np.complex128)
                    p = phi.copy()
                    p[~real_mask] = 0.0
                    nrm = complex(c_product(p, p))
                    if nrm != 0:
                        phi = phi / np.sqrt(nrm)
                    phi = _align_phase(phi, prev)
                    prev = phi
                    last = (complex(E_pole), E_pole.real - v0R, phi)
                    window = (
                        E_pole.real - re_half_width,
                        E_pole.real + re_half_width,
                        E_pole.imag - im_half_width,
                        E_pole.imag + im_half_width,
                    )
                    energies[0, idx] = E_pole
                    if with_states:
                        states[0, :, idx] = phi
                    continue
            broken = True
        if last is None:
            raise ConvergenceError(
                "resonance_curve: the pole finder failed at the seed edge "
                f"(R = {R:.4f}); widen seed_window or move the outer nuclear node"
            )
        energies[0, idx] = v0R + last[1] + 1j * last[0].imag
        if with_states:
            states[0, :, idx] = last[2]

    if tail.size:
        assert last is not None
        # Analytic continuation, exactly as `lcp._assemble_lcp` does for V_d:
        # the ASYMPTOTIC electronic shift laid on v0(z), and the outermost real
        # eigenvector frozen. Recorded at walk[0] -- the outermost real node.
        outer = walk[0]
        shift_inf = float(energies[0, outer].real) - float(
            np.real(model.v0(np.asarray(float(pts[outer].real))))
        )
        energies[0, tail] = model.v0(pts[tail]) + shift_inf
        if with_states:
            states[0][:, tail] = states[0][:, outer][:, None]

    return ElectronicCurves(energies=energies, states=states)


@dataclass(frozen=True)
class BoState:
    """One Born-Oppenheimer product state, with the identity that labels it."""

    psi: npt.NDArray[np.complex128]  # flat (C order over (r, R)), c-normalized
    # The vibrational level's own energy (Hartree). REAL for a bound electronic
    # curve; genuinely complex (`E_v - i Gamma_v/2`) for a resonance curve,
    # where the level is quasi-bound and its width is part of the answer.
    energy: complex
    curve: int  # j -- the electronic curve index (the SUPERscript)
    vib: int  # v -- the vibrational quantum number (the SUBscript)

    @property
    def label(self) -> str:
        r"""`omega^j_v` as a LaTeX string, the published labelling convention."""
        return rf"$\omega^{{{self.curve}}}_{{{self.vib}}}$"


@dataclass(frozen=True)
class BoBasis:
    """A `(curve, vib)`-keyed family of BO product states plus its level table.

    `energies` is the rectangular level table (`NaN` where a curve supports
    fewer than `n_vib` levels); `states` is empty when the basis was built from
    curves carrying no eigenvectors, which is the level-table-only path.
    """

    energies: npt.NDArray[np.complex128]  # (n_curves, n_vib), NaN-padded
    states: dict[tuple[int, int], BoState]

    def __contains__(self, key: object) -> bool:
        """`(curve, v) in basis` membership test."""
        return key in self.states

    def __getitem__(self, key: tuple[int, int]) -> BoState:
        """The `BoState` stored under `(curve, v)`."""
        return self.states[key]

    def __len__(self) -> int:
        """Number of BO product states in the basis."""
        return len(self.states)

    def items(self) -> list[tuple[tuple[int, int], BoState]]:
        """All `((curve, v), BoState)` pairs, dict-style."""
        return list(self.states.items())

    @property
    def has_states(self) -> bool:
        """True when the basis holds at least one state."""
        return bool(self.states)

    def levels(self) -> list[tuple[int, int]]:
        """Every `(curve, vib)` whose level energy is finite, curve-major."""
        j_idx, v_idx = np.nonzero(np.isfinite(self.energies))
        return [(int(j), int(v)) for j, v in zip(j_idx, v_idx, strict=True)]

    def flat(self) -> tuple[npt.NDArray[np.complex128], list[tuple[int, int]]]:
        """`(energy, key)` for the finite levels, ascending in `Re E`."""
        keys = self.levels()
        e = np.array([self.energies[j, v] for j, v in keys], dtype=np.complex128)
        order = np.argsort(e.real)
        return e[order], [keys[i] for i in order]


def bo_basis(
    curves: ElectronicCurves,
    g_R: FemDvrEcsGrid,
    mu: float,
    *,
    n_vib: int,
    allow_partial: bool = False,
) -> BoBasis:
    """Vibrational ladders in `curves`, and the product states they define.

    Each curve is fed to `qscat.core.vibrational_states` as the nuclear
    potential. That solver builds `T_nuc(mu) + diag(v0(grid.points))`, calling
    `v0` with EXACTLY `grid.points`, so a closure returning the curve already
    tabulated on those same points is an exact lookup rather than an
    interpolation.

    When `curves` carries eigenvectors, the product `phi_j(r; R) chi_v(R)` is
    formed and c-product-normalized (the bilinear, non-conjugated ECS pairing --
    a conjugated norm would weight the exponentially growing ECS tail instead of
    cancelling it). Without eigenvectors only the level table is produced, which
    is the cheap path when no overlap is going to be taken.

    Parameters
    ----------
    curves : ElectronicCurves
        From `electronic_curves` or `resonance_curve`. Its `energies` must be
        tabulated on `g_R`.
    g_R : FemDvrEcsGrid
        The nuclear grid the curves live on.
    mu : float
        Nuclear reduced mass (atomic units) -- `model.mu`.
    n_vib : int
        Vibrational levels REQUESTED per curve; see `allow_partial`.
    allow_partial : bool, optional
        Ask each curve for as many of `n_vib` as it supports, padding the row
        with `NaN`, instead of raising when one cannot supply all of them.

    Returns
    -------
    BoBasis

    Raises
    ------
    GridError
        If `curves.energies` is not tabulated on `g_R`, or `n_vib < 1`.
    """
    if n_vib < 1:
        raise GridError(f"n_vib must be >= 1, got {n_vib}")
    if curves.energies.shape[1] != g_R.n:
        raise GridError(
            f"curves are tabulated on {curves.energies.shape[1]} nuclear points "
            f"but g_R has {g_R.n} -- they must be the same grid"
        )

    # A RESONANCE curve carries a width, so `T + diag(curve)` has genuinely
    # complex eigenvalues and `vibrational_states`' bound-state gate rejects
    # every one of them -- correctly, since they are not bound levels. Say so
    # here rather than letting that gate's message ("requested more states than
    # there are bound levels") send a caller looking for a grid problem.
    real_nodes = g_R.points.imag == 0.0
    if np.max(np.abs(curves.energies[:, real_nodes].imag)) > _REAL_CURVE_TOL:
        raise GridError(
            "bo_basis: this curve is complex in the real nuclear region, i.e. it "
            "carries a width, so its levels are quasi-bound rather than bound. "
            "Build them with qscat.core.lcp.resonance_levels and combine via "
            "bo_basis_from_levels instead."
        )

    n_curves = curves.n_curves
    energies = np.full((n_curves, n_vib), np.nan, dtype=np.complex128)
    states: dict[tuple[int, int], BoState] = {}

    for j in range(n_curves):

        def v_n(
            _R: npt.ArrayLike, _c: npt.NDArray[np.complex128] = curves.energies[j]
        ) -> npt.NDArray[np.complex128]:
            # Exact lookup on grid.points, not interpolation -- see docstring.
            return np.asarray(_c, dtype=np.complex128)

        basis = None
        if allow_partial:
            # Walk down to the deepest count this curve supports.
            # `vibrational_states` selects the `n` lowest-Re(E) eigenvalues and
            # rejects the batch if ANY is quasi-continuum, so a smaller `n` is a
            # strictly cleaner subset: there is no `n` that raises where `n-1`
            # contains a level `n` did not.
            for count in range(n_vib, 0, -1):
                try:
                    basis = vibrational_states(g_R, mu, count, v_n)
                except GridError:
                    continue
                break
            if basis is None:
                continue  # this curve supports no clean bound level at all
        else:
            basis = vibrational_states(g_R, mu, n_vib, v_n)

        eps = np.asarray(basis.eps, dtype=np.float64)
        energies[j, : eps.size] = eps.astype(np.complex128)
        if not curves.has_states:
            continue
        for v in range(eps.size):
            psi = (curves.states[j] * basis.chi[v][None, :]).ravel()
            nrm = complex(c_product(psi, psi))
            if nrm == 0:
                continue
            states[(j, v)] = BoState(
                psi=np.asarray(psi / np.sqrt(nrm), dtype=np.complex128),
                energy=complex(eps[v]),
                curve=j,
                vib=v,
            )
    return BoBasis(energies=energies, states=states)


def bo_basis_from_levels(
    curves: ElectronicCurves,
    level_energies: npt.ArrayLike,
    level_states: npt.NDArray[np.complex128],
    *,
    curve: int = 0,
) -> BoBasis:
    """Product states from an ALREADY-SOLVED nuclear problem -- the neutral path.

    For a resonance curve the nuclear levels are quasi-bound, and
    `qscat.core.lcp.resonance_levels` already solves for them properly: it
    diagonalizes `T(mu) + diag(V_d - i Gamma/2)` on two nuclear grids and keeps
    what is angle-stable, which is a strictly harder problem than
    `vibrational_states`' bound-state solve. Rather than reimplement it, this
    takes its output and forms the products::

        levels = resonance_levels(model, nu_a, nu_b, el_a, el_b, ...)
        cur = resonance_curve(model, el_a, el_b, nu_a, seed_window)
        basis = bo_basis_from_levels(cur, levels.energies, levels.states)

    Parameters
    ----------
    curves : ElectronicCurves
        Must carry states. Only index `curve` is used.
    level_energies : array_like of complex, shape (n_levels,)
        `E_v - i Gamma_v / 2`.
    level_states : ndarray of complex, shape (n_levels, n_R)
        Nuclear eigenvectors, one row per level, on the SAME nuclear grid the
        curve is tabulated on. `ResonanceLevels.states` has exactly this shape
        and normalization.
    curve : int, optional
        Which electronic curve to multiply in (default 0, the only one a
        `resonance_curve` family has).

    Returns
    -------
    BoBasis
        Keyed `(curve, v)`, with `energies` a `(1, n_levels)` row.

    Raises
    ------
    GridError
        On a curve/level grid mismatch, or a missing curve index or states.
    """
    if not curves.has_states:
        raise GridError("bo_basis_from_levels needs curves built with with_states=True")
    if not 0 <= curve < curves.n_curves:
        raise GridError(f"curve={curve} is outside the {curves.n_curves} curves given")
    chi = np.atleast_2d(np.asarray(level_states, dtype=np.complex128))
    eps = np.atleast_1d(np.asarray(level_energies, dtype=np.complex128))
    if chi.shape[0] != eps.size:
        raise GridError(f"{eps.size} level energies but {chi.shape[0]} level states")
    if chi.shape[1] != curves.energies.shape[1]:
        raise GridError(
            f"levels live on {chi.shape[1]} nuclear points but the curve is "
            f"tabulated on {curves.energies.shape[1]} -- they must be the same grid"
        )

    energies = np.full((1, eps.size), np.nan, dtype=np.complex128)
    states: dict[tuple[int, int], BoState] = {}
    phi = curves.states[curve]
    for v in range(eps.size):
        energies[0, v] = eps[v]
        psi = (phi * chi[v][None, :]).ravel()
        nrm = complex(c_product(psi, psi))
        if nrm == 0:
            continue
        states[(curve, v)] = BoState(
            psi=np.asarray(psi / np.sqrt(nrm), dtype=np.complex128),
            energy=complex(eps[v]),
            curve=curve,
            vib=v,
        )
    return BoBasis(energies=energies, states=states)


def n_eff(e_tot: float, thresholds: npt.ArrayLike) -> float:
    """Hydrogenic effective quantum number of a state at `e_tot`.

    `n_eff = 1/sqrt(2 * binding)`, with `binding` measured to the nearest
    threshold ABOVE `e_tot` -- NOT to the lowest one. The per-level threshold is
    the physically relevant one: a Rydberg state is bound against the ion level
    it sits below, not against the incident channel.

    Raises `ValueError` if `e_tot` lies above every threshold, where the
    quantity is undefined rather than large.
    """
    thr = np.sort(np.asarray(thresholds, dtype=np.float64))
    above = thr[thr > e_tot]
    if above.size == 0:
        raise ValueError(f"e_tot={e_tot} sits above every threshold given")
    return float(1.0 / np.sqrt(2.0 * float(above[0] - e_tot)))


def admissible_levels(
    e_tot: float, thresholds: npt.ArrayLike, *, n_eff_max: float | None = None
) -> list[tuple[int, int]]:
    """The `(curve, vib)` levels that can exist AT this energy, from energy alone.

    A Rydberg series is attached to a CLOSED channel: above `eps[v]` that
    vibrational channel is open and states attached to it are continuum rather
    than bound. Only thresholds above `e_tot` therefore contribute, and each
    contributes exactly one index::

        binding = eps[v] - e_tot
        n_eff   = 1/sqrt(2 * binding)
        curve   ~ n_eff - 1

    (the last from the measured series, where `Ry_j` has `n_eff ~ j+1`).

    The consequence is a strong constraint and the reason this function exists:
    at fixed energy a HIGHER vibrational level needs a LARGER binding and so a
    LOWER Rydberg index. The admissible set is finite and small, so a basis can
    be checked for covering it -- which turns "is this pole spurious, or is its
    partner merely missing?" from a judgement call into a computation. See
    `basis_covers`, and `qscat.core.assignment.pair_by_overlap` which uses it.

    The set is finite only away from an accumulation region: as `e_tot ->
    eps[v]` the binding tends to zero and the admissible index diverges. That
    happens in the last ~1 mHa below each threshold. Pass `n_eff_max` to cut the
    series there explicitly rather than returning an index no basis will hold.
    """
    thr = np.sort(np.asarray(thresholds, dtype=np.float64))
    out: list[tuple[int, int]] = []
    for v, t in enumerate(thr):
        if t <= e_tot:
            continue
        n = 1.0 / np.sqrt(2.0 * float(t - e_tot))
        if n_eff_max is not None and n > n_eff_max:
            continue
        j = round(n - 1.0)
        if j >= 0:
            out.append((j, v))
    return out


def basis_covers(
    e_tot: float,
    thresholds: npt.ArrayLike,
    basis: BoBasis,
    *,
    curve_tol: int = 1,
    n_eff_max: float | None = None,
) -> bool:
    """Does `basis` contain every level energetically admissible at `e_tot`?

    When it does, a low overlap means the state has no BO partner at all --
    spurious. When it does not, a low overlap is uninformative. Without this
    distinction the two are indistinguishable, and conflating them once nearly
    cost eight genuine `Ry_12..Ry_16` states on H2+.

    `curve_tol` accepts a curve index within +/-1 of the predicted one because
    the `n_eff ~ j+1` mapping is the asymptotic hydrogenic relation and the low
    curves depart from it.
    """
    for j, v in admissible_levels(e_tot, thresholds, n_eff_max=n_eff_max):
        if not any((jj, v) in basis for jj in range(j - curve_tol, j + curve_tol + 1)):
            return False
    return True
