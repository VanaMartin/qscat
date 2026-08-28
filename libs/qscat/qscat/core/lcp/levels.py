"""Born-Oppenheimer nuclear eigenproblem in the LCP curve.

`ResonanceLevels`, `lcp_resonance_levels`, `resonance_levels`.
"""

from __future__ import annotations

import os
import warnings
from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal, overload

import numpy as np
import numpy.typing as npt

from qscat.dvr import FemDvrEcsGrid, eigen, kinetic
from qscat.ecs import match_angle_stable
from qscat.linalg import c_product

from .._archive import load_dataclass_npz, save_dataclass_npz
from ..grids import assert_shared_real_nodes
from .curve import _assemble_lcp, _walk_from_anion_seed

if TYPE_CHECKING:
    from qscat.model import ResonanceModel

# `_levels_from`: smallest |sum_i c_i^2| a state may have and still be normalized
# by the bilinear c-product. `eigen` returns Euclidean-unit vectors, so this
# quantity is <= 1 and equals 1 for a real vector; a value near 0 means the state
# is (numerically) self-orthogonal and `c / sqrt(sum c^2)` blows it up into noise.
_C_NORM_TOL = 1e-8

# `resonance_levels`: largest `Gamma(R)` tolerated, without warning, in the
# region where the anion curve lies BELOW the neutral one and autodetachment
# is energetically closed (Vana 2017, Sec. 1.5) -- there Gamma must vanish, so
# anything above this is pole-finder noise leaking into a closed region.
_CLOSED_REGION_GAMMA_TOL = 1e-6


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
    - `widths`: `Gamma_v = -2 Im E_v` (Hartree), UNCLAMPED -- after the
      `Im E <= atol` physicality filter the most negative representable value
      is `-2*atol`; a small negative width is a round-off diagnostic, and
      hiding it behind a clamp is how it goes unnoticed. Same convention as
      `ExactResonanceStates.widths`. A level below the anion dissociation
      limit carries only the ELECTRONIC autodetachment width; one above it
      also carries a NUCLEAR (dissociative) width. Both come out of the one
      diagonalization.
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

    def save(self, path: str | os.PathLike[str]) -> None:
        """Write to an `.npz` archive under the dataclass's own field names.

        Shares its mechanism with `ExactResonanceStates.save` -- see
        `qscat.core._archive` for why the pair round-trips through the
        dataclass's own field names rather than a hand-rolled cache.
        """
        save_dataclass_npz(self, path)

    @classmethod
    def load(cls, path: str | os.PathLike[str]) -> ResonanceLevels:
        """Read back a `save()` file, checking every field is present."""
        return cls(**load_dataclass_npz(cls, path))


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

    widths = -2.0 * energies.imag
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
