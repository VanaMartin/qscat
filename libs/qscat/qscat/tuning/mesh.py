"""Equidistribution mesh generator + h/p quadrature sweep.

`equidistribution_elements` lays out real-region FEM-DVR element boundaries
so each element carries a ~constant de Broglie phase `phase_per_element =
integral of k dx` over the element. In classically forbidden stretches (`k ~
0`, `kappa > 0`) phase does not accumulate, so element length there is
instead capped by the local `kappa`-decay length. `optimal_real_mesh` sweeps
a handful of DVR orders and picks the `(mesh, order)` combination that gives
the fewest total DVR points for a given target accuracy -- the h/p optimum.

`order_for_wavenumber` and `refine_elements_in_window` are the two levers
the resonance-aware nuclear mesh (`qscat.tuning.propose`,
`channel="dissociation"`) uses instead of a worst-case profile merge: a
fixed high DVR order sized to resolve the fast dissociation exit wave, plus
a LOCAL super-refinement of the narrow resonance-crossing region (overriding
`min_len` only there) -- see docs/physics/discretisation-tuning.md.
"""

from __future__ import annotations

import math

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import cumulative_trapezoid

from .analyze import PotentialProfile

FloatArray = NDArray[np.float64]

# Phase-per-(order-1) coefficient `C` in `phase_per_element = C * (order -
# 1)`. Calibrated against the eMoScat N2/NO/F2 decks
# (`validation.tuning.calibrate`): the smallest C making F2's nuclear grid
# reproduce-and-beat the eMoScat DA deck (probe_channel_representation on the
# K~78 dissociative-attachment wave converges at rtol=1e-3 with 609 points vs
# the deck's 974) -- see docs/physics/discretisation-tuning.md for the full
# sweep and the N2/NO findings (their nuclear grids cost more points than
# their decks at this C, traced to `_NUCLEAR_X_MAX_DEFAULT` exceeding their
# per-molecule real-region extent -- a limitation of the a-priori adapter's
# fixed extent default, not a miscalibration of C).
_PHASE_COEFF = 0.10

# How many kappa-decay-lengths (1/kappa) a forbidden-region element may
# span, before it is clamped by max_len anyway.
_DECAY_LENGTHS_PER_ELEMENT = 1.0

# Turning points / singularities are "near" a boundary if within this
# fraction of the local element length -- triggers the halving refinement.
_REFINE_FRACTION = 0.5


def equidistribution_elements(
    profile: PotentialProfile,
    order: int,
    *,
    phase_per_element: float,
    min_len: float,
    max_len: float,
) -> list[float]:
    """Return real-region element lengths equidistributing de Broglie phase.

    Boundaries are placed at `Phi(x) = j * phase_per_element`, where
    `Phi(x) = cumulative_trapezoid(profile.k, profile.x)`. Where `Phi` fails
    to advance (classically forbidden: `k ~ 0`, `kappa > 0`), the local
    element length is instead capped by a `kappa`-based decay length. All
    lengths are then brought into `[min_len, max_len]` -- oversized elements
    subdivided, undersized ones merged forward -- and elements adjacent to a
    turning point or singularity are halved as a post-pass refinement.
    `max_len` is a HARD cap (never exceeded); `min_len` is a SOFT floor (met
    except where doing so would breach `max_len`). Every step only
    regroups/splits existing length, so the returned lengths always sum to
    the real-region domain span `profile.x[-1] - profile.x[0]`.
    """
    x = profile.x
    k = profile.k
    total_span = float(x[-1] - x[0])
    if total_span <= 0.0:
        return []

    phi = cumulative_trapezoid(k, x, initial=0)
    phi_total = float(phi[-1])

    boundary_list: list[float]
    if phi_total > 0.0:
        n_elem = max(int(np.ceil(phi_total / phase_per_element)), 1)
        targets = np.arange(1, n_elem) * phase_per_element
        targets = targets[targets < phi_total]
        # Interior boundaries only ever come from targets < phi_total, so
        # np.interp always lands strictly inside (x[0], x[-1]); the clip is
        # just a defensive guard against float round-off at the phi_total
        # edge. The final boundary is *always* appended explicitly below,
        # so the mesh covers the whole domain regardless of where the last
        # phase target landed -- there is no separate "missing last
        # boundary" gap to close.
        interior = np.clip(np.interp(targets, phi, x), x[0], x[-1])
        boundary_list = [float(x[0]), *interior.tolist(), float(x[-1])]
    else:
        boundary_list = [float(x[0]), float(x[-1])]

    boundary_list = _dedupe_sorted(boundary_list)
    boundary_list = _subdivide_forbidden_gaps(boundary_list, x, profile.kappa, max_len)
    lengths = np.diff(np.asarray(boundary_list, dtype=np.float64))

    # Clamp to [min_len, max_len] WITHOUT changing the total span: an
    # oversized element is subdivided into equal pieces (mirrors
    # `_subdivide_forbidden_gaps`), an undersized one is merged forward
    # into its neighbors until the accumulated length clears the floor.
    # A blind `np.clip` here would silently drop or inflate the domain
    # whenever it actually fired.
    lengths = _clamp_lengths_span_preserving(lengths, min_len, max_len)
    boundaries: FloatArray = np.concatenate(
        [[boundary_list[0]], boundary_list[0] + np.cumsum(lengths)]
    )

    boundaries = _refine_near_features(
        boundaries, profile.turning_points, profile.singularities, min_len
    )

    result: list[float] = np.diff(boundaries).tolist()
    return result


def optimal_real_mesh(
    profile: PotentialProfile,
    *,
    orders: tuple[int, ...] = (6, 8, 10, 14),
    phase_coeff: float = _PHASE_COEFF,
    min_len: float,
    max_len: float,
) -> tuple[list[float], int]:
    """Sweep DVR `orders`, returning the `(mesh, order)` with fewest points.

    For each candidate `order`, `phase_per_element = phase_coeff * (order -
    1)` sets the per-element phase budget; the resulting mesh's DVR point
    count is estimated as `len(elements) * (order - 1)`. The h/p optimum is
    the combination minimizing that count.
    """
    best_mesh: list[float] | None = None
    best_order = orders[0]
    best_points = None

    for order in orders:
        phase_per_element = phase_coeff * (order - 1)
        mesh = equidistribution_elements(
            profile,
            order,
            phase_per_element=phase_per_element,
            min_len=min_len,
            max_len=max_len,
        )
        n_points = len(mesh) * (order - 1)
        if best_points is None or n_points < best_points:
            best_points = n_points
            best_mesh = mesh
            best_order = order

    assert best_mesh is not None
    return best_mesh, best_order


def order_for_wavenumber(
    k: float,
    element_len: float,
    *,
    target_ppw: float = 6.0,
    orders: tuple[int, ...] = (6, 8, 10, 14),
) -> int:
    """Smallest DVR `order` in `orders` resolving wavenumber `k` at a fixed
    real element length `element_len`, in points-per-wavelength.

    Points-per-wavelength at a given `order` is `order * (2*pi/k) /
    element_len` -- `order` DVR points span each element, and `2*pi/k` is
    the de Broglie wavelength of the wave being resolved. Returns the
    smallest `order` (searched in the given, ascending order) with
    `ppw >= target_ppw`; falls back to `max(orders)` if none qualifies (a
    caller asking to resolve an implausibly fast wave at a fixed element
    length still gets the best available order rather than an exception).

    `k <= 0` (no wave to resolve) returns `orders[0]` -- the cheapest order,
    since there is nothing to resolve.
    """
    if k <= 0.0:
        return orders[0]
    wavelength = 2.0 * math.pi / k
    for order in orders:
        ppw = order * wavelength / element_len
        if ppw >= target_ppw:
            return order
    return max(orders)


def refine_elements_in_window(
    real_lengths: list[float],
    x_min: float,
    R_lo: float,
    R_hi: float,
    target_len: float,
) -> list[float]:
    """Span-preservingly subdivide every real element overlapping `[R_lo,
    R_hi]` until each piece is `<= target_len`; elements entirely outside
    the window are returned untouched.

    This is a LOCAL override of `min_len` -- unlike `optimal_real_mesh`'s
    global equidistribution sweep, this targets one narrow feature (e.g. the
    resonance-crossing region a resonant nuclear mesh must super-refine)
    without perturbing element lengths anywhere else. An element only
    counts as "overlapping" if it shares more than a single boundary point
    with `[R_lo, R_hi]` (`hi > R_lo and lo < R_hi`), so a neighbor that
    merely touches the window at one endpoint is left alone. Every input
    length is either kept whole or subdivided into equal pieces, so
    `sum(result) == sum(real_lengths)` exactly (up to float round-off) --
    nothing outside the window is coarsened, and nothing inside it is left
    above `target_len`.
    """
    lengths = np.asarray(real_lengths, dtype=np.float64)
    boundaries = np.concatenate([[float(x_min)], float(x_min) + np.cumsum(lengths)])

    out: list[float] = []
    for lo, hi, length in zip(boundaries[:-1], boundaries[1:], lengths, strict=True):
        overlaps = hi > R_lo and lo < R_hi
        if overlaps and length > target_len:
            n_sub = max(int(np.ceil(length / target_len)), 1)
            out.extend([float(length) / n_sub] * n_sub)
        else:
            out.append(float(length))
    return out


def _clamp_lengths_span_preserving(
    lengths: FloatArray, min_len: float, max_len: float
) -> FloatArray:
    """Clamp element lengths to `[min_len, max_len]` preserving their sum.

    `max_len` is a HARD cap: no emitted element ever exceeds it. `min_len` is
    a SOFT floor: met whenever achievable, but not at the cost of breaching
    `max_len` -- a small element wedged next to a near-`max_len` neighbor may
    end up emitted below `min_len` rather than glued onto that neighbor.

    Oversized elements are subdivided into equal pieces (each `<= max_len`).
    The (now all `<= max_len`) pieces are then merged forward: consecutive
    pieces accumulate until the running total would exceed `max_len` (emit
    what's accumulated so far and start fresh) or clears `min_len` (emit and
    start fresh). A trailing sub-`min_len` remainder is folded into the
    previous emitted element only if that keeps it `<= max_len`; otherwise
    it is emitted on its own (below `min_len`, but the cap is inviolable).
    Every input length is consumed exactly once by either subdivision or
    merging, so `sum(result) == sum(lengths)` exactly (up to float
    round-off) -- unlike `np.clip`, which changes individual lengths without
    redistributing and so silently shrinks or inflates the total domain span
    whenever it actually fires.
    """
    subdivided: list[float] = []
    for length in lengths:
        length = float(length)
        if length > max_len:
            n_sub = max(int(np.ceil(length / max_len)), 1)
            subdivided.extend([length / n_sub] * n_sub)
        else:
            subdivided.append(length)

    merged: list[float] = []
    acc = 0.0
    for length in subdivided:
        if acc > 0.0 and acc + length > max_len:
            merged.append(acc)  # adding would exceed the hard cap -> emit what we have
            acc = 0.0
        acc += length
        if acc >= min_len:
            merged.append(acc)  # reached the floor -> emit
            acc = 0.0
    if acc > 0.0:  # leftover < min_len
        if merged and merged[-1] + acc <= max_len:
            merged[-1] += acc  # fold into previous only if it stays <= max_len
        else:
            merged.append(acc)  # else emit as-is (below min_len: best-effort floor)

    return np.asarray(merged, dtype=np.float64)


def _dedupe_sorted(boundaries: list[float], tol: float = 1e-12) -> list[float]:
    out = [boundaries[0]]
    for b in boundaries[1:]:
        if b - out[-1] > tol:
            out.append(b)
    return out


def _subdivide_forbidden_gaps(
    boundaries: list[float],
    x: FloatArray,
    kappa: FloatArray,
    max_len: float,
) -> list[float]:
    """Insert extra boundaries inside any gap whose local kappa demands it.

    Within a forbidden stretch, cap the element length at
    `_DECAY_LENGTHS_PER_ELEMENT / kappa_local` (the local decay length),
    where `kappa_local` is the mean `kappa` sampled over that gap.
    """
    out = [boundaries[0]]
    for lo, hi in zip(boundaries[:-1], boundaries[1:], strict=True):
        gap = hi - lo
        mask = (x >= lo) & (x <= hi)
        kappa_local = float(np.mean(kappa[mask])) if np.any(mask) else 0.0
        if kappa_local > 0.0:
            decay_len = _DECAY_LENGTHS_PER_ELEMENT / kappa_local
            cap = min(decay_len, max_len)
            if gap > cap:
                n_sub = max(int(np.ceil(gap / cap)), 1)
                sub_bounds = np.linspace(lo, hi, n_sub + 1)[1:].tolist()
                out.extend(sub_bounds)
                continue
        out.append(hi)
    return out


def _refine_near_features(
    boundaries: FloatArray,
    turning_points: FloatArray,
    singularities: FloatArray,
    min_len: float,
) -> FloatArray:
    """Halve the element adjacent to any turning point / singularity."""
    features = np.concatenate([turning_points, singularities])
    if features.size == 0:
        return boundaries

    out = boundaries.tolist()
    changed = True
    # Bounded number of passes: refine can only add finitely many boundaries
    # before hitting min_len, but guard against pathological loops anyway.
    max_passes = 20
    passes = 0
    while changed and passes < max_passes:
        changed = False
        passes += 1
        new_out = [out[0]]
        for lo, hi in zip(out[:-1], out[1:], strict=True):
            length = hi - lo
            near = np.any(
                (features >= lo - _REFINE_FRACTION * length)
                & (features <= hi + _REFINE_FRACTION * length)
            )
            half = length / 2.0
            if near and half >= min_len:
                mid = lo + half
                new_out.append(mid)
                new_out.append(hi)
                changed = True
            else:
                new_out.append(hi)
        out = new_out
    return np.asarray(out, dtype=np.float64)
