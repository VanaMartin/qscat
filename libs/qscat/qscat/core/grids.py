"""Model-independent FEM-DVR-ECS radial grid builders.

Promoted verbatim from `projects/n2_2d_cross_section/electronic_grid.py`
(`n2_electronic_grid`) and `projects/n2_ti_cross_section/nuclear_grid.py`
(`n2_nuclear_grid`) -- see
`docs/superpowers/specs/2026-07-27-diatomic-ve-scattering-library-design.md`.
Only the element LAYOUT (segment boundaries/lengths, ECS tail construction)
moved here; nothing molecule-specific -- extents/orders are caller-supplied
parameters, not baked-in N2 constants.

`electronic_grid`: layout follows eMoScat's
`input/experimental/N2-model.json` `grids.electronic` -- a finely-resolved
region near the origin where the interaction `-lambda(R) exp(-alpha_c r^2)`
lives, coarsening outward, then an ECS tail of exponentially growing
elements at a single angle. EVERY parameter is exposed because the 2-D
convergence studies vary all of them.

`nuclear_grid`: builds a `qscat.dvr.FemDvrEcsGrid` covering the nuclear
coordinate R -- a real region segmented per `N2.json`'s nuclear layout
(`reference/eMoScat/input/experimental/N2.json`'s `grids.nuclear.real`), then
a `n_complex`-element ECS tail at `angle_deg` out to `r_max`, giving the
outgoing dissociative-attachment boundary condition (see
`docs/physics/diatomic-ve-cross-sections.md`). The default
segment boundaries/lengths (1.5/3.0/4.0/12.0 bohr with 0.5/0.15/0.5/1.0
element lengths) are an exact division (3, 10, 2, 8 elements respectively --
23 real elements total), so no rounding/truncation occurs at the N2
defaults; a caller passing non-default `r_max` still divides exactly since
`n_complex` uniform elements tile the remaining [R0, r_max] span exactly by
construction.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from qscat.dvr import ElementSpec, FemDvrEcsGrid, GridSpec
from qscat.exceptions import GridError

__all__ = ["electronic_grid", "fem_grid_exp_tail", "nuclear_grid", "segmented_grid"]

# --- electronic_grid layout -------------------------------------------------

# (segment_end_bohr, element_length_bohr); the final segment runs to r_max.
_INNER_SEGMENTS: tuple[tuple[float, float], ...] = ((1.0, 0.2), (5.0, 1.0), (7.0, 2.0), (10.0, 3.0))
_OUTER_LENGTH = 4.0


def _ecs_tail(base: float, n: int, *, skip: int, alpha: float) -> list[float]:
    """eMoScat `uniform_increment`/`exp`: `skip` elements at `base`, then growing."""
    return [base if i < skip else base * float(np.exp(alpha * (i - skip + 1))) for i in range(n)]


def electronic_grid(
    *,
    r_max: float = 30.0,
    angle_deg: float = 35.0,
    order: int = 8,
    n_complex: int = 8,
    tail_alpha: float = 0.2,
    tail_skip: int = 2,
) -> FemDvrEcsGrid:
    """Electronic radial grid: real region [0, r_max] + an ECS tail at `angle_deg`.

    The ECS pivot is ``R0 == r_max`` by construction.

    Parameters
    ----------
    r_max : float, optional
        Outer edge of the real region (bohr); the ECS pivot. Must exceed the
        fixed inner segments.
    angle_deg : float, optional
        ECS rotation angle of the complex tail, in degrees.
    order : int, optional
        DVR quadrature order (points per element).
    n_complex : int, optional
        Number of complex (ECS-tail) elements.
    tail_alpha : float, optional
        Geometric growth factor of the tail element lengths.
    tail_skip : int, optional
        Number of leading tail elements before growth begins.

    Returns
    -------
    FemDvrEcsGrid
        The assembled electronic radial grid.

    Raises
    ------
    GridError
        If ``r_max`` does not exceed the fixed inner segments.
    """
    if r_max <= _INNER_SEGMENTS[-1][0]:
        raise GridError(f"r_max must exceed {_INNER_SEGMENTS[-1][0]} bohr, got {r_max}")

    elements: list[ElementSpec] = []
    start = 0.0
    for end, length in _INNER_SEGMENTS:
        k = round((end - start) / length)
        elements += [ElementSpec((end - start) / k) for _ in range(k)]
        start = end

    k_out = max(1, round((r_max - start) / _OUTER_LENGTH))
    elements += [ElementSpec((r_max - start) / k_out) for _ in range(k_out)]

    base = (r_max - start) / k_out
    elements += [
        ElementSpec(h, angle_deg)
        for h in _ecs_tail(base, n_complex, skip=tail_skip, alpha=tail_alpha)
    ]
    return FemDvrEcsGrid(GridSpec(quadrature=order, elements=elements, x_min=0.0))


# --- nuclear_grid layout -----------------------------------------------------

# (segment_end_bohr, element_length_bohr) pairs, per N2.json's nuclear grid
# real region: start=0, lengths=[0.5, 0.15, 0.5, 1.0], points=[1.5, 3.0, 4.0, 12.0].
_REAL_SEGMENTS: tuple[tuple[float, float], ...] = (
    (1.5, 0.5),
    (3.0, 0.15),
    (4.0, 0.5),
    (12.0, 1.0),
)


def nuclear_grid(
    *,
    angle_deg: float = 35.0,
    r_max: float = 40.0,
    n_complex: int = 10,
    quadrature: int = 14,
) -> FemDvrEcsGrid:
    """Build the nuclear FEM-DVR-ECS grid: real region [0, 12] + ECS tail.

    Real elements tile [0, 12] bohr in four segments (see `_REAL_SEGMENTS`),
    finely resolved (0.15 bohr) around the equilibrium bond length so the low
    vibrational states converge. `n_complex` uniform complex elements at
    `angle_deg` tile [12, r_max]. The ECS pivot `R0 == 12.0` bohr by
    construction (GridSpec computes R0 as x_min plus the sum of real element
    lengths).

    Parameters
    ----------
    angle_deg : float, optional
        ECS rotation angle of the complex tail, in degrees.
    r_max : float, optional
        Outer edge of the complex tail (bohr).
    n_complex : int, optional
        Number of complex (ECS-tail) elements tiling ``[12, r_max]``.
    quadrature : int, optional
        DVR quadrature order (points per element).

    Returns
    -------
    FemDvrEcsGrid
        The assembled nuclear radial grid.
    """
    elements: list[ElementSpec] = []
    start = 0.0
    for end, length in _REAL_SEGMENTS:
        span = end - start
        n_seg = round(span / length)
        if abs(n_seg * length - span) > 1e-9:
            raise GridError(f"segment [{start}, {end}] is not an exact multiple of length {length}")
        elements += [ElementSpec(length) for _ in range(n_seg)]
        start = end

    r_pivot = start  # == 12.0
    complex_length = (r_max - r_pivot) / n_complex
    elements += [ElementSpec(complex_length, angle_deg) for _ in range(n_complex)]

    spec = GridSpec(quadrature=quadrature, elements=elements, x_min=0.0)
    return FemDvrEcsGrid(spec)


# --- segmented_grid layout ---------------------------------------------------


def segmented_grid(
    real_segments: Sequence[tuple[int, float]],
    complex_segments: Sequence[tuple[int, float]],
    *,
    angle_deg: float,
    quadrature: int,
    x_min: float = 0.0,
) -> FemDvrEcsGrid:
    """A FEM-DVR-ECS grid from eMoScat's `grids.txt` segment format.

    `real_segments` / `complex_segments` are `(n_elements, endpoint)` pairs:
    from `x_min`, each segment tiles `n` uniform elements up to `endpoint`.
    The complex part is an ECS tail at `angle_deg`; the ECS pivot `R0` is the
    last real endpoint. `complex_segments` may be empty (a pure real grid).
    This is the per-molecule discretisation route -- see
    docs/physics/diatomic-ve-cross-sections.md (DA nuclear grids).
    """
    if quadrature < 2:
        raise GridError(f"quadrature must be >= 2, got {quadrature}")
    elements: list[ElementSpec] = []
    start = x_min
    for label, segs, angle in (
        ("real", real_segments, None),
        ("complex", complex_segments, angle_deg),
    ):
        for n, end in segs:
            if n < 1:
                raise GridError(f"{label} segment ({n}, {end}) has n_elements < 1")
            if end <= start:
                raise GridError(f"{label} endpoint {end} must exceed previous {start}")
            h = (end - start) / n
            elements += [
                ElementSpec(h) if angle is None else ElementSpec(h, angle) for _ in range(n)
            ]
            start = end
    return FemDvrEcsGrid(GridSpec(quadrature=quadrature, elements=elements, x_min=x_min))


# --- fem_grid_exp_tail layout -------------------------------------------------


def fem_grid_exp_tail(
    real_segments: Sequence[tuple[int, float]],
    *,
    angle_deg: float,
    quadrature: int,
    tail_n: int,
    tail_alpha: float = 0.2,
    tail_skip: int = 2,
    x_min: float = 0.0,
) -> FemDvrEcsGrid:
    """A FEM-DVR-ECS grid whose real region is uniform-per-segment (like
    `segmented_grid`) but whose ECS tail grows EXPONENTIALLY rather than
    tiling uniformly -- the H2+ electronic/nuclear deck layout (eMoScat's
    `uniform_increment`/`exp` element spec), which `segmented_grid` cannot
    express.

    `real_segments` are `(n_elements, endpoint)` pairs from `x_min`, tiled
    exactly like `segmented_grid`'s real region. The tail appends `tail_n`
    complex `ElementSpec`s at `angle_deg` whose lengths are `_ecs_tail(base,
    tail_n, skip=tail_skip, alpha=tail_alpha)`, with `base` the length of the
    LAST real element. The ECS pivot `R0` is the last real endpoint.
    """
    if quadrature < 2:
        raise GridError(f"quadrature must be >= 2, got {quadrature}")
    if not real_segments:
        raise GridError("real_segments must be non-empty")

    elements: list[ElementSpec] = []
    start = x_min
    for n, end in real_segments:
        if n < 1:
            raise GridError(f"real segment ({n}, {end}) has n_elements < 1")
        if end <= start:
            raise GridError(f"real endpoint {end} must exceed previous {start}")
        h = (end - start) / n
        elements += [ElementSpec(h) for _ in range(n)]
        start = end

    base = h
    elements += [
        ElementSpec(length, angle_deg)
        for length in _ecs_tail(base, tail_n, skip=tail_skip, alpha=tail_alpha)
    ]
    return FemDvrEcsGrid(GridSpec(quadrature=quadrature, elements=elements, x_min=x_min))
