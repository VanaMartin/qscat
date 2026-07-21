"""Hand-built FEM-DVR-ECS nuclear grid factory for the N2 vibrational /
time-independent cross-section problem (sub-project #3, Task 1).

Builds a `qscat.dvr.FemDvrEcsGrid` covering the nuclear coordinate R:
  - a real region [0, 12] bohr, segmented per `N2.json`'s nuclear layout
    (`reference/eMoScat/input/experimental/N2.json`'s `grids.nuclear.real`):
    boundaries at R = 1.5, 3.0, 4.0, 12.0 bohr with element lengths 0.5,
    0.15, 0.5, 1.0 respectively -- i.e. finely resolved (0.15 bohr elements)
    around the N2 equilibrium bond length R0 = 2.01943 bohr, coarser
    elsewhere. This mirrors `projects/n2_resonance/grid_n2.py`'s pattern
    (hand-built `GridSpec`/`ElementSpec` list) but tiles multiple segments of
    differing element size instead of one uniform real region;
  - a 35 degree ECS tail of `n_complex` uniform elements out to `r_max`,
    giving the outgoing dissociative-attachment boundary condition (see
    `.superpowers/sdd/ti-cross-section-extraction.md` section 7).

The segment boundaries/lengths are an exact division (3, 10, 2, 8 elements
respectively -- 23 real elements total), so no rounding/truncation occurs.
"""

from __future__ import annotations

from qscat.dvr import ElementSpec, FemDvrEcsGrid, GridSpec

__all__ = ["n2_nuclear_grid"]

# (segment_end_bohr, element_length_bohr) pairs, per N2.json's nuclear grid
# real region: start=0, lengths=[0.5, 0.15, 0.5, 1.0], points=[1.5, 3.0, 4.0, 12.0].
_REAL_SEGMENTS: tuple[tuple[float, float], ...] = (
    (1.5, 0.5),
    (3.0, 0.15),
    (4.0, 0.5),
    (12.0, 1.0),
)


def n2_nuclear_grid(
    *,
    angle_deg: float = 35.0,
    r_max: float = 40.0,
    n_complex: int = 10,
    quadrature: int = 14,
) -> FemDvrEcsGrid:
    """Build the N2 nuclear FEM-DVR-ECS grid: real region [0, 12] + ECS tail.

    Real elements tile [0, 12] bohr in four segments (see `_REAL_SEGMENTS`),
    finely resolved (0.15 bohr) around the equilibrium bond length so the low
    vibrational states converge. `n_complex` uniform complex elements at
    `angle_deg` tile [12, r_max]. The ECS pivot `R0 == 12.0` bohr by
    construction (GridSpec computes R0 as x_min plus the sum of real element
    lengths).
    """
    elements: list[ElementSpec] = []
    start = 0.0
    for end, length in _REAL_SEGMENTS:
        span = end - start
        n_seg = round(span / length)
        if abs(n_seg * length - span) > 1e-9:
            raise ValueError(
                f"segment [{start}, {end}] is not an exact multiple of length {length}"
            )
        elements += [ElementSpec(length) for _ in range(n_seg)]
        start = end

    r_pivot = start  # == 12.0
    complex_length = (r_max - r_pivot) / n_complex
    elements += [ElementSpec(complex_length, angle_deg) for _ in range(n_complex)]

    spec = GridSpec(quadrature=quadrature, elements=elements, x_min=0.0)
    return FemDvrEcsGrid(spec)
