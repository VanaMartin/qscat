"""Hand-built FEM-DVR-ECS grid factory for the N2 electronic problem
(sub-project #2, Task 1).

Builds a `qscat.dvr.FemDvrEcsGrid` with:
  - `n_real` uniform real elements tiling [0, r_pivot] (the well + centrifugal
    region, where V_eff_el lives);
  - `n_complex` uniform complex (ECS) elements tiling [r_pivot, r_max] at a
    single shared `angle_deg`, capturing the outgoing resonance wave.

Two grids built at different angles but identical (r_pivot, n_real, r_max,
n_complex) MUST share identical `real_points` -- only the tail rotation
differs. This holds by construction: `FemDvrEcsGrid.real_points` is computed
purely from element `length`s and `x_min` (never from `angle_deg`), so as
long as the element length sequence is unchanged, the real (unscaled)
bookkeeping points are bit-identical across angles. The pole finder (a later
task) relies on this to compare Hamiltonians built at different ECS angles.
"""

from __future__ import annotations

from qscat.dvr import ElementSpec, FemDvrEcsGrid, GridSpec


def n2_electronic_grid(
    angle_deg: float,
    *,
    r_pivot: float = 10.0,
    n_real: int = 8,
    r_max: float = 30.0,
    n_complex: int = 8,
    quadrature: int = 8,
) -> FemDvrEcsGrid:
    """Build the N2 electronic FEM-DVR-ECS grid: real region + ECS tail.

    `n_real` uniform real elements tile [0, r_pivot]; `n_complex` uniform
    complex elements at `angle_deg` tile [r_pivot, r_max]. The ECS pivot
    `R0 == r_pivot` by construction (GridSpec computes R0 as x_min plus the
    sum of real element lengths).
    """
    real_length = r_pivot / n_real
    complex_length = (r_max - r_pivot) / n_complex

    elements = [ElementSpec(real_length) for _ in range(n_real)]
    elements += [ElementSpec(complex_length, angle_deg) for _ in range(n_complex)]

    spec = GridSpec(quadrature=quadrature, elements=elements, x_min=0.0)
    return FemDvrEcsGrid(spec)
