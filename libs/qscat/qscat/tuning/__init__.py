"""Automatic FEM-DVR-ECS discretisation tuner.

Public API:
  - `analyze_potential` -- pure potential analysis (no models): samples a
    plain callable `V(x)` on a dense real grid and returns the local-
    wavenumber / forbidden-region-decay profile, classical turning points,
    and boundary singularities (e.g. the `-1/r` origin). This is the sole
    input the mesh/ECS generators consume.
  - `PotentialProfile` -- the frozen dataclass `analyze_potential` returns.
  - `equidistribution_elements` -- adaptive real-region element lengths
    equidistributing de Broglie phase per element (capped by kappa-decay
    length in classically forbidden runs, refined near turning points /
    singularities).
  - `optimal_real_mesh` -- h/p sweep over DVR orders, returning the
    `(mesh, order)` combination with the fewest DVR points.
  - `order_for_wavenumber` -- the smallest DVR order (from a fixed
    candidate set) resolving a given wavenumber at a fixed element length,
    in points-per-wavelength -- see `qscat.tuning.mesh`.
  - `refine_elements_in_window` -- span-preservingly subdivide every real
    element overlapping an `[R_lo, R_hi]` window until each piece is below a
    target length; elements outside the window are untouched -- the LOCAL,
    `min_len`-overriding refinement the resonance-aware nuclear mesh uses to
    super-refine a narrow resonance crossing -- see `qscat.tuning.mesh`.
  - `max_stable_angle` -- the largest ECS rotation angle (capped at the
    double-ECS bound, ~35 deg) for which a potential `V` stays bounded on
    the rotated tail contour.
  - `tune_ecs_tail` -- exp-growth ECS-tail element lengths sized to absorb
    a given outgoing wavenumber `K` down to a target decay.
  - `ProbeResult`, `refine`, `probe_nuclear`, `probe_electronic`,
    `probe_channel_representation` -- the decoupled 1-D convergence probes:
    empirical validators that tell the tuner whether a candidate grid
    resolves the physics (nuclear vibrational levels, the electronic
    bound-state energy, and the cheap/diagnostic channel-representation
    check that catches an unresolved fast outgoing wave) -- see
    `qscat.tuning.probes`.
  - `grid_cost`, `tensor_cost` -- the cost model: exact DVR point counts,
    plus ROUGH anchored estimates (nnz, factor memory/time) of the tensor-
    product problem's sparse-LU cost, for RELATIVE ranking of candidate
    grids -- see `qscat.tuning.metrics`.
  - `propose_grid` -- the one-shot a-priori grid assembler: model adapter ->
    `analyze_potential` -> `optimal_real_mesh` -> `max_stable_angle` +
    `tune_ecs_tail` -> a complete `FemDvrEcsGrid`, the a-priori half of the
    hybrid tuner -- see `qscat.tuning.propose`.
  - `IncidentSpec`, `required_extent`, `tw_analysis` -- incident-state /
    test-function placement: the TD Gaussian-wavepacket spec, the real-
    region extent it forces, and the (best-effort) Tannor-Weeks auto-tune
    that places it for a target energy range -- see `qscat.tuning.incident`.
  - `interaction_region` -- the R-window where the electron-molecule
    interaction `V_int(r, R)` is non-negligible, from `model.v_int` alone
    -- see `qscat.tuning.resonance`.
  - `resonance_curve_arrays` -- the efficient adiabatic resonance-curve
    sampler `(R, V_d(R), Gamma(R))`: dense inside `interaction_region`, a
    single far point at the asymptote -- see `qscat.tuning.resonance`.
  - `refine_to_2d_convergence` -- the general, model-agnostic 2-D-convergence
    FALLBACK: given any scalar observable closing over a real cross-section
    (or a synthetic test function), iteratively `refine`s whichever of the
    electronic/nuclear grids moves the observable more, until the larger
    relative move is under `rtol` or `max_iter` adopted steps is hit -- the
    step-6 spot-check's supervised loop, generalized -- see
    `qscat.tuning.refine2d`.
"""

from __future__ import annotations

from .analyze import PotentialProfile, analyze_potential
from .ecs import max_stable_angle, tune_ecs_tail
from .incident import IncidentSpec, required_extent, tw_analysis
from .mesh import (
    equidistribution_elements,
    optimal_real_mesh,
    order_for_wavenumber,
    refine_elements_in_window,
)
from .metrics import grid_cost, tensor_cost
from .probes import (
    ProbeResult,
    probe_channel_representation,
    probe_electronic,
    probe_nuclear,
    refine,
)
from .propose import propose_grid
from .refine2d import refine_to_2d_convergence
from .resonance import interaction_region, resonance_curve_arrays

__all__ = [
    "IncidentSpec",
    "PotentialProfile",
    "ProbeResult",
    "analyze_potential",
    "equidistribution_elements",
    "grid_cost",
    "interaction_region",
    "max_stable_angle",
    "optimal_real_mesh",
    "order_for_wavenumber",
    "probe_channel_representation",
    "probe_electronic",
    "probe_nuclear",
    "propose_grid",
    "refine",
    "refine_elements_in_window",
    "refine_to_2d_convergence",
    "required_extent",
    "resonance_curve_arrays",
    "tensor_cost",
    "tune_ecs_tail",
    "tw_analysis",
]
