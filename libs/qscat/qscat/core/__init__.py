"""The model-independent VE-scattering engine.

`qscat.core` holds everything the electron-diatomic vibrational-excitation
(VE) solver stack needs that does NOT depend on which molecule is being
solved: the FEM-DVR-ECS grid builders (`grids`), the neutral-molecule
vibrational-states solver (`vibrational`), the asymptotic channel functions
(`channels`), the exact TI driven-equation VE cross section (`driven`),
the incident/outgoing wavepacket construction (`wavepacket`), the
Tannor-Weeks deconvolution factors (`correlation`), the time-dependent
Pade-propagation VE cross section (`time_dependent`), and generic sigma(E)
plotting (`plot`), promoted from the N2 projects (sub-project #3, Tasks 4-5)
-- see
`docs/superpowers/specs/2026-07-27-diatomic-ve-scattering-library-design.md`.

**Hard boundary: `qscat.core` must never import `qscat.model` (nor any
`projects.*`) at runtime.** Anything molecule-specific (a potential-energy
surface, a parameter set) is passed in by the caller -- e.g.
`vibrational_states` takes `v0` as a callable, and `driven.ve_cross_section`/
`time_dependent.td_ve_cross_section` take a `model: qscat.model.ResonanceModel`
(imported only under `TYPE_CHECKING`) rather than importing a hardcoded
potential/Hamiltonian. This keeps `core` reusable for NO, F2, and any future
model the `ResonanceModel` protocol admits.

Public API:
  - `electronic_grid`, `nuclear_grid` -- FEM-DVR-ECS radial grid builders,
    parameterized (extents/orders are config, not baked in).
  - `vibrational_states` -- the `n` lowest bound eigenpairs of
    `T_nuc(mu) + diag(v0(R))` on a nuclear grid.
  - `channel_vector` -- DVR coefficients of the asymptotic channel function
    `F_{E,l}(r) chi_v(R)`, masked to the unscaled region.
  - `ve_cross_section` -- the exact TI driven Lippmann-Schwinger VE cross
    section, `sigma_{v_init->v'}(E)`, for any `model`.
  - `gaussian_coeffs`, `initial_state`, `outgoing_channel` -- the incident
    Gaussian electron wavepacket and the 2-D initial/outgoing states.
  - `eta_incident`, `eta_outgoing` -- Tannor-Weeks deconvolution factors.
  - `propagate`, `sigma_from_correlations`, `td_ve_cross_section` -- the
    time-dependent Pade-propagation route to the same VE cross section.
  - `plot_cross_sections` -- generic sigma(E) plotting (no physics baked in).
"""

from __future__ import annotations

from .channels import channel_vector
from .correlation import eta_incident, eta_outgoing, outgoing_channel
from .driven import ve_cross_section
from .grids import electronic_grid, nuclear_grid
from .plot import plot_cross_sections
from .time_dependent import propagate, sigma_from_correlations, td_ve_cross_section
from .vibrational import vibrational_states
from .wavepacket import gaussian_coeffs, initial_state

__all__ = [
    "electronic_grid",
    "nuclear_grid",
    "vibrational_states",
    "channel_vector",
    "ve_cross_section",
    "gaussian_coeffs",
    "initial_state",
    "outgoing_channel",
    "eta_incident",
    "eta_outgoing",
    "propagate",
    "sigma_from_correlations",
    "td_ve_cross_section",
    "plot_cross_sections",
]
