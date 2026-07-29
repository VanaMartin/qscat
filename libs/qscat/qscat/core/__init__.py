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
  - `electronic_grid`, `nuclear_grid`, `segmented_grid`, `fem_grid_exp_tail`
    -- FEM-DVR-ECS radial grid builders; `segmented_grid` takes eMoScat's
    per-molecule `grids.txt` format (segment-based, uniform-per-segment
    element layout); `fem_grid_exp_tail` is the same real-region layout with
    an EXPONENTIALLY growing ECS tail (H2+'s deck). Parameterized
    (extents/orders are config, not baked in).
  - `vibrational_states` -- the `n` lowest bound eigenpairs of
    `T_nuc(mu) + diag(v0(R))` on a nuclear grid.
  - `channel_vector` -- DVR coefficients of the asymptotic channel function
    `F_{E,l}(r) chi_v(R)`, masked to the unscaled region.
  - `ve_cross_section` -- the exact TI driven Lippmann-Schwinger VE cross
    section, `sigma_{v_init->v'}(E)`, for any `model`.
  - `anion_electronic_states`, `v_dr_diag`, `da_cross_section` -- the anion's
    bound electronic state(s) at `R_inf`, the `V_DR` rearrangement
    interaction, and the exact TI driven-equation dissociative-attachment
    (DA) cross section `sigma_DA(E)`, for any `model`.
  - `dr_cross_section` -- `da_cross_section` generalized for a CHARGED target
    (H2+): a Coulomb incident channel and a loop over `n_channels` Rydberg
    electronic exit states, `sigma_DR(E)`.
  - `local_complex_potential` -- the LCP reduction `(V_d(R), Gamma(R))` of the
    fixed-R electronic resonance to a single complex number per R.
  - `lcp_da_cross_section` -- the LCP dissociative-attachment cross section
    `sigma_DA(E)` via the TI resolvent (1-D nuclear doorway, boundary-value
    outgoing flux; the approximation under test vs. `da_cross_section`).
  - `gaussian_coeffs`, `initial_state`, `outgoing_channel` -- the incident
    Gaussian electron wavepacket and the 2-D initial/outgoing states.
  - `eta_incident`, `eta_outgoing` -- Tannor-Weeks deconvolution factors.
  - `hankel_point_value` -- the outgoing-Hankel-half VALUE at a single
    electronic coordinate (the delta extractor's deconvolution factor).
  - `Extractor`, `propagate`, `sigma_from_correlations`, `td_ve_cross_section`
    -- the time-dependent Pade-propagation route to the same VE cross
    section: `propagate` drives a LIST of `Extractor`s (recorder+transform
    pairs) from one shared trajectory; `TannorWeeks` (below) is the current
    one, wired in via `td_ve_cross_section(method="tw")` (the default).
  - `TannorWeeks` -- the Tannor-Weeks `Extractor` (eta deconvolution +
    elastic free-reference subtraction).
  - `Dirac` -- the delta-distribution `Extractor` (eMoScat
    `DiracTestFunction2d`): TW with a delta test function instead of the
    Gaussian test packet; a flow extractor is a pending sub-project task.
  - `plot_cross_sections` -- generic sigma(E) plotting (no physics baked in).
"""

from __future__ import annotations

from .channels import channel_vector
from .correlation import eta_incident, eta_outgoing, hankel_point_value, outgoing_channel
from .dissociation import (
    anion_electronic_states,
    da_cross_section,
    dr_cross_section,
    v_dr_diag,
)
from .driven import ve_cross_section
from .grids import electronic_grid, fem_grid_exp_tail, nuclear_grid, segmented_grid
from .lcp import lcp_da_cross_section, local_complex_potential
from .plot import plot_cross_sections
from .td_extractors import Dirac, TannorWeeks
from .time_dependent import Extractor, propagate, sigma_from_correlations, td_ve_cross_section
from .vibrational import vibrational_states
from .wavepacket import gaussian_coeffs, initial_state

__all__ = [
    "electronic_grid",
    "fem_grid_exp_tail",
    "nuclear_grid",
    "segmented_grid",
    "vibrational_states",
    "channel_vector",
    "ve_cross_section",
    "anion_electronic_states",
    "v_dr_diag",
    "da_cross_section",
    "dr_cross_section",
    "local_complex_potential",
    "lcp_da_cross_section",
    "gaussian_coeffs",
    "initial_state",
    "outgoing_channel",
    "eta_incident",
    "eta_outgoing",
    "hankel_point_value",
    "Extractor",
    "propagate",
    "sigma_from_correlations",
    "td_ve_cross_section",
    "TannorWeeks",
    "Dirac",
    "plot_cross_sections",
]
