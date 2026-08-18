"""The model-independent VE-scattering engine.

`qscat.core` holds everything the electron-diatomic vibrational-excitation
(VE) solver stack needs that does NOT depend on which molecule is being
solved: the FEM-DVR-ECS grid builders (`grids`), the neutral-molecule
vibrational-states solver (`vibrational`), the asymptotic channel functions
(`channels`), the exact TI driven-equation VE cross section (`driven`),
the incident/outgoing wavepacket construction (`wavepacket`), the
Tannor-Weeks deconvolution factors (`correlation`), the time-dependent
Pade-propagation VE cross section (`time_dependent`), and generic sigma(E)
plotting (`plot`), promoted from the N2 projects
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
  - `outgoing_channel_nuclear` -- the NUCLEAR-axis transpose of
    `outgoing_channel` (`phi_c(r) g_out(R)`, an anion electronic bound state
    times a nuclear outgoing test packet): the nuclear `TannorWeeks`
    dissociative-attachment extractor's channel function.
  - `eta_incident`, `eta_outgoing` -- Tannor-Weeks deconvolution factors;
    `eta_outgoing` takes a `mass` keyword (default `1.0`, electronic,
    byte-identical) so the SAME function serves the nuclear-axis
    `TannorWeeks` DA extractor (`mass=model.mu`, `l=0`).
  - `hankel_point_value` -- the outgoing-Hankel-half VALUE at a single
    electronic coordinate (the delta extractor's deconvolution factor).
  - `outgoing_surface_wave` -- the outgoing-Hankel-half VALUE *and* its
    spatial derivative at a single electronic coordinate (the flow/flux
    extractor's per-channel deconvolution pair).
  - `Extractor`, `propagate`, `sigma_from_correlations`, `td_ve_cross_section`,
    `td_ve_cross_sections_all` -- the time-dependent Pade-propagation route
    to the same VE cross section: `propagate` drives a LIST of `Extractor`s
    (recorder+transform pairs) from one shared trajectory;
    `td_ve_cross_section(method="tw"|"delta"|"flow")` selects one extractor
    (`TannorWeeks`/`Dirac`/`Flux`, below; `"tw"` is the default);
    `td_ve_cross_sections_all` runs ONE propagation driving all three and
    returns `{"tw":, "delta":, "flow":}` sigma(E) -- the honest, identical-
    dynamics three-way comparison.
  - `TannorWeeks` -- the Tannor-Weeks `Extractor` (eta deconvolution +
    elastic free-reference subtraction).
  - `Dirac` -- the delta-distribution `Extractor` (eMoScat
    `DiracTestFunction2d`): TW with a delta test function instead of the
    Gaussian test packet.
  - `Flux` -- the flow `Extractor` (eMoScat `FluxTestFunction2d`): the
    time-energy Fourier transform of the probability flux (value AND
    electronic-coordinate derivative) projected onto the outgoing channel at
    a fixed electronic surface -- the Wronskian-like transform built on the
    new `qscat.dvr.dvr_first_derivative_at_node` primitive. All three
    extractors also implement `axis="nuclear"` -- the DISSOCIATIVE
    ATTACHMENT (DA) generalization: the outgoing side moves to the nuclear
    coordinate R, projecting onto `n_channels` anion electronic bound states
    (`anion_electronic_states`) instead of neutral vibrational levels; no
    elastic free-reference subtraction (DA is a pure rearrangement channel).
  - `td_da_cross_section(method="flow"|"delta"|"tw")` -- the DA sibling of
    `td_ve_cross_section`: builds the requested `axis="nuclear"` extractor,
    propagates once, returns `sigma_DA(E)`. `td_da_cross_sections_all` runs
    ONE propagation driving all three nuclear extractors and returns
    `{"flow":, "delta":, "tw":}` -- see `docs/physics/td-da.md`.
  - `plot_cross_sections` -- generic sigma(E) plotting (no physics baked in).
  - `plot_resonance_levels` -- generic complex-level plotting (position vs
    width, plus per-level differences against a chosen baseline series).
"""

from __future__ import annotations

from .assignment import (
    OverlapPair,
    PeakAlignment,
    overlap,
    pair_by_overlap,
    pair_one_to_one,
    peak_alignment,
    peak_positions,
)
from .bo import (
    BoBasis,
    BoState,
    ElectronicCurves,
    admissible_levels,
    basis_covers,
    bo_basis,
    bo_basis_from_levels,
    electronic_curves,
    n_eff,
    resonance_curve,
)
from .channels import channel_vector
from .correlation import (
    eta_incident,
    eta_outgoing,
    hankel_point_value,
    outgoing_channel,
    outgoing_channel_nuclear,
    outgoing_surface_wave,
)
from .dissociation import (
    anion_electronic_states,
    da_cross_section,
    dr_cross_section,
    v_dr_diag,
)
from .driven import ve_cross_section
from .grids import (
    assert_shared_real_nodes,
    ecs_angle_family,
    electronic_grid,
    fem_grid_exp_tail,
    nuclear_grid,
    segmented_grid,
)
from .lcp import (
    ResonanceLevels,
    lcp_da_cross_section,
    lcp_resonance_levels,
    local_complex_potential,
    resonance_levels,
)
from .plot import plot_cross_sections, plot_resonance_levels
from .problem import ScatteringProblem
from .resonance import ExactResonanceStates, exact_resonance_states
from .td_extractors import Dirac, Flux, TannorWeeks
from .time_dependent import (
    Extractor,
    propagate,
    sigma_from_correlations,
    td_da_cross_section,
    td_da_cross_sections_all,
    td_ve_cross_section,
    td_ve_cross_sections_all,
)
from .vibrational import VibrationalBasis, vibrational_states
from .wavepacket import gaussian_coeffs, initial_state

__all__ = [
    "assert_shared_real_nodes",
    "ecs_angle_family",
    "electronic_grid",
    "fem_grid_exp_tail",
    "nuclear_grid",
    "segmented_grid",
    "vibrational_states",
    "VibrationalBasis",
    "ScatteringProblem",
    "channel_vector",
    "ve_cross_section",
    "anion_electronic_states",
    "v_dr_diag",
    "da_cross_section",
    "dr_cross_section",
    "local_complex_potential",
    "lcp_da_cross_section",
    "ResonanceLevels",
    "lcp_resonance_levels",
    "resonance_levels",
    "ExactResonanceStates",
    "exact_resonance_states",
    "ElectronicCurves",
    "electronic_curves",
    "resonance_curve",
    "BoState",
    "BoBasis",
    "bo_basis",
    "bo_basis_from_levels",
    "n_eff",
    "admissible_levels",
    "basis_covers",
    "OverlapPair",
    "overlap",
    "pair_by_overlap",
    "pair_one_to_one",
    "PeakAlignment",
    "peak_positions",
    "peak_alignment",
    "gaussian_coeffs",
    "initial_state",
    "outgoing_channel",
    "outgoing_channel_nuclear",
    "eta_incident",
    "eta_outgoing",
    "hankel_point_value",
    "outgoing_surface_wave",
    "Extractor",
    "propagate",
    "sigma_from_correlations",
    "td_ve_cross_section",
    "td_ve_cross_sections_all",
    "td_da_cross_section",
    "td_da_cross_sections_all",
    "TannorWeeks",
    "Dirac",
    "Flux",
    "plot_cross_sections",
    "plot_resonance_levels",
]
