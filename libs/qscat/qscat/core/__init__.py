"""The model-independent VE-scattering engine.

`qscat.core` holds everything the electron-diatomic vibrational-excitation
(VE) solver stack needs that does NOT depend on which molecule is being
solved: the FEM-DVR-ECS grid builders (`grids`), the neutral-molecule
vibrational-states solver (`vibrational`), the asymptotic channel functions
(`channels`), the exact TI driven-equation VE cross section (`driven`),
the incident/outgoing wavepacket construction (`wavepacket`), the
Tannor-Weeks deconvolution factors (`correlation`), the time-dependent
Pade-propagation VE cross section (`time_dependent`), and generic sigma(E)
plotting (`plot`) -- promoted from the N2 projects. See
docs/physics/qscat-core-scattering.md.

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
    an EXPONENTIALLY growing ECS tail (H2+'s deck). Extents/orders are
    parameters; the real-region element layout is not (see `grids`).
  - `ecs_angle_family`, `assert_shared_real_nodes` -- the three-grid family
    (base grid plus one partner per rotated ECS angle) that
    `exact_resonance_states` needs, and the check that two grids share every
    real node while differing only in ECS angle.
  - `vibrational_states`, `VibrationalBasis` -- the `n` lowest bound
    eigenpairs of `T_nuc(mu) + diag(v0(R))` on a nuclear grid, returned as a
    named `(eps, chi)` tuple.
  - `ScatteringProblem` -- the high-level object API bundling
    `(grid, model, n_vib)` once and exposing each observable below as a
    method. The recommended entry point; the functional solvers are the
    low-level layer it delegates to (ADR 0004).
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
  - `resonance_levels`, `lcp_resonance_levels`, `ResonanceLevels` -- the
    BORN-OPPENHEIMER quasi-bound levels `E_v - i*Gamma_v/2`: the nuclear
    eigenvalue problem in the complex curve `V_d - i*Gamma/2`, with the
    golden-rule comparator riding along. NOT Siegert pseudostates -- see
    docs/physics/lcp-resonance-levels.md.
  - `exact_resonance_states`, `ExactResonanceStates` -- the
    approximation-free counterpart of those levels: poles of the FULL 2-D
    S-matrix, accepted only if stable under BOTH ECS angles moved
    independently (docs/physics/exact-2d-resonances.md). Seeds are passed
    in, so it never calls the approximation it measures.
  - `electronic_curves`, `resonance_curve`, `ElectronicCurves`, `bo_basis`,
    `bo_basis_from_levels`, `BoState`, `BoBasis` -- the Born-Oppenheimer
    reference states `phi_j(r;R) chi_v(R)` a pole is verified against:
    `electronic_curves` for a BOUND (ion/Rydberg) electronic curve,
    `resonance_curve` for a RESONANT (neutral/anion) one, and the two
    `bo_basis*` builders putting a vibrational ladder in either, phase-
    aligned across R.
  - `n_eff`, `admissible_levels`, `basis_covers` -- the closed-channel energy
    constraint that separates a SPURIOUS pole from a merely BASIS-LIMITED
    one: `n_eff = 1/sqrt(2*binding)` to the nearest threshold above, so a
    higher vibrational level admits only a lower Rydberg index.
  - `overlap`, `pair_by_overlap`, `pair_one_to_one`, `real_weight`,
    `OverlapPair` -- pairing a pole to a BO level BY OVERLAP (the c-product,
    which is bilinear, so values above 1 are legitimate), the Hungarian
    bijection cross-check, and the box-containment check the overlap cannot
    make. Together they return one of seven verdicts -- see
    docs/physics/h2plus-resonance-states.md.
  - `peak_positions`, `peak_alignment`, `PeakAlignment` -- the distance from
    a level to an observed cross-section peak IN UNITS OF A RESONANCE WIDTH,
    the only scale on which "lands on the peak" means anything.
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
  - `td_da_cross_section` (`method="flow"|"delta"|"tw"`) -- the DA sibling of
    `td_ve_cross_section`: builds the requested `axis="nuclear"` extractor,
    propagates once, returns `sigma_DA(E)`. `td_da_cross_sections_all` runs
    ONE propagation driving all three nuclear extractors and returns
    `{"flow":, "delta":, "tw":}` -- see `docs/physics/td-da.md`.
  - `plot_cross_sections` -- generic sigma(E) plotting (no physics baked in).
  - `plot_route_comparison`, `ComparisonPanel` -- generic LINEAR-scale grid
    of panels, one curve per METHOD (`plot_cross_sections`' sibling, which is
    log-scale with one curve per CHANNEL); `ComparisonPanel` is one panel's
    labelled series plus optional pinned axis limits.
  - `plot_resonance_levels` -- generic complex-level plotting (position vs
    width, plus per-level differences against a chosen baseline series).

`qscat.core.nrm` (a separate subpackage, NOT re-exported here) is the
Feshbach-projection nonlocal resonance model (NRM) for dissociative
attachment (Houfek, Rescigno & McCurdy, PRA 77, 012710 (2008)) -- it sits
between `local_complex_potential`'s LCP reduction and the exact 2-D solvers
above. It is left out of this module's own imports/`__all__` on purpose:
importing it here would defeat the hard boundary above the moment `nrm` ever
grew a runtime `qscat.model` import, since `import qscat.core` would then pull
it in transitively. Import it explicitly:
`from qscat.core.nrm import nrm_da_cross_section`. See
`docs/physics/nonlocal-resonance-model.md`.
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
    real_weight,
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
from .plot import (
    ComparisonPanel,
    plot_cross_sections,
    plot_resonance_levels,
    plot_route_comparison,
)
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
    "BoBasis",
    "BoState",
    "ComparisonPanel",
    "Dirac",
    "ElectronicCurves",
    "ExactResonanceStates",
    "Extractor",
    "Flux",
    "OverlapPair",
    "PeakAlignment",
    "ResonanceLevels",
    "ScatteringProblem",
    "TannorWeeks",
    "VibrationalBasis",
    "admissible_levels",
    "anion_electronic_states",
    "assert_shared_real_nodes",
    "basis_covers",
    "bo_basis",
    "bo_basis_from_levels",
    "channel_vector",
    "da_cross_section",
    "dr_cross_section",
    "ecs_angle_family",
    "electronic_curves",
    "electronic_grid",
    "eta_incident",
    "eta_outgoing",
    "exact_resonance_states",
    "fem_grid_exp_tail",
    "gaussian_coeffs",
    "hankel_point_value",
    "initial_state",
    "lcp_da_cross_section",
    "lcp_resonance_levels",
    "local_complex_potential",
    "n_eff",
    "nuclear_grid",
    "outgoing_channel",
    "outgoing_channel_nuclear",
    "outgoing_surface_wave",
    "overlap",
    "pair_by_overlap",
    "pair_one_to_one",
    "peak_alignment",
    "peak_positions",
    "plot_cross_sections",
    "plot_resonance_levels",
    "plot_route_comparison",
    "propagate",
    "real_weight",
    "resonance_curve",
    "resonance_levels",
    "segmented_grid",
    "sigma_from_correlations",
    "td_da_cross_section",
    "td_da_cross_sections_all",
    "td_ve_cross_section",
    "td_ve_cross_sections_all",
    "v_dr_diag",
    "ve_cross_section",
    "vibrational_states",
]
