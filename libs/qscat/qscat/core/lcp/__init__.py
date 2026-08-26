"""Model-independent local complex potential (LCP) approximation.

The local-complex-potential (LCP) approximation of dissociative attachment
(and vibrational excitation) reduces the full 2-D (electronic r x nuclear R)
resonance problem to a 1-D nuclear problem by replacing the fixed-R
electronic resonance with a single complex number at each R -- eMoScat's
`ModelLCP` (`v0(R) + E_res(R)` real part, width `-2*Im(E_res(R))`). This is
the RESEARCH-PROGRAM "approximation under test": the exact 2-D solver
(`qscat.core.dissociation`/`driven`) is the oracle, and this package's
reduction is the thing being measured against it -- not a description of the
"real" physics.

Three submodules, one capability each:

- `.curve`: the fixed-`R` electronic pole machinery (`local_complex_potential`,
  `resonance_pole_walk`, `resonance_eigenstate`,
  `resonance_eigenstate_at_peak_width`) that produces the `(V_d(R), Gamma(R))`
  curve.
- `.cross_section`: the 1-D TI resolvent solvers on that curve
  (`lcp_da_cross_section`, `lcp_ve_cross_section`).
- `.levels`: the Born-Oppenheimer nuclear eigenproblem in that curve
  (`ResonanceLevels`, `lcp_resonance_levels`, `resonance_levels`).

See `docs/physics/lcp-resonance-levels.md` and
`docs/physics/diatomic-ve-cross-sections.md`.
"""

from .cross_section import lcp_da_cross_section, lcp_ve_cross_section
from .curve import (
    local_complex_potential,
    resonance_eigenstate,
    resonance_eigenstate_at_peak_width,
    resonance_pole_walk,
)
from .levels import ResonanceLevels, lcp_resonance_levels, resonance_levels

__all__ = [
    "ResonanceLevels",
    "lcp_da_cross_section",
    "lcp_resonance_levels",
    "lcp_ve_cross_section",
    "local_complex_potential",
    "resonance_eigenstate",
    "resonance_eigenstate_at_peak_width",
    "resonance_levels",
    "resonance_pole_walk",
]
