"""Born-Oppenheimer quasi-bound levels in the H2+ Rydberg curves.

The neutral's bound electronic states at fixed R form the Rydberg series
`E_Ryn(R)`; the vibrational levels those curves support are the quasi-bound
states a DR cross-section peak is conventionally assigned to. They carry NO
width here: a Rydberg state is bound in its own curve, and the width it
really has comes from the coupling to the dissociative channel the BO
picture drops. This is `docs/physics/lcp-resonance-levels.md`'s
angle-stability caveat applying in reverse: these levels are the reference
BEFORE that coupling is added back in (see `exact_poles.py`), not a claim
they are the final, physical widths.

Index convention (MEASURED against the thesis's own text, not assumed):
Fig. 4.3's caption (Vana 2017, p. 64) states the convention explicitly --
"the vibrational levels in the electron energy potential
`E_Ry_j` [are] labeled with `omega_i^j`, where `i` is the vibrational level
and `j` stands for the corresponding Rydberg state number." That is the
OPPOSITE of what an earlier draft of this docstring assumed. In this
module's arrays the FIRST index of `curves`/`energies` is the Rydberg curve
(built by looping `n_curves` electronic eigenvalues -- the thesis's `j`;
this construction reproduces Table 4.1's Ry_0/Ry_1/Ry_2 asymptotic electron
energies to ~1e-4 on this module's reduced-electronic-grid deck: measured
-1.384927 / -0.124999 / -0.054810 Ha vs published -1.38492776 / -0.12499996
/ -0.05481037 Ha). The SECOND index is the vibrational quantum number
within that curve -- the thesis's `i`. So
`energies[n, v]` is the thesis's `omega_v^n`, NOT `omega_n^v`. Fig. 4.7's
three energy-window panels are bounded by consecutive cation vibrational
thresholds, and each panel's dominant dashed lines share one fixed (low)
`i` across several curves `j` -- the windowing tracks the vibrational index
(this module's SECOND axis), not the curve index (this module's FIRST
axis); `energies[n, :]` sweeping across all three published windows as `v`
increases (measured directly for curve `Ry_2`: its levels run
-0.034 / -0.025 / -0.016 / -0.008 / -0.001 / +0.006 / +0.012 / +0.018 /
+0.022 / +0.026 Ha for `v` = 0..9, crossing all three windows) is
exactly this pattern, not a bug.

Grid limitation (measured): on `proxy_grid()`'s electronic grid (~60 bohr
real region, 377 pts -- the reduced deck, NOT the full 1300-bohr/1406-pt
one), curves `n >= 5` support ZERO numerically clean bound vibrational levels
(`vibrational_states` already raises at the first, `n_vib=1`, probe) when
combined with `full_grid()`'s real nuclear grid -- the truncated electronic
box is too small to contain those more diffuse Rydberg orbitals. Curves 0-4
are unaffected (Ry_0/Ry_1/Ry_2's asymptotes reproduce Table 4.1 to ~5e-5,
see `test_curve_asymptotes_match_table_4_1` in `test_rydberg_levels.py`).
A caller needing curves beyond
`Ry_4` must call this function with the FULL electronic grid
(`full_grid().grids[0]`), not `proxy_grid()`'s.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from qscat.core.bo import bo_basis, electronic_curves
from qscat.dvr import FemDvrEcsGrid
from qscat.model import ResonanceModel

__all__ = ["RydbergLevels", "rydberg_levels"]


@dataclass(frozen=True)
class RydbergLevels:
    curves: npt.NDArray[np.complex128]  # (n_curves, n_R): E_Ryn(R) on g_R.points
    energies: npt.NDArray[np.float64]  # (n_curves, n_vib): real, ascending per curve


def rydberg_levels(
    model: ResonanceModel,
    g_r: FemDvrEcsGrid,
    g_R: FemDvrEcsGrid,
    *,
    n_curves: int,
    n_vib: int,
    allow_partial: bool = False,
) -> RydbergLevels:
    """`E_Ryn(R)` and the vibrational levels each curve supports.

    A thin adapter over `qscat.core.bo`: `electronic_curves` builds the Rydberg
    series over the nuclear grid and `bo_basis` puts a vibrational ladder in
    each. Both were promoted out of this module -- nothing here is specific to
    H2+ except the published index convention documented above, which is why
    this shim survives while the numerics moved.

    Levels are returned REAL. These curves are genuinely bound, so their levels
    carry no width; the library keeps them complex because the same function
    also serves resonance curves, where the width is part of the answer.

    `n_vib` is a per-curve REQUEST -- see `bo_basis` on `allow_partial`.
    """
    curves = electronic_curves(model, g_r, g_R, n_curves=n_curves, with_states=False)
    basis = bo_basis(curves, g_R, model.mu, n_vib=n_vib, allow_partial=allow_partial)
    return RydbergLevels(
        curves=curves.energies,
        energies=np.asarray(basis.energies.real, dtype=np.float64),
    )
