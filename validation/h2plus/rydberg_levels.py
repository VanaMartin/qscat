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
from qscat.core import vibrational_states
from qscat.dvr import FemDvrEcsGrid, eigen, kinetic
from qscat.exceptions import GridError
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

    For each nuclear grid point `R` (including the complex ECS-tail points,
    NOT just the real region), this diagonalizes the frozen-nucleus
    electronic problem `-1/2 d^2/dr^2 + surface(r, R_inf=R)` (electron mass
    1) on `g_r` and takes the `n_curves` lowest-Re(E) eigenvalues -- the
    Rydberg series `E_Ry0(R) < E_Ry1(R) < ...` at that R. Stacking those
    across `g_R.points` gives each curve tabulated on the nuclear grid.

    This reproduces `qscat.core.anion_electronic_states` bit-for-bit
    (verified to ~1e-14) everywhere that function's own "genuinely bound"
    gate (`|Im(E)| < 1e-6` AND `Re(E) < v0(R_inf)`) succeeds -- the real
    nuclear region, where `R_inf` is real. It does NOT reuse that gate,
    because it fails on the ECS tail: feeding a COMPLEX `R_inf` into
    `surface` makes `v0(R_inf)`/`v_int(r, R_inf)` complex, so every
    eigenvalue there -- including genuinely Rydberg-like ones -- picks up an
    O(Im(v0(R_inf))) imaginary shift (measured ~1.8e-6 Ha at the first tail
    point, i.e. already past the library's 1e-6 tolerance a hair's width
    into the tail). That shift is a curve-parametrization artifact of the
    nuclear ECS rotation -- exactly like `model.v0(R)` itself going complex
    on the tail everywhere else in this codebase -- NOT a resonance width;
    gating on it would raise spuriously and is why this function composes
    `kinetic`/`eigen` directly instead of calling the library gate.

    `n_vib` is a per-curve count, and the curves do NOT all support the same
    one: measured on `proxy_grid()`'s electronic grid with `full_grid()`'s
    real nuclear grid, curve 0 supports 5 clean bound levels while curves 1-4
    each support at least 12. A single `n_vib` that fits the shallowest curve
    therefore truncates the others, and one that fits the others makes
    `vibrational_states` raise on the shallowest. `allow_partial=True` asks
    each curve for as many of `n_vib` as it can supply, padding the rest of
    that row with `NaN`; callers must skip non-finite entries. The default
    keeps the strict behaviour -- every curve supplies all `n_vib` or the
    call raises -- so a caller that assumes a full rectangular table still
    finds out when it is wrong.

    Each curve is then fed to `vibrational_states` as the nuclear potential.
    `vibrational_states` builds `T_nuc(mu) + diag(v0(grid.points))`, calling
    `v0` with EXACTLY `grid.points` (`qscat/core/vibrational.py`) -- so a
    closure that ignores its argument and returns the curve already
    tabulated on those same points is an exact lookup, not an
    approximation/interpolation.
    """
    pts = g_R.points
    curves = np.empty((n_curves, pts.size), dtype=np.complex128)
    for j, R in enumerate(pts):
        H_el = kinetic(g_r, 1.0) + np.diag(model.surface(g_r.points, complex(R)))
        E, _ = eigen(H_el)  # ascending Re(E), same formula as anion_electronic_states
        curves[:, j] = E[:n_curves]

    energies = np.full((n_curves, n_vib), np.nan, dtype=np.float64)
    for n in range(n_curves):
        curve = curves[n]

        def v_n(
            _R: npt.ArrayLike, _curve: npt.NDArray[np.complex128] = curve
        ) -> npt.NDArray[np.complex128]:
            # Exact lookup, not interpolation -- see the docstring above.
            return _curve

        if not allow_partial:
            basis = vibrational_states(g_R, model.mu, n_vib, v_n)
            energies[n] = np.asarray(basis.eps, dtype=np.float64)
            continue

        # Walk down from n_vib to the deepest count this curve actually
        # supports. `vibrational_states` selects the n lowest-Re(E)
        # eigenvalues and rejects the batch if ANY is quasi-continuum, so a
        # smaller n is a strictly cleaner subset -- there is no n where it
        # raises and n-1 contains a level n did not.
        for count in range(n_vib, 0, -1):
            try:
                basis = vibrational_states(g_R, model.mu, count, v_n)
            except GridError:
                continue
            energies[n, :count] = np.asarray(basis.eps, dtype=np.float64)
            break

    return RydbergLevels(curves=curves, energies=energies)
