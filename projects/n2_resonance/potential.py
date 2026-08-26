"""N2 electronic potentials for the ²Π_g resonance pole search (sub-project #2).

Same physics as `qscat.model.N2`: closed-form LCP potentials extracted from
reference/eMoScat, verified there against `libs/qscat/tests/test_model.py`.
This module re-implements the same formulas (not shared code) but builds its
`PARAMS` from the same `qscat.model.N2` deck constants. Consistency between
the two independent implementations is NOT guaranteed by construction -- the
formulas could still silently drift apart on edit -- it is guaranteed BY
TEST: `test_potential.py::test_matches_library_model_to_1e_12` cross-checks
this module against `qscat.model.N2` to 1e-12 on every change.

E_res(R)/Gamma(R) are NOT closed form (ECS eigenvalue pole) and are out of
scope here -- that's the pole finder built on top of this potential + the
grid factory in `grid_n2.py`.

NOTE -- model vs. reality (read before comparing v0(R) to N2 spectroscopy):
eMoScat's neutral N2 Morse curve (`v0`, params D_0=0.75102 Ha =~ 20.4 eV,
alpha_0=1.1535) is a MODEL potential built for the resonance study, NOT a
spectroscopic fit to real N2. D_0 is =~2x real N2's actual dissociation
energy (=~9.8 eV), so the resulting neutral vibrational spacing
(omega_e =~ 0.0124 Ha) is =~16% larger than real N2's (0.01074 Ha /
2358 cm^-1). The model's resonance parameters E_res(R0)/Gamma(R0)
(=~2.44 eV / 0.46 eV) DO match real N2 electron-scattering data -- it is
only the *neutral* vibrational ladder that departs from real N2
spectroscopy. This gap is inherited by, and folded into, the
LCP-vs-Houfek-2D differences seen in the cross-section benchmark
(`projects/n2_ti_cross_section/`). See `docs/physics/n2-resonance.md`'s
"Model caveat" section and the model-caveat analysis.
"""

from __future__ import annotations

import numpy as np
from qscat.model import N2

# The N2 deck constants, read from the layer-neutral single source
# `qscat.model.N2` (locked to the eMoScat deck by
# libs/qscat/tests/test_model.py). The dict shape mirrors the historical
# validation/n2/config.json layout so existing PARAMS consumers are
# unaffected.
PARAMS: dict = {
    "reduced_mass": N2.mu,
    "impulsemomentum": N2.ell,
    "potential": {
        "D_0": N2.D0,
        "alpha_0": N2.alpha0,
        "R_0": N2.R0,
        "lambda_inf": N2.lambda_inf,
        "lambda_1": N2.lambda_1,
        "R_lambda": N2.R_lambda,
        "lambda_c": N2.lambda_c,
        "R_c": N2.R_c,
        "alpha_c": N2.alpha_c,
    },
}


def v0(R):
    """Neutral N2 Morse potential (Hartree). Minimum -D_0 at R_0."""
    p = PARAMS["potential"]
    a, R0, D0 = p["alpha_0"], p["R_0"], p["D_0"]
    return D0 * (np.exp(-2 * a * (R - R0)) - 2 * np.exp(-a * (R - R0)))


def lam(R):
    """Interaction strength lambda(R); lambda(R_c) == lambda_c."""
    p = PARAMS["potential"]
    li, l1, Rl, lc, Rc = (
        p["lambda_inf"],
        p["lambda_1"],
        p["R_lambda"],
        p["lambda_c"],
        p["R_c"],
    )
    lam0 = (lc - li) * (1 + np.exp(l1 * (Rc - Rl)))
    return li + lam0 / (1 + np.exp(l1 * (R - Rl)))


def v_int(r, R):
    """Electron-molecule interaction potential (Hartree)."""
    return -lam(R) * np.exp(-PARAMS["potential"]["alpha_c"] * np.asarray(r) ** 2)


def v_eff_el(r, R):
    """Fixed-R electronic effective potential incl. l(l+1)/2r^2 centrifugal term.

    `r` may be complex (ECS-rotated tail points): both `v_int` and the
    centrifugal term are analytic in `r`, so this must NOT coerce to
    `dtype=float` -- doing so silently discards Im(r) and corrupts the
    analytic continuation the exterior-complex-scaling method relies on.
    """
    ell = PARAMS["impulsemomentum"]
    r = np.asarray(r)
    return v_int(r, R) + ell * (ell + 1) / (2 * r**2)
