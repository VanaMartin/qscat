"""N2 electronic potentials for the ²Π_g resonance pole search (sub-project #2).

`v0`/`lam`/`v_int` below ARE `qscat.model.N2`'s bound methods (the single
runtime source) -- not a lockstep copy -- so this module cannot drift from
the library: `test_potential.py::test_potential_is_the_library_model` fails
immediately if a local copy is reintroduced. `PARAMS` still mirrors
`qscat.model.N2`'s deck constants (see below) for existing consumers of the
historical dict shape.

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


v0 = N2.v0
lam = N2.lam
v_int = N2.v_int


def v_eff_el(r, R):
    """Fixed-R electronic effective potential incl. l(l+1)/2r^2 centrifugal term.

    `N2.surface` includes `v0(R)`; this deliberately does NOT, so it is
    `N2.v_int` plus the centrifugal term. `r` may be complex (ECS-rotated
    tail points): both `v_int` and the centrifugal term are analytic in `r`,
    so this must use `dtype=complex128` -- coercing to `dtype=float` would
    silently discard Im(r) and corrupt the analytic continuation the
    exterior-complex-scaling method relies on.
    """
    rc = np.asarray(r, dtype=np.complex128)
    return N2.v_int(rc, R) + N2.ell * (N2.ell + 1) / (2 * rc**2)
