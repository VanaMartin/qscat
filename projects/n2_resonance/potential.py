"""N2 electronic potentials for the ²Π_g resonance pole search (sub-project #2).

Same physics as the already-validated `validation/n2/model.py`: closed-form
LCP potentials extracted from reference/eMoScat, verified there against
`model_checks()`. This module re-implements the same formulas, loading
parameters from the same `validation/n2/config.json`, so the two stay in
lockstep by construction rather than by copy-paste drift.

E_res(R)/Gamma(R) are NOT closed form (ECS eigenvalue pole) and are out of
scope here -- that's the pole finder built on top of this potential + the
grid factory in `grid_n2.py`.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

PARAMS: dict = json.loads(
    (Path(__file__).resolve().parents[2] / "validation" / "n2" / "config.json").read_text()
)


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
