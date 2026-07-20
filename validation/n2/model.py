"""Closed-form N₂ LCP potentials (extracted from reference/eMoScat, verified).

E_res(R)/Γ(R) are NOT closed form (ECS eigenvalue pole) and are out of scope here.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

PARAMS: dict = json.loads((Path(__file__).parent / "config.json").read_text())


def v0(R):
    """Neutral N₂ Morse potential (Hartree). Minimum -D_0 at R_0."""
    p = PARAMS["potential"]
    a, R0, D0 = p["alpha_0"], p["R_0"], p["D_0"]
    return D0 * (np.exp(-2 * a * (R - R0)) - 2 * np.exp(-a * (R - R0)))


def lam(R):
    """Interaction strength λ(R); λ(R_c) == λ_c."""
    p = PARAMS["potential"]
    li, l1, Rl, lc, Rc = (p["lambda_inf"], p["lambda_1"], p["R_lambda"],
                          p["lambda_c"], p["R_c"])
    lam0 = (lc - li) * (1 + np.exp(l1 * (Rc - Rl)))
    return li + lam0 / (1 + np.exp(l1 * (R - Rl)))


def v_int(r, R):
    """Electron–molecule interaction potential (Hartree)."""
    return -lam(R) * np.exp(-PARAMS["potential"]["alpha_c"] * np.asarray(r) ** 2)


def v_eff_el(r, R):
    """Fixed-R electronic effective potential incl. l(l+1)/2r² centrifugal term."""
    l = PARAMS["impulsemomentum"]
    r = np.asarray(r, dtype=float)
    return v_int(r, R) + l * (l + 1) / (2 * r**2)


def model_checks() -> list[tuple[str, bool, str]]:
    p = PARAMS["potential"]
    R0, D0, Rc, lc = p["R_0"], p["D_0"], p["R_c"], p["lambda_c"]
    out: list[tuple[str, bool, str]] = []
    out.append(("A1 V0(R0) == -D_0", abs(float(v0(R0)) + D0) < 1e-12, f"{float(v0(R0)):.6f} Ha"))
    Rg = np.linspace(1.0, 6.0, 200001)
    out.append(("A2 Morse minimum at R0", abs(Rg[np.argmin(v0(Rg))] - R0) < 1e-3, f"R0={R0}"))
    out.append(("A3 V0(inf) -> 0", abs(float(v0(20.0))) < 1e-6, f"{float(v0(20.0)):.2e}"))
    out.append(("A4 lambda(Rc) == lambda_c", abs(float(lam(Rc)) - lc) < 1e-12, f"{float(lam(Rc)):.6f}"))
    out.append(("A5 V_int negative well", float(v_int(1.0, R0)) < 0.0, f"{float(v_int(1.0, R0)):.4f} Ha"))
    r_t, ell = 2.0, PARAMS["impulsemomentum"]
    centrifugal_ok = abs(float(v_eff_el(r_t, R0)) - (float(v_int(r_t, R0)) + ell * (ell + 1) / (2 * r_t**2))) < 1e-12
    decays = abs(float(v_int(10.0, R0))) < abs(float(v_int(1.0, R0)))
    out.append(("A6 V_eff_el l=2 centrifugal + V_int r-decay", centrifugal_ok and decays,
                f"l(l+1)/2r^2 at r={r_t}, l={ell}"))
    return out
