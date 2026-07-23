"""Closed-form N₂ LCP potentials (extracted from reference/eMoScat, verified).

E_res(R)/Γ(R) are NOT closed form (ECS eigenvalue pole) and are out of scope here.

NOTE -- model vs. reality (read before comparing v0(R) to N2 spectroscopy):
eMoScat's neutral N₂ Morse curve (`v0`, params `D_0=0.75102` Ha =~ 20.4 eV,
`alpha_0=1.1535`) is a MODEL potential built for the resonance study, NOT a
spectroscopic fit to real N₂. `D_0` is =~2x real N₂'s actual dissociation
energy (=~9.8 eV), so the resulting neutral vibrational spacing (FEM-DVR
`eps1-eps0 =~ 0.0124` Ha; analytic Morse `omega_e =~ 0.0125` Ha) is =~16%
larger than real N₂'s (0.01074 Ha / 2358 cm⁻¹). The model's resonance
parameters `E_res(R0)`/`Γ(R0)` (=~2.44 eV / 0.46 eV) DO match real N₂
electron-scattering data -- it is only the *neutral* vibrational ladder
that departs from real N₂ spectroscopy. This model-vs-reality gap is
inherited by, and folded into,
the LCP-vs-Houfek-2D differences seen in the cross-section benchmark
(`projects/n2_ti_cross_section/`). See
`.superpowers/sdd/task1fix-report.md` and
`docs/physics/n2-resonance.md`'s "Model caveat" section for the full
analysis; validated by
`projects/n2_ti_cross_section/test_vibrational.py`, which checks the
FEM-DVR vibrational solver against the ANALYTIC Morse spectrum of THIS
potential (model-consistent), not against real N₂.
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
    """Fixed-R electronic effective potential incl. l(l+1)/2r² centrifugal term.

    `r` may be complex (ECS-rotated tail points): both `v_int` and the
    centrifugal term are analytic in `r`, so this must NOT coerce to
    `dtype=float` -- doing so silently discards Im(r) and corrupts the
    analytic continuation the exterior-complex-scaling method relies on
    (see `projects/n2_resonance/potential.v_eff_el`, the lockstep copy this
    is the single source for).
    """
    l = PARAMS["impulsemomentum"]
    r = np.asarray(r)
    return v_int(r, R) + l * (l + 1) / (2 * r**2)


def model_checks() -> list[tuple[str, bool, str]]:
    p = PARAMS["potential"]
    R0, D0, Rc, lc = p["R_0"], p["D_0"], p["R_c"], p["lambda_c"]
    out: list[tuple[str, bool, str]] = []
    out.append(("A1 V0(R0) == -D_0", abs(float(v0(R0)) + D0) < 1e-12, f"{float(v0(R0)):.6f} Ha"))
    Rg = np.linspace(1.0, 6.0, 200001)
    out.append(("A2 Morse minimum at R0", abs(Rg[np.argmin(v0(Rg))] - R0) < 1e-3, f"R0={R0}"))
    out.append(("A3 V0(inf) -> 0", abs(float(v0(20.0))) < 1e-6, f"{float(v0(20.0)):.2e}"))
    out.append(
        ("A4 lambda(Rc) == lambda_c", abs(float(lam(Rc)) - lc) < 1e-12, f"{float(lam(Rc)):.6f}")
    )
    out.append(
        ("A5 V_int negative well", float(v_int(1.0, R0)) < 0.0, f"{float(v_int(1.0, R0)):.4f} Ha")
    )
    r_t, ell = 2.0, PARAMS["impulsemomentum"]
    centrifugal_ok = (
        abs(float(v_eff_el(r_t, R0)) - (float(v_int(r_t, R0)) + ell * (ell + 1) / (2 * r_t**2)))
        < 1e-12
    )
    decays = abs(float(v_int(10.0, R0))) < abs(float(v_int(1.0, R0)))
    out.append(("A6 V_eff_el l=2 centrifugal + V_int r-decay", centrifugal_ok and decays,
                f"l(l+1)/2r^2 at r={r_t}, l={ell}"))
    return out
