"""The O2 Target for the potential factory.

Data: Alt & Houfek, Phys. Rev. A 103, 032829 (2021) — Fig. 2 (p. 032829-3) for
the curves (vector-extracted, see `extract_fig2.py`), Table I (p. 032829-3) for
the electron affinities and D_0, Table II (p. 032829-4) with Eqs. (24)-(27) for
the nonlocal model's energy-dependent width. Spectroscopic constants of 16O2
(R_e, omega_e) from Huber & Herzberg. Nothing here is an observable.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np

from projects.potential_factory.ansatz import FlexibleDiatomicModel, SmoothR, TailR, with_params
from projects.potential_factory.report import FitReport
from projects.potential_factory.target import (
    CouplingTarget,
    Curve,
    NeutralTarget,
    Provenance,
    ResonanceTarget,
    Target,
    polarisation_tail,
)
from validation.factory.targets.o2_data import EV, load_o2, load_so_split

__all__ = ["ALPHA_D_O", "O2_MU", "O2_R_INF", "o2_model_from_report", "o2_seed", "o2_target"]

_PAPER = "Alt & Houfek, Phys. Rev. A 103, 032829 (2021)"
O2_MU = 15.99491461956 * 1822.888486 / 2.0  # m(16O)/2 in electron masses = 14578.4
_EA_O = 1.4611 / EV  # Table I (expt.), p. 032829-3
_D0 = 5.165 / EV  # Table I (expt.), p. 032829-3
_OMEGA_E = 1580.19 / 219474.63  # cm^-1 -> Ha, Huber & Herzberg
_R_E = 2.2819  # bohr (1.20752 A), Huber & Herzberg
# The asymptotic form of the anion curve is THEORY, not the figure: O + O^-
# at -EA(O), approached through the ion--atom polarisation -alpha_d/(2 R^4)
# with the neutral atom's static dipole polarisability alpha_d(O) = 5.3(2)
# a.u. (Schwerdtfeger & Nagle, Mol. Phys. 117, 1200 (2019), Table 1, Z = 8;
# reference/literature/schwerdtfeger-nagle-2019-molphys117-1200.md).
ALPHA_D_O = 5.3
# Where the anion curve is at that form: the remaining binding at finite R
# is the O^- <-> O charge resonance, decaying as exp(-kappa R) with
# kappa = sqrt(2 EA(O)) = 0.33/bohr. Fig. 2 still shows 0.15 eV of it at
# 6 bohr beyond the polarisation term, so it is under the 0.02 eV extraction
# floor from ~12 bohr and under 5 meV from ~16; 14 bohr is the operator's
# choice here -- a per-molecule judgement, not a library constant.
O2_R_INF = 14.0


def o2_target(R_range: tuple[float, float] = (1.85, 6.0), so: int = 0) -> Target:
    """The O2 target; `so = -1 / +1` selects the 2Pi_{1/2} / 2Pi_{3/2}
    spin-orbit component, whose anion curve lies at `V_ion + so * Delta_SO(R)/2`
    with the paper's Fig. 1 splitting (Sec. III A: the two curves lie
    symmetrically around 2Pi_g; the width is the same for both, and the
    atom + ion asymptote moves with the curve, by `so * Delta_SO(inf)/2`).
    `so = 0` is the unsplit 2Pi_g curve (statistical factor 2/3; each
    component carries 1/3)."""
    c = load_o2()
    v_ion = c.v_ion
    ea = _EA_O
    label = ""
    if so:
        if so not in (-1, 1):
            raise ValueError(f"so must be -1, 0 or +1, got {so}")
        R_so, d_so = load_so_split()
        v_ion = v_ion + so * 0.5 * np.interp(c.R, R_so, d_so)
        ea = _EA_O - so * 0.5 * float(d_so[-1])  # the curve's own asymptote
        label = f", 2Pi_{'1/2' if so < 0 else '3/2'} component (Fig. 1 splitting)"
    return Target(
        name=f"O2 (Alt & Houfek 2021, Fig. 2 vector extraction){label}",
        mu=O2_MU,
        ell=2,
        charge=0,
        coordinates=("R",),
        neutral=NeutralTarget(
            Curve.from_table(c.R, c.v0),
            # omega_e deliberately NOT passed: the target is the paper's MRCI curve, whose
            # own ladder (G(1)-G(0) ~ 1607 cm^-1 as extracted) is ~2% stiffer than the
            # spectroscopic 1551 cm^-1 -- the T0 check compares against the curve's ladder.
            {"R_e": _R_E, "D_e": _D0 + 0.5 * _OMEGA_E},
            R_range,
        ),
        resonance=ResonanceTarget(
            Curve.from_table(c.R, v_ion),
            Curve.from_table(c.R, c.gamma),
            ea,
            R_range,
            R_inf=O2_R_INF,
            tail=polarisation_tail(ALPHA_D_O),
        ),
        coupling=CouplingTarget.from_alt_houfek(
            a0=13.836690,
            a1=0.892095,
            a2=-0.935987,
            b0=3.015014,
            b1=0.718160,
            alpha=2.5,
            R_range=(1.85, 2.25),
            eps_window=(0.002, 0.22),
        ),
        provenance={
            "neutral": Provenance(
                _PAPER, "Fig. 2 p. 032829-3 (vector-extracted); Table I p. 032829-3"
            ),
            "resonance": Provenance(_PAPER, "Fig. 2 p. 032829-3 (vector-extracted); EA(O) Table I"),
            "coupling": Provenance(_PAPER, "Table II p. 032829-4; Eqs. (24)-(27) p. 032829-4"),
        },
    )


def o2_seed() -> FlexibleDiatomicModel:
    """A d-wave well of N2-like depth on O2's frame; the fit moves everything.

    `lam(R)` is the long-range-correct `TailR` form (`q = 4`, the polarisation
    power), NOT a sigmoid: O2's anion curve has to rise from the resonant
    wall, peak on the bound branch and settle onto `-EA - alpha_d/(2R^4)`,
    and the sigmoid x polynomial form could hold the table or the asymptote
    but not both (measured: pinning the asymptote cost T1 20 meV and 30 %
    of Gamma, and still missed by 3.5 mHa). Five polynomial coefficients
    are carried from the start (`fit_o2.py` keeps them: `lam_coeffs=5`).
    """
    return FlexibleDiatomicModel(
        mu=O2_MU,
        ell=2,
        D_e=0.19,
        R_e=_R_E,
        betas=(1.4,),
        p=3,
        lam=TailR(f_inf=5.3, coeffs=(1.0, 1.5, 0.0, 0.0, 0.0), R_e=_R_E, p=3, q=4),
        alpha=SmoothR(f_inf=0.45, f_0=0.0, f_1=1.0, R_f=0.0, R_e=_R_E),
        shell=None,
        alpha_b=2.0,
        r_b=3.0,
    )


def o2_model_from_report(report: FitReport) -> FlexibleDiatomicModel:
    """The fitted O2 model, rebuilt from a `FitReport`'s flat parameters on
    the seed's frame: the polynomial lengths of `lam`/`alpha` and the presence
    of a shell are read off the parameter names, then every value is set."""
    p = report.parameters
    m = o2_seed()
    n_lam = sum(1 for k in p if k.startswith("lam.c"))
    n_alpha = sum(1 for k in p if k.startswith("alpha.c"))
    n_beta = sum(1 for k in p if k.startswith("beta"))
    m = replace(
        m,
        betas=tuple(0.0 for _ in range(n_beta)),
        lam=replace(m.lam, coeffs=tuple(0.0 for _ in range(n_lam))),
        alpha=replace(m.alpha, coeffs=tuple(0.0 for _ in range(n_alpha))),
    )
    if "shell.f_inf" in p:
        m = m.with_shell(
            SmoothR(f_inf=0.0, f_0=0.0, f_1=1.0, R_f=m.R_e, R_e=m.R_e), p["alpha_b"], p["r_b"]
        )
    return with_params(m, p)
