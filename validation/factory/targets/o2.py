"""The O2 Target for the potential factory.

Data: Alt & Houfek, Phys. Rev. A 103, 032829 (2021) — Fig. 2 (p. 032829-3) for
the curves (vector-extracted, see `extract_fig2.py`), Table I (p. 032829-3) for
the electron affinities and D_0, Table II (p. 032829-4) with Eqs. (24)-(27) for
the nonlocal model's energy-dependent width. Spectroscopic constants of 16O2
(R_e, omega_e) from Huber & Herzberg. Nothing here is an observable.
"""

from __future__ import annotations

from projects.potential_factory.ansatz import FlexibleDiatomicModel, SmoothR
from projects.potential_factory.target import (
    CouplingTarget,
    Curve,
    NeutralTarget,
    Provenance,
    ResonanceTarget,
    Target,
)
from validation.factory.targets.o2_data import EV, load_o2

__all__ = ["O2_MU", "o2_target", "o2_seed"]

_PAPER = "Alt & Houfek, Phys. Rev. A 103, 032829 (2021)"
O2_MU = 15.99491461956 * 1822.888486 / 2.0  # m(16O)/2 in electron masses = 14578.4
_EA_O = 1.4611 / EV  # Table I (expt.), p. 032829-3
_D0 = 5.165 / EV  # Table I (expt.), p. 032829-3
_OMEGA_E = 1580.19 / 219474.63  # cm^-1 -> Ha, Huber & Herzberg
_R_E = 2.2819  # bohr (1.20752 A), Huber & Herzberg


def o2_target(R_range: tuple[float, float] = (1.85, 6.0)) -> Target:
    c = load_o2()
    return Target(
        name="O2 (Alt & Houfek 2021, Fig. 2 vector extraction)",
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
            Curve.from_table(c.R, c.v_ion), Curve.from_table(c.R, c.gamma), _EA_O, R_range
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
    """A d-wave well of N2-like depth on O2's frame; the fit moves everything."""
    return FlexibleDiatomicModel(
        mu=O2_MU,
        ell=2,
        D_e=0.19,
        R_e=_R_E,
        betas=(1.4,),
        p=3,
        lam=SmoothR(f_inf=6.0, f_0=1.5, f_1=5.0, R_f=2.3, R_e=_R_E),
        alpha=SmoothR(f_inf=0.45, f_0=0.0, f_1=1.0, R_f=0.0, R_e=_R_E),
        shell=None,
        alpha_b=2.0,
        r_b=3.0,
    )
