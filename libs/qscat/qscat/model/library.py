"""Per-molecule `DiatomicResonanceModel` parameter registry.

These constants are PUBLISHED, and the eMoScat decks
(`reference/eMoScat/input/{N2,NO,F2}/2D_model.txt`) agree with them:

- N2 and NO: Houfek, Rescigno & McCurdy, Phys. Rev. A 73, 032721 (2006),
  Table I -- the paper that introduced this 2-D model.
- F2: Houfek, Rescigno & McCurdy, Phys. Rev. A 77, 012710 (2008), Table I.
- H2P: Hvizdos, master's thesis (Charles Univ., 2016) Table 1.1; Vana 2017
  thesis Table 1.2; Hvizdos et al., Phys. Rev. A 97, 022704 (2018), Sec. II A
  -- all three give the reduced mass 918.076 = m_p/2.

(The GRIDS remain eMoScat-deck provenance; only the model constants are
published.) `N2` reproduces `validation/n2/config.json`'s
`reduced_mass`/`impulsemomentum`/`potential` block, and
`projects/n2_2d_cross_section/hamiltonian2d.py`'s `MU`/`ELL`
(`libs/qscat/tests/test_model.py` gates this). This is the single place a
new molecule is added to the registry.
"""

from __future__ import annotations

from .diatomic import DiatomicResonanceModel
from .flexible import FlexibleDiatomicModel, SmoothR, TailR
from .ionic import IonicResonanceModel

__all__ = ["F2", "H2P", "N2", "NO", "O2", "O2_SO12", "O2_SO32"]

N2 = DiatomicResonanceModel(
    mu=12766.36,
    ell=2,
    D0=0.75102,
    alpha0=1.1535,
    R0=2.01943,
    lambda_inf=6.21066,
    lambda_1=1.05708,
    R_lambda=-27.9833,
    lambda_c=5.38022,
    R_c=2.405,
    alpha_c=0.4,
)

NO = DiatomicResonanceModel(
    mu=13614.16,
    ell=1,
    D0=0.2363,
    alpha0=1.5710,
    R0=2.1570,
    lambda_inf=6.3670,
    lambda_1=5.0000,
    R_lambda=2.0843,
    lambda_c=6.0500,
    R_c=2.2850,
    alpha_c=1.0,
)

F2 = DiatomicResonanceModel(
    mu=17315.99,
    ell=1,
    D0=0.05980,
    alpha0=1.51610,
    R0=2.69060,
    lambda_inf=18.8490,
    lambda_1=3.21300,
    R_lambda=1.8320,
    lambda_c=18.1450,
    R_c=2.5950,
    alpha_c=3.0,
)

H2P = IonicResonanceModel(
    # m_p/2 for the modern proton mass (1836.15267/2), as given by Vana 2017
    # Table 1.2 and Hvizdos et al., Phys. Rev. A 97, 022704 (2018) Sec. II A.
    # eMoScat's JSON deck carries 918.25, which both publications contradict.
    mu=918.076,
    ell=1,
    charge=-1,
    V0=0.1027,
    R0=2.0,
    alpha=0.69,
    a1=1.6435,
    a2=6.2,
    a3=0.0125,
    a4=1.15,
)

# O2: NOT a published parameter set but the potential factory's FIT to the
# curves of Alt & Houfek, Phys. Rev. A 103, 032829 (2021), Fig. 2 (vector-
# extracted, ~0.02 eV) -- the 2-D model whose fixed-R resonance curve
# E_res(R), width Gamma(R) and anion curve reproduce the paper's to the
# extraction floor (E_res rms 20 meV, Gamma 8 %/14 %, crossing 2.289 bohr),
# with the anion ending at O + O^- (-EA(O) = -1.4611 eV) through the
# polarisation tail -alpha_d(O)/(2R^4) from R_inf = 14 bohr. Every constant
# below is the committed report `validation/factory/results/o2-fit-report.json`
# verbatim (`validation/factory/test_o2_report.py` locks them); the fit is
# described in docs/physics/potential-factory.md. The 2Pi_g resonance is a
# d-wave (ell = 2) like N2's. Nothing here is fitted to an observable.
_O2_R_E = 2.268012257109915  # the FITTED equilibrium (EMO R_e)
# The y_p reference radius inside lam/alpha/shell is a FRAME constant, not a
# fitted parameter: it stays at the seed's Huber & Herzberg R_e = 2.2819 bohr
# (`validation/factory/targets/o2.py::o2_seed`), and the report does not
# record it -- the surface-equality lock in test_o2_report.py is what holds it.
_O2_Y_REF = 2.2819
O2 = FlexibleDiatomicModel(
    mu=15.99491461956 * 1822.888486 / 2.0,  # m(16O)/2 in electron masses = 14578.47
    ell=2,
    D_e=0.19331866564928865,
    R_e=_O2_R_E,
    betas=(
        1.4523575847133374,
        0.10381040322431209,
        -0.004015128162685024,
        0.3910799122325903,
        1.9183648370285529,
    ),
    p=3,
    lam=TailR(
        f_inf=5.064307074513824,
        coeffs=(
            -0.20821474894387013,
            0.8116574159293211,
            0.10252797159229596,
            0.5444402680526038,
            -1.3026360633646719,
            3.3083561788901363,
            32.788738338929384,
            -73.14541076473161,
            37.40650902876745,
        ),
        R_e=_O2_Y_REF,
        p=3,
        q=4,
    ),
    alpha=SmoothR(
        f_inf=0.3612778078765902, f_0=0.0, f_1=1.0, R_f=0.0, coeffs=(0.0, 0.0, 0.0), R_e=_O2_Y_REF
    ),
    # T3's shell fit left a zero-strength shell (1e-9); it is carried so the
    # parameter set matches the report key for key, and contributes nothing.
    shell=SmoothR(f_inf=1e-09, f_0=1e-09, f_1=1.0, R_f=2.05, R_e=_O2_Y_REF),
    alpha_b=2.0,
    r_b=3.0,
)

# The two spin-orbit components of O2: the same fit with the anion curve moved
# by -/+ Delta_SO(R)/2 (Alt & Houfek Fig. 1, Sec. III A -- the 2Pi_1/2 and 2Pi_3/2
# curves lie symmetrically around 2Pi_g, same width, asymptotes at
# -EA(O) +/- Delta_SO(inf)/2), obtained by polishing O2's lam(R)/alpha(R) only
# (`validation/factory/fit_o2_so.py`; reports o2-so12-/o2-so32-fit-report.json,
# locked by test_o2_report.py). Each carries the statistical factor 1/3; their
# sum is the spin-orbit-resolved VE cross section (p. 032829-4).
O2_SO12 = FlexibleDiatomicModel(
    mu=15.99491461956 * 1822.888486 / 2.0,
    ell=2,
    D_e=0.19331866564928865,
    R_e=2.268012257109915,
    betas=(
        1.4523575847133374,
        0.10381040322431209,
        -0.004015128162685024,
        0.3910799122325903,
        1.9183648370285529,
    ),
    p=3,
    lam=TailR(
        f_inf=5.025820694299848,
        coeffs=(
            -0.20746494863175372,
            0.811406107431042,
            0.10120916277113931,
            0.5537703046708998,
            -1.308995528994924,
            3.2579374412303177,
            32.862701258238495,
            -73.13157294390486,
            37.36870749341204,
        ),
        R_e=_O2_Y_REF,
        p=3,
        q=4,
    ),
    alpha=SmoothR(
        f_inf=0.35836239142537873,
        f_0=0.0,
        f_1=1.0,
        R_f=0.0,
        coeffs=(0.0, 0.0, 0.0),
        R_e=_O2_Y_REF,
    ),
    shell=SmoothR(
        f_inf=1e-09,
        f_0=1e-09,
        f_1=1.0,
        R_f=2.05,
        R_e=_O2_Y_REF,
    ),
    alpha_b=2.0,
    r_b=3.0,
)

O2_SO32 = FlexibleDiatomicModel(
    mu=15.99491461956 * 1822.888486 / 2.0,
    ell=2,
    D_e=0.19331866564928865,
    R_e=2.268012257109915,
    betas=(
        1.4523575847133374,
        0.10381040322431209,
        -0.004015128162685024,
        0.3910799122325903,
        1.9183648370285529,
    ),
    p=3,
    lam=TailR(
        f_inf=5.102921705841919,
        coeffs=(
            -0.2089645731643405,
            0.8118910038742814,
            0.10394620036618253,
            0.535378201217016,
            -1.2986782171831563,
            3.3593321281929445,
            32.728595606603704,
            -73.18273419451207,
            37.455649430870515,
        ),
        R_e=_O2_Y_REF,
        p=3,
        q=4,
    ),
    alpha=SmoothR(
        f_inf=0.3642029360689995,
        f_0=0.0,
        f_1=1.0,
        R_f=0.0,
        coeffs=(0.0, 0.0, 0.0),
        R_e=_O2_Y_REF,
    ),
    shell=SmoothR(
        f_inf=1e-09,
        f_0=1e-09,
        f_1=1.0,
        R_f=2.05,
        R_e=_O2_Y_REF,
    ),
    alpha_b=2.0,
    r_b=3.0,
)
