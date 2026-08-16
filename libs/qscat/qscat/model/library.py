"""Per-molecule `DiatomicResonanceModel` parameter registry.

Verbatim eMoScat deck parameters (`reference/eMoScat/input/{N2,NO,F2}/
2D_model.txt`), transcribed exactly (no re-derivation) -- see
`docs/superpowers/specs/2026-07-27-diatomic-ve-scattering-library-design.md`'s
parameter table. `N2` reproduces `validation/n2/config.json`'s
`reduced_mass`/`impulsemomentum`/`potential` block, and
`projects/n2_2d_cross_section/hamiltonian2d.py`'s `MU`/`ELL`
(`libs/qscat/tests/test_model.py` gates this). This is the single place a
new molecule is added to the registry.
"""

from __future__ import annotations

from .diatomic import DiatomicResonanceModel
from .ionic import IonicResonanceModel

__all__ = ["N2", "NO", "F2", "H2P"]

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
