"""Static-typing fixture: `ScatteringProblem`'s detail flags narrow by literal.

This module is NEVER executed. `test_static_types.py` runs a type checker over
it; `typing.assert_type` fails there if an inferred type drifts. Narrowing is
not observable at run time -- `return_wavefunction=True` and a runtime `bool`
that happens to be `True` produce the identical object -- so only a type
checker can tell the two apart, and only these fixtures make it do so.

Three call shapes are exercised at every flagged surface:

  * a literal `True`   -> the detailed tuple,
  * a literal `False`, and the flag omitted -> the bare cross section,
  * a non-literal `bool` -> the honest union, which the caller must
    discriminate. That last one is the correctness check on the overloads: a
    set of overloads that narrowed a runtime `bool` would be lying.

The wavefunction half of every tuple stays OPTIONAL (`_Psi` includes `None`),
because a channel that is closed at an energy has no wavefunction there. The
overloads narrow which SHAPE comes back, never whether the physics produced a
wavefunction.
"""

from __future__ import annotations

from typing import assert_type

import numpy as np
import numpy.typing as npt
from qscat.core import ResonanceLevels, ScatteringProblem
from qscat.dvr import FemDvrEcsGrid

# The conventions `qscat.core.problem` names privately, restated here so the
# assertions below read as the documented return contract rather than as a
# reference to a private alias.
Sigma = npt.NDArray[np.float64]
Psi = npt.NDArray[np.complex128] | None
PsiOut = Psi | list[Psi]
Amp = npt.NDArray[np.complex128]
Curve = tuple[ResonanceLevels, npt.NDArray[np.complex128], npt.NDArray[np.float64]]


def check_ve_cross_section(prob: ScatteringProblem, E: float, flag: bool) -> None:
    assert_type(prob.ve_cross_section([0, 1], E), Sigma)
    assert_type(prob.ve_cross_section([0, 1], E, return_wavefunction=False), Sigma)
    assert_type(prob.ve_cross_section([0, 1], E, return_wavefunction=True), tuple[Sigma, PsiOut])
    assert_type(
        prob.ve_cross_section([0, 1], E, return_wavefunction=flag),
        Sigma | tuple[Sigma, PsiOut],
    )

    # Plain annotated bindings: these are only accepted because the literal
    # narrowed the return. Before the overloads, both were union-typed and
    # neither assignment was legal.
    sigma: Sigma = prob.ve_cross_section([0, 1], E)
    pair: tuple[Sigma, PsiOut] = prob.ve_cross_section([0, 1], E, return_wavefunction=True)
    del sigma, pair

    # ...and the runtime flag must NOT narrow. The suppression IS the
    # assertion: strict mode reports an UNUSED ignore, failing this fixture,
    # if the bool catch-all ever starts claiming the bare array.
    from_bool: Sigma = prob.ve_cross_section(  # type: ignore[assignment]
        [0, 1], E, return_wavefunction=flag
    )
    del from_bool


def check_da_cross_section(prob: ScatteringProblem, E: float, flag: bool) -> None:
    assert_type(prob.da_cross_section(E), Sigma)
    assert_type(prob.da_cross_section(E, return_wavefunction=False), Sigma)
    assert_type(prob.da_cross_section(E, return_wavefunction=True), tuple[Sigma, PsiOut])
    assert_type(prob.da_cross_section(E, return_wavefunction=flag), Sigma | tuple[Sigma, PsiOut])

    from_bool: Sigma = prob.da_cross_section(E, return_wavefunction=flag)  # type: ignore[assignment]
    del from_bool


def check_dr_cross_section(prob: ScatteringProblem, E: float, flag: bool) -> None:
    # Two independent flags, so four literal combinations. This method is the
    # only DR route that takes them at all: `qscat.core.dr_cross_section` is
    # sigma-only, and `dr_solve` returns one `DrResult` whose `psi` and
    # `amplitude` are Optional whatever was asked for.
    assert_type(prob.dr_cross_section(E), Sigma)
    assert_type(prob.dr_cross_section(E, return_wavefunction=True), tuple[Sigma, PsiOut])
    assert_type(prob.dr_cross_section(E, return_amplitude=True), tuple[Sigma, Amp])
    assert_type(
        prob.dr_cross_section(E, return_wavefunction=True, return_amplitude=True),
        tuple[Sigma, PsiOut, Amp],
    )
    assert_type(
        prob.dr_cross_section(E, return_wavefunction=flag, return_amplitude=flag),
        Sigma | tuple[Sigma, PsiOut] | tuple[Sigma, Amp] | tuple[Sigma, PsiOut, Amp],
    )

    # The amplitude is non-Optional under a literal flag -- the solver stores
    # it exactly when asked -- while the wavefunction stays Optional, since a
    # closed channel has none.
    _sigma, amplitude = prob.dr_cross_section(E, return_amplitude=True)
    assert_type(amplitude, Amp)

    from_bool: Sigma = prob.dr_cross_section(E, return_amplitude=flag)  # type: ignore[assignment]
    del from_bool


def check_lcp_cross_sections(
    prob: ScatteringProblem,
    E: float,
    Vd: npt.NDArray[np.complex128],
    Gamma: npt.NDArray[np.float64],
    flag: bool,
) -> None:
    assert_type(prob.lcp_da_cross_section(E, Vd=Vd, Gamma=Gamma), Sigma)
    assert_type(
        prob.lcp_da_cross_section(E, Vd=Vd, Gamma=Gamma, return_wavefunction=True),
        tuple[Sigma, PsiOut],
    )
    assert_type(
        prob.lcp_da_cross_section(E, Vd=Vd, Gamma=Gamma, return_wavefunction=flag),
        Sigma | tuple[Sigma, PsiOut],
    )

    assert_type(prob.lcp_ve_cross_section([0, 1], E, Vd=Vd, Gamma=Gamma), Sigma)
    assert_type(
        prob.lcp_ve_cross_section([0, 1], E, Vd=Vd, Gamma=Gamma, return_wavefunction=False),
        Sigma,
    )
    assert_type(
        prob.lcp_ve_cross_section([0, 1], E, Vd=Vd, Gamma=Gamma, return_wavefunction=True),
        tuple[Sigma, PsiOut],
    )
    assert_type(
        prob.lcp_ve_cross_section([0, 1], E, Vd=Vd, Gamma=Gamma, return_wavefunction=flag),
        Sigma | tuple[Sigma, PsiOut],
    )

    from_bool: Sigma = prob.lcp_ve_cross_section(  # type: ignore[assignment]
        [0, 1], E, Vd=Vd, Gamma=Gamma, return_wavefunction=flag
    )
    del from_bool


def check_resonance_levels(
    prob: ScatteringProblem,
    nuclear_grid_b: FemDvrEcsGrid,
    elec_grid_b: FemDvrEcsGrid,
    flag: bool,
) -> None:
    assert_type(prob.resonance_levels(nuclear_grid_b, elec_grid_b), ResonanceLevels)
    assert_type(
        prob.resonance_levels(nuclear_grid_b, elec_grid_b, return_curve=False),
        ResonanceLevels,
    )
    assert_type(prob.resonance_levels(nuclear_grid_b, elec_grid_b, return_curve=True), Curve)
    assert_type(
        prob.resonance_levels(nuclear_grid_b, elec_grid_b, return_curve=flag),
        ResonanceLevels | Curve,
    )

    # `return_curve=True` is how a caller gets the `(Vd, Gamma)` pair the LCP
    # cross sections take, so this unpacking has to be legal without a cast.
    levels, Vd, Gamma = prob.resonance_levels(nuclear_grid_b, elec_grid_b, return_curve=True)
    assert_type(levels, ResonanceLevels)
    assert_type(Vd, npt.NDArray[np.complex128])
    assert_type(Gamma, npt.NDArray[np.float64])

    from_bool: ResonanceLevels = prob.resonance_levels(  # type: ignore[assignment]
        nuclear_grid_b, elec_grid_b, return_curve=flag
    )
    del from_bool
