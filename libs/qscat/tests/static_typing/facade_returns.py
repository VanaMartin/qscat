"""What a checker infers from `ScatteringProblem`'s declared return types.

Never imported or executed -- `test_static_types.py` runs mypy over this
directory. A runtime test cannot stand in: the object a call returns is the
same whatever the annotation says, so only a checker sees the declared union,
and only a checker can tell whether the exported aliases are usable in a
caller's own annotations. That second half is what shipping `py.typed`
promises, and nothing else in the suite exercises it.
"""

from __future__ import annotations

from typing import assert_type

from qscat.core.problem import (
    Amplitude,
    CrossSection,
    DrCrossSection,
    ScatteringProblem,
    Wavefunction,
)


def _returns(prob: ScatteringProblem) -> None:
    """One signature per method, returning the union of the shapes its flag selects."""
    assert_type(prob.ve_cross_section([1], 0.1), CrossSection | tuple[CrossSection, Wavefunction])
    assert_type(prob.da_cross_section(0.1), CrossSection | tuple[CrossSection, Wavefunction])
    assert_type(prob.dr_cross_section(0.1), DrCrossSection)


def _annotatable(sigma: CrossSection, psi: Wavefunction, amp: Amplitude) -> DrCrossSection:
    """The aliases work in a caller's own annotations, which is why they are exported."""
    return sigma, psi, amp
