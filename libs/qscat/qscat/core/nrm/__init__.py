"""Nonlocal resonance model (NRM) -- the Feshbach projection-operator
approximation to the 2-D electron-diatomic model.

The method is Houfek, Rescigno & McCurdy, Phys. Rev. A 77, 012710 (2008)
(`reference/literature/houfek-2008-pra77-012710.md`); every equation number in
this package refers to that paper. It sits between `qscat.core.lcp` (the local
complex potential, which discards the energy dependence and the nonlocality)
and `qscat.core.driven`/`qscat.core.dissociation` (the exact 2-D solver, the
oracle). See `docs/physics/nonlocal-resonance-model.md`.

NAMING: this package's `v_d_discrete` is the paper's discrete-state potential
`V_d(R) = V_0(R) + <phi_d|H_el|phi_d>` (Eq. 20). It is NOT `qscat.core.lcp`'s
`Vd`, which is the real part of the LCP curve. The paper states the two only
"almost coincide", and only for the physical discrete state (p. 012710-8).

`reference/eMoScat`'s `module_NRM.cpp` is NOT a reference for this package. The
NRM was never delivered as a working capability there; nothing in it is
treated as correct and no code here is derived from it.
"""

from __future__ import annotations

from .discrete_state import (
    AsymptoticDiscreteState,
    DiscreteState,
    PhysicalDiscreteState,
)
from .dissociation import nrm_da_cross_section
from .ingredients import NrmIngredients, nrm_ingredients

__all__ = [
    "AsymptoticDiscreteState",
    "DiscreteState",
    "NrmIngredients",
    "PhysicalDiscreteState",
    "nrm_da_cross_section",
    "nrm_ingredients",
]
