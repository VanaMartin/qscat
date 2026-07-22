"""Crank-Nicolson time propagator -- thin re-export of the promoted primitive.

The general Crank-Nicolson stepper (Cayley-form time propagator for the
time-dependent Schrodinger equation, valid for general non-Hermitian `H`)
has no N2-specific or FEM-DVR-ECS-specific structure, so it was promoted to
`qscat.evolution.make_cn_stepper` (Task 3 of this sub-project). This module
keeps the project's original import path (`projects.n2_td_cross_section.
propagator.make_cn_stepper`) working as a thin re-export so existing
callers/tests are undisturbed.
"""

from __future__ import annotations

from qscat.evolution import make_cn_stepper

__all__ = ["make_cn_stepper"]
