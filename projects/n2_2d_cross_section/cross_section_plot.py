"""Re-export of the promoted, generic `qscat.core.plot` (sub-project #A,
Task 5) -- see there for the full docstring (`plot_cross_sections`; no
physics, no Houfek, no N2 baked in).

Kept as a module (not deleted) so existing callers/imports in this project
(and `validation/n2/ti_curve.py`) are unaffected by the move; no new physics
or numerics live here.
"""

from __future__ import annotations

from qscat.core.plot import plot_cross_sections

__all__ = ["plot_cross_sections"]
