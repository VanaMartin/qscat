"""Re-export of the promoted, model-agnostic `qscat.core.wavepacket`
(sub-project #A, Task 5) -- see there for the full docstring
(`gaussian_coeffs`/`initial_state` formulas and conventions).

Kept as a module (not deleted) so existing callers/imports in this project
(and its tests) are unaffected by the move; no new physics or numerics live
here.
"""

from __future__ import annotations

from qscat.core.wavepacket import gaussian_coeffs, initial_state

__all__ = ["gaussian_coeffs", "initial_state"]
