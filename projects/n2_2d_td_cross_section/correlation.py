"""Re-export of the promoted, model-agnostic `qscat.core.correlation`
(sub-project #A, Task 5) -- see there for the full docstring (the
Tannor-Weeks eta deconvolution factors, `outgoing_channel`, and why `F_out`
is the outgoing Hankel half rather than the regular free function).

Kept as a module (not deleted) so existing callers/imports in this project
(and its tests) are unaffected by the move; no new physics or numerics live
here.
"""

from __future__ import annotations

from qscat.core.correlation import eta_incident, eta_outgoing, outgoing_channel

__all__ = ["outgoing_channel", "eta_incident", "eta_outgoing"]
