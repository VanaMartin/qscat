"""Asymptotic channel functions for the exact 2-D e-N2 scattering problem.

The entrance/exit channel of a VE transition is a free electron of momentum
`k` in partial wave `l`, times a neutral vibrational state. The electronic
factor is the ENERGY-NORMALIZED regular free radial solution
`riccati_bessel_en`, promoted to `qscat.special` -- see its docstring there
for the full derivation/eMoScat provenance. Re-exported here so existing
callers in this project are unaffected by the move.
"""

from __future__ import annotations

from qscat.special import riccati_bessel_en

__all__ = ["riccati_bessel_en"]
