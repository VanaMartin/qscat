"""`Extractor` implementations for the shared propagate-once TD engine
(`qscat.core.time_dependent.propagate`).

`TannorWeeks` is the current (and, as of this task, only) energy-extraction
route: it records `c_{v'}(t_n) = c_product(Phi_{v'}, Psi(t_n))` every step
(the SAME per-step bookkeeping `time_dependent.propagate`'s legacy
`out_channels` path performs) and its `sigma(E)` reproduces the existing
Tannor-Weeks transform (`time_dependent.sigma_from_correlations` /
`_sigma_one_energy`) VERBATIM -- eta_out/eta_in deconvolution, the elastic
free-reference subtraction, `sigma = pi*|S - S_ref|^2 / (2E)` -- by building
an ad hoc `PropagationResult` from its own recorded series and delegating to
`sigma_from_correlations` (least churn: no transform logic is duplicated
here). `td_ve_cross_section(method="tw")` builds a `TannorWeeks` extractor
for the full propagation and, when `subtract_free_reference` applies, a
second one for the `V_int=0` free-reference propagation, then calls
`full.sigma(E, free=free_extractor)` -- see `time_dependent.py`'s module
docstring for the full physics writeup.

Delta/flow extractors (sub-project tasks 2-3) will live alongside
`TannorWeeks` here, each implementing the same `Extractor` protocol
(`record`/`sigma`) so `propagate` can drive any combination of them from one
shared trajectory.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

from qscat.dvr import TensorGrid
from qscat.linalg import c_product

from .correlation import outgoing_channel
from .time_dependent import PropagationResult, sigma_from_correlations

if TYPE_CHECKING:
    from qscat.model import ResonanceModel

__all__ = ["TannorWeeks"]

_WpIn = dict[str, float]
_WpOut = dict[str, float]


class TannorWeeks:
    """The Tannor-Weeks `Extractor`: records `c_{v'}(t)`, transforms via the
    existing eta-deconvolution + elastic free-reference logic.

    Construct one instance per propagation run (the full run, and -- when
    `subtract_free_reference` applies -- a second instance for the `V_int=0`
    free-reference run); pass it to `time_dependent.propagate(...,
    extractors=[tw])`. `sigma(E, free=...)` accepts the companion free-
    reference `TannorWeeks` (built from ITS OWN, separately propagated,
    run) so the elastic channel can subtract `S_free(E)` exactly as
    `_sigma_one_energy` documents -- this is the "second TannorWeeks" the
    `time_dependent.propagate` design notes describe.
    """

    def __init__(
        self,
        tgrid: TensorGrid,
        model: ResonanceModel,
        eps: npt.NDArray[np.float64],
        chi: npt.NDArray[np.complex128],
        v_init: int,
        vprimes: list[int],
        wp_out: _WpOut,
        *,
        wp_in: _WpIn,
        dt: float,
    ) -> None:
        self._tgrid = tgrid
        self._model = model
        self._eps = eps
        self._v_init = v_init
        self._vprimes = vprimes
        self._wp_in = wp_in
        self._wp_out = wp_out
        self._dt = dt
        self._out_channels = [outgoing_channel(tgrid, chi[vp], **wp_out) for vp in vprimes]
        self._rows: list[npt.NDArray[np.complex128]] = []

    def record(self, psi: npt.NDArray[np.complex128]) -> None:
        """Append this step's `c_{v'}(t_n) = c_product(Phi_{v'}, psi)` row.

        Same per-channel c-product, same channel order, as the legacy
        `out_channels` loop in `time_dependent.propagate` -- so a run driven
        purely through this extractor reproduces `PropagationResult.c` to
        machine precision.
        """
        row = np.empty(len(self._out_channels), dtype=np.complex128)
        for k, ch in enumerate(self._out_channels):
            row[k] = c_product(ch, psi)
        self._rows.append(row)

    @property
    def result(self) -> PropagationResult:
        """The recorded series as a `PropagationResult` (`norm`/`snapshots`
        are not tracked by this extractor and are left empty/zero -- only
        `t`/`c` feed `sigma_from_correlations`)."""
        n_t = len(self._rows)
        n_ch = len(self._out_channels)
        c = np.stack(self._rows) if self._rows else np.zeros((0, n_ch), dtype=np.complex128)
        t = np.arange(n_t, dtype=np.float64) * self._dt
        norm = np.zeros(n_t, dtype=np.float64)
        return PropagationResult(t=t, c=c, norm=norm, snapshots=[])

    def sigma(
        self, E: float | npt.ArrayLike, *, free: TannorWeeks | None = None
    ) -> npt.NDArray[np.float64]:
        """`sigma_{v_init->v'}(E)` (bohr^2) via the Tannor-Weeks transform.

        `free`, when given, is a companion `TannorWeeks` extractor recorded
        from a SEPARATE `V_int=0` propagation (same wavepacket/grid); its
        series supplies `sigma_from_correlations`'s `free_result` -- the
        elastic (`v'==v_init`) channel then subtracts the free-particle
        `S_free(E)` instead of a literal 1 (see
        `time_dependent._sigma_one_energy`). Leave `None` to reproduce the
        literal-1 fallback.
        """
        free_result = free.result if free is not None else None
        return sigma_from_correlations(
            self._tgrid,
            self._model,
            self.result,
            self._eps,
            self._v_init,
            self._vprimes,
            E,
            dt=self._dt,
            wp_in=self._wp_in,
            wp_out=self._wp_out,
            free_result=free_result,
        )
