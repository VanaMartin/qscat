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

`Dirac` (this task) is the delta-distribution sibling: eMoScat's
`DiracTestFunction2d` -- "Tannor-Weeks with a delta-distribution test
function instead of the Gaussian test packet". Its `record` is a LINE
PROJECTION of `psi` onto `chi_{v'}` at a FIXED electronic DVR index
`position` (no outgoing wavepacket, no propagation against `Phi_{v'}` --
the projection is onto the ALREADY-propagated `psi` directly), converted
from a coefficient to a wavefunction VALUE via `/sqrt(w_r[position])` (the
same coefficient<->value convention the LCP boundary-flux code uses, see
`qscat.core.lcp`); its `sigma` is TW's transform with `eta_out_i ->
hankel_point_value(...)` (a delta test function's `F_out` is unintegrated,
so the deconvolution factor is the raw outgoing-Hankel-half VALUE at
`position`, not an overlap with a Gaussian) -- a small self-contained
transform (not routed through `sigma_from_correlations`, which is TW-
specific) so this task cannot perturb TW's byte-identical golden regression
(`test_td_extractors.py::test_tw_method_matches_prerefactor_golden_*`).

`Flux` (this task) is the flow (flux) sibling: eMoScat's `FluxTestFunction2d`
-- "the time-energy Fourier transform of the probability flux projected to
the outgoing state". Its `record` appends BOTH the value `b_{v'}(t) =
<chi_{v'}|psi(surface,.)>` (`Dirac`'s line projection, at a FIXED electronic
surface -- an element border past the interaction) AND its electronic-
coordinate derivative `d_{v'}(t) = <chi_{v'}| d/dr psi(surface,.)>`, using
the new `qscat.dvr.dvr_first_derivative_at_node` primitive applied along the
electronic axis before the nuclear c-product projection. Its `sigma` is the
Wronskian-like flux transform (module docstring of `correlation.
outgoing_surface_wave`): `S_i = -i/(2*mu_e*ifc_i) * sum_j w_j *
(conj(phi_out_i)*d_{v'}(t_j) - b_{v'}(t_j)*conj(dphi_out_i)) *
exp(i*E_tot*t_j)*dt`, `mu_e = 1` (electronic reduced mass, a.u.) -- again a
small self-contained transform (like `Dirac`'s), not routed through
`sigma_from_correlations`, so this task cannot perturb TW's byte-identical
golden regression either.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

from qscat.dvr import FemDvrEcsGrid, TensorGrid, dvr_first_derivative_at_node
from qscat.linalg import c_product

from .correlation import eta_incident, hankel_point_value, outgoing_channel, outgoing_surface_wave
from .time_dependent import PropagationResult, _quadrature_weights, sigma_from_correlations

if TYPE_CHECKING:
    from qscat.model import ResonanceModel

__all__ = ["TannorWeeks", "Dirac", "Flux"]

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


def _dirac_s_vector_one_energy(
    grid: FemDvrEcsGrid,
    model: ResonanceModel,
    result: PropagationResult,
    eps: npt.NDArray[np.float64],
    v_init: int,
    vprimes: list[int],
    z_position: float,
    E: float,
    dt: float,
    wp_in: _WpIn,
) -> npt.NDArray[np.complex128]:
    """The raw delta-transform S-matrix, `TannorWeeks._s_vector_one_energy`'s
    twin with `eta_out_i -> hankel_point_value(...)` (module docstring)."""
    S = np.zeros(len(vprimes), dtype=np.complex128)
    if E <= 0.0:
        return S
    weights = _quadrature_weights(result.t.size)
    e_tot = E + eps[v_init]
    k = float(np.sqrt(2.0 * E))
    eta_in = eta_incident(grid, k, model.ell, **wp_in)
    phase = np.exp(1j * e_tot * result.t)
    for j, vp in enumerate(vprimes):
        excess = e_tot - eps[vp]
        if excess <= 0.0:
            continue  # closed channel
        kp = float(np.sqrt(2.0 * excess))
        f_i = hankel_point_value(grid, z_position, kp, model.ell, model.charge)
        s_raw = np.sum(weights * phase * result.c[:, j]) * dt
        S[j] = s_raw / (2.0 * np.pi * np.conj(f_i) * eta_in)
    return S


def _dirac_sigma_one_energy(
    grid: FemDvrEcsGrid,
    model: ResonanceModel,
    result: PropagationResult,
    eps: npt.NDArray[np.float64],
    v_init: int,
    vprimes: list[int],
    z_position: float,
    E: float,
    dt: float,
    wp_in: _WpIn,
    free_result: PropagationResult | None,
) -> npt.NDArray[np.float64]:
    """`sigma_{v_init->v'}(E)` (bohr^2) via the delta transform, one energy.

    Same elastic free-reference pattern as `time_dependent._sigma_one_energy`
    (see that function's docstring): `S_free(E)` from a companion `V_int=0`
    propagation subtracts on the diagonal (`v'==v_init`) channel instead of a
    literal 1.
    """
    sigma = np.zeros(len(vprimes), dtype=np.float64)
    if E <= 0.0:
        return sigma
    s_full = _dirac_s_vector_one_energy(
        grid, model, result, eps, v_init, vprimes, z_position, E, dt, wp_in
    )
    s_free = None
    if free_result is not None:
        s_free = _dirac_s_vector_one_energy(
            grid, model, free_result, eps, v_init, vprimes, z_position, E, dt, wp_in
        )
    e_tot = E + eps[v_init]
    for j, vp in enumerate(vprimes):
        if e_tot - eps[vp] <= 0.0:
            continue  # closed channel
        if vp == v_init:
            ref = complex(s_free[j]) if s_free is not None else 1.0 + 0.0j
        else:
            ref = 0.0 + 0.0j
        sigma[j] = np.pi * abs(s_full[j] - ref) ** 2 / (2.0 * E)
    return sigma


class Dirac:
    """The delta (Dirac) `Extractor`: eMoScat's `DiracTestFunction2d` --
    Tannor-Weeks with a delta-distribution test function instead of the
    Gaussian test packet.

    `record` projects `psi` onto `chi_{v'}` in the nuclear coordinate at a
    FIXED electronic DVR index `position` -- no propagated outgoing test
    function, unlike `TannorWeeks`. `position` must land in the real
    (unscaled) electronic region (`grid.real_points[position] <= grid.R0`),
    typically an element end well past the interaction range (mirroring
    `TannorWeeks`'s `wp_out` standoff).

    Construct one instance per propagation run (the full run, and -- when an
    elastic free reference is wanted -- a second instance for a SEPARATE
    `V_int=0` run); pass it to `time_dependent.propagate(..., extractors=
    [dirac])`. `sigma(E, free=...)` accepts the companion free-reference
    `Dirac` (built from ITS OWN, separately propagated, run), matching
    `TannorWeeks.sigma`'s `free` contract.
    """

    def __init__(
        self,
        tgrid: TensorGrid,
        model: ResonanceModel,
        eps: npt.NDArray[np.float64],
        chi: npt.NDArray[np.complex128],
        v_init: int,
        vprimes: list[int],
        position: int,
        *,
        wp_in: _WpIn,
        dt: float,
    ) -> None:
        grid = tgrid.grids[0]
        if not (0 <= position < grid.n):
            raise ValueError(f"position {position} out of range for grid of size {grid.n}")
        if grid.real_points[position] > grid.R0:
            raise ValueError(
                f"position {position} (r={grid.real_points[position]}) is not in the real "
                f"(unscaled) electronic region (R0={grid.R0}) -- pick an index with "
                "real_points[position] <= R0"
            )
        self._tgrid = tgrid
        self._model = model
        self._eps = eps
        self._v_init = v_init
        self._vprimes = vprimes
        self._wp_in = wp_in
        self._dt = dt
        self._position = position
        self._z_position = float(grid.real_points[position])
        self._inv_sqrt_w = 1.0 / np.sqrt(np.asarray(grid.weights[position], dtype=np.complex128))
        self._chi = [np.asarray(chi[vp], dtype=np.complex128) for vp in vprimes]
        self._rows: list[npt.NDArray[np.complex128]] = []

    def record(self, psi: npt.NDArray[np.complex128]) -> None:
        """Append this step's `b_{v'}(t_n) = <chi_{v'}|psi(position, .)> /
        sqrt(w_r[position])` row -- the DVR-coefficient-to-VALUE conversion
        `qscat.core.lcp`'s boundary flux also relies on."""
        block = psi.reshape(self._tgrid.shape)
        psi_row = block[self._position, :]
        row = np.empty(len(self._vprimes), dtype=np.complex128)
        for k, chi_vp in enumerate(self._chi):
            row[k] = c_product(chi_vp, psi_row) * self._inv_sqrt_w
        self._rows.append(row)

    @property
    def result(self) -> PropagationResult:
        """The recorded series as a `PropagationResult` (same shape/role as
        `TannorWeeks.result`; `norm`/`snapshots` are left empty/zero)."""
        n_t = len(self._rows)
        n_ch = len(self._vprimes)
        c = np.stack(self._rows) if self._rows else np.zeros((0, n_ch), dtype=np.complex128)
        t = np.arange(n_t, dtype=np.float64) * self._dt
        norm = np.zeros(n_t, dtype=np.float64)
        return PropagationResult(t=t, c=c, norm=norm, snapshots=[])

    def sigma(
        self, E: float | npt.ArrayLike, *, free: Dirac | None = None
    ) -> npt.NDArray[np.float64]:
        """`sigma_{v_init->v'}(E)` (bohr^2) via the delta transform (module
        docstring): TW's transform with `eta_out_i -> hankel_point_value(...)`
        at `position`.

        `free`, when given, is a companion `Dirac` extractor recorded from a
        SEPARATE `V_int=0` propagation -- same elastic free-reference
        contract as `TannorWeeks.sigma`.
        """
        free_result = free.result if free is not None else None
        e_arr = np.atleast_1d(np.asarray(E, dtype=np.float64))
        result = self.result
        grid = self._tgrid.grids[0]
        out = np.stack(
            [
                _dirac_sigma_one_energy(
                    grid,
                    self._model,
                    result,
                    self._eps,
                    self._v_init,
                    self._vprimes,
                    self._z_position,
                    float(e),
                    self._dt,
                    self._wp_in,
                    free_result,
                )
                for e in e_arr
            ]
        )
        scalar = np.isscalar(E) or (isinstance(E, np.ndarray) and np.ndim(E) == 0)
        if scalar:
            return np.asarray(out[0], dtype=np.float64)
        return out


def _flux_s_vector_one_energy(
    grid: FemDvrEcsGrid,
    model: ResonanceModel,
    t: npt.NDArray[np.float64],
    b: npt.NDArray[np.complex128],
    d: npt.NDArray[np.complex128],
    eps: npt.NDArray[np.float64],
    v_init: int,
    vprimes: list[int],
    z_surface: float,
    E: float,
    dt: float,
    wp_in: _WpIn,
) -> npt.NDArray[np.complex128]:
    """The raw flux-transform S-matrix (module docstring's Wronskian formula):

        S_i = -i/(2*mu_e*ifc_i) * sum_j w_j *
              (conj(phi_out_i)*d_j - b_j*conj(dphi_out_i)) * exp(i*E_tot*t_j) * dt

    `mu_e = 1.0` (electronic reduced mass, a.u.); `phi_out_i`/`dphi_out_i` are
    `outgoing_surface_wave`'s pair at the channel's outgoing momentum `k' =
    sqrt(2*(E_tot - eps[v']))` -- the SAME per-channel momentum `TannorWeeks`/
    `Dirac` use for their own outgoing deconvolution factor.
    """
    S = np.zeros(len(vprimes), dtype=np.complex128)
    if E <= 0.0:
        return S
    mu_e = 1.0
    weights = _quadrature_weights(t.size)
    e_tot = E + eps[v_init]
    k = float(np.sqrt(2.0 * E))
    eta_in = eta_incident(grid, k, model.ell, **wp_in)
    phase = np.exp(1j * e_tot * t)
    for j, vp in enumerate(vprimes):
        excess = e_tot - eps[vp]
        if excess <= 0.0:
            continue  # closed channel
        kp = float(np.sqrt(2.0 * excess))
        phi_out, dphi_out = outgoing_surface_wave(grid, z_surface, kp, model.ell, model.charge)
        wronskian = np.conj(phi_out) * d[:, j] - b[:, j] * np.conj(dphi_out)
        s_raw = np.sum(weights * wronskian * phase) * dt
        S[j] = (-1j / (2.0 * mu_e * eta_in)) * s_raw
    return S


def _flux_sigma_one_energy(
    grid: FemDvrEcsGrid,
    model: ResonanceModel,
    t: npt.NDArray[np.float64],
    b: npt.NDArray[np.complex128],
    d: npt.NDArray[np.complex128],
    eps: npt.NDArray[np.float64],
    v_init: int,
    vprimes: list[int],
    z_surface: float,
    E: float,
    dt: float,
    wp_in: _WpIn,
    free: tuple[npt.NDArray[np.float64], npt.NDArray[np.complex128], npt.NDArray[np.complex128]]
    | None,
) -> npt.NDArray[np.float64]:
    """`sigma_{v_init->v'}(E)` (bohr^2) via the flux transform, one energy --
    same elastic free-reference pattern as `time_dependent._sigma_one_energy`
    / `_dirac_sigma_one_energy`."""
    sigma = np.zeros(len(vprimes), dtype=np.float64)
    if E <= 0.0:
        return sigma
    s_full = _flux_s_vector_one_energy(
        grid, model, t, b, d, eps, v_init, vprimes, z_surface, E, dt, wp_in
    )
    s_free = None
    if free is not None:
        t_free, b_free, d_free = free
        s_free = _flux_s_vector_one_energy(
            grid, model, t_free, b_free, d_free, eps, v_init, vprimes, z_surface, E, dt, wp_in
        )
    e_tot = E + eps[v_init]
    for j, vp in enumerate(vprimes):
        if e_tot - eps[vp] <= 0.0:
            continue  # closed channel
        if vp == v_init:
            ref = complex(s_free[j]) if s_free is not None else 1.0 + 0.0j
        else:
            ref = 0.0 + 0.0j
        sigma[j] = np.pi * abs(s_full[j] - ref) ** 2 / (2.0 * E)
    return sigma


class Flux:
    """The flow (flux) `Extractor`: eMoScat's `FluxTestFunction2d` -- the
    time-energy Fourier transform of the probability flux projected onto the
    outgoing channel, at a FIXED electronic surface.

    `record` appends, per `v'`, BOTH the value `b_{v'}(t) = <chi_{v'}|
    psi(surface,.)>` (`Dirac`'s line projection) AND its electronic-
    coordinate derivative `d_{v'}(t) = <chi_{v'}| d/dr psi(surface,.)>`
    (via `qscat.dvr.dvr_first_derivative_at_node` applied along the
    electronic axis, then a nuclear c-product onto `chi_{v'}` -- no extra
    `/sqrt(w)` needed there, unlike `b_{v'}`: `dvr_first_derivative_at_node`
    already converts coefficient->value internally, see its docstring).
    `sigma(E, free=...)` is the Wronskian-like flux transform (module
    docstring); same elastic free-reference contract as
    `TannorWeeks.sigma`/`Dirac.sigma`.

    `surface` must be a real (unscaled) electronic DVR index in the
    asymptotic region (past the interaction), same requirement as `Dirac`'s
    `position`.
    """

    def __init__(
        self,
        tgrid: TensorGrid,
        model: ResonanceModel,
        eps: npt.NDArray[np.float64],
        chi: npt.NDArray[np.complex128],
        v_init: int,
        vprimes: list[int],
        surface: int,
        *,
        wp_in: _WpIn,
        dt: float,
    ) -> None:
        grid = tgrid.grids[0]
        if not (0 <= surface < grid.n):
            raise ValueError(f"surface {surface} out of range for grid of size {grid.n}")
        if grid.real_points[surface] > grid.R0:
            raise ValueError(
                f"surface {surface} (r={grid.real_points[surface]}) is not in the real "
                f"(unscaled) electronic region (R0={grid.R0}) -- pick an index with "
                "real_points[surface] <= R0"
            )
        self._tgrid = tgrid
        self._model = model
        self._eps = eps
        self._v_init = v_init
        self._vprimes = vprimes
        self._wp_in = wp_in
        self._dt = dt
        self._surface = surface
        self._z_surface = float(grid.real_points[surface])
        self._inv_sqrt_w = 1.0 / np.sqrt(np.asarray(grid.weights[surface], dtype=np.complex128))
        self._deriv_row = dvr_first_derivative_at_node(grid, surface)
        self._chi = [np.asarray(chi[vp], dtype=np.complex128) for vp in vprimes]
        self._b_rows: list[npt.NDArray[np.complex128]] = []
        self._d_rows: list[npt.NDArray[np.complex128]] = []

    def record(self, psi: npt.NDArray[np.complex128]) -> None:
        """Append this step's `b_{v'}(t_n)` and `d_{v'}(t_n)` rows (docstring)."""
        block = psi.reshape(self._tgrid.shape)
        psi_row = block[self._surface, :]
        dpsi_row = self._deriv_row @ block  # (n_nuclear,): d/dr psi(surface, R), coeff-in-R
        row_b = np.empty(len(self._vprimes), dtype=np.complex128)
        row_d = np.empty(len(self._vprimes), dtype=np.complex128)
        for k, chi_vp in enumerate(self._chi):
            row_b[k] = c_product(chi_vp, psi_row) * self._inv_sqrt_w
            row_d[k] = c_product(chi_vp, dpsi_row)
        self._b_rows.append(row_b)
        self._d_rows.append(row_d)

    @property
    def _arrays(
        self,
    ) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.complex128], npt.NDArray[np.complex128]]:
        n_t = len(self._b_rows)
        n_ch = len(self._vprimes)
        b = np.stack(self._b_rows) if self._b_rows else np.zeros((0, n_ch), dtype=np.complex128)
        d = np.stack(self._d_rows) if self._d_rows else np.zeros((0, n_ch), dtype=np.complex128)
        t = np.arange(n_t, dtype=np.float64) * self._dt
        return t, b, d

    def sigma(
        self, E: float | npt.ArrayLike, *, free: Flux | None = None
    ) -> npt.NDArray[np.float64]:
        """`sigma_{v_init->v'}(E)` (bohr^2) via the flux transform (module
        docstring). `free`, when given, is a companion `Flux` extractor
        recorded from a SEPARATE `V_int=0` propagation -- same elastic
        free-reference contract as `TannorWeeks.sigma`/`Dirac.sigma`."""
        free_arrays = free._arrays if free is not None else None
        e_arr = np.atleast_1d(np.asarray(E, dtype=np.float64))
        t, b, d = self._arrays
        grid = self._tgrid.grids[0]
        out = np.stack(
            [
                _flux_sigma_one_energy(
                    grid,
                    self._model,
                    t,
                    b,
                    d,
                    self._eps,
                    self._v_init,
                    self._vprimes,
                    self._z_surface,
                    float(e),
                    self._dt,
                    self._wp_in,
                    free_arrays,
                )
                for e in e_arr
            ]
        )
        scalar = np.isscalar(E) or (isinstance(E, np.ndarray) and np.ndim(E) == 0)
        if scalar:
            return np.asarray(out[0], dtype=np.float64)
        return out
