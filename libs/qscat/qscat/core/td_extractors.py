"""`Extractor` implementations for the shared propagate-once TD engine
(`qscat.core.time_dependent.propagate`).

`TannorWeeks` (Tannor & Weeks, J. Chem. Phys. 98, 3884 (1993)) records
`c_{v'}(t_n) = c_product(Phi_{v'}, Psi(t_n))` every step
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

`Dirac` is the delta-distribution sibling: eMoScat's
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
specific) so nothing here can perturb TW's byte-identical golden regression
(`test_td_extractors.py::test_tw_method_matches_prerefactor_golden_*`).

`Flux` is the flow (flux) sibling: eMoScat's `FluxTestFunction2d`
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
`sigma_from_correlations`, so it cannot perturb TW's byte-identical golden
regression either.

`Flux(axis="nuclear")` is the SAME flow
extractor with the electronic/nuclear roles of eMoScat's `FluxTestFunction2d`
`axis_=='y'` branch: the surface is a NUCLEAR node (`R = R_surface`, real
region), and the "bound state" projected onto is not a single nuclear
vibrational level but one of `n_channels` ANION ELECTRONIC states at the
dissociation limit `R_inf = tgrid.grids[1].R0`
(`qscat.core.dissociation.anion_electronic_states`) -- the DIRECT nuclear-
axis analog of `b_{v'}(t)`/`d_{v'}(t)`, transposed: `b_c(t) = <phi_c|
psi(.,R=surface)>` (a c-product on the ELECTRONIC axis, at fixed nuclear
node `surface`) and `d_c(t) = <phi_c| d/dR psi(.,R=surface)>` (via
`dvr_first_derivative_at_node` applied along the NUCLEAR axis, contracted
into the electronic-axis vector BEFORE projecting onto `phi_c` -- bilinear,
non-conjugated `c_product` makes this order interchangeable with "project
then differentiate", eMoScat's literal order). The outgoing wave is the mass-`mu_R`
Hankel half (`correlation.outgoing_surface_wave(..., mass=model.mu)`,
`l=0`) at `K_R = sqrt(2 mu_R (E_tot - eps_e_c))`, `mu_R = model.mu`
(eMoScat's `reduced_mass()` for `axis_=='y'`: `mu_y_ = mass`, NOT 1) -- the
Wronskian flux prefactor is `1/(2*mu_R*ifc)`, `ifc = eta_incident` STILL on
the electronic incident axis (`tgrid.grids[0]`, the incident electron --
unchanged: only the OUTGOING side moves to the nuclear axis). `sigma_c(E) =
C_DA * |S_c|^2 / (2E)` per anion channel `c` (no elastic free-reference
subtraction -- DA is a pure rearrangement channel, there is no `v'==v_init`
diagonal to subtract a literal 1 from). `C_DA = pi`: the SAME constant the
electronic-axis `Flux`/`TannorWeeks` use for an INELASTIC channel (`sigma =
pi*|S|^2/2E`), not the TI DA oracle's literal `4 pi^3` prefactor -- because
the TI oracle's `4 pi^3 |T|^2/2E` (`qscat.core.dissociation.da_cross_section`,
`qscat.core.driven.ve_cross_section`'s own docstring) is ALREADY `pi*|S|^2/
2E` written in terms of the driven-equation T-matix via the general
partial-wave identity `S = 1 - 2 pi i T` (`|S|^2 = 4 pi^2 |T|^2` for an
off-diagonal/inelastic element, so `pi|S|^2/2E == 4 pi^3|T|^2/2E`): the flux
Wronskian transform is a channel-agnostic way to extract "the" S-matrix
element for ANY exit channel (electronic or nuclear), so it lands on the
SAME `pi*|S|^2/2E` convention regardless of which axis carries the outgoing
flux -- confirmed empirically by the TI-convergence gate (see
`docs/physics/td-da.md`), not merely asserted.

`Dirac(axis="nuclear")` is the delta (point-
projection) sibling of `Flux(axis="nuclear")` -- the SAME anion-channel
scaffolding (`n_channels`, `anion_electronic_states` at `R_inf =
tgrid.grids[1].R0`), but `record` keeps only the point VALUE, no derivative:
`b_c(t) = <phi_c|psi(.,R=position)> / sqrt(w_R[position])` (a c-product on
the ELECTRONIC axis at fixed nuclear node `position`) -- the direct nuclear-
axis transpose of the electronic `Dirac`'s `b_{v'}(t) = <chi_{v'}|
psi(position,.)> / sqrt(w_r[position])`. Its `sigma` is `_dirac_s_vector_
one_energy`'s nuclear-axis twin: `eta_out_c -> hankel_point_value(g_nuc,
R_position, K_R, l=0, model.charge, mass=mu_R)` (the point-VALUE half-Hankel,
not `outgoing_surface_wave`'s Wronskian pair -- a delta test function has no
derivative side), `K_R = sqrt(2 mu_R (E_tot - eps_e_c))`, `mu_R = model.mu`,
`eta_in = eta_incident` STILL on the electronic incident axis (unchanged --
only the outgoing side moves to the nuclear axis, exactly as `Flux(axis=
"nuclear")`'s docstring explains). `sigma_DA,c(E) = C_DA * |S_c|^2 / (2E)`,
the SAME `C_DA = pi` (no elastic free-reference subtraction: DA has no
`v'==v_init` diagonal, `free != None` raises `ValueError`, matching
`Flux(axis="nuclear")`).

`TannorWeeks(axis="nuclear")` is the NUCLEAR-
axis TannorWeeks: the DA sibling of the electronic `TannorWeeks` above, the
same `n_channels` anion-electronic-state scaffolding as `Flux(axis=
"nuclear")`/`Dirac(axis="nuclear")`, but keeping TW's own defining trait --
a PROPAGATED outgoing Gaussian test packet, not a fixed point/surface. The
outgoing test function moves to the NUCLEAR coordinate (`correlation.
outgoing_channel_nuclear`, the transpose of `outgoing_channel`):
`Phi_c = phi_c(r) g_out(R)`, `phi_c` one of the anion electronic bound
states (`R_inf = tgrid.grids[1].R0`) instead of a nuclear vibrational
level. `record` keeps the IDENTICAL per-channel c-product loop the
electronic path already uses (`c_c(t) = c_product(Phi_c, psi)`) -- only the
channel vectors differ, so this needed no change. `sigma(E)` is the DA
analog of `_s_vector_one_energy`/`_sigma_one_energy` (this module's own
`_tw_da_s_vector_one_energy`/`_tw_da_sigma_one_energy`, defined below, NOT
routed through `sigma_from_correlations`: that helper is hardwired to
electronic vibrational levels/axis, see its own docstring in
`time_dependent.py`): the outgoing deconvolution factor is `eta_outgoing`
on the NUCLEAR axis (`tgrid.grids[1]`, mass `mu_R = model.mu`, `l=0`,
`correlation.eta_outgoing`'s new `mass` keyword), evaluated against the
propagated `c_c(t)` series itself (unlike `Flux`'s Wronskian pair or
`Dirac`'s point value); the incident deconvolution `eta_in` stays on the
ELECTRONIC incident axis, unchanged. `sigma_DA,c(E) = C_DA * |S_c|^2 /
(2E)`, the SAME `C_DA = pi` convention, no elastic free-reference
subtraction (DA has no `v'==v_init` diagonal; `free != None` raises
`ValueError`, matching `Flux`/`Dirac(axis="nuclear")`).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

from qscat.dvr import FemDvrEcsGrid, TensorGrid, dvr_first_derivative_at_node
from qscat.linalg import c_product

from .correlation import (
    eta_incident,
    eta_outgoing,
    hankel_point_value,
    outgoing_channel,
    outgoing_channel_nuclear,
    outgoing_surface_wave,
)
from .dissociation import anion_electronic_states
from .time_dependent import (
    Extractor,
    PropagationResult,
    quadrature_weights,
    sigma_from_correlations,
)

if TYPE_CHECKING:
    from qscat.model import ResonanceModel

__all__ = ["Dirac", "Flux", "TannorWeeks"]

_WpIn = dict[str, float]
_WpOut = dict[str, float]

_AXES = ("electronic", "nuclear")

# The DA sigma prefactor: `pi`, NOT the TI oracle's literal `4 pi^3` -- see
# the module docstring's `Flux(axis="nuclear")` section for the `S = 1 -
# 2 pi i T` derivation, and `docs/physics/td-da.md` for the empirical
# TI-convergence confirmation.
_C_DA = np.pi


def _check_axis(axis: str, cls_name: str) -> None:
    """Validate `axis`. All three extractors implement both axes; anything
    outside `_AXES` is a `ValueError`."""
    if axis not in _AXES:
        raise ValueError(f"{cls_name}: axis must be one of {_AXES}, got {axis!r}")


def _axis_grid_index(axis: str) -> int:
    """`TensorGrid.grids` index for `axis` (`0` electronic, `1` nuclear)."""
    return 0 if axis == "electronic" else 1


def _tw_da_s_vector_one_energy(
    g_elec: FemDvrEcsGrid,
    g_nuc: FemDvrEcsGrid,
    model: ResonanceModel,
    mu_r: float,
    result: PropagationResult,
    eps: npt.NDArray[np.float64],
    v_init: int,
    eps_e: npt.NDArray[np.float64],
    wp_out: _WpOut,
    E: float,
    dt: float,
    wp_in: _WpIn,
) -> npt.NDArray[np.complex128]:
    """The raw DA Tannor-Weeks S-matrix, per anion dissociation channel --
    `_s_vector_one_energy`'s (`time_dependent.py`) nuclear-axis twin (module
    docstring's `TannorWeeks(axis="nuclear")` section): the outgoing
    deconvolution factor is `eta_outgoing` evaluated on the NUCLEAR axis
    (mass `mu_r`, `l=0`) against the PROPAGATED test-packet correlation
    `result.c[:, c]` itself (unlike `Flux`'s Wronskian pair or `Dirac`'s
    point value); the incident deconvolution `eta_in` stays on the
    ELECTRONIC incident axis.
    """
    n_channels = len(eps_e)
    S = np.zeros(n_channels, dtype=np.complex128)
    if E <= 0.0:
        return S
    weights = quadrature_weights(result.t.size)
    e_tot = E + eps[v_init]
    k = float(np.sqrt(2.0 * E))
    eta_in = eta_incident(g_elec, k, model.ell, **wp_in)
    phase = np.exp(1j * e_tot * result.t)
    for c in range(n_channels):
        e_dr = e_tot - eps_e[c]
        if e_dr <= 0.0:
            continue  # closed dissociation channel
        k_r = float(np.sqrt(2.0 * mu_r * e_dr))
        eta_out = eta_outgoing(g_nuc, k_r, 0, mass=mu_r, **wp_out)
        s_raw = np.sum(weights * phase * result.c[:, c]) * dt
        S[c] = s_raw / (2.0 * np.pi * np.conj(eta_out) * eta_in)
    return S


def _tw_da_sigma_one_energy(
    g_elec: FemDvrEcsGrid,
    g_nuc: FemDvrEcsGrid,
    model: ResonanceModel,
    mu_r: float,
    result: PropagationResult,
    eps: npt.NDArray[np.float64],
    v_init: int,
    eps_e: npt.NDArray[np.float64],
    wp_out: _WpOut,
    E: float,
    dt: float,
    wp_in: _WpIn,
) -> npt.NDArray[np.float64]:
    """`sigma_DA,c(E)` (bohr^2) per anion dissociation channel `c`, via the
    nuclear-axis Tannor-Weeks transform -- no elastic free-reference
    subtraction (DA is a pure rearrangement channel, no `v'==v_init`
    diagonal), the SAME `_C_DA = pi` convention `_flux_da_sigma_one_energy`/
    `_dirac_da_sigma_one_energy` use."""
    n_channels = len(eps_e)
    sigma = np.zeros(n_channels, dtype=np.float64)
    if E <= 0.0:
        return sigma
    s_full = _tw_da_s_vector_one_energy(
        g_elec, g_nuc, model, mu_r, result, eps, v_init, eps_e, wp_out, E, dt, wp_in
    )
    e_tot = E + eps[v_init]
    for c in range(n_channels):
        if e_tot - eps_e[c] <= 0.0:
            continue  # closed dissociation channel
        sigma[c] = _C_DA * abs(s_full[c]) ** 2 / (2.0 * E)
    return sigma


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

    Nuclear (`axis="nuclear"`): the DISSOCIATIVE ATTACHMENT (DA) sibling of
    `Flux(axis="nuclear")`/`Dirac(axis="nuclear")` (module docstring's
    `TannorWeeks(axis="nuclear")` section) -- `wp_out` is now the NUCLEAR
    outgoing test packet's parameters (`r0_out`/`p0_out`/`sigma_out` in `R`),
    `n_channels` selects how many anion electronic bound states
    (`qscat.core.dissociation.anion_electronic_states`, at `R_inf =
    tgrid.grids[1].R0`) are tracked as exit channels; `vprimes` is unused
    (pass `[]`). `record` is UNCHANGED (the same per-channel c-product loop
    -- only the channel test functions differ). `sigma(E)` returns the DA
    transform's `sigma_DA,c(E)` per channel, shape `(n_channels,)` for
    scalar `E` (matching `Flux`/`Dirac(axis="nuclear")`'s contract). `free`
    is not supported (DA has no elastic diagonal to subtract a reference
    from; passing it raises `ValueError`).
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
        axis: str = "electronic",
        n_channels: int = 1,
    ) -> None:
        _check_axis(axis, "TannorWeeks")
        self._axis = axis
        self._tgrid = tgrid
        self._model = model
        self._eps = eps
        self._v_init = v_init
        self._vprimes = vprimes
        self._wp_in = wp_in
        self._wp_out = wp_out
        self._dt = dt
        if axis == "electronic":
            self._out_channels = [outgoing_channel(tgrid, chi[vp], **wp_out) for vp in vprimes]
        else:  # nuclear
            eps_e, phi = anion_electronic_states(
                tgrid.grids[0], model, R_inf=tgrid.grids[1].R0, n_states=n_channels
            )
            self._n_channels = n_channels
            self._eps_e = eps_e
            self._out_channels = [
                outgoing_channel_nuclear(tgrid, phi[c], **wp_out) for c in range(n_channels)
            ]
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

    def _result(self, n_steps: int | None = None) -> PropagationResult:
        """`.result`, optionally truncated to the FIRST `n_steps` recorded
        samples (`self._rows[:n_steps]`). `n_steps=None` (the default) is the
        full series -- byte-identical to `.result`/pre-truncation behavior.
        The moment-resolved `cross_section_vs_time` runner artifact calls
        this (via `sigma`'s own `n_steps`) with `n_steps = round(t_i / dt) +
        1` to read sigma(E) as of an earlier time `t_i` without re-running
        the propagation.
        """
        rows = self._rows if n_steps is None else self._rows[:n_steps]
        n_t = len(rows)
        n_ch = len(self._out_channels)
        c = np.stack(rows) if rows else np.zeros((0, n_ch), dtype=np.complex128)
        t = np.arange(n_t, dtype=np.float64) * self._dt
        norm = np.zeros(n_t, dtype=np.float64)
        return PropagationResult(t=t, c=c, norm=norm, snapshots=[])

    @property
    def result(self) -> PropagationResult:
        """The recorded series as a `PropagationResult` (`norm`/`snapshots`
        are not tracked by this extractor and are left empty/zero -- only
        `t`/`c` feed `sigma_from_correlations` (electronic axis) or
        `_tw_da_sigma_one_energy` (nuclear axis)). The full series
        (`self._result(None)`); see `sigma`'s `n_steps` for a truncated read.
        """
        return self._result()

    def sigma(
        self,
        E: float | npt.ArrayLike,
        *,
        free: Extractor | None = None,
        n_steps: int | None = None,
    ) -> npt.NDArray[np.float64]:
        """`sigma_{v_init->v'}(E)` (bohr^2, electronic axis) or `sigma_DA,c(E)`
        (bohr^2 per anion channel `c`, nuclear axis) via the Tannor-Weeks
        transform (class docstring).

        Electronic: `free`, when given, must be a companion `TannorWeeks`
        extractor recorded from a SEPARATE `V_int=0` propagation (same
        wavepacket/grid); its series supplies `sigma_from_correlations`'s
        `free_result` -- the elastic (`v'==v_init`) channel then subtracts
        the free-particle `S_free(E)` instead of a literal 1 (see
        `time_dependent._sigma_one_energy`). Leave `None` to reproduce the
        literal-1 fallback. `free` is typed as the general `Extractor` protocol
        so all three extractors share one `sigma` signature, but only a
        `TannorWeeks` is meaningful here; anything else raises `TypeError`.

        Nuclear: DA has no elastic diagonal to subtract a reference from --
        `free` must be `None` (`ValueError` otherwise).

        `n_steps` (keyword-only, default `None`): transform only the FIRST
        `n_steps` recorded samples instead of the full series -- `None`
        reproduces today's behavior byte-identically (`self._rows[:None] ==
        self._rows[:]`). When `free` is also given, ITS series is truncated
        to the SAME `n_steps` (both runs share the same step schedule).
        """
        if self._axis == "nuclear":
            if free is not None:
                raise ValueError(
                    "TannorWeeks.sigma(axis='nuclear'): DA has no elastic free-reference "
                    "subtraction -- free must be None"
                )
            e_arr = np.atleast_1d(np.asarray(E, dtype=np.float64))
            result = self._result(n_steps)
            g_elec = self._tgrid.grids[0]
            g_nuc = self._tgrid.grids[1]
            mu_r = self._model.mu
            out = np.stack(
                [
                    _tw_da_sigma_one_energy(
                        g_elec,
                        g_nuc,
                        self._model,
                        mu_r,
                        result,
                        self._eps,
                        self._v_init,
                        self._eps_e,
                        self._wp_out,
                        float(e),
                        self._dt,
                        self._wp_in,
                    )
                    for e in e_arr
                ]
            )
            scalar = np.isscalar(E) or (isinstance(E, np.ndarray) and np.ndim(E) == 0)
            if scalar:
                return np.asarray(out[0], dtype=np.float64)
            return out

        free_result: PropagationResult | None = None
        if free is not None:
            if not isinstance(free, TannorWeeks):
                raise TypeError(f"TannorWeeks.sigma: free must be a TannorWeeks, got {type(free)}")
            free_result = free._result(n_steps)
        return sigma_from_correlations(
            self._tgrid,
            self._model,
            self._result(n_steps),
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
    weights = quadrature_weights(result.t.size)
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


def _dirac_da_s_vector_one_energy(
    g_elec: FemDvrEcsGrid,
    g_nuc: FemDvrEcsGrid,
    model: ResonanceModel,
    mu_r: float,
    result: PropagationResult,
    eps: npt.NDArray[np.float64],
    v_init: int,
    eps_e: npt.NDArray[np.float64],
    R_position: float,
    E: float,
    dt: float,
    wp_in: _WpIn,
) -> npt.NDArray[np.complex128]:
    """The raw DA delta-transform S-matrix, per anion dissociation channel --
    `_dirac_s_vector_one_energy`'s nuclear-axis twin (module docstring's
    `Flux(axis="nuclear")` section, adapted to a point VALUE rather than a
    Wronskian): the outgoing test function moves to the nuclear axis at mass
    `mu_r` (`l=0`, `hankel_point_value` instead of `outgoing_surface_wave`),
    the incident deconvolution `eta_in` stays on the ELECTRONIC incident axis.
    """
    n_channels = len(eps_e)
    S = np.zeros(n_channels, dtype=np.complex128)
    if E <= 0.0:
        return S
    weights = quadrature_weights(result.t.size)
    e_tot = E + eps[v_init]
    k = float(np.sqrt(2.0 * E))
    eta_in = eta_incident(g_elec, k, model.ell, **wp_in)
    phase = np.exp(1j * e_tot * result.t)
    for c in range(n_channels):
        e_dr = e_tot - eps_e[c]
        if e_dr <= 0.0:
            continue  # closed dissociation channel
        k_r = float(np.sqrt(2.0 * mu_r * e_dr))
        f_c = hankel_point_value(g_nuc, R_position, k_r, 0, model.charge, mass=mu_r)
        s_raw = np.sum(weights * phase * result.c[:, c]) * dt
        S[c] = s_raw / (2.0 * np.pi * np.conj(f_c) * eta_in)
    return S


def _dirac_da_sigma_one_energy(
    g_elec: FemDvrEcsGrid,
    g_nuc: FemDvrEcsGrid,
    model: ResonanceModel,
    mu_r: float,
    result: PropagationResult,
    eps: npt.NDArray[np.float64],
    v_init: int,
    eps_e: npt.NDArray[np.float64],
    R_position: float,
    E: float,
    dt: float,
    wp_in: _WpIn,
) -> npt.NDArray[np.float64]:
    """`sigma_DA,c(E)` (bohr^2) per anion dissociation channel `c`, via the
    nuclear-axis delta transform -- no elastic free-reference subtraction
    (DA is a pure rearrangement channel, no `v'==v_init` diagonal), the SAME
    `_C_DA = pi` convention `_flux_da_sigma_one_energy` uses."""
    n_channels = len(eps_e)
    sigma = np.zeros(n_channels, dtype=np.float64)
    if E <= 0.0:
        return sigma
    s_full = _dirac_da_s_vector_one_energy(
        g_elec, g_nuc, model, mu_r, result, eps, v_init, eps_e, R_position, E, dt, wp_in
    )
    e_tot = E + eps[v_init]
    for c in range(n_channels):
        if e_tot - eps_e[c] <= 0.0:
            continue  # closed dissociation channel
        sigma[c] = _C_DA * abs(s_full[c]) ** 2 / (2.0 * E)
    return sigma


class Dirac:
    """The delta (Dirac) `Extractor`: eMoScat's `DiracTestFunction2d` --
    Tannor-Weeks with a delta-distribution test function instead of the
    Gaussian test packet.

    Electronic (`axis="electronic"`, the default): `record` projects `psi`
    onto `chi_{v'}` in the nuclear coordinate at a FIXED electronic DVR index
    `position` -- no propagated outgoing test function, unlike `TannorWeeks`.
    `position` must land in the real (unscaled) electronic region
    (`grid.real_points[position] <= grid.R0`), typically an element end well
    past the interaction range (mirroring `TannorWeeks`'s `wp_out` standoff).
    `sigma(E, free=...)` accepts the companion free-reference `Dirac` (built
    from ITS OWN, separately propagated, run), matching `TannorWeeks.sigma`'s
    `free` contract.

    Nuclear (`axis="nuclear"`): the DISSOCIATIVE ATTACHMENT (DA) sibling of
    `Flux(axis="nuclear")` (module docstring's `Flux(axis="nuclear")`
    section) -- `position` is a NUCLEAR DVR index, `n_channels` selects how
    many anion electronic bound states (`qscat.core.dissociation.
    anion_electronic_states`, at `R_inf = tgrid.grids[1].R0`) are tracked as
    exit channels; `vprimes` is unused (pass `[]`). `record` projects `psi`
    onto `phi_c(r)` (the electronic anion state) at the FIXED nuclear node
    `position` -- the SAME point-VALUE projection as the electronic path,
    transposed to the other axis (no derivative, unlike `Flux`: a delta test
    function needs only the point value). `sigma(E)` returns the delta
    transform's `sigma_DA,c(E)` per channel (module docstring), `eta_out ->
    hankel_point_value(..., mass=model.mu)` at `position` -- shape
    `(n_channels,)` for scalar `E`, matching `Flux(axis="nuclear")`'s
    contract. `free` is not supported (DA has no elastic diagonal to
    subtract a reference from; passing it raises `ValueError`).

    Construct one instance per propagation run (the full run, and -- when an
    elastic free reference is wanted, electronic axis only -- a second
    instance for a SEPARATE `V_int=0` run); pass it to `time_dependent.
    propagate(..., extractors=[dirac])`.

    `position` must be a real (unscaled) DVR index, on the axis-appropriate
    grid, in the asymptotic region (past the interaction) -- same
    requirement as `Flux`'s `surface`.
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
        axis: str = "electronic",
        n_channels: int = 1,
    ) -> None:
        _check_axis(axis, "Dirac")
        self._axis = axis
        grid = tgrid.grids[_axis_grid_index(axis)]
        region = "electronic" if axis == "electronic" else "nuclear"
        if not (0 <= position < grid.n):
            raise ValueError(f"position {position} out of range for grid of size {grid.n}")
        if grid.real_points[position] > grid.R0:
            raise ValueError(
                f"position {position} (r={grid.real_points[position]}) is not in the real "
                f"(unscaled) {region} region (R0={grid.R0}) -- pick an index with "
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
        if axis == "electronic":
            self._n_channels = len(vprimes)
            self._chi = [np.asarray(chi[vp], dtype=np.complex128) for vp in vprimes]
        else:  # nuclear
            eps_e, phi = anion_electronic_states(
                tgrid.grids[0], model, R_inf=tgrid.grids[1].R0, n_states=n_channels
            )
            self._n_channels = n_channels
            self._eps_e = eps_e
            self._phi = [np.asarray(phi[c], dtype=np.complex128) for c in range(n_channels)]
        self._rows: list[npt.NDArray[np.complex128]] = []

    def record(self, psi: npt.NDArray[np.complex128]) -> None:
        """Append this step's point-value projection row (class docstring):
        `b_{v'}(t_n) = <chi_{v'}|psi(position, .)> / sqrt(w_r[position])`
        (electronic axis) or `b_c(t_n) = <phi_c|psi(., R=position)> /
        sqrt(w_R[position])` (nuclear axis) -- the DVR-coefficient-to-VALUE
        conversion `qscat.core.lcp`'s boundary flux also relies on."""
        block = psi.reshape(self._tgrid.shape)
        if self._axis == "electronic":
            psi_row = block[self._position, :]
            row = np.empty(len(self._vprimes), dtype=np.complex128)
            for k, chi_vp in enumerate(self._chi):
                row[k] = c_product(chi_vp, psi_row) * self._inv_sqrt_w
        else:  # nuclear
            psi_col = block[:, self._position]
            row = np.empty(self._n_channels, dtype=np.complex128)
            for k, phi_c in enumerate(self._phi):
                row[k] = c_product(phi_c, psi_col) * self._inv_sqrt_w
        self._rows.append(row)

    def _result(self, n_steps: int | None = None) -> PropagationResult:
        """`.result`, optionally truncated to the FIRST `n_steps` recorded
        samples -- see `TannorWeeks._result` for the exact contract
        (`n_steps=None` is byte-identical to the pre-truncation `.result`)."""
        rows = self._rows if n_steps is None else self._rows[:n_steps]
        n_t = len(rows)
        n_ch = len(self._vprimes) if self._axis == "electronic" else self._n_channels
        c = np.stack(rows) if rows else np.zeros((0, n_ch), dtype=np.complex128)
        t = np.arange(n_t, dtype=np.float64) * self._dt
        norm = np.zeros(n_t, dtype=np.float64)
        return PropagationResult(t=t, c=c, norm=norm, snapshots=[])

    @property
    def result(self) -> PropagationResult:
        """The recorded series as a `PropagationResult` (same shape/role as
        `TannorWeeks.result`; `norm`/`snapshots` are left empty/zero). The
        full series (`self._result(None)`); see `sigma`'s `n_steps` for a
        truncated read."""
        return self._result()

    def sigma(
        self,
        E: float | npt.ArrayLike,
        *,
        free: Extractor | None = None,
        n_steps: int | None = None,
    ) -> npt.NDArray[np.float64]:
        """`sigma_{v_init->v'}(E)` (bohr^2, electronic axis) or `sigma_DA,c(E)`
        (bohr^2 per anion channel `c`, nuclear axis) via the delta transform
        (module docstring / class docstring).

        Electronic: `free`, when given, must be a companion `Dirac` extractor
        recorded from a SEPARATE `V_int=0` propagation -- same elastic
        free-reference contract as `TannorWeeks.sigma`. `free` is typed as the
        general `Extractor` protocol (one shared `sigma` signature) but only
        a `Dirac` is meaningful here; anything else raises `TypeError`.

        Nuclear: DA has no elastic diagonal to subtract a reference from --
        `free` must be `None` (`ValueError` otherwise).

        `n_steps` (keyword-only, default `None`): transform only the FIRST
        `n_steps` recorded samples -- see `TannorWeeks.sigma`'s `n_steps` for
        the exact contract (byte-identical to today when `None`; `free`'s
        series is truncated to the SAME `n_steps` when both are given).
        """
        if self._axis == "nuclear":
            if free is not None:
                raise ValueError(
                    "Dirac.sigma(axis='nuclear'): DA has no elastic free-reference "
                    "subtraction -- free must be None"
                )
            e_arr = np.atleast_1d(np.asarray(E, dtype=np.float64))
            result = self._result(n_steps)
            g_elec = self._tgrid.grids[0]
            g_nuc = self._tgrid.grids[1]
            mu_r = self._model.mu
            out = np.stack(
                [
                    _dirac_da_sigma_one_energy(
                        g_elec,
                        g_nuc,
                        self._model,
                        mu_r,
                        result,
                        self._eps,
                        self._v_init,
                        self._eps_e,
                        self._z_position,
                        float(e),
                        self._dt,
                        self._wp_in,
                    )
                    for e in e_arr
                ]
            )
            scalar = np.isscalar(E) or (isinstance(E, np.ndarray) and np.ndim(E) == 0)
            if scalar:
                return np.asarray(out[0], dtype=np.float64)
            return out

        free_result: PropagationResult | None = None
        if free is not None:
            if not isinstance(free, Dirac):
                raise TypeError(f"Dirac.sigma: free must be a Dirac, got {type(free)}")
            free_result = free._result(n_steps)
        e_arr = np.atleast_1d(np.asarray(E, dtype=np.float64))
        result = self._result(n_steps)
        grid = self._tgrid.grids[_axis_grid_index(self._axis)]
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
    weights = quadrature_weights(t.size)
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


def _flux_da_s_vector_one_energy(
    g_elec: FemDvrEcsGrid,
    g_nuc: FemDvrEcsGrid,
    model: ResonanceModel,
    mu_r: float,
    t: npt.NDArray[np.float64],
    b: npt.NDArray[np.complex128],
    d: npt.NDArray[np.complex128],
    eps: npt.NDArray[np.float64],
    v_init: int,
    eps_e: npt.NDArray[np.float64],
    R_surface: float,
    E: float,
    dt: float,
    wp_in: _WpIn,
) -> npt.NDArray[np.complex128]:
    """The raw DA flux-transform S-matrix, per anion dissociation channel --
    `_flux_s_vector_one_energy`'s nuclear-axis twin (module docstring): the
    outgoing wave moves to the nuclear axis at mass `mu_r` (`l=0`), the
    incident deconvolution `eta_in` stays on the ELECTRONIC incident axis.
    """
    n_channels = len(eps_e)
    S = np.zeros(n_channels, dtype=np.complex128)
    if E <= 0.0:
        return S
    weights = quadrature_weights(t.size)
    e_tot = E + eps[v_init]
    k = float(np.sqrt(2.0 * E))
    eta_in = eta_incident(g_elec, k, model.ell, **wp_in)
    phase = np.exp(1j * e_tot * t)
    for c in range(n_channels):
        e_dr = e_tot - eps_e[c]
        if e_dr <= 0.0:
            continue  # closed dissociation channel
        k_r = float(np.sqrt(2.0 * mu_r * e_dr))
        phi_out, dphi_out = outgoing_surface_wave(g_nuc, R_surface, k_r, 0, model.charge, mass=mu_r)
        wronskian = np.conj(phi_out) * d[:, c] - b[:, c] * np.conj(dphi_out)
        s_raw = np.sum(weights * wronskian * phase) * dt
        S[c] = (-1j / (2.0 * mu_r * eta_in)) * s_raw
    return S


def _flux_da_sigma_one_energy(
    g_elec: FemDvrEcsGrid,
    g_nuc: FemDvrEcsGrid,
    model: ResonanceModel,
    mu_r: float,
    t: npt.NDArray[np.float64],
    b: npt.NDArray[np.complex128],
    d: npt.NDArray[np.complex128],
    eps: npt.NDArray[np.float64],
    v_init: int,
    eps_e: npt.NDArray[np.float64],
    R_surface: float,
    E: float,
    dt: float,
    wp_in: _WpIn,
) -> npt.NDArray[np.float64]:
    """`sigma_DA,c(E)` (bohr^2) per anion dissociation channel `c`, via the
    nuclear-axis flux transform -- no elastic free-reference subtraction
    (DA is a pure rearrangement channel, no `v'==v_init` diagonal)."""
    n_channels = len(eps_e)
    sigma = np.zeros(n_channels, dtype=np.float64)
    if E <= 0.0:
        return sigma
    s_full = _flux_da_s_vector_one_energy(
        g_elec, g_nuc, model, mu_r, t, b, d, eps, v_init, eps_e, R_surface, E, dt, wp_in
    )
    e_tot = E + eps[v_init]
    for c in range(n_channels):
        if e_tot - eps_e[c] <= 0.0:
            continue  # closed dissociation channel
        sigma[c] = _C_DA * abs(s_full[c]) ** 2 / (2.0 * E)
    return sigma


class Flux:
    """The flow (flux) `Extractor`: eMoScat's `FluxTestFunction2d` -- the
    time-energy Fourier transform of the probability flux projected onto the
    outgoing channel, at a FIXED surface.

    Electronic (`axis="electronic"`, the default): `record` appends, per
    `v'`, BOTH the value `b_{v'}(t) = <chi_{v'}|psi(surface,.)>` (`Dirac`'s
    line projection) AND its electronic-coordinate derivative `d_{v'}(t) =
    <chi_{v'}| d/dr psi(surface,.)>` (via `qscat.dvr.
    dvr_first_derivative_at_node` applied along the electronic axis, then a
    nuclear c-product onto `chi_{v'}` -- no extra `/sqrt(w)` needed there,
    unlike `b_{v'}`: `dvr_first_derivative_at_node` already converts
    coefficient->value internally, see its docstring). `sigma(E, free=...)`
    is the Wronskian-like flux transform (module docstring); same elastic
    free-reference contract as `TannorWeeks.sigma`/`Dirac.sigma`. `vprimes`
    is the VE exit-vibrational-level list (electronic-axis param only).

    Nuclear (`axis="nuclear"`): the DISSOCIATIVE ATTACHMENT (DA) extractor
    (module docstring's `Flux(axis="nuclear")` section) -- `surface` is a
    NUCLEAR DVR index, `n_channels` selects how many anion electronic bound
    states (`qscat.core.dissociation.anion_electronic_states`, at `R_inf =
    tgrid.grids[1].R0`) are tracked as exit channels; `vprimes` is unused
    (pass `[]`). `record` appends `b_c(t)`/`d_c(t)` per anion channel `c`.
    `sigma(E)` returns `sigma_DA,c(E)` per channel, shape `(n_channels,)` for
    scalar `E` (matching `dissociation.da_cross_section`'s per-channel
    return contract) -- `free` is not supported (DA has no elastic diagonal
    to subtract a reference from; passing it raises `ValueError`).

    `surface` must be a real (unscaled) DVR index, on the axis-appropriate
    grid, in the asymptotic region (past the interaction) -- same
    requirement as `Dirac`'s `position`.
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
        axis: str = "electronic",
        n_channels: int = 1,
    ) -> None:
        _check_axis(axis, "Flux")
        self._axis = axis
        grid = tgrid.grids[_axis_grid_index(axis)]
        region = "electronic" if axis == "electronic" else "nuclear"
        if not (0 <= surface < grid.n):
            raise ValueError(f"surface {surface} out of range for grid of size {grid.n}")
        if grid.real_points[surface] > grid.R0:
            raise ValueError(
                f"surface {surface} (r={grid.real_points[surface]}) is not in the real "
                f"(unscaled) {region} region (R0={grid.R0}) -- pick an index with "
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
        if axis == "electronic":
            self._n_channels = len(vprimes)
            self._chi = [np.asarray(chi[vp], dtype=np.complex128) for vp in vprimes]
        else:  # nuclear
            eps_e, phi = anion_electronic_states(
                tgrid.grids[0], model, R_inf=tgrid.grids[1].R0, n_states=n_channels
            )
            self._n_channels = n_channels
            self._eps_e = eps_e
            self._phi = [np.asarray(phi[c], dtype=np.complex128) for c in range(n_channels)]
        self._b_rows: list[npt.NDArray[np.complex128]] = []
        self._d_rows: list[npt.NDArray[np.complex128]] = []

    def record(self, psi: npt.NDArray[np.complex128]) -> None:
        """Append this step's `b(t_n)`/`d(t_n)` rows (docstring): per `v'`
        (electronic axis) or per anion channel `c` (nuclear axis)."""
        block = psi.reshape(self._tgrid.shape)
        row_b = np.empty(self._n_channels, dtype=np.complex128)
        row_d = np.empty(self._n_channels, dtype=np.complex128)
        if self._axis == "electronic":
            psi_row = block[self._surface, :]
            dpsi_row = self._deriv_row @ block  # (n_nuclear,): d/dr psi(surface, R), coeff-in-R
            for k, chi_vp in enumerate(self._chi):
                row_b[k] = c_product(chi_vp, psi_row) * self._inv_sqrt_w
                row_d[k] = c_product(chi_vp, dpsi_row)
        else:  # nuclear
            psi_col = block[:, self._surface]
            dpsi_col = block @ self._deriv_row  # (n_electronic,): d/dR psi(r, surface), coeff-in-r
            for k, phi_c in enumerate(self._phi):
                row_b[k] = c_product(phi_c, psi_col) * self._inv_sqrt_w
                row_d[k] = c_product(phi_c, dpsi_col)
        self._b_rows.append(row_b)
        self._d_rows.append(row_d)

    def _arrays_n(
        self, n_steps: int | None = None
    ) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.complex128], npt.NDArray[np.complex128]]:
        """`.series`, optionally truncated to the FIRST `n_steps` recorded
        samples -- see `TannorWeeks._result` for the exact contract
        (`n_steps=None` is byte-identical to the full `.series` read)."""
        b_rows = self._b_rows if n_steps is None else self._b_rows[:n_steps]
        d_rows = self._d_rows if n_steps is None else self._d_rows[:n_steps]
        n_t = len(b_rows)
        n_ch = self._n_channels
        b = np.stack(b_rows) if b_rows else np.zeros((0, n_ch), dtype=np.complex128)
        d = np.stack(d_rows) if d_rows else np.zeros((0, n_ch), dtype=np.complex128)
        t = np.arange(n_t, dtype=np.float64) * self._dt
        return t, b, d

    @property
    def series(
        self,
    ) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.complex128], npt.NDArray[np.complex128]]:
        """The full recorded flux series `(t, boundary_value, derivative)` --
        the raw material of the flux transform, the `Flux` analogue of
        `TannorWeeks`/`Dirac`'s `.result` (which carry `.t`/`.c` instead).
        Public so a caller can persist the raw series (e.g. the `qscat-run`
        CLI's optional `correlations` artifact); see `sigma`'s `n_steps` for a
        truncated read."""
        return self._arrays_n()

    def sigma(
        self,
        E: float | npt.ArrayLike,
        *,
        free: Extractor | None = None,
        n_steps: int | None = None,
    ) -> npt.NDArray[np.float64]:
        """`sigma_{v_init->v'}(E)` (bohr^2, electronic axis) or `sigma_DA,c(E)`
        (bohr^2 per anion channel `c`, nuclear axis) via the flux transform
        (module docstring / class docstring).

        Electronic: `free`, when given, must be a companion `Flux` extractor
        recorded from a SEPARATE `V_int=0` propagation -- same elastic
        free-reference contract as `TannorWeeks.sigma`/`Dirac.sigma`. `free`
        is typed as the general `Extractor` protocol (one shared `sigma`
        signature) but only a `Flux` is meaningful here; anything else
        raises `TypeError`.

        Nuclear: DA has no elastic diagonal to subtract a reference from --
        `free` must be `None` (`ValueError` otherwise).

        `n_steps` (keyword-only, default `None`): transform only the FIRST
        `n_steps` recorded samples -- see `TannorWeeks.sigma`'s `n_steps` for
        the exact contract (byte-identical to today when `None`; `free`'s
        series is truncated to the SAME `n_steps` when both are given).
        """
        if self._axis == "nuclear":
            if free is not None:
                raise ValueError(
                    "Flux.sigma(axis='nuclear'): DA has no elastic free-reference "
                    "subtraction -- free must be None"
                )
            e_arr = np.atleast_1d(np.asarray(E, dtype=np.float64))
            t, b, d = self._arrays_n(n_steps)
            g_elec = self._tgrid.grids[0]
            g_nuc = self._tgrid.grids[1]
            mu_r = self._model.mu
            out = np.stack(
                [
                    _flux_da_sigma_one_energy(
                        g_elec,
                        g_nuc,
                        self._model,
                        mu_r,
                        t,
                        b,
                        d,
                        self._eps,
                        self._v_init,
                        self._eps_e,
                        self._z_surface,
                        float(e),
                        self._dt,
                        self._wp_in,
                    )
                    for e in e_arr
                ]
            )
            scalar = np.isscalar(E) or (isinstance(E, np.ndarray) and np.ndim(E) == 0)
            if scalar:
                return np.asarray(out[0], dtype=np.float64)
            return out

        free_arrays: (
            tuple[npt.NDArray[np.float64], npt.NDArray[np.complex128], npt.NDArray[np.complex128]]
            | None
        ) = None
        if free is not None:
            if not isinstance(free, Flux):
                raise TypeError(f"Flux.sigma: free must be a Flux, got {type(free)}")
            free_arrays = free._arrays_n(n_steps)
        e_arr = np.atleast_1d(np.asarray(E, dtype=np.float64))
        t, b, d = self._arrays_n(n_steps)
        grid = self._tgrid.grids[_axis_grid_index(self._axis)]
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
