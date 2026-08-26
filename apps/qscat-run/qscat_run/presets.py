"""The single source of the validated per-molecule PRODUCTION decks.

Consolidates what today lives in `validation/diatomic/config.py` (F2/NO) and
`validation/h2plus/config.py` (H2P) -- as literal DATA, copied by value, not
imported: `qscat_run` must never depend on `validation.*`/`projects.*` (see
`tests/test_no_validation_import.py`). N2's deck has no `validation/`
counterpart of its own (its numbers live in `projects.n2_2d_cross_section
.convergence.WORKING_GRID`/`projects.n2_2d_td_cross_section.convergence
.TD_WORKING_GRID`) and are transcribed here the same way.

Two grid roles per molecule, per the CLI design spec:
  - `ti_grid()` -- the driven-solve deck (small electronic box; F2/NO reuse
    their `da_grid()`-style fine per-molecule nuclear deck so DA is resolved
    too).
  - `td_grid()` -- the LAUNCH-BOX deck: a larger electronic `r_max` so an
    incident wavepacket launches, interacts, and (for DA/DR) its outgoing
    test function sits, entirely inside the real region. For F2/NO this is
    documented as a REASONABLE, not yet independently convergence-tested,
    choice (r_max=30, within the spec's suggested 25-40 bohr range) --
    unlike N2 (r_max=50) and H2+ (r_max=1300/60), which come from an actual
    Task-4-style convergence study or an eMoScat deck. See the module's
    `CONCERNS` note in the Task 1 report.

`PRESETS` keys are `"{molecule}:{variant}"` (e.g. `"F2:emoscat"`,
`"H2P:proxy"`) -- a molecule maps to one variant (N2/NO/F2, laptop-feasible
already) or two (H2P: `emoscat` the ~1.15M-unknown Docker/MUMPS deck,
`proxy` a laptop-sized reduction, both from `validation/h2plus/config.py`).
`DEFAULT_PRESET` (`"emoscat"`) is what an omitted `grid.preset` resolves to.

Kept pure: data + builder functions only, no CLI/YAML concerns.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, replace

from qscat.core.grids import electronic_grid, fem_grid_exp_tail, nuclear_grid, segmented_grid
from qscat.dvr import FemDvrEcsGrid, TensorGrid
from qscat.model import F2, H2P, N2, NO, O2, O2_SO12, O2_SO32, ResonanceModel

from qscat_run.config import (
    ConfigError,
    EnergySpec,
    ExperimentConfig,
    IncidentSpec,
    NrmSpec,
    Observable,
    SegmentSpec,
    TdSpec,
    TestFunctionSpec,
)

__all__ = [
    "DEFAULT_PRESET",
    "MODELS",
    "PRESETS",
    "VALIDITY",
    "WARN_OBSERVABLES",
    "MoleculePreset",
    "available_presets",
    "nuclear_angle_b",
    "nuclear_grid_at_angle",
    "resolve_defaults",
    "resolve_grid",
    "resolve_lcp_grids",
    "resolve_nrm_grids",
    "resolve_surface_r",
    "resolve_test_function",
]

MODELS: dict[str, ResonanceModel] = {
    "N2": N2,
    "NO": NO,
    "F2": F2,
    "H2P": H2P,
    "O2": O2,
    "O2_SO12": O2_SO12,
    "O2_SO32": O2_SO32,
}

# The (molecule, observable.kind) validity matrix from the design spec.
# N2 "da" is CLOSED-IN-RANGE -- allowed, not rejected, but flagged
# (`WARN_OBSERVABLES`) since the channel is only weakly open near the top of
# the studied energy range.
VALIDITY: dict[str, frozenset[str]] = {
    "N2": frozenset({"ve", "da"}),
    "NO": frozenset({"ve", "da", "resonance_levels"}),
    "F2": frozenset({"ve", "da", "resonance_levels"}),
    "H2P": frozenset({"dr"}),
    # O2 (the fitted model): DA is closed until 3.7 eV, above the whole
    # window the model was fitted for (0-2.7 eV) -- VE only.
    "O2": frozenset({"ve"}),
    "O2_SO12": frozenset({"ve"}),
    "O2_SO32": frozenset({"ve"}),
}
WARN_OBSERVABLES: dict[str, frozenset[str]] = {"N2": frozenset({"da"})}

DEFAULT_PRESET = "emoscat"


@dataclass(frozen=True)
class MoleculePreset:
    """One molecule/variant's numerical deck: the TI and TD grid builders,
    the default energy sweep, the default TD incident + PER-KIND outgoing
    test-function parameters, the observables this molecule supports, and
    the vibrational basis size to diagonalize.

    `ve_test_function` is the ELECTRONIC outgoing packet (in `r`) `td.
    ve.TannorWeeks`/`Dirac`/`Flux` use for a `ve` observable -- also the
    fixed electronic analysis point (`_electronic_index_near`'s `r_value`),
    unchanged from before this fix. `da_test_function`/`dr_test_function`
    are the NUCLEAR outgoing packet (in `R`) for `da`/`dr` -- a physically
    DIFFERENT scale from the electronic packet (this is the bug this fix
    corrects: previously the single electronic packet doubled as the
    nuclear one too). `da_surface_R`/`dr_surface_R` are the fixed nuclear
    analysis point `Dirac`/`Flux` read at for that kind -- NOT necessarily
    the same coordinate as the nuclear packet's own `r0_out` (e.g. F2's
    validated DA deck: the outgoing packet is centered at R=8 but the
    surface/point extractors read at R=6, see `libs/qscat/tests/
    test_td_extractors.py::test_nuclear_flux_da_converges_to_ti_oracle`).
    `None` for a kind this molecule does not (or does not yet validatedly)
    support in TD.
    """

    molecule: str
    variant: str
    ti_grid: Callable[[], TensorGrid]
    td_grid: Callable[[], TensorGrid]
    default_energies: EnergySpec
    default_incident: IncidentSpec
    valid_observables: frozenset[str]
    n_vib: int
    ve_test_function: TestFunctionSpec | None = None
    da_test_function: TestFunctionSpec | None = None
    da_surface_R: float | None = None
    dr_test_function: TestFunctionSpec | None = None
    dr_surface_R: float | None = None
    # `(nuclear, elec_a, elec_b)` grid builder for the LCP method (the two
    # ECS-angle electronic decks the resonance-pole match needs + the fine
    # nuclear deck). `None` for molecules with no LCP path (H2P -- DR, not
    # DA/VE-LCP).
    lcp_grids: Callable[[], tuple[FemDvrEcsGrid, FemDvrEcsGrid, FemDvrEcsGrid]] | None = None
    # This molecule's electronic deck rebuilt at the NRM's SECOND ECS angle
    # (choice A's two-angle resonance-pole walk). The NRM's other two grids are
    # `ti_grid()`'s own factors, so `nrm` and `ti` in one run are computed on
    # one discretisation and a ratio across method prefixes means something.
    # `None` for molecules with no NRM path (H2P -- DR, not VE/DA).
    nrm_elec_b: Callable[[], FemDvrEcsGrid] | None = None


# --- N2 -----------------------------------------------------------------
# `projects.n2_2d_cross_section.convergence.WORKING_GRID` (TI) /
# `projects.n2_2d_td_cross_section.convergence.TD_WORKING_GRID` (TD).


# N2's electronic deck, named because the NRM has to rebuild it at a second ECS
# angle and must not drift from the one the exact `ti` solve runs on.
_N2_ELEC: dict[str, float | int] = {
    "r_max": 16.0,
    "angle_deg": 35.0,
    "order": 7,
    "n_complex": 5,
}


def _n2_ti_grid() -> TensorGrid:
    return TensorGrid(
        [
            electronic_grid(**_N2_ELEC),  # type: ignore[arg-type]
            nuclear_grid(angle_deg=35.0, r_max=20.0, n_complex=5, quadrature=10),
        ]
    )


def _n2_td_grid() -> TensorGrid:
    return TensorGrid(
        [
            electronic_grid(r_max=50.0, angle_deg=35.0, order=8, n_complex=6),
            nuclear_grid(angle_deg=35.0, r_max=22.0, n_complex=5, quadrature=10),
        ]
    )


# --- NO / F2 --------------------------------------------------------------
# `validation.diatomic.config.CONFIGS["NO"/"F2"]`: electronic_grid (VE deck)
# x the eMoScat per-molecule nuclear deck (`da_grid()`, via `segmented_grid`).

_NO_NUC_REAL = ((1, 1.0), (1, 1.6), (37, 9.0))
_NO_NUC_COMPLEX = ((1, 9.25), (1, 10.0), (1, 12.0), (4, 42.0))
_NO_NUC_ANGLE, _NO_NUC_QUAD = 45.0, 14

_F2_NUC_REAL = ((9, 1.8), (1, 2.0), (5, 2.5), (4, 2.596908), (4, 2.7), (40, 10.7))
_F2_NUC_COMPLEX = (
    (1, 10.8),
    (1, 11.0),
    (1, 11.5),
    (1, 12.5),
    (1, 14.0),
    (1, 18.0),
    (4, 30.0),
    (2, 101.0),
)
_F2_NUC_ANGLE, _F2_NUC_QUAD = 35.0, 14


def _no_nuc_grid(angle_deg: float = _NO_NUC_ANGLE) -> FemDvrEcsGrid:
    return segmented_grid(
        _NO_NUC_REAL, _NO_NUC_COMPLEX, angle_deg=angle_deg, quadrature=_NO_NUC_QUAD
    )


def _no_ti_grid() -> TensorGrid:
    return TensorGrid([electronic_grid(r_max=16.0, order=8, n_complex=6), _no_nuc_grid()])


def _no_td_grid() -> TensorGrid:
    # Launch-box: same fine nuclear deck, a larger electronic r_max (see
    # module docstring's caveat -- not independently convergence-tested).
    return TensorGrid([electronic_grid(r_max=30.0, order=8, n_complex=6), _no_nuc_grid()])


def _f2_nuc_grid(angle_deg: float = _F2_NUC_ANGLE) -> FemDvrEcsGrid:
    return segmented_grid(
        _F2_NUC_REAL, _F2_NUC_COMPLEX, angle_deg=angle_deg, quadrature=_F2_NUC_QUAD
    )


def _f2_ti_grid() -> TensorGrid:
    return TensorGrid([electronic_grid(r_max=16.0, order=8, n_complex=6), _f2_nuc_grid()])


def _f2_td_grid() -> TensorGrid:
    return TensorGrid([electronic_grid(r_max=30.0, order=8, n_complex=6), _f2_nuc_grid()])


# --- O2 (fitted model; discretisation-tuner deck) --------------------------
# Both grids are `qscat.tuning.propose_grid(O2, ..., (0.002, 0.10))` -- the
# tuner's a-priori mesh for the 0-2.7 eV VE window, NOT an eMoScat deck (there
# is none: O2 is the factory's fit, docs/physics/potential-factory.md). The
# nuclear real region is the tuner's, cut at 8 bohr with the tuner's own ECS
# tail re-attached there: the VE path's fixed 18-bohr extent spent 62
# elements on empty space (DA closed until 3.7 eV; the anion's outer turning
# point at 2.3 eV is 4.0 bohr) -- and then h-REFINED ONCE (every real
# element halved): the 2-D spot check found one nuclear refinement moving
# sigma(0->1) at 1.36 eV by 69 %, after which the pair is converged (< 2 %).
# A comb of 1-8 meV peaks needs its levels far tighter than the 1-D probe's
# 1e-3. `validation/factory/o2_grids.py` regenerates and probes both, and
# `validation/factory/test_o2_grids.py` locks these numbers to it.
# 324 x 549 = 178k unknowns: MUMPS territory.
_O2_ELEC_REAL = (
    (26, 1.335887),
    (1, 1.713557),
    (1, 2.200777),
    (1, 2.477971),
    (1, 2.616568),
    (4, 2.893762),
    (1, 3.032359),
    (1, 3.309552),
    (3, 4.972715),
    (1, 5.249909),
    (1, 5.388506),
    (4, 5.665699),
    (1, 5.804296),
    (1, 6.081490),
    (1, 6.635877),
    (1, 8.298146),
    (1, 9.709861),
    (1, 11.027766),
    (1, 12.294768),
    (1, 13.529715),
    (1, 14.742692),
    (1, 15.939749),
    (1, 17.124803),
    (1, 18.300528),
    (1, 19.468831),
    (1, 19.99998),
)
_O2_ELEC_COMPLEX = (
    (2, 46.929626),
    (1, 63.375598),
    (1, 83.462753),
    (1, 107.997261),
    (1, 137.963775),
)
_O2_NUC_REAL = (
    (22, 1.650755),
    (2, 1.804038),
    (2, 1.954796),
    (2, 2.109889),
    (2, 2.262361),
    (2, 2.414378),
    (2, 2.569938),
    (2, 2.725830),
    (2, 2.877165),
    (2, 3.030356),
    (2, 3.185877),
    (2, 3.336426),
    (2, 3.488046),
    (2, 3.639623),
    (2, 3.798035),
    (2, 3.953462),
    (2, 4.103964),
    (2, 4.257053),
    (2, 4.411877),
    (2, 4.567795),
    (2, 4.724370),
    (2, 4.881324),
    (2, 5.038491),
    (2, 5.195774),
    (2, 5.353120),
    (2, 5.510500),
    (2, 5.667897),
    (2, 5.825303),
    (2, 5.982714),
    (2, 6.140127),
    (2, 6.297541),
    (2, 6.454956),
    (2, 6.612371),
    (2, 6.769787),
    (2, 6.927202),
    (2, 7.084618),
    (2, 7.242034),
    (2, 7.399449),
    (2, 7.556865),
    (2, 7.714281),
    (4, 8.029112),
)
_O2_NUC_COMPLEX = (
    (2, 8.252148),
    (1, 8.388356),
    (1, 8.554721),
    (1, 8.757920),
    (1, 9.006108),
)
_O2_ANGLE, _O2_QUAD = 35.0, 6


def _o2_elec_grid() -> FemDvrEcsGrid:
    return segmented_grid(_O2_ELEC_REAL, _O2_ELEC_COMPLEX, angle_deg=_O2_ANGLE, quadrature=_O2_QUAD)


def _o2_nuc_grid() -> FemDvrEcsGrid:
    return segmented_grid(_O2_NUC_REAL, _O2_NUC_COMPLEX, angle_deg=_O2_ANGLE, quadrature=_O2_QUAD)


def _o2_ti_grid() -> TensorGrid:
    return TensorGrid([_o2_elec_grid(), _o2_nuc_grid()])


# LCP (local-complex-potential) pole-matching uses TWO fixed-R electronic grids
# at distinct ECS angles + the fine per-molecule nuclear deck (the same deck DA
# uses) -- the exact recipe validation/diatomic/config.py's MoleculeConfig
# encodes (lcp_angle_a=35, lcp_angle_b=44).
_LCP_ANGLE_A, _LCP_ANGLE_B = 35.0, 44.0

# The NRM's second electronic ECS angle. It is NOT the LCP's 44 deg: the NRM
# runs its electronic Hamiltonian on the SAME deck the exact 2-D `ti` solve
# uses, and `PhysicalDiscreteState`'s two-angle pole walk needs a second angle
# near it. 40/35 is the pairing `validation/diatomic/nrm.py` and
# `validation/diatomic/ve_nrm.py` measured every recorded number on, so the app
# and the validation drivers agree by construction.
_NRM_ANGLE_B = 40.0


def _n2_lcp_grids() -> tuple[FemDvrEcsGrid, FemDvrEcsGrid, FemDvrEcsGrid]:
    """(nuclear, elec_a, elec_b) for N2 LCP -- the TI deck's own factors plus
    the electronic deck rebuilt at the eMoScat LCP partner angle. 35/44 is
    the pairing F2/NO use here AND the pairing the projects N2 pole walk
    itself uses (`projects/n2_ti_cross_section/vres.py`: _ANGLE_A_DEG=35,
    _ANGLE_B_DEG=44). N2's LCP observable is VE (its DA is
    closed-in-range); the partner angle only gates two-angle pole
    STABILITY -- accepted pole values come from grid a's spectrum."""
    return (
        nuclear_grid(angle_deg=35.0, r_max=20.0, n_complex=5, quadrature=10),
        electronic_grid(**_N2_ELEC),  # type: ignore[arg-type]
        electronic_grid(**{**_N2_ELEC, "angle_deg": _LCP_ANGLE_B}),  # type: ignore[arg-type]
    )


def _lcp_elec(angle_deg: float) -> FemDvrEcsGrid:
    return electronic_grid(r_max=16.0, order=8, n_complex=6, angle_deg=angle_deg)


def _diatomic_nrm_elec_b() -> FemDvrEcsGrid:
    """NO/F2's electronic deck at the NRM's second ECS angle."""
    return _lcp_elec(_NRM_ANGLE_B)


def _n2_nrm_elec_b() -> FemDvrEcsGrid:
    """N2's electronic deck at the NRM's second ECS angle."""
    return electronic_grid(**{**_N2_ELEC, "angle_deg": _NRM_ANGLE_B})  # type: ignore[arg-type]


def _f2_lcp_grids() -> tuple[FemDvrEcsGrid, FemDvrEcsGrid, FemDvrEcsGrid]:
    """(nuclear, elec_a, elec_b) for F2 LCP -- fine DA nuclear deck + two ECS angles."""
    return _f2_nuc_grid(), _lcp_elec(_LCP_ANGLE_A), _lcp_elec(_LCP_ANGLE_B)


def _no_lcp_grids() -> tuple[FemDvrEcsGrid, FemDvrEcsGrid, FemDvrEcsGrid]:
    """(nuclear, elec_a, elec_b) for NO LCP -- fine DA nuclear deck + two ECS angles."""
    return _no_nuc_grid(), _lcp_elec(_LCP_ANGLE_A), _lcp_elec(_LCP_ANGLE_B)


# --- H2+ -------------------------------------------------------------------
# `validation.h2plus.config.full_grid`/`proxy_grid`, transcribed verbatim.
# TI and TD share the SAME grid here (no separate launch box): the full
# deck's electronic real region already runs to 1300 bohr, comfortably
# holding both the driven-solve deck and eMoScat's r0=800 incident.


def _h2p_full_grid() -> TensorGrid:
    electronic = fem_grid_exp_tail(
        ((10, 1.0), (10, 4.0), (16, 20.0), (20, 100.0), (120, 1300.0)),
        angle_deg=5.0,
        quadrature=8,
        tail_n=25,
    )
    nuclear = fem_grid_exp_tail(
        ((5, 1.0), (20, 4.0), (67, 14.0)),
        angle_deg=22.0,
        quadrature=8,
        tail_n=25,
    )
    return TensorGrid([electronic, nuclear])


def _h2p_proxy_grid() -> TensorGrid:
    electronic = fem_grid_exp_tail(
        ((10, 1.0), (10, 4.0), (16, 20.0), (10, 60.0)),
        angle_deg=5.0,
        quadrature=8,
        tail_n=8,
    )
    nuclear = fem_grid_exp_tail(
        ((5, 1.0), (20, 4.0), (40, 14.0)),
        angle_deg=22.0,
        quadrature=8,
        tail_n=8,
    )
    return TensorGrid([electronic, nuclear])


# The VALIDATED F2 nuclear (DA) outgoing test packet + surface: the outgoing
# packet is centered INWARD of the packet's own r0_out at the fixed
# analysis point -- see `libs/qscat/tests/test_td_extractors.py`'s
# `test_nuclear_flux_da_converges_to_ti_oracle`/`test_nuclear_tw_da_converges
# _to_ti_oracle` docstrings (SP2 validation) for the full rationale. NO has
# no validated TD DA experiment of its own (see CLAUDE.md's diatomic note);
# its nuclear real region is a comparable scale (segmented to R=9.0 vs F2's
# R=10.7), so it reuses F2's validated packet SHAPE verbatim as a reasonable,
# UNVALIDATED placeholder rather than leaving `da` TD unsupported for NO.
# N2's own nuclear grid (`nuclear_grid()`, real region to R=12.0) also
# comfortably fits this packet; N2's `da` channel is separately flagged
# closed-in-range (`WARN_OBSERVABLES`), so this too is a documented,
# UNVALIDATED placeholder, not a tuned N2 deck.
_F2_DA_TEST_FUNCTION = TestFunctionSpec(r0_out=8.0, p0_out=72.0, sigma_out=0.07)
_F2_DA_SURFACE_R = 6.0

# The VALIDATED H2+ nuclear (DR) outgoing test packet + surface. These values
# come from the TD-DR driver that lived under `validation/h2plus/`, which was
# removed in the qscat-run consolidation; they were copied by VALUE and no
# longer have a source to be checked against -- there is no lock test tying
# them to anything, so treat this block as the definition, not as a copy.
_H2P_DR_TEST_FUNCTION = TestFunctionSpec(r0_out=12.0, p0_out=15.0, sigma_out=0.4)
_H2P_DR_SURFACE_R = 12.0

PRESETS: dict[str, MoleculePreset] = {
    "N2:emoscat": MoleculePreset(
        molecule="N2",
        variant="emoscat",
        ti_grid=_n2_ti_grid,
        td_grid=_n2_td_grid,
        # Approximates the dense curve's np.linspace(0.005, 0.2, 60)
        # (validation/n2/ti_curve.py) as a min/max/step sweep.
        default_energies=EnergySpec(min=0.005, max=0.20, step=0.005),
        default_incident=IncidentSpec(r0=25.0, p0=-0.5, sigma=5.0),
        valid_observables=VALIDITY["N2"],
        n_vib=6,
        ve_test_function=TestFunctionSpec(r0_out=35.0, p0_out=0.5, sigma_out=4.0),
        da_test_function=_F2_DA_TEST_FUNCTION,
        da_surface_R=_F2_DA_SURFACE_R,
        lcp_grids=_n2_lcp_grids,
        nrm_elec_b=_n2_nrm_elec_b,
    ),
    "NO:emoscat": MoleculePreset(
        molecule="NO",
        variant="emoscat",
        ti_grid=_no_ti_grid,
        td_grid=_no_td_grid,
        default_energies=EnergySpec(min=0.004, max=0.120, step=0.004),
        # NO has no validated TD experiment yet (see CLAUDE.md's diatomic
        # note) -- these incident/test-function defaults are a reasonable,
        # UNVALIDATED scaling of N2's (see module docstring's caveat).
        default_incident=IncidentSpec(r0=20.0, p0=-0.4, sigma=4.0),
        valid_observables=VALIDITY["NO"],
        n_vib=4,
        ve_test_function=TestFunctionSpec(r0_out=24.0, p0_out=0.4, sigma_out=3.0),
        da_test_function=_F2_DA_TEST_FUNCTION,
        da_surface_R=_F2_DA_SURFACE_R,
        lcp_grids=_no_lcp_grids,
        nrm_elec_b=_diatomic_nrm_elec_b,
    ),
    "F2:emoscat": MoleculePreset(
        molecule="F2",
        variant="emoscat",
        ti_grid=_f2_ti_grid,
        td_grid=_f2_td_grid,
        default_energies=EnergySpec(min=0.004, max=0.100, step=0.004),
        # Same caveat as NO: unvalidated TD VE defaults, scaled from N2's
        # (the DA packet below IS validated -- see the module note above).
        default_incident=IncidentSpec(r0=20.0, p0=-0.4, sigma=4.0),
        valid_observables=VALIDITY["F2"],
        n_vib=4,
        ve_test_function=TestFunctionSpec(r0_out=24.0, p0_out=0.4, sigma_out=3.0),
        da_test_function=_F2_DA_TEST_FUNCTION,
        da_surface_R=_F2_DA_SURFACE_R,
        lcp_grids=_f2_lcp_grids,
        nrm_elec_b=_diatomic_nrm_elec_b,
    ),
    "O2:tuner": MoleculePreset(
        molecule="O2",
        variant="tuner",
        ti_grid=_o2_ti_grid,
        # TD on O2 is NOT validated; the TI deck stands in and the incident /
        # test-function packets are N2's scaled to the 20-bohr electronic box.
        td_grid=_o2_ti_grid,
        # A uniform sweep is a placeholder for O2: its VE peaks are 0.01-8
        # meV wide, so a real run uses the level-aware `energies: {values}`
        # list `validation/factory/o2_ve_energies.py` writes.
        default_energies=EnergySpec(min=0.002, max=0.100, step=0.001),
        default_incident=IncidentSpec(r0=12.0, p0=-0.5, sigma=3.0),
        valid_observables=VALIDITY["O2"],
        n_vib=12,
        ve_test_function=TestFunctionSpec(r0_out=14.0, p0_out=0.5, sigma_out=3.0),
    ),
    # The spin-orbit components share O2's deck: their anion curves differ
    # from O2's by +-10 meV, far under anything the discretisation resolves
    # differently. The level-aware mesh is per component (their levels sit
    # +-Delta_SO/2 apart).
    "O2_SO12:tuner": MoleculePreset(
        molecule="O2_SO12",
        variant="tuner",
        ti_grid=_o2_ti_grid,
        td_grid=_o2_ti_grid,
        default_energies=EnergySpec(min=0.002, max=0.100, step=0.001),
        default_incident=IncidentSpec(r0=12.0, p0=-0.5, sigma=3.0),
        valid_observables=VALIDITY["O2_SO12"],
        n_vib=12,
        ve_test_function=TestFunctionSpec(r0_out=14.0, p0_out=0.5, sigma_out=3.0),
    ),
    "O2_SO32:tuner": MoleculePreset(
        molecule="O2_SO32",
        variant="tuner",
        ti_grid=_o2_ti_grid,
        td_grid=_o2_ti_grid,
        default_energies=EnergySpec(min=0.002, max=0.100, step=0.001),
        default_incident=IncidentSpec(r0=12.0, p0=-0.5, sigma=3.0),
        valid_observables=VALIDITY["O2_SO32"],
        n_vib=12,
        ve_test_function=TestFunctionSpec(r0_out=14.0, p0_out=0.5, sigma_out=3.0),
    ),
    "H2P:emoscat": MoleculePreset(
        molecule="H2P",
        variant="emoscat",
        ti_grid=_h2p_full_grid,
        td_grid=_h2p_full_grid,
        default_energies=EnergySpec(min=0.001, max=0.050, step=0.001),
        default_incident=IncidentSpec(r0=800.0, p0=-0.25, sigma=8.0),
        valid_observables=VALIDITY["H2P"],
        n_vib=4,
        dr_test_function=_H2P_DR_TEST_FUNCTION,
        dr_surface_R=_H2P_DR_SURFACE_R,
    ),
    "H2P:proxy": MoleculePreset(
        molecule="H2P",
        variant="proxy",
        ti_grid=_h2p_proxy_grid,
        td_grid=_h2p_proxy_grid,
        default_energies=EnergySpec(min=0.001, max=0.050, step=0.001),
        # r0 scaled DOWN from the full deck's 800 to fit inside the proxy
        # grid's ~60-bohr electronic real region (an off-box incident lands
        # in the ECS tail and diverges -- a finding of the retired
        # `validation/h2plus/` TD-DR driver, removed in the qscat-run
        # consolidation; the reasoning is recorded in
        # docs/physics/td-da.md, the code is not in the tree).
        default_incident=IncidentSpec(r0=40.0, p0=-0.25, sigma=8.0),
        valid_observables=VALIDITY["H2P"],
        n_vib=4,
        dr_test_function=_H2P_DR_TEST_FUNCTION,
        dr_surface_R=_H2P_DR_SURFACE_R,
    ),
}


def available_presets(molecule: str) -> frozenset[str]:
    """The preset variant names (e.g. `{"emoscat", "proxy"}`) defined for `molecule`."""
    prefix = f"{molecule}:"
    return frozenset(key[len(prefix) :] for key in PRESETS if key.startswith(prefix))


def _build_explicit_segment(seg: SegmentSpec) -> FemDvrEcsGrid:
    """One explicit-grid axis -> a `FemDvrEcsGrid`, via `fem_grid_exp_tail`
    (the ECS spec gives `angle`/`elements`/`quadrature` but no explicit tail
    endpoint, so the exponential-tail builder -- not `segmented_grid`, which
    needs an endpoint for every segment -- is the one that fits)."""
    return fem_grid_exp_tail(
        seg.real,
        angle_deg=seg.ecs.angle,
        quadrature=seg.ecs.quadrature,
        tail_n=seg.ecs.elements,
    )


def resolve_grid(cfg: ExperimentConfig, method: str) -> TensorGrid:
    """The `TensorGrid` for `method` (`"ti"` or `"td"`): an explicit grid
    (if `cfg.grid.electronic`/`nuclear` are both given) used for BOTH
    methods, else the named preset's `ti_grid()`/`td_grid()`."""
    grid = cfg.grid
    if grid.electronic is not None and grid.nuclear is not None:
        return TensorGrid(
            [_build_explicit_segment(grid.electronic), _build_explicit_segment(grid.nuclear)]
        )

    variant = grid.preset or DEFAULT_PRESET
    key = f"{cfg.molecule}:{variant}"
    preset = PRESETS.get(key)
    if preset is None:
        raise ConfigError(
            f"no preset {variant!r} for molecule {cfg.molecule!r}; "
            f"available: {sorted(available_presets(cfg.molecule))}"
        )
    if method == "ti":
        return preset.ti_grid()
    if method == "td":
        return preset.td_grid()
    raise ConfigError(f"unknown method {method!r}; choose 'ti' or 'td'")


def resolve_lcp_grids(
    cfg: ExperimentConfig,
) -> tuple[FemDvrEcsGrid, FemDvrEcsGrid, FemDvrEcsGrid]:
    """The `(nuclear, elec_a, elec_b)` LCP grid triple for `cfg`'s molecule.

    LCP has no explicit-grid path (its two ECS-angle electronic decks + fine
    nuclear deck cannot be expressed by the single electronic/nuclear explicit
    schema), so it always resolves through the named preset's `lcp_grids`.
    Raises `ConfigError` if the molecule/variant has no LCP path.
    """
    variant = cfg.grid.preset or DEFAULT_PRESET
    key = f"{cfg.molecule}:{variant}"
    preset = PRESETS.get(key)
    if preset is None or preset.lcp_grids is None:
        raise ConfigError(
            f"the 'lcp' method is not available for {cfg.molecule} ({variant}); "
            "LCP is wired for the neutral diatomics (N2, F2, NO); H2P's observable is DR"
        )
    return preset.lcp_grids()


def resolve_nrm_grids(
    cfg: ExperimentConfig,
) -> tuple[FemDvrEcsGrid, FemDvrEcsGrid, FemDvrEcsGrid]:
    """The `(nuclear, elec_a, elec_b)` NRM grid triple for `cfg`'s molecule.

    `elec_a` and `nuclear` are `ti_grid()`'s OWN factors, so `nrm` and `ti` in
    one run are computed on one discretisation -- a `methods: [ti, nrm]` ratio
    then measures the model reduction rather than two discretisations. (For
    F2/NO those factors are also exactly the LCP path's decks: `_lcp_elec(
    _LCP_ANGLE_A)` is `_f2_ti_grid`'s electronic factor and `_f2_nuc_grid()`
    its nuclear one, so all three methods land on the same grid.) `elec_b` is
    the same electronic deck rebuilt at `_NRM_ANGLE_B`, for choice A's
    two-angle pole walk.

    Raises `ConfigError` if the molecule has no NRM deck.
    """
    variant = cfg.grid.preset or DEFAULT_PRESET
    preset = PRESETS.get(f"{cfg.molecule}:{variant}")
    if preset is None or preset.nrm_elec_b is None:
        raise ConfigError(
            f"the 'nrm' method is not available for {cfg.molecule} ({variant}); the "
            "nonlocal resonance model is wired for the neutral diatomics (N2, NO, F2)"
        )
    elec_a, nuc = preset.ti_grid().grids
    return nuc, elec_a, preset.nrm_elec_b()


_NUC_GRID_BUILDERS = {"NO": _no_nuc_grid, "F2": _f2_nuc_grid}

# How far below grid a's tail angle grid b sits, when not set explicitly.
# eMoScat's electronic LCP decks pair 44/35 and 40/30 -- about ten degrees.
_ANGLE_B_OFFSET = 10.0


def nuclear_angle_b(cfg: ExperimentConfig) -> float:
    """The second nuclear grid's tail angle: explicit, or `angle_a - 10` deg.

    Always moves DOWNWARD, which is unconditionally safe against the ionic
    model's `max_nuclear_ecs_angle_deg` divergence bound.
    """
    if cfg.grid.nuclear_angle_b is not None:
        return float(cfg.grid.nuclear_angle_b)
    g_a, _ea, _eb = resolve_lcp_grids(cfg)
    return max(el.angle_deg for el in g_a.spec.elements) - _ANGLE_B_OFFSET


def nuclear_grid_at_angle(cfg: ExperimentConfig, angle_deg: float) -> FemDvrEcsGrid:
    """This molecule's nuclear deck rebuilt at a different ECS tail angle.

    Same real segments and quadrature -- only the tail rotates -- so every real
    node is shared with `resolve_lcp_grids`'s nuclear grid, which is what
    `qscat.core.lcp.lcp_resonance_levels` requires of its two grids.
    """
    builder = _NUC_GRID_BUILDERS.get(cfg.molecule)
    if builder is None:
        raise ConfigError(
            f"no LCP nuclear deck for {cfg.molecule}; available: {sorted(_NUC_GRID_BUILDERS)}"
        )
    return builder(angle_deg)


def _default_channels(kind: str, preset: MoleculePreset) -> int | tuple[int, ...] | None:
    """A reasonable channel default when the config omits `channels`: for
    `ve`, the first few excited vibrational levels (bounded by `n_vib`);
    for `da`/`dr`, a single channel (the usual case); for `resonance_levels`,
    `None` -- there is no natural count, so report EVERY angle-stable level
    inside the default window (`qscat.core.lcp.lcp_resonance_levels`'s own
    `n_levels=None`), which is what the design spec and the shipped example
    both promise for an omitted `channels`."""
    if kind == "ve":
        return tuple(range(1, min(preset.n_vib, 4)))
    if kind == "resonance_levels":
        return None
    return 1


def _kinds_in(observables: tuple[Observable, ...]) -> set[str]:
    """The distinct `ve`/`da`/`dr` observable kinds a run actually requests
    (ignores any other kind -- unreachable in practice since `validate_config`
    already rejects unknown kinds)."""
    return {obs.kind for obs in observables if obs.kind in ("ve", "da", "dr")}


def _preset_test_function(preset: MoleculePreset, kind: str) -> TestFunctionSpec | None:
    return {
        "ve": preset.ve_test_function,
        "da": preset.da_test_function,
        "dr": preset.dr_test_function,
    }.get(kind)


def _test_function_for_kind(
    td: TdSpec, preset: MoleculePreset | None, kind: str, kinds: set[str]
) -> TestFunctionSpec | None:
    """The TD outgoing test-packet for one observable `kind`, given the full
    set of `kinds` this run requests -- the resolution order both
    `resolve_test_function` (below) and `resolve_defaults` share:

    1. An explicit per-kind override (`td.test_functions[kind]`).
    2. The flat back-compat `td.test_function`, but ONLY when it is
       unambiguous: either this run requests just ONE kind (there's nothing
       else it could mean), or `kind == "ve"` (the flat block's historical,
       sole meaning before this fix -- the electronic outgoing packet). A
       MIXED run's `da`/`dr` kind is NEVER filled from a flat block: that is
       precisely the bug this fix corrects (an electronic-scale packet
       silently reused as the nuclear one).
    3. The resolved preset's per-kind default
       (`ve_test_function`/`da_test_function`/`dr_test_function`).

    Returns `None` if none of the three resolves it (e.g. no preset AND no
    matching override) -- callers decide whether that is fatal.
    """
    if td.test_functions is not None and kind in td.test_functions:
        return td.test_functions[kind]
    if td.test_function is not None and (len(kinds) == 1 or kind == "ve"):
        return td.test_function
    if preset is not None:
        return _preset_test_function(preset, kind)
    return None


def resolve_test_function(cfg: ExperimentConfig, kind: str) -> TestFunctionSpec:
    """The resolved TD outgoing test-packet for observable `kind` (`"ve"`
    electronic, `"da"`/`"dr"` nuclear) -- see `_test_function_for_kind` for
    the resolution order. Works whether or not `resolve_defaults` has already
    run (it re-derives the preset itself), so `_run_td` can call it directly
    even for an explicit-grid config with no matching preset (as long as the
    user supplied enough of `td.test_function`/`test_functions` to cover
    `kind`). Raises `ConfigError` with an actionable message if `kind`
    remains unresolved.
    """
    if cfg.td is None:
        raise ConfigError("no 'td' block resolved for this config")
    kinds = _kinds_in(cfg.observables)
    variant = cfg.grid.preset or DEFAULT_PRESET
    preset = PRESETS.get(f"{cfg.molecule}:{variant}")
    tf = _test_function_for_kind(cfg.td, preset, kind, kinds)
    if tf is None:
        raise ConfigError(
            f"no TD test-function resolved for observable kind {kind!r} "
            f"(molecule={cfg.molecule!r}); supply 'td.test_function' (a flat "
            "{r0_out, p0_out, sigma_out} block, valid for a single-observable-kind "
            "run) or a per-kind 'td.test_function: {ve: {...}, da: {...}, dr: {...}}' "
            "mapping, or use a named grid preset that provides a default for this kind"
        )
    return tf


def resolve_surface_r(cfg: ExperimentConfig, kind: str) -> float:
    """The fixed analysis-surface coordinate (bohr) the `Dirac`/`Flux`
    extractors read at for observable `kind`: for `ve`, the resolved
    electronic test function's OWN `r0_out` (unchanged -- the electronic
    packet's center doubles as its analysis point, as before this fix). For
    `da`/`dr`, the resolved preset's `da_surface_R`/`dr_surface_R` -- a
    DIFFERENT coordinate from the nuclear packet's `r0_out` in general (see
    `MoleculePreset`'s docstring) -- when a preset resolves; otherwise (an
    explicit custom grid with no matching preset) falls back to the resolved
    nuclear test function's own `r0_out`, mirroring `ve`'s convention (the
    user built that grid and picked `r0_out` for it already, so it is a
    reasonable analysis point absent a preset-specific one)."""
    if kind == "ve":
        return resolve_test_function(cfg, "ve").r0_out
    variant = cfg.grid.preset or DEFAULT_PRESET
    preset = PRESETS.get(f"{cfg.molecule}:{variant}")
    if preset is not None:
        surface = preset.da_surface_R if kind == "da" else preset.dr_surface_R
        if surface is not None:
            return surface
    return resolve_test_function(cfg, kind).r0_out


def resolve_defaults(cfg: ExperimentConfig) -> ExperimentConfig:
    """Fill omitted `energies`/`td.incident`/`td.test_functions`/
    `Observable.channels` from the resolved preset. A no-op (returns `cfg`
    unchanged) if no preset matches (an explicit grid with no matching
    preset key, or an unknown molecule/variant) -- `validate_config` is
    responsible for rejecting that case with an actionable message.

    `td.test_functions` is filled to a COMPLETE per-kind mapping (covering
    every kind this run's `observables` request) via `_test_function_for_kind`
    -- purely for API convenience/inspection (e.g. `qscat-run validate`'s
    resolved-config introspection); `_run_td` does not depend on this having
    run first, since `resolve_test_function`/`resolve_surface_r` re-derive
    the same resolution directly from `cfg` (so an explicit-grid config with
    no matching preset, where this function no-ops, still resolves correctly
    at point of use provided the user supplied enough of `td.test_function`/
    `test_functions`).
    """
    variant = cfg.grid.preset or DEFAULT_PRESET
    preset = PRESETS.get(f"{cfg.molecule}:{variant}")
    if preset is None:
        return cfg

    # A levels-only run has no energy sweep to fill: `resonance_levels` needs
    # a molecule and two nuclear grids, nothing else. Leave `energies` None so
    # the runner can tell "not requested" from "requested but unresolved".
    needs_energies = any(obs.kind != "resonance_levels" for obs in cfg.observables)
    energies = cfg.energies
    if energies is None and needs_energies:
        energies = preset.default_energies
    observables = tuple(
        obs
        if obs.channels is not None
        else replace(obs, channels=_default_channels(obs.kind, preset))
        for obs in cfg.observables
    )
    td = cfg.td
    if td is not None:
        incident = td.incident if td.incident is not None else preset.default_incident
        kinds = _kinds_in(observables)
        test_functions = dict(td.test_functions or {})
        for kind in kinds:
            tf = test_functions.get(kind) or _test_function_for_kind(td, preset, kind, kinds)
            if tf is not None:
                test_functions[kind] = tf
        td = replace(td, incident=incident, test_functions=test_functions or None)

    # `nrm` has no per-molecule preset data -- its defaults are the measured
    # ones on `NrmSpec` itself. Materializing the block here is what puts the
    # discrete-state choice and `n_states` actually used into
    # `config.resolved.yaml`, rather than leaving a bare `nrm: null`.
    nrm = cfg.nrm
    if nrm is None and "nrm" in cfg.methods:
        nrm = NrmSpec()

    return replace(cfg, energies=energies, observables=observables, td=td, nrm=nrm)
