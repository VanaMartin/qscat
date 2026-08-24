"""The committed figure set for the TIME-DEPENDENT nonlocal resonance model.

    uv run --no-sync python -m validation.diatomic.td_nrm_figures            # all six
    uv run --no-sync python -m validation.diatomic.td_nrm_figures launch-rank  # one
    uv run --no-sync python -m validation.diatomic.td_nrm_figures --smoke      # cheap

`qscat.core.nrm`'s time-dependent route (`extended.py`, `propagation.py`,
`td_cross_section.py`) resums Gertitschke & Domcke, Phys. Rev. A 47, 1031
(1993) Eq. (2.1)'s memory integral into a time-LOCAL propagation under a
sparse arrow block Hamiltonian, and half-Fourier-transforms the propagated
packet back to the SAME `Psi_d(R;E)` the time-independent route
(Houfek/Rescigno/McCurdy, PRA 77, 012710 (2008) Eq. 52) solves for per
energy. That identity is gated in `libs/qscat/tests/`; these six figures are
what it LOOKS like, which is a different deliverable -- a reader deciding by
eye whether the two routes agree and whether the dynamics are physical.

Each figure writes `docs/physics/figures/<name>.png` and a sibling
`<name>.npz` holding the arrays behind it, exactly as
`validation/diatomic/ve_nrm_figure.py` does.

    td-vs-ti      f2-da-nrm-td-vs-ti.png     + f2-da-nrm-td-packet.png
    convergence   nrm-td-convergence.png
    launch-rank   nrm-td-launch-rank.png
    truncation    nrm-td-truncation-diverges.png
    vector        n2-nrm-td-vs-ti-vector.png

`td-vs-ti` renders TWO figures from ONE propagation on purpose: the cross
section and the packet diagnostics are the same run seen two ways, and
running it twice would cost another ~40 minutes for nothing.

**`td-vs-ti` HAS NEVER BEEN RUN TO COMPLETION, so its two PNGs are absent
from `docs/physics/figures/` rather than stale.** Three attempts were killed
for memory (see COST below). Do not judge the code by that: the path is
exercised end to end by `--smoke`, and every other figure here is a real
measurement. It also cannot be cheapened into fitting -- at `T = 40` instead
of `T = 12000` the transform is so far from converged that `sigma_TD` reads
~1e-32 against a `sigma_TI` of ~1, i.e. a shortened run does not produce a
noisy figure, it produces an empty one.

COST (12-core dev machine, `OMP_NUM_THREADS=8`): `td-vs-ti` ~50 min, of
which 39 min is the fine-deck propagation (6000 steps at 0.39 s on a
53570-square `H_ext`) and ~7 min the coarse-deck overlay; `vector` ~3 min;
`launch-rank` ~1 min; `truncation` ~3 min; `convergence` is instant (it
plots recorded numbers, see `_N2_BUDGET`/`_F2_SIGMA`). Run them ONE AT A
TIME with a pinned thread count -- concurrent unpinned sweeps cost ~300x
here (`docs/physics/optimization-targets.md`).

MEMORY IS THE BINDING CONSTRAINT, NOT TIME, AND `td-vs-ti` DOES NOT FIT ON A
LAPTOP. What dominates is the SPARSE LU OF THE PADE DENOMINATORS: an order-3
diagonal Pade stepper factors three shifted copies of `H_ext` and holds all
three factorizations for the whole propagation. `H_ext` is
`(1 + n_states) * N_R` square, so F2's DA deck is 53570 with the reduced
55-point electronic grid used here and 128568 with the production 132-point
one -- and SuperLU's fill-in on these complex-symmetric ECS matrices is what
sets peak RSS, not the 9e5 stored nonzeros. Measured 2026-08-24: three
attempts at the 53570 case were killed by the OS on a 12-core laptop, one of
them while a 6210-square dense `eigvals` (~0.6 GB) ran alongside it. Treat
`td-vs-ti` as a batch job for a machine with room, run it ALONE, and prefer
the MUMPS backend (`qscat[mumps]`, ~9x lower peak RSS at 143k unknowns --
`docs/physics/mumps-sparse-backend.md`) over the SuperLU fallback a bare Mac
gives you. `truncation`'s dense `numpy.linalg.eigvals` is the other memory
item and is why its fixture is kept small.

EVERY PROPAGATION USES THE COMPLETE ARM SET (`n_states=None`). Truncating it
leaves `H_ext` with growing eigenmodes, the transform's premise fails, and
`psi_d` comes back exponentially wrong rather than under-converged. The one
exception is `truncation`, whose whole subject is that failure.

`validation/` may import `qscat` and `projects`; the reverse is forbidden.
"""

from __future__ import annotations

import sys
import time
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt
from qscat.core.grids import electronic_grid, nuclear_grid, segmented_grid
from qscat.core.nrm import (
    AsymptoticDiscreteState,
    NrmIngredients,
    PhysicalDiscreteState,
    nrm_da_cross_section,
    nrm_ingredients,
)
from qscat.core.nrm.discrete_state import DiscreteState

# `da_sigma_from_psi` and `_boundary_node` are how `td_nrm_da_cross_section`
# itself turns a propagated packet into `sigma_DA` (Eq. 54). `td-vs-ti` needs
# the `TdNrmResult` as well as the cross section -- the packet diagnostics ARE
# figure 2 -- so it drives the same four calls that entry point makes rather
# than calling it and then propagating a second time for the diagnostics.
# Importing the private `_boundary_node` keeps `eps_e` byte-identical to both
# routes' own; recomputing it here is exactly the kind of second definition
# that would make the comparison quietly meaningless.
from qscat.core.nrm.dissociation import _boundary_node, da_sigma_from_psi
from qscat.core.nrm.extended import LaunchBasis, extended_hamiltonian, initial_packet
from qscat.core.nrm.nonlocal_potential import continue_to_tail
from qscat.core.nrm.propagation import TdNrmResult, propagate_nrm
from qscat.core.vibrational import vibrational_states
from qscat.dvr import FemDvrEcsGrid, dvr_interpolation_matrix
from qscat.model import F2, N2, ResonanceModel

__all__ = [
    "FIGURE_DIR",
    "Deck",
    "main",
    "render_convergence",
    "render_launch_rank",
    "render_td_vs_ti",
    "render_truncation",
    "render_vector",
]

_REPO_ROOT = Path(__file__).resolve().parents[2]
FIGURE_DIR = _REPO_ROOT / "docs" / "physics" / "figures"

_V_INIT = 0

# --- decks ------------------------------------------------------------------
#
# The propagation decks are the shipped test fixtures', transcribed here so the
# figures show the same runs the gates measure. `validation/` may import
# `qscat`, but `libs/qscat/tests` is not an importable package, so a literal
# copy is the only way to share them; each block names its source.

# `libs/qscat/tests/test_nrm_td_cross_section.py::f2_deck` -- the eMoScat F2
# NUCLEAR deck (`reference/eMoScat/input/F2/grids.txt`, 2nd declaration), 974
# points, 65 points/bohr over R = 3-10.7. It has to be this fine: F2's
# dissociating wave carries `K_R = 55.6-64.2` over E = 0.02-0.05 Ha
# (wavelength 0.098-0.113 bohr), and a coarse deck cannot represent it.
_F2_FINE_REAL = ((9, 1.8), (1, 2.0), (5, 2.5), (4, 2.596908), (4, 2.7), (40, 10.7))
_F2_FINE_COMPLEX = (
    (1, 10.8),
    (1, 11.0),
    (1, 11.5),
    (1, 12.5),
    (1, 14.0),
    (1, 18.0),
    (4, 30.0),
    (2, 101.0),
)

# The deliberately COARSE F2 nuclear deck of figure 2's dashed overlay: N2-style
# element counts (~15 points/bohr over the dissociation region) stretched onto
# F2's own R range, so `R0`, the electronic grid, and the discrete state are
# unchanged and the ONLY difference from `_F2_FINE_*` is nuclear resolution.
# Built here rather than borrowed: no shipped deck is this coarse for F2,
# because no shipped calculation may use one.
_F2_COARSE_REAL = ((3, 1.8), (8, 3.0), (2, 4.0), (6, 10.7))
_F2_COARSE_COMPLEX = ((1, 11.0), (1, 12.5), (1, 14.0), (3, 30.0))

# `libs/qscat/tests/test_nrm_td_cross_section.py::small_deck` -- small enough
# for a DENSE `numpy.linalg.eigvals` of the complete `H_ext`, which is what
# figure 5 needs and ARPACK cannot supply (it under-reports `max Im` by ~6x on
# these strongly non-normal matrices).
#
# THE ECS TAIL IS REFINED against that fixture's own `((2, 20.0),)`, and it has
# to be. Measured here 2026-08-24: with two tail elements spanning 6 -> 20 bohr
# the COMPLETE arm set's `max Im E` reads +6.98e-05 -- four orders above the
# +2.4e-12 the production ingredients give, and enough to make figure 5's own
# point look false. Refining the tail to four elements drops it to +4.15e-08
# (six, to +4.71e-08 -- converged), i.e. it was a discretisation artefact of an
# under-resolved rotated continuum, not a property of the complete sum. The
# offending eigenvalue sits at `Re E ~ +666 Ha`, a high-energy tail mode with
# nothing to do with the packet. What does NOT move is the truncated set:
# `n_states=3` gives +3.305e-04 at `Re E = -0.0754 Ha` on ALL THREE tails, in
# the energy range the packet actually occupies. That contrast -- artefact
# moves under refinement, real growing mode does not -- is why the figure can
# be trusted at all.
_F2_SMALL_REAL = ((3, 2.5), (4, 6.0))
_F2_SMALL_COMPLEX = ((1, 7.0), (1, 9.0), (1, 14.0), (1, 20.0))

# `libs/qscat/tests/test_nrm_propagation.py::n2_deck` -- the fixture the
# vector-to-vector gate (rel = 1.7264e-04) runs on.
_N2_GATE_REAL = ((3, 1.5), (8, 3.0), (2, 4.0), (4, 8.0))
_N2_GATE_COMPLEX = ((3, 20.0),)

# `libs/qscat/tests/test_nrm_extended.py`'s F2 fixture -- the deck the launch
# ranks of `extended.py`'s module docstring were measured on.
_F2_RANK_REAL = ((9, 1.8), (1, 2.0), (5, 2.5), (4, 2.7), (20, 10.7))
_F2_RANK_COMPLEX = ((1, 11.0), (1, 12.5), (1, 14.0), (3, 30.0))


@dataclass(frozen=True)
class Deck:
    """The five things every propagation in this module needs."""

    label: str
    model: ResonanceModel
    nuc: FemDvrEcsGrid
    elec: FemDvrEcsGrid
    phi_d: DiscreteState
    ing: NrmIngredients
    eps: npt.NDArray[np.float64]
    chi: npt.NDArray[np.complex128]

    @property
    def eps_e(self) -> float:
        """The dissociation-channel asymptote `V_d(X)` at the flux surface --
        `nrm_da_cross_section`'s own, read at the same node."""
        v_d_full = continue_to_tail(self.ing.v_d_discrete, self.ing.R, self.nuc)
        return float(v_d_full[_boundary_node(self.nuc)].real)


def _deck(
    label: str,
    model: ResonanceModel,
    nuc: FemDvrEcsGrid,
    elec: FemDvrEcsGrid,
    n_vib: int = 4,
    phi_d: DiscreteState | None = None,
) -> Deck:
    """Assemble a `Deck` -- choice B (`AsymptoticDiscreteState`) unless a
    discrete state is passed in.

    `R_inf` is a NUCLEAR coordinate (`nuc.R0`), never the electronic grid's.
    `nrm_ingredients` requires strictly DESCENDING `R`.
    """
    ds = AsymptoticDiscreteState(elec, model, R_inf=nuc.R0) if phi_d is None else phi_d
    r_desc = np.sort(nuc.points[nuc.points.imag == 0.0].real)[::-1]
    t0 = time.time()
    ing = nrm_ingredients(elec, model, ds, r_desc)
    eps, chi = vibrational_states(nuc, model.mu, n_vib, model.v0)
    print(
        f"  deck {label}: N_R={nuc.n} N_elec={elec.n} states={ing.E_n.shape[1]} "
        f"H_ext={(1 + ing.E_n.shape[1]) * nuc.n} ({time.time() - t0:.1f} s)",
        flush=True,
    )
    return Deck(label, model, nuc, elec, ds, ing, eps, chi)


def f2_fine_deck() -> Deck:
    """F2 on the production nuclear deck with the REDUCED (55-point)
    electronic grid of the shipped TD gate.

    The reduction is legitimate because this comparison is DIFFERENTIAL --
    both routes run on the same ingredients on the same grids -- but the
    absolute `sigma_DA` it produces is NOT the converged F2 cross section and
    must not be quoted as one (`validation/diatomic` owns that; see
    `test_nrm_td_cross_section.py::f2_deck`).
    """
    nuc = segmented_grid(_F2_FINE_REAL, _F2_FINE_COMPLEX, angle_deg=35.0, quadrature=14)
    return _deck("F2 fine", F2, nuc, electronic_grid(r_max=13.0, order=5, n_complex=2))


def f2_coarse_deck() -> Deck:
    """The same molecule, discrete state and electronic grid as
    `f2_fine_deck`, on a nuclear grid that cannot resolve the exit wave."""
    nuc = segmented_grid(_F2_COARSE_REAL, _F2_COARSE_COMPLEX, angle_deg=35.0, quadrature=10)
    return _deck("F2 coarse", F2, nuc, electronic_grid(r_max=13.0, order=5, n_complex=2))


def f2_small_deck() -> Deck:
    """A small F2 fixture -- nothing here is physically converged; it exists
    so a dense eigensolver can see the whole spectrum of `H_ext`."""
    nuc = segmented_grid(_F2_SMALL_REAL, _F2_SMALL_COMPLEX, angle_deg=35.0, quadrature=8)
    return _deck("F2 small", F2, nuc, electronic_grid(r_max=11.0, order=6, n_complex=2), n_vib=3)


def n2_gate_deck() -> Deck:
    """The N2 fixture of the vector-to-vector propagation gate."""
    nuc = segmented_grid(_N2_GATE_REAL, _N2_GATE_COMPLEX, angle_deg=35.0, quadrature=10)
    return _deck("N2 gate", N2, nuc, electronic_grid(r_max=11.0, order=6, n_complex=3))


# --- plotting helpers -------------------------------------------------------


def _plt() -> Any:
    """matplotlib with a non-interactive backend and readable defaults.

    No `style.use` of anything that could reach the network, and no font
    package assumptions -- only rcParams matplotlib ships with.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.size": 10.0,
            "axes.titlesize": 11.0,
            "axes.labelsize": 10.5,
            "legend.fontsize": 9.0,
            "figure.titlesize": 12.5,
            "axes.grid": True,
            "grid.alpha": 0.25,
            "lines.linewidth": 1.4,
        }
    )
    return plt


def _save(fig: Any, name: str, **arrays: Any) -> Path:
    """Write `<name>.png` and `<name>.npz` into `FIGURE_DIR`."""
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    png = FIGURE_DIR / f"{name}.png"
    fig.savefig(png, dpi=150, bbox_inches="tight")
    np.savez(png.with_suffix(".npz"), **arrays)
    print(f"  wrote {png} and {png.with_suffix('.npz').name}", flush=True)
    return png


def _report(warned: list[warnings.WarningMessage]) -> None:
    """Print the warnings a run raised instead of swallowing them.

    Several runs here are EXPECTED to warn -- F2's `S(t)` plateaus above
    `unabsorbed_tol` for spectral reasons, and figure 5's truncated run is
    supposed to diverge -- so they are recorded and printed rather than
    filtered away, which would also hide an unexpected one.
    """
    for w in warned:
        print(f"    [warning] {str(w.message).splitlines()[0][:160]}", flush=True)


# --- figures 1 and 2: the cross section, and the packet that produced it ----

# Nine energies spanning the F2 DA window. Explicit rather than `linspace` so
# the three shipped gate anchors (0.02 / 0.03 / 0.05) are among them and the
# figure can be read against `test_nrm_td_cross_section.py`'s recorded ratios.
# Nine cost almost nothing over three: `H_ext` is energy-independent, so the
# whole batch is ONE propagation of the launch basis's few singular-vector
# columns (rank 4 here) and only the launch SVD and the extraction scale with
# the energy count.
TD_VS_TI_ENERGIES = np.array([0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040, 0.050, 0.060])

#: The energy figure 2 shows the dynamics of.
PACKET_ENERGY = 0.030

#: `dt = 2`, `T = 12000` -- the shipped gate's settings, and both are measured
#: rather than chosen. `T = 8000` would need a ~1.5e-1 tolerance and `T = 4000`
#: a ~3e-1 one; beyond `T = 12000` the ratios stop improving (residual
#: oscillation, not truncation). At order 3 the `dt` error falls as `dt^6`, and
#: a `dt = 1` rerun reproduces the `dt = 2` ratios to 3-4 significant figures.
TD_DT, TD_STEPS = 2.0, 6000


def _propagate(
    deck: Deck, energies: npt.NDArray[np.float64], dt: float, n_steps: int
) -> tuple[TdNrmResult, LaunchBasis, npt.NDArray[np.float64]]:
    """Propagate `deck`'s launch basis and extract `sigma_DA` per energy.

    This is `td_nrm_da_cross_section`'s body with the `TdNrmResult` kept:
    same `initial_packet`, same `extended_hamiltonian`, same `propagate_nrm`,
    same `da_sigma_from_psi` at the same `eps_e`.

    The one thing it does NOT reproduce is that function's `unabsorbed_tol`
    warning, because on F2 it fires on every converged run -- `S(t)` plateaus
    at 0.006-0.009 for spectral reasons, above the 1e-2 default, which is
    documented behaviour and not a signal here. The fraction still left in the
    real region is printed per energy instead, and lands in the `.npz`.
    """
    t0 = time.time()
    launch = initial_packet(
        deck.nuc, deck.elec, deck.model, deck.phi_d, deck.ing, deck.eps, deck.chi, _V_INIT, energies
    )
    h_ext = extended_hamiltonian(deck.ing, deck.nuc, deck.model)
    print(
        f"  {deck.label}: launch rank {launch.rank} (truncation "
        f"{launch.truncation_error:.2e}), H_ext {h_ext.shape[0]}, "
        f"T={dt * n_steps:g} at dt={dt:g}",
        flush=True,
    )
    with warnings.catch_warnings(record=True) as warned:
        warnings.simplefilter("always")
        res = propagate_nrm(h_ext, launch, deck.nuc, dt=dt, n_steps=n_steps, order=3)
    _report(warned)

    eps_e = deck.eps_e
    sigma = np.array(
        [
            da_sigma_from_psi(
                deck.nuc,
                deck.model.mu,
                res.psi_d[:, j],
                float(launch.e_total[j]),
                eps_e,
                float(energies[j]),
            )
            for j in range(energies.size)
        ]
    )
    print(f"  {deck.label}: propagated in {time.time() - t0:.0f} s", flush=True)
    return res, launch, sigma


def render_td_vs_ti(smoke: bool = False) -> dict[str, Any]:
    """Figures 1 and 2 -- `sigma_DA(E)` two ways, and the packet behind it.

    One propagation on the fine deck serves both: figure 1 reads its
    transform (the cross section at nine energies), figure 2 reads its
    diagnostics at `PACKET_ENERGY`. A second, cheap propagation on
    `f2_coarse_deck` supplies figure 2's dashed overlay -- the same physics
    on a nuclear grid that cannot represent the `K_R ~ 58` exit wave.
    """
    plt = _plt()
    dt, n_steps = (TD_DT, 20) if smoke else (TD_DT, TD_STEPS)
    fine = f2_fine_deck()
    res, launch, sigma_td = _propagate(fine, TD_VS_TI_ENERGIES, dt, n_steps)

    t0 = time.time()
    sigma_ti = np.asarray(
        nrm_da_cross_section(
            fine.nuc,
            fine.elec,
            F2,
            fine.phi_d,
            fine.eps,
            fine.chi,
            _V_INIT,
            TD_VS_TI_ENERGIES,
            ingredients=fine.ing,
        ),
        dtype=np.float64,
    )
    print(f"  TI oracle: {time.time() - t0:.0f} s", flush=True)
    ratio = sigma_td / sigma_ti

    print("\n  E (Ha)   sigma_TI (a0^2)   sigma_TD (a0^2)   TD/TI    S(T)/S(0)")
    left = res.unabsorbed / res.survival[0]
    for j, e in enumerate(TD_VS_TI_ENERGIES):
        print(
            f"  {e:6.3f}   {sigma_ti[j]:15.6e}   {sigma_td[j]:15.6e}   "
            f"{ratio[j]:6.4f}   {left[j]:.3f}",
            flush=True,
        )
    print(f"  worst |TD/TI - 1| = {np.max(np.abs(ratio - 1.0)):.4f}\n", flush=True)

    fig, (ax0, ax1) = plt.subplots(
        2, 1, figsize=(7.2, 6.6), sharex=True, gridspec_kw={"height_ratios": [2.0, 1.0]}
    )
    ax0.semilogy(
        TD_VS_TI_ENERGIES, sigma_ti, "-", color="black", label=r"TI  $-$ per-energy Eq. (52) solve"
    )
    ax0.semilogy(
        TD_VS_TI_ENERGIES,
        sigma_td,
        "o",
        color="tab:red",
        markerfacecolor="none",
        markeredgewidth=1.4,
        markersize=7.0,
        label=rf"TD  $-$ one propagation, $T={dt * n_steps:g}$, $\Delta t={dt:g}$ (a.u. of time)",
    )
    ax0.set_ylabel(r"$\sigma_{\mathrm{DA}}$  ($a_0^2$)")
    ax0.legend(loc="lower left")
    ax0.set_title(
        rf"one propagation, {TD_VS_TI_ENERGIES.size} energies "
        rf"(launch rank {launch.rank} of {TD_VS_TI_ENERGIES.size})"
    )

    ax1.axhspan(0.95, 1.05, color="tab:blue", alpha=0.12, label=r"$\pm 5\%$")
    ax1.axhline(1.0, color="black", linewidth=0.9)
    ax1.plot(TD_VS_TI_ENERGIES, ratio, "o-", color="tab:red", markersize=5.0)
    ax1.set_ylabel(r"$\sigma_{\mathrm{TD}} / \sigma_{\mathrm{TI}}$")
    ax1.set_xlabel(r"incident electron energy  $E$  (hartree)")
    # At least +/-10%, but never clip a point: a ratio outside the window
    # would otherwise vanish and the panel would read as agreement.
    ax1.set_ylim(min(0.90, float(ratio.min()) * 0.98), max(1.10, float(ratio.max()) * 1.02))
    ax1.legend(loc="upper right")

    fig.suptitle(
        "F$_2$ dissociative attachment, nonlocal resonance model\n"
        r"discrete state B (asymptotic, $R$-independent), $v_i = 0$, "
        "complete arm set: time-dependent vs time-independent"
    )
    fig.tight_layout()
    _save(
        fig,
        "f2-da-nrm-td-vs-ti",
        energies=TD_VS_TI_ENERGIES,
        sigma_ti=sigma_ti,
        sigma_td=sigma_td,
        ratio=ratio,
        unabsorbed_fraction=left,
        dt=np.asarray(dt),
        n_steps=np.asarray(n_steps),
        launch_rank=np.asarray(launch.rank),
        launch_truncation_error=np.asarray(launch.truncation_error),
    )
    plt.close(fig)

    coarse = f2_coarse_deck()
    res_c, _, sigma_c = _propagate(coarse, TD_VS_TI_ENERGIES, dt, n_steps)
    _render_packet(plt, fine, res, coarse, res_c, dt * n_steps)
    return {
        "energies": TD_VS_TI_ENERGIES,
        "sigma_ti": sigma_ti,
        "sigma_td": sigma_td,
        "ratio": ratio,
        "sigma_td_coarse": sigma_c,
    }


def _k_dissociation(deck: Deck, e_kin: float) -> float:
    """`K_R = sqrt(2 mu (E_tot - eps_e))` -- the momentum the exit wave
    carries once it is past the interaction region, i.e. what `<P>_t` must
    rise to if the packet is genuinely dissociating."""
    e_total = float(e_kin) + float(deck.eps[_V_INIT])
    return float(np.sqrt(2.0 * deck.model.mu * (e_total - deck.eps_e)))


def _render_packet(
    plt: Any, fine: Deck, res: TdNrmResult, coarse: Deck, res_c: TdNrmResult, t_max: float
) -> None:
    """Figure 2 -- survival, centroid and momentum, fine deck vs coarse."""
    j = int(np.argmin(np.abs(TD_VS_TI_ENERGIES - PACKET_ENERGY)))
    e_kin = float(TD_VS_TI_ENERGIES[j])
    k_r = _k_dissociation(fine, e_kin)

    t, t_c = res.time, res_c.time
    s = res.survival[:, j] / res.survival[0, j]
    s_c = res_c.survival[:, j] / res_c.survival[0, j]
    r_t, r_c = res.centroid[:, j], res_c.centroid[:, j]
    p_t, p_c = res.momentum[:, j], res_c.momentum[:, j]

    print(
        f"\n  packet at E = {e_kin:.3f} Ha:  K_R = {k_r:.1f}\n"
        f"    fine   <R>: {r_t[0]:.2f} -> {np.nanmax(r_t):.2f} bohr,  "
        f"<P> max {np.nanmax(p_t):.1f},  S(T)/S(0) = {s[-1]:.3g}\n"
        f"    coarse <R>: {r_c[0]:.2f} -> {np.nanmax(r_c):.2f} bohr,  "
        f"<P> max {np.nanmax(p_c):.1f},  S(T)/S(0) = {s_c[-1]:.3g}\n",
        flush=True,
    )

    fig, axes = plt.subplots(3, 1, figsize=(7.2, 8.4), sharex=True)
    fine_kw = {"color": "tab:blue", "linestyle": "-"}
    coarse_kw = {"color": "tab:orange", "linestyle": "--"}
    fine_label = f"fine nuclear deck ({fine.nuc.n} pts, 65 pts/bohr)"
    coarse_label = f"coarse nuclear deck ({coarse.nuc.n} pts, ~15 pts/bohr)"

    axes[0].semilogy(t, s, label=fine_label, **fine_kw)
    axes[0].semilogy(t_c, s_c, label=coarse_label, **coarse_kw)
    axes[0].set_ylabel(r"$S(t)\,/\,S(0)$")
    axes[0].set_title("survival in the real nuclear region (PRA 47 Eq. 4.4)")
    axes[0].legend(loc="lower left")

    axes[1].plot(t, r_t, **fine_kw)
    axes[1].plot(t_c, r_c, **coarse_kw)
    axes[1].axhline(
        fine.nuc.R0,
        color="black",
        linestyle=":",
        linewidth=1.2,
        label=rf"ECS boundary $R_0 = {fine.nuc.R0:g}$ bohr",
    )
    axes[1].set_ylabel(r"$\langle R \rangle_t$  (bohr)")
    axes[1].set_title("centroid (Eq. 4.5): the fine deck dissociates, the coarse deck does not")
    axes[1].legend(loc="lower right")

    axes[2].plot(t, p_t, **fine_kw)
    axes[2].plot(t_c, p_c, **coarse_kw)
    axes[2].axhline(
        k_r,
        color="black",
        linestyle=":",
        linewidth=1.2,
        label=rf"$K_R = {k_r:.1f}$ (a.u. of momentum)",
    )
    axes[2].set_ylabel(r"$\langle P \rangle_t$  (a.u.)")
    axes[2].set_xlabel("time  $t$  (a.u. of time)")
    axes[2].set_title("momentum (Eq. 4.6): the exit wave's own wavenumber, or nothing")
    axes[2].legend(loc="lower right")
    axes[2].set_xlim(0.0, t_max)

    fig.suptitle(
        rf"F$_2$ dissociative attachment, nonlocal resonance model, $E = {e_kin:.3f}$ Ha"
        "\n"
        r"discrete state B (asymptotic), $v_i = 0$, complete arm set: "
        "the packet the TD route propagates"
    )
    fig.tight_layout()
    _save(
        fig,
        "f2-da-nrm-td-packet",
        energy=np.asarray(e_kin),
        k_dissociation=np.asarray(k_r),
        ecs_boundary=np.asarray(fine.nuc.R0),
        time_fine=t,
        survival_fine=s,
        centroid_fine=r_t,
        momentum_fine=p_t,
        time_coarse=t_c,
        survival_coarse=s_c,
        centroid_coarse=r_c,
        momentum_coarse=p_c,
        n_nuclear_fine=np.asarray(fine.nuc.n),
        n_nuclear_coarse=np.asarray(coarse.nuc.n),
    )
    plt.close(fig)


# --- figure 3: convergence, from RECORDED measurements ----------------------
#
# Everything below is TRANSCRIBED, not recomputed. Re-measuring it would cost
# hours (the F2 column alone is six propagations of the 39-minute run) and add
# nothing: these are the numbers the shipped gates were set from.
#
# N2 vector-gate `rel` and `S(T)/S(0)` vs `T` at `dt = 1`, on
# `n2_gate_deck()`: `.superpowers/sdd/2026-08-19-nrm-td/task-4-report.md`,
# "Refinement tables on the converging fixture" (measured 2026-08-19; the
# T = 4000 row is the shipped gate, 1.7264e-04).
_N2_T = np.array([500.0, 1000.0, 2000.0, 3000.0, 4000.0, 5000.0])
_N2_REL = np.array([2.8955e-01, 4.0175e-02, 3.9963e-03, 6.4580e-04, 1.7264e-04, 1.4390e-04])
_N2_SURVIVAL = np.array([1.08e-01, 4.96e-03, 1.05e-04, 2.63e-06, 6.68e-08, 1.69e-09])
# The error budget fitted to those rows in the same report: the truncation and
# propagation errors are separable and add IN QUADRATURE.
_N2_BUDGET = (0.40, 1.43e-4)  # truncation(T) = 0.40*sqrt(S(T)/S(0)); propagation(dt=1)

# F2 DA worst-of-three-energies `max |sigma_TD/sigma_TI - 1|` and `S(T)/S(0)`
# vs `T` at `dt = 2`, on `f2_fine_deck()`:
# `.superpowers/sdd/2026-08-19-nrm-td/task-5-report.md`, "Convergence in T"
# (measured 2026-08-19; T = 12000 is the shipped gate's setting).
_F2_T = np.array([4000.0, 6000.0, 8000.0, 10000.0, 12000.0, 14000.0])
_F2_SIGMA = np.array([0.29, 0.13, 0.065, 0.024, 0.014, 0.022])
_F2_SURVIVAL = np.array([0.25, 0.045, 0.018, 0.011, 0.009, 0.008])


def render_convergence(smoke: bool = False) -> dict[str, Any]:
    """Figure 3 -- what `T` buys, on the two observables that were measured.

    LEFT, N2: the vector-to-vector error against the time-independent solve
    falls four decades and then stops, because it hits the `dt = 1`
    propagation floor. The dashed curve is not a fit to this panel's points
    -- it is the two-parameter budget of `_N2_BUDGET` evaluated at the
    measured `S(T)/S(0)`, and it reproduces every row.

    RIGHT, F2 DA: the cross section converges while the norm does NOT. `S(t)`
    plateaus at 0.006-0.009 because F2's `V_d(R)` well holds >= 24 near-real
    modes, so no absolute survival floor is ever met on this molecule --
    convergence has to be read off `sigma_DA` being stationary in `T`, which
    is what the left-hand axis shows.
    """
    plt = _plt()
    c, floor = _N2_BUDGET
    truncation = c * np.sqrt(_N2_SURVIVAL)
    budget = np.sqrt(truncation**2 + floor**2)

    print("\n  N2 vector gate (recorded, dt=1):  T / rel / budget")
    for tt, rel, bud in zip(_N2_T, _N2_REL, budget, strict=True):
        print(f"    {tt:7.0f}  {rel:.4e}  {bud:.4e}", flush=True)
    print("  F2 DA (recorded, dt=2):  T / max|ratio-1| / S(T)/S(0)")
    for tt, sg, sv in zip(_F2_T, _F2_SIGMA, _F2_SURVIVAL, strict=True):
        print(f"    {tt:7.0f}  {sg:.3f}  {sv:.3f}", flush=True)

    fig, (axl, axr) = plt.subplots(1, 2, figsize=(11.0, 4.8))

    axl.loglog(_N2_T, _N2_REL, "o", color="tab:blue", markersize=7.0, label="measured")
    axl.loglog(
        _N2_T,
        budget,
        "--",
        color="tab:red",
        label=rf"$\sqrt{{\mathrm{{trunc}}^2 + \mathrm{{prop}}^2}}$,  "
        rf"trunc $= {c:.2f}\sqrt{{S(T)/S(0)}}$",
    )
    axl.axhline(
        floor,
        color="black",
        linestyle=":",
        linewidth=1.2,
        label=rf"propagation floor at $\Delta t = 1$: ${floor:.2e}$",
    )
    axl.set_xlabel("propagation time  $T$  (a.u. of time)")
    axl.set_ylabel(r"$\|\Psi_d^{\mathrm{TD}} - \Psi_d^{\mathrm{TI}}\| / \|\Psi_d^{\mathrm{TI}}\|$")
    axl.set_title(r"N$_2$ VE: vector-to-vector error, $\Delta t = 1$")
    axl.legend(loc="lower left")
    # The default log locator crowds 6e2/1e3/2e3/3e3/4e3 into one another;
    # label the measured points themselves instead.
    ticks = [500.0, 1000.0, 2000.0, 5000.0]  # 3000/4000 would collide with 5000
    axl.set_xticks(ticks)
    axl.set_xticklabels([f"{t:.0f}" for t in ticks])
    axl.set_xticks([], minor=True)
    # The budget is an ASYMPTOTIC fit and the figure should not hide where it
    # stops holding: it reproduces T >= 2000 to ~3% but under-predicts the two
    # shortest runs (2.2x at T=500, 1.4x at T=1000), where most of the packet
    # is still inside the box and `truncation ~ sqrt(S)` is not yet in its
    # regime. Reported as measured; see this module's report.
    axl.annotate(
        "budget under-predicts\nwhile the packet is\nstill in the box",
        xy=(float(_N2_T[0]), float(_N2_REL[0])),
        xytext=(0.30, 0.86),
        textcoords="axes fraction",
        fontsize=8.5,
        ha="left",
        arrowprops={"arrowstyle": "->", "lw": 0.9, "color": "0.35"},
    )

    axr.semilogy(
        _F2_T,
        _F2_SIGMA,
        "o-",
        color="tab:blue",
        markersize=7.0,
        label=r"$\max|\sigma_{TD}/\sigma_{TI} - 1|$",
    )
    axr.set_xlabel("propagation time  $T$  (a.u. of time)")
    axr.set_ylabel(r"$\max_E |\sigma_{\mathrm{TD}}/\sigma_{\mathrm{TI}} - 1|$", color="tab:blue")
    axr.tick_params(axis="y", labelcolor="tab:blue")
    axr.set_title(r"F$_2$ DA: the cross section converges, the norm plateaus")

    axr2 = axr.twinx()
    axr2.semilogy(
        _F2_T, _F2_SURVIVAL, "s--", color="tab:orange", markersize=6.0, label=r"$S(T)/S(0)$"
    )
    axr2.set_ylabel(r"$S(T)/S(0)$  (unabsorbed norm)", color="tab:orange")
    axr2.tick_params(axis="y", labelcolor="tab:orange")
    axr2.grid(False)
    lines = axr.get_lines() + axr2.get_lines()
    axr.legend(lines, [ln.get_label() for ln in lines], loc="upper right")

    fig.suptitle(
        "Time-dependent nonlocal resonance model: convergence in the propagation time\n"
        "(recorded measurements, 2026-08-19 -- see this module's source for provenance)"
    )
    fig.tight_layout()
    _save(
        fig,
        "nrm-td-convergence",
        n2_time=_N2_T,
        n2_rel=_N2_REL,
        n2_survival=_N2_SURVIVAL,
        n2_budget=budget,
        n2_budget_coefficient=np.asarray(c),
        n2_propagation_floor=np.asarray(floor),
        f2_time=_F2_T,
        f2_sigma_error=_F2_SIGMA,
        f2_survival=_F2_SURVIVAL,
    )
    plt.close(fig)
    return {"n2_budget": budget}


# --- figure 4: the launch basis is low-rank ---------------------------------

#: The two published energy windows: F2's DA panel (PRA 77 Fig. 9) and N2's VE
#: panels (Figs. 4 and 8).
RANK_WINDOWS = {
    "F$_2$ DA, 0.010-0.050 Ha": np.linspace(0.010, 0.050, 9),
    "N$_2$ VE, 0.060-0.160 Ha": np.linspace(0.060, 0.160, 9),
}

#: `initial_packet`'s own default -- the line figure 4 marks.
RANK_TOL = 1e-6


def _singular_values(deck: Deck, energies: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """`sigma_j / sigma_1` of the launch matrix `M[R, j] = V_dk(R;E_j) chi_v(R)`.

    Read back off an EXACT (`rank_tol = 0.0`) `LaunchBasis` rather than by
    rebuilding `M`: there `coeffs[m] = sigma_m * vh[m]` with `vh` rows of unit
    norm, so `sigma_m = ||coeffs[m]||`. Using `initial_packet` itself means
    the figure measures the shipped factorization, not a re-implementation of
    it. `n_states=0` keeps the embedding one block wide -- the SVD is over the
    nuclear axis and does not depend on how many arm blocks follow it.
    """
    launch = initial_packet(
        deck.nuc,
        deck.elec,
        deck.model,
        deck.phi_d,
        deck.ing,
        deck.eps,
        deck.chi,
        _V_INIT,
        energies,
        n_states=0,
        rank_tol=0.0,
    )
    sv = np.linalg.norm(launch.coeffs, axis=1)
    return sv / sv[0]


def render_launch_rank(smoke: bool = False) -> dict[str, Any]:
    """Figure 4 -- why one propagation can serve a whole energy sweep.

    PRA 47 Eq. (2.17) removes the launch state's energy dependence
    analytically, but only under Eq. (2.16)'s separability, which these models
    satisfy numerically rather than exactly. The figure measures how nearly:
    choice B's spectrum collapses within two or three singular values (Eq.
    2.17's own claim, near enough), while choice A's -- whose `phi_d` is
    R-dependent, so its launch matrix is further from separable -- decays
    several times more slowly and needs about twice the columns at the same
    tolerance.
    """
    plt = _plt()
    elec_f2 = electronic_grid(r_max=16.0, order=8, n_complex=6)
    elec_f2_b = electronic_grid(r_max=16.0, order=8, n_complex=6, angle_deg=40.0)
    nuc_f2 = segmented_grid(_F2_RANK_REAL, _F2_RANK_COMPLEX, angle_deg=45.0, quadrature=14)
    # `qscat.core.grids.nuclear_grid`'s N2 default -- the grid the recorded
    # ranks were measured on (task-2-report.md's caveat: it is not a
    # crossing-refined deck, so `nrm_ingredients`' adiabatic-tracking warning
    # fires on it; that is pre-existing and does not touch the SVD, which is
    # exact given whatever `M` it was handed).
    nuc_n2 = nuclear_grid(r_max=22.0, n_complex=6, quadrature=10)
    elec_n2 = electronic_grid(r_max=16.0, order=7, n_complex=5)
    elec_n2_b = electronic_grid(r_max=16.0, order=7, n_complex=5, angle_deg=40.0)

    series: dict[str, npt.NDArray[np.float64]] = {}
    for window, (model, nuc, elec, elec_b) in zip(
        RANK_WINDOWS,
        ((F2, nuc_f2, elec_f2, elec_f2_b), (N2, nuc_n2, elec_n2, elec_n2_b)),
        strict=True,
    ):
        energies = RANK_WINDOWS[window]
        r_desc = np.sort(nuc.points[nuc.points.imag == 0.0].real)[::-1]
        for choice, phi_d in (
            ("B (asymptotic)", AsymptoticDiscreteState(elec, model, R_inf=nuc.R0)),
            ("A (physical)", PhysicalDiscreteState(elec, model, r_desc, elec_b)),
        ):
            with warnings.catch_warnings(record=True) as warned:
                warnings.simplefilter("always")
                deck = _deck(f"{window} {choice}", model, nuc, elec, phi_d=phi_d)
                sv = _singular_values(deck, energies)
            _report(warned)
            series[f"{window}, choice {choice}"] = sv

    print("\n  launch-matrix singular values (sigma_j / sigma_1), 9 energies:")
    ranks: dict[str, int] = {}
    for label, sv in series.items():
        ranks[label] = int(np.sum(sv > RANK_TOL))
        print(f"    {label}: rank@1e-6 = {ranks[label]}", flush=True)
        print("      " + ", ".join(f"{v:.2e}" for v in sv[:6]), flush=True)

    fig, ax = plt.subplots(figsize=(7.6, 5.4))
    styles = {
        "F$_2$ DA, 0.010-0.050 Ha, choice B (asymptotic)": ("tab:blue", "o", "-"),
        "F$_2$ DA, 0.010-0.050 Ha, choice A (physical)": ("tab:blue", "s", "--"),
        "N$_2$ VE, 0.060-0.160 Ha, choice B (asymptotic)": ("tab:red", "o", "-"),
        "N$_2$ VE, 0.060-0.160 Ha, choice A (physical)": ("tab:red", "s", "--"),
    }
    for label, sv in series.items():
        color, marker, ls = styles[label]
        idx = np.arange(1, sv.size + 1)
        # A singular value that underflowed to exactly zero cannot be drawn on
        # a log axis; matplotlib would drop it silently, so clip it to the
        # bottom of the plotted range and say so in the axis limit.
        ax.semilogy(
            idx,
            np.maximum(sv, 1e-17),
            marker=marker,
            linestyle=ls,
            color=color,
            markersize=6.0,
            label=f"{label} (rank {ranks[label]})",
        )
    ax.axhline(
        RANK_TOL,
        color="black",
        linestyle=":",
        linewidth=1.3,
        label=rf"truncation tolerance $= {RANK_TOL:.0e}$ (`initial_packet`'s default)",
    )
    ax.set_xlabel("singular-value index  $j$")
    ax.set_ylabel(r"$\sigma_j / \sigma_1$")
    ax.set_ylim(1e-17, 3.0)
    ax.set_xticks(np.arange(1, 10))
    ax.legend(loc="lower left")
    ax.set_title(
        "Launch matrix  $M[R, E_j] = V_{dk}(R; E_j)\\,\\chi_{v=0}(R)$:\n"
        "nine energies, a handful of columns"
    )
    fig.tight_layout()
    _save(
        fig,
        "nrm-td-launch-rank",
        rank_tol=np.asarray(RANK_TOL),
        labels=np.asarray(list(series)),
        singular_values=np.vstack([series[k] for k in series]),
        ranks=np.asarray([ranks[k] for k in series]),
        f2_energies=RANK_WINDOWS["F$_2$ DA, 0.010-0.050 Ha"],
        n2_energies=RANK_WINDOWS["N$_2$ VE, 0.060-0.160 Ha"],
    )
    plt.close(fig)
    return {"ranks": ranks, "singular_values": series}


# --- figure 5: truncating the arm set breaks the transform ------------------

#: `n_states` values figure 5's left panel takes a dense spectrum at. The last
#: entry is replaced by the COMPLETE count at runtime.
TRUNCATION_LADDER = (1, 2, 3, 6, 12, 24)

#: The propagation the right panel's two runs share.
TRUNCATION_DT, TRUNCATION_STEPS = 2.0, 2000


def render_truncation(smoke: bool = False) -> dict[str, Any]:
    """Figure 5 -- the one place `n_states` is not `None`, and why.

    `V_dn` and `E_n` are complex (the electronic grid is exterior-complex-
    scaled), so `H_ext`'s anti-Hermitian part is INDEFINITE: a truncated arm
    set can leave eigenvalues in the upper half-plane even though every
    diagonal block on its own is dissipative. The half-Fourier transform
    assumes they all decay, so when one does not, `psi_d` is exponentially
    wrong rather than under-converged, and nothing in the returned object
    says so except the runtime guard.

    The spectra are taken with a DENSE `numpy.linalg.eigvals`, not ARPACK:
    `eigs(which="LI")` under-reports `max Im` by ~6x on these strongly
    non-normal matrices, which would turn the figure's own point into a
    rounding artefact. That is what limits the fixture size.
    """
    plt = _plt()
    deck = f2_small_deck()
    n_avail = int(deck.ing.E_n.shape[1])
    ladder = [n for n in TRUNCATION_LADDER if n < n_avail]
    if smoke:
        ladder = ladder[:2]
    else:
        ladder.append(n_avail)

    max_im: list[float] = []
    for n in ladder:
        t0 = time.time()
        h = extended_hamiltonian(deck.ing, deck.nuc, deck.model, n_states=n)
        vals = np.linalg.eigvals(np.asarray(h.todense()))
        max_im.append(float(vals.imag.max()))
        print(
            f"    n_states={n:3d}  H_ext={h.shape[0]:5d}  max Im(E) = "
            f"{max_im[-1]:+.3e}  ({time.time() - t0:.0f} s)",
            flush=True,
        )

    # The right panel's truncated run is CHOSEN by the left panel, not assumed:
    # whichever truncation the dense spectrum says is worst. The failure is not
    # monotone in `n_states` -- a fixture can lose its growing modes again as
    # arms are added -- so hard-coding a value would risk a panel that shows a
    # truncation which happens to be benign here and calls it the general case.
    worst = int(ladder[int(np.argmax(max_im))])
    dt, n_steps = (TRUNCATION_DT, 50) if smoke else (TRUNCATION_DT, TRUNCATION_STEPS)
    energies = np.array([PACKET_ENERGY])
    runs: dict[str, TdNrmResult] = {}
    truncated_label = f"truncated: n_states = {worst}"
    complete_label = f"complete: n_states = None ({n_avail} arms)"
    for label, n in ((truncated_label, worst), (complete_label, None)):
        launch = initial_packet(
            deck.nuc,
            deck.elec,
            deck.model,
            deck.phi_d,
            deck.ing,
            deck.eps,
            deck.chi,
            _V_INIT,
            energies,
            n_states=n,
        )
        h = extended_hamiltonian(deck.ing, deck.nuc, deck.model, n_states=n)
        with warnings.catch_warnings(record=True) as warned:
            warnings.simplefilter("always")
            runs[label] = propagate_nrm(h, launch, deck.nuc, dt=dt, n_steps=n_steps, order=3)
        print(
            f"    {label}: S(T)/S(0) = {runs[label].unabsorbed[0] / runs[label].survival[0, 0]:.3e}"
        )
        _report(warned)

    fig, (axl, axr) = plt.subplots(1, 2, figsize=(11.0, 4.8))
    # NOT a bare sign test. A mode with `Im E = g` multiplies the packet by
    # `exp(2 g T)` over the propagation, so what matters is whether `g` is big
    # enough to matter over THIS `T` -- and a discretisation as small as this
    # fixture's leaves a residual `+4e-08` even on the complete sum (see
    # `_F2_SMALL_COMPLEX`), which a `> 0` test would paint as divergence. The
    # line is drawn where a mode would grow the packet by 1% over the
    # propagation; the truncations sit four orders above it, the complete sum
    # four orders below.
    t_total = dt * n_steps
    growth_tol = float(np.log(1.01) / (2.0 * t_total))
    im = np.array(max_im)
    hot = im > growth_tol
    x = np.arange(len(ladder))
    axl.axhline(0.0, color="black", linewidth=1.0)
    axl.plot(
        x[hot],
        im[hot],
        "o",
        color="tab:red",
        markersize=9.0,
        label="propagation diverges",
    )
    axl.plot(
        x[~hot],
        im[~hot],
        "o",
        color="tab:green",
        markersize=9.0,
        label="dissipative to solver tolerance",
    )
    axl.axhline(
        growth_tol,
        color="tab:red",
        linestyle=":",
        linewidth=1.2,
        label=rf"1% growth over $T = {t_total:g}$  (${growth_tol:.1e}$)",
    )
    axl.set_yscale("symlog", linthresh=1e-12)
    axl.axhspan(growth_tol, max(1e-2, float(im.max()) * 3.0), color="tab:red", alpha=0.08)
    labels = [str(n) for n in ladder]
    if not smoke:
        labels[-1] = f"all\n({n_avail})"  # the complete sum, the only safe value
    axl.set_xticks(x)
    axl.set_xticklabels(labels)
    axl.set_xlabel("number of projected electronic arms  `n_states`")
    axl.set_ylabel(r"$\max \mathrm{Im}\, E$  of  $H_{\mathrm{ext}}$  (hartree)")
    axl.set_title("dense spectrum of the arrow Hamiltonian")
    axl.set_ylim(-1e-12, max(1e-2, float(im.max()) * 3.0))
    axl.legend(loc="lower left")

    for label, res in runs.items():
        color = "tab:red" if label.startswith("truncated") else "tab:green"
        axr.semilogy(
            res.time,
            res.survival[:, 0] / res.survival[0, 0],
            color=color,
            linestyle="--" if color == "tab:red" else "-",
            label=label,
        )
    axr.set_xlabel("time  $t$  (a.u. of time)")
    axr.set_ylabel(r"$S(t)\,/\,S(0)$")
    axr.set_title(rf"propagated survival, F$_2$ at $E = {PACKET_ENERGY:.3f}$ Ha")
    axr.legend(loc="upper right")

    fig.suptitle(
        "Why the time-dependent route needs the COMPLETE arm set\n"
        r"(F$_2$, small fixture, discrete state B; $H_{\mathrm{ext}}$ is "
        "complex symmetric, not Hermitian)"
    )
    fig.tight_layout()
    _save(
        fig,
        "nrm-td-truncation-diverges",
        n_states=np.asarray(ladder),
        max_imaginary_part=np.asarray(max_im),
        growth_tolerance=np.asarray(growth_tol),
        n_states_available=np.asarray(n_avail),
        time=next(iter(runs.values())).time,
        survival_truncated=runs[truncated_label].survival[:, 0],
        survival_complete=runs[complete_label].survival[:, 0],
        truncated_n_states=np.asarray(worst),
        energy=np.asarray(PACKET_ENERGY),
    )
    plt.close(fig)
    return {"n_states": ladder, "max_im": max_im, "truncated": worst}


# --- figure 6: the gate, as a picture ---------------------------------------

#: `test_nrm_propagation.py::test_propagated_psi_d_reproduces_the_time_independent_solution`
VECTOR_ENERGY = 0.10
VECTOR_DT, VECTOR_STEPS = 1.0, 4000


def render_vector(smoke: bool = False) -> dict[str, Any]:
    """Figure 6 -- `Psi_d^TD(R;E)` against `Psi_d^TI(R;E)`, node by node.

    The shipped gate asserts one number (relative vector error 1.73e-4); this
    is that number drawn. It is a VECTOR comparison on purpose: agreeing on a
    cross section can hide a compensating error in the extraction, agreeing on
    the whole nuclear wavefunction cannot.
    """
    plt = _plt()
    # `_psi_d_for_energy` is the time-independent route's own per-energy solve,
    # the thing the propagation has to reproduce. The gate imports it the same
    # way; there is no public entry point returning `Psi_d` rather than a cross
    # section.
    from qscat.core.nrm.vibrational_excitation import _psi_d_for_energy

    deck = n2_gate_deck()
    dt, n_steps = (VECTOR_DT, 50) if smoke else (VECTOR_DT, VECTOR_STEPS)

    t0 = time.time()
    h = extended_hamiltonian(deck.ing, deck.nuc, deck.model)
    launch = initial_packet(
        deck.nuc,
        deck.elec,
        deck.model,
        deck.phi_d,
        deck.ing,
        deck.eps,
        deck.chi,
        _V_INIT,
        np.array([VECTOR_ENERGY]),
        rank_tol=1e-10,
    )
    with warnings.catch_warnings(record=True) as warned:
        warnings.simplefilter("always")
        res = propagate_nrm(h, launch, deck.nuc, dt=dt, n_steps=n_steps, order=3)
    _report(warned)
    got = res.psi_d[:, 0]
    want = _psi_d_for_energy(
        deck.nuc,
        deck.elec,
        deck.model,
        deck.phi_d,
        deck.eps,
        deck.chi,
        _V_INIT,
        VECTOR_ENERGY,
        deck.ing,
        None,
    )
    rel = float(np.linalg.norm(got - want) / np.linalg.norm(want))
    print(
        f"\n  N2 vector gate: rel = {rel:.4e}, unabsorbed/S(0) = "
        f"{res.unabsorbed[0] / res.survival[0, 0]:.2e}, {time.time() - t0:.0f} s\n",
        flush=True,
    )

    real = deck.nuc.points.imag == 0.0
    r = deck.nuc.points[real].real
    order = np.argsort(r)
    r = r[order]
    ti = np.asarray(want)[real][order]
    td = got[real][order]
    diff = np.abs(td - ti)

    # `Psi_d` is a resonance wavefunction: on this grid it is confined to
    # R ~ 1.7-3.2 bohr and is flat zero over the rest of the 0-8 bohr real
    # region. Drawn full width, the comparison is a spike and the markers pile
    # up in it, so the panels are windowed to where there is a wavefunction to
    # compare -- the `rel` in the title is over the WHOLE vector regardless,
    # and the window is reported in the `.npz`.
    support = np.flatnonzero(np.abs(ti) > 1e-4 * np.abs(ti).max())
    lo = float(r[max(0, support[0] - 2)])
    hi = float(r[min(r.size - 1, support[-1] + 2)])
    window = (r >= lo) & (r <= hi)
    every = max(1, int(window.sum()) // 30)
    mark = np.flatnonzero(window)[::every]

    # `ti`/`td` above are raw DVR COEFFICIENTS, c_j = psi(x_j) * sqrt(w_j), not
    # wavefunction values -- correct for the difference panel below (the gate's
    # rel = 1.7264e-4 is a coefficient-space norm), wrong to plot directly
    # against R. On this nuclear deck the quadrature weight w_j varies 39x
    # across this window (spacing ratio 21.9 grid-wide), so sqrt(w_j) alone
    # manufactures visible cusps at the FEM element bridges even though
    # `Psi_d` itself is smooth (confirmed by a leave-one-out neighbour-
    # interpolation check and by the curvature of the projection below, ~10 --
    # a genuine cusp reads in the hundreds). The fix is to evaluate the DVR
    # INTERPOLANT between nodes rather than rescale node values:
    # `dvr_interpolation_matrix` is exactly that operator and, being built for
    # the sqrt(w)-scaled coefficient vector `qscat.viz` states are stored as,
    # consumes `want`/`got` directly. The TD curve stays at the propagation's
    # own checkpoints (values c_j / sqrt(w_j) at the nodes) rather than being
    # interpolated too, so the overlay shows one smooth reference (TI) against
    # the actual discrete samples (TD) instead of two interpolants that could
    # hide a disagreement between them.
    axis = np.linspace(lo, hi, 400)
    ti_proj = np.asarray(dvr_interpolation_matrix(deck.nuc, axis) @ want)
    w_nodes = deck.nuc.weights[real][order]
    td_value = td / np.sqrt(w_nodes)

    fig, (ax0, ax1) = plt.subplots(
        2, 1, figsize=(7.6, 6.6), sharex=True, gridspec_kw={"height_ratios": [2.0, 1.0]}
    )
    ax0.plot(axis, ti_proj.real, "-", color="tab:blue", label=r"TI  $\mathrm{Re}\,\Psi_d$")
    ax0.plot(axis, ti_proj.imag, "-", color="tab:orange", label=r"TI  $\mathrm{Im}\,\Psi_d$")
    ax0.plot(
        r[mark],
        td_value.real[mark],
        "o",
        color="tab:blue",
        markerfacecolor="none",
        markersize=6.0,
        label=r"TD  $\mathrm{Re}\,\Psi_d$ (nodal values)",
    )
    ax0.plot(
        r[mark],
        td_value.imag[mark],
        "s",
        color="tab:orange",
        markerfacecolor="none",
        markersize=6.0,
        label=r"TD  $\mathrm{Im}\,\Psi_d$ (nodal values)",
    )
    ax0.set_ylabel(r"$\Psi_d(R; E)$  (bohr$^{-1/2}$)")
    ax0.legend(loc="upper right", ncol=2)
    ax0.set_title(rf"$E = {VECTOR_ENERGY:g}$ Ha,  $T = {dt * n_steps:g}$,  $\Delta t = {dt:g}$")

    ax1.semilogy(r, np.maximum(diff, 1e-30), color="tab:red")
    ax1.set_ylabel(r"$|\Psi_d^{\mathrm{TD}} - \Psi_d^{\mathrm{TI}}|$  (coeffs)")
    ax1.set_xlabel(r"internuclear distance  $R$  (bohr)")
    ax1.set_title(
        f"pointwise difference, DVR coefficients -- relative vector error {rel:.3e},\n"
        "the same coefficient-space norm the gate uses (whole vector)"
    )
    ax1.set_xlim(lo, hi)
    top = float(diff[window].max())
    ax1.set_ylim(top * 1e-4, top * 5.0)

    fig.suptitle(
        "N$_2$ vibrational excitation, nonlocal resonance model: the propagated packet's\n"
        "half-Fourier transform IS the time-independent solution "
        r"(discrete state B, $v_i = 0$, complete arm set)"
    )
    fig.tight_layout()
    _save(
        fig,
        "n2-nrm-td-vs-ti-vector",
        R=r,
        psi_ti=ti,
        psi_td=td,
        abs_difference=diff,
        plot_window=np.asarray([lo, hi]),
        relative_error=np.asarray(rel),
        energy=np.asarray(VECTOR_ENERGY),
        dt=np.asarray(dt),
        n_steps=np.asarray(n_steps),
        R_projected=axis,
        psi_ti_projected=ti_proj,
        psi_td_nodal_values=td_value,
    )
    plt.close(fig)
    return {"relative_error": rel}


# --- driver -----------------------------------------------------------------

_FIGURES = {
    "td-vs-ti": render_td_vs_ti,
    "convergence": render_convergence,
    "launch-rank": render_launch_rank,
    "truncation": render_truncation,
    "vector": render_vector,
}


def main(argv: list[str] | None = None) -> None:
    """Render the requested figures (all of them by default).

    `--smoke` shrinks every propagation to a few dozen steps: it exercises
    the whole path -- grids, ingredients, launch, propagation, extraction,
    plotting, `.npz` -- in a couple of minutes, and produces figures that are
    numerically meaningless. Use it to check a change, never to publish.
    """
    args = list(sys.argv[1:] if argv is None else argv)
    smoke = "--smoke" in args
    names = [a for a in args if not a.startswith("-")] or list(_FIGURES)
    unknown = [n for n in names if n not in _FIGURES]
    if unknown:
        raise SystemExit(f"unknown figure(s) {unknown}; choose from {list(_FIGURES)}")
    for name in names:
        print(f"\n=== {name}{' [SMOKE]' if smoke else ''} ===", flush=True)
        t0 = time.time()
        _FIGURES[name](smoke=smoke)
        print(f"=== {name} done in {time.time() - t0:.0f} s ===", flush=True)


if __name__ == "__main__":  # pragma: no cover - a driver, not a test
    main()
