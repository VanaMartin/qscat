"""Differential oracle for the 2026-08-25 kernel consolidation
(docs/superpowers/plans/2026-08-25-kernel-consolidation.md).

Pins the outputs of every code path the consolidation touches -- the six TD
S-matrix/sigma extractor transforms (VE and DA x tw/delta/flow), the TI
`ve/da/dr_cross_section` solvers, the correlation deconvolution factors, and
the four radial special functions -- against a committed golden file,
through PUBLIC entry points only, so the same tests run unchanged before and
after every refactor task.

Two tolerance layers (see the plan's Global Constraints for the measured
justification):

* TRANSFORM layer -- extractor `record`/`sigma` arithmetic driven by
  closed-form synthetic recorded series (no propagation, no SparseLU, no
  threaded BLAS). Compared at `_RTOL`, default 1e-7 (portable across the
  CI architecture); the refactor gate runs it at 1e-12 via
  `QSCAT_KERNEL_ORACLE_RTOL=1e-12`, which holds because nothing on these
  paths is nondeterministic.
* PROPAGATED layer -- end-to-end paths through `make_pade_stepper` /
  `SparseLU`. Compared at `_RTOL_PROPAGATED = max(_RTOL, 1e-6)`: the
  sparse-LU solves route through a multi-threaded BLAS with a measured
  ~4e-9 cross-process drift (see `test_core_td.py`'s module docstring), so
  1e-12 is physically unattainable here and 1e-6 is the established bound.

Regenerating the golden (only legitimate BEFORE the refactor starts, or
when a deliberate, documented physics change lands):

    QSCAT_KERNEL_ORACLE_WRITE=1 uv run --no-sync pytest \
        libs/qscat/tests/test_kernel_consolidation_oracle.py -q

Run WITHOUT `-m "not slow"` so the slow `dr_cross_section` key is written
too, and WITHOUT `-n` (the write-merge fixture is per-process). Write mode
merges into an existing file, so a fast-only rewrite preserves the dr key.

This module stays after the consolidation lands, as a permanent regression
net at the portable default tolerances.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import numpy.typing as npt
import pytest
from qscat.core.correlation import (
    eta_incident,
    eta_outgoing,
    hankel_point_value,
    outgoing_surface_wave,
)
from qscat.core.dissociation import da_cross_section, dr_cross_section
from qscat.core.driven import ve_cross_section
from qscat.core.grids import electronic_grid, nuclear_grid
from qscat.core.td_extractors import Dirac, Flux, TannorWeeks
from qscat.core.time_dependent import (
    PropagationResult,
    sigma_from_correlations,
    td_da_cross_section,
    td_ve_cross_section,
)
from qscat.core.vibrational import vibrational_states
from qscat.dvr import TensorGrid
from qscat.model import H2P, N2
from qscat.special import (
    riccati_bessel_en,
    riccati_bessel_en_mass,
    riccati_hankel_en,
    riccati_hankel_en_mass,
)

GOLDEN_PATH = Path(__file__).with_name("kernel_consolidation_golden.npz")
_WRITE = os.environ.get("QSCAT_KERNEL_ORACLE_WRITE") == "1"
_RTOL = float(os.environ.get("QSCAT_KERNEL_ORACLE_RTOL", "1e-7"))
_RTOL_PROPAGATED = max(_RTOL, 1e-6)

_STAGED: dict[str, npt.NDArray[np.complex128] | npt.NDArray[np.float64]] = {}
_GOLDEN_CACHE: dict[str, np.ndarray] | None = None


def _golden() -> dict[str, np.ndarray]:
    global _GOLDEN_CACHE
    if _GOLDEN_CACHE is None:
        if not GOLDEN_PATH.exists():
            pytest.fail(
                f"{GOLDEN_PATH.name} missing -- regenerate with "
                "QSCAT_KERNEL_ORACLE_WRITE=1 (module docstring)"
            )
        with np.load(GOLDEN_PATH) as z:
            _GOLDEN_CACHE = {name: z[name] for name in z.files}
    return _GOLDEN_CACHE


def _check(key: str, value: npt.ArrayLike, *, rtol: float) -> None:
    arr = np.asarray(value)
    if _WRITE:
        _STAGED[key] = arr
        return
    golden = _golden()
    assert key in golden, f"{key!r} missing from {GOLDEN_PATH.name} -- regenerate"
    np.testing.assert_allclose(arr, golden[key], rtol=rtol, atol=0.0)


@pytest.fixture(scope="module", autouse=True)
def _write_golden_after_module() -> object:
    yield None
    if _WRITE and _STAGED:
        existing: dict[str, np.ndarray] = {}
        if GOLDEN_PATH.exists():
            with np.load(GOLDEN_PATH) as z:
                existing = {name: z[name] for name in z.files}
        existing.update(_STAGED)
        np.savez(GOLDEN_PATH, **existing)


# --- The tiny/fast N2 deck: copied from test_core_td.py / test_td_extractors.py
# (deliberately unconverged -- a regression pin, not physics).

TG = TensorGrid(
    [
        electronic_grid(r_max=12.0, order=5, n_complex=3),
        nuclear_grid(quadrature=6, r_max=14.0, n_complex=3),
    ]
)
EPS, CHI = vibrational_states(TG.grids[1], N2.mu, 4, N2.v0)

V_INIT = 0
VPRIMES = [0, 1]  # includes the elastic (v'=v_init) channel
WP_IN = {"r0": 4.0, "p0": -0.5, "sigma": 1.2}
WP_OUT = {"r0_out": 6.0, "p0_out": 0.5, "sigma_out": 1.0}
NUCLEAR_WP_OUT = {"r0_out": 7.0, "p0_out": 5.0, "sigma_out": 1.0}
DT = 0.2
N_STEPS = 5
POSITION = 37  # electronic DVR index in the real region (test_td_extractors.py)
NUCLEAR_SURFACE = 90  # nuclear DVR index, R=7.12 bohr (test_td_extractors.py)

# E grids chosen to pin every branch: E<=0 early return, a closed v'=1
# channel (0.001), open channels (0.10/0.15); DA: closed (0.10, N2's DA
# threshold sits above 0.3 Ha) and open (0.60).
E_VE = np.array([-0.05, 0.001, 0.10, 0.15])
E_DA = np.array([0.10, 0.60])

N_T = 7  # odd -> Simpson weights; the n_steps=6 truncation pins the trapezoid fallback


def _synthetic_states(n_t: int, seed_phase: float) -> list[npt.NDArray[np.complex128]]:
    """Deterministic, closed-form fake trajectory (size TG.size) -- elementwise
    numpy only, so the transform layer built on it is reproducible to the
    last bit within one machine (no RNG, no BLAS reductions, no LU)."""
    j = np.arange(TG.size, dtype=np.float64)
    envelope = np.exp(-(((j / TG.size) - 0.35) ** 2) / 0.02)
    states = []
    for n in range(n_t):
        re = envelope * np.cos(0.037 * j + seed_phase + 0.61 * n)
        im = 0.8 * envelope * np.sin(0.023 * j + 1.3 * seed_phase + 0.4 * n)
        states.append((re + 1j * im).astype(np.complex128))
    return states


def _recorded(ext, states):  # type: ignore[no-untyped-def]
    for psi in states:
        ext.record(psi)
    return ext


# --- Transform layer: the six extractor paths on synthetic series ------------


def test_ve_tw_transform_golden() -> None:
    tw = _recorded(
        TannorWeeks(TG, N2, EPS, CHI, V_INIT, VPRIMES, WP_OUT, wp_in=WP_IN, dt=DT),
        _synthetic_states(N_T, 0.0),
    )
    tw_free = _recorded(
        TannorWeeks(TG, N2, EPS, CHI, V_INIT, VPRIMES, WP_OUT, wp_in=WP_IN, dt=DT),
        _synthetic_states(N_T, 0.9),
    )
    _check("ve_tw_sigma_nofree", tw.sigma(E_VE), rtol=_RTOL)
    _check("ve_tw_sigma_free", tw.sigma(E_VE, free=tw_free), rtol=_RTOL)
    _check("ve_tw_sigma_trunc6", tw.sigma(E_VE, free=tw_free, n_steps=6), rtol=_RTOL)


def test_ve_delta_transform_golden() -> None:
    di = _recorded(
        Dirac(TG, N2, EPS, CHI, V_INIT, VPRIMES, POSITION, wp_in=WP_IN, dt=DT),
        _synthetic_states(N_T, 0.2),
    )
    di_free = _recorded(
        Dirac(TG, N2, EPS, CHI, V_INIT, VPRIMES, POSITION, wp_in=WP_IN, dt=DT),
        _synthetic_states(N_T, 1.1),
    )
    _check("ve_delta_sigma_nofree", di.sigma(E_VE), rtol=_RTOL)
    _check("ve_delta_sigma_free", di.sigma(E_VE, free=di_free), rtol=_RTOL)


def test_ve_flow_transform_golden() -> None:
    fl = _recorded(
        Flux(TG, N2, EPS, CHI, V_INIT, VPRIMES, POSITION, wp_in=WP_IN, dt=DT),
        _synthetic_states(N_T, 0.4),
    )
    fl_free = _recorded(
        Flux(TG, N2, EPS, CHI, V_INIT, VPRIMES, POSITION, wp_in=WP_IN, dt=DT),
        _synthetic_states(N_T, 1.3),
    )
    _check("ve_flow_sigma_nofree", fl.sigma(E_VE), rtol=_RTOL)
    _check("ve_flow_sigma_free", fl.sigma(E_VE, free=fl_free), rtol=_RTOL)


def test_da_tw_transform_golden() -> None:
    tw = _recorded(
        TannorWeeks(
            TG,
            N2,
            EPS,
            CHI,
            V_INIT,
            [],
            NUCLEAR_WP_OUT,
            wp_in=WP_IN,
            dt=DT,
            axis="nuclear",
            n_channels=1,
        ),
        _synthetic_states(N_T, 0.6),
    )
    _check("da_tw_sigma", tw.sigma(E_DA), rtol=_RTOL)


def test_da_delta_transform_golden() -> None:
    di = _recorded(
        Dirac(
            TG,
            N2,
            EPS,
            CHI,
            V_INIT,
            [],
            NUCLEAR_SURFACE,
            wp_in=WP_IN,
            dt=DT,
            axis="nuclear",
            n_channels=1,
        ),
        _synthetic_states(N_T, 0.7),
    )
    _check("da_delta_sigma", di.sigma(E_DA), rtol=_RTOL)


def test_da_flow_transform_golden() -> None:
    fl = _recorded(
        Flux(
            TG,
            N2,
            EPS,
            CHI,
            V_INIT,
            [],
            NUCLEAR_SURFACE,
            wp_in=WP_IN,
            dt=DT,
            axis="nuclear",
            n_channels=1,
        ),
        _synthetic_states(N_T, 0.8),
    )
    _check("da_flow_sigma", fl.sigma(E_DA), rtol=_RTOL)


def test_sigma_from_correlations_golden() -> None:
    """The public batch transform on a hand-built PropagationResult -- pins
    the TW-VE kernel without going through an extractor instance at all."""
    t = np.arange(N_T, dtype=np.float64) * DT
    n = np.arange(N_T, dtype=np.float64)[:, None]
    ch = np.arange(len(VPRIMES), dtype=np.float64)[None, :]
    c = (np.cos(0.5 * n + ch) + 1j * np.sin(0.3 * n - 0.2 * ch)).astype(np.complex128)
    result = PropagationResult(t=t, c=c, norm=np.zeros(N_T), snapshots=[])
    free = PropagationResult(t=t, c=(0.7 * c + (0.1 + 0.05j)), norm=np.zeros(N_T), snapshots=[])
    sigma = sigma_from_correlations(
        TG,
        N2,
        result,
        EPS,
        V_INIT,
        VPRIMES,
        E_VE,
        dt=DT,
        wp_in=WP_IN,
        wp_out=WP_OUT,
        free_result=free,
    )
    _check("sigma_from_correlations", sigma, rtol=_RTOL)


# --- Transform layer: the correlation deconvolution factors (m5/m6) ----------


def test_correlation_factor_golden() -> None:
    g_e, g_n = TG.grids[0], TG.grids[1]
    vals = np.array(
        [
            eta_incident(g_e, 0.45, N2.ell, **WP_IN),
            eta_outgoing(g_e, 0.45, N2.ell, **WP_OUT),
            eta_outgoing(g_n, 1.7, 0, mass=N2.mu, **NUCLEAR_WP_OUT),
            hankel_point_value(g_e, 7.5, 0.45, N2.ell, 0),
            hankel_point_value(g_n, 7.12, 1.7, 0, 0, mass=N2.mu),
            *outgoing_surface_wave(g_e, 7.5, 0.45, N2.ell, 0),
            *outgoing_surface_wave(g_n, 7.12, 1.7, 0, 0, mass=N2.mu),
            *outgoing_surface_wave(g_e, 7.5, 0.45, 0, -1, mass=1.0),  # charged (mpmath)
            hankel_point_value(g_e, 7.5, 0.45, 0, -1, mass=1.0),  # charged (mpmath)
        ],
        dtype=np.complex128,
    )
    _check("correlation_factors", vals, rtol=_RTOL)


# --- Transform layer: the radial special functions (M17) ----------------------


def test_radial_golden_and_mu1_identity() -> None:
    r = np.linspace(0.1, 12.0, 40)
    blocks = []
    for k, l in ((0.7, 0), (1.3, 2)):
        blocks.append(riccati_bessel_en(r, k, l).astype(np.complex128))
        blocks.append(riccati_hankel_en(r, k, l))
        blocks.append(riccati_bessel_en_mass(r, k, l, 2.5).astype(np.complex128))
        blocks.append(riccati_hankel_en_mass(r, k, l, 2.5))
        # The mu=1 reduction is BYTE-identical (2.0*1.0 == 2.0 exactly) --
        # the invariant M17's merge must preserve.
        assert np.array_equal(riccati_bessel_en_mass(r, k, l, 1.0), riccati_bessel_en(r, k, l))
        assert np.array_equal(riccati_hankel_en_mass(r, k, l, 1.0), riccati_hankel_en(r, k, l))
    _check("radial_functions", np.stack(blocks), rtol=_RTOL)


# --- Propagated/LU layer ------------------------------------------------------


def test_ti_ve_golden() -> None:
    sigma = ve_cross_section(TG, N2, EPS, CHI, V_INIT, VPRIMES, [0.10, 0.15])
    _check("ti_ve_sigma", sigma, rtol=_RTOL_PROPAGATED)


def test_ti_da_golden() -> None:
    sigma = da_cross_section(TG, N2, EPS, CHI, V_INIT, E_DA, n_channels=1)
    _check("ti_da_sigma", sigma, rtol=_RTOL_PROPAGATED)


@pytest.mark.parametrize("method", ["tw", "delta", "flow"])
def test_td_ve_end_to_end_golden(method: str) -> None:
    sigma = td_ve_cross_section(
        TG,
        N2,
        EPS,
        CHI,
        V_INIT,
        VPRIMES,
        [0.10, 0.15],
        dt=DT,
        n_steps=N_STEPS,
        wp_in=WP_IN,
        wp_out=WP_OUT,
        method=method,
        position=POSITION,
        surface=POSITION,
    )
    _check(f"td_ve_{method}", sigma, rtol=_RTOL_PROPAGATED)


@pytest.mark.parametrize("method", ["flow", "delta", "tw"])
def test_td_da_end_to_end_golden(method: str) -> None:
    sigma = td_da_cross_section(
        TG,
        N2,
        EPS,
        CHI,
        V_INIT,
        E_DA,
        dt=DT,
        n_steps=N_STEPS,
        wp_in=WP_IN,
        method=method,
        surface=NUCLEAR_SURFACE,
        position=NUCLEAR_SURFACE,
        wp_out=NUCLEAR_WP_OUT,
        n_channels=1,
    )
    _check(f"td_da_{method}", sigma, rtol=_RTOL_PROPAGATED)


@pytest.mark.slow
def test_dr_golden() -> None:
    """dr_cross_section on the laptop-fast H2+ proxy deck (the same shape
    test_dissociation.py's slow dr tests use) -- gates the lib-M2 rewrite.
    Slow tier: the mpmath Coulomb channel functions dominate the runtime."""
    tg = TensorGrid(
        [
            electronic_grid(r_max=60.0, order=8, n_complex=6),
            nuclear_grid(r_max=22.0, n_complex=6, quadrature=10),
        ]
    )
    eps, chi = vibrational_states(tg.grids[1], H2P.mu, 3, H2P.v0)
    sigma, amp = dr_cross_section(
        tg, H2P, eps, chi, 0, np.array([0.01, 0.03]), n_channels=2, return_amplitude=True
    )
    _check("dr_sigma", sigma, rtol=_RTOL_PROPAGATED)
    _check("dr_amp", amp, rtol=_RTOL_PROPAGATED)
