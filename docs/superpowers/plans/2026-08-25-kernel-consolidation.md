# Kernel Consolidation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Collapse the six copy-pasted TD S-matrix transform kernels (and their five sibling duplications: the three-way `td_ve_cross_section` method blocks, `dr_cross_section`'s inlined driven sweep, the `correlation.py` coefficient/value helpers, and the four half-duplicated radial functions) into single parameterized kernels, with a committed differential oracle proving every output unchanged.

**Architecture:** One shared S-vector skeleton (`s_vector_transform`) plus one shared sigma-with-elastic-reference kernel (`sigma_from_s`) live in `qscat.core.time_dependent`; the per-method differences (Tannor-Weeks eta, Dirac point value, Flux Wronskian pair — each on either the electronic or nuclear exit axis) become small strategy closures built by two factories (`correlation_channel_s`, `flux_channel_s`). Everything else stays where it is: the six existing per-energy sigma helpers become thin assemblies over the shared kernels, keeping their names and signatures because `projects/n2_2d_td_cross_section/td_cross_section.py` and `libs/qscat/tests/test_core_td.py` import two of them.

**Tech Stack:** Python 3.12, numpy/scipy, pytest (`uv run --no-sync pytest`), mypy --strict, ruff. No new dependencies, no Rust changes, no Docker changes.

**Spec:** the "Findings addressed" section below (self-contained; from the 2026-08-25 release review).

## Global Constraints

PyPI release is DEFERRED until the peer-reviewed article publishes — repo-only distribution; no release/publishing tasks. After every task: `uv run --no-sync pytest -m "not slow" -n auto --dist loadfile` green, `uv run --no-sync mypy libs/qscat/qscat apps/qscat-run/qscat_run` clean, `uv run --no-sync ruff check .` and `ruff format --check .` clean. The consolidation must be numerically IDENTITY-PRESERVING: the differential oracle (Task 1) pins current outputs at rtol=1e-12 and every later task re-runs it. Solver changes need a `validate:core` + `validate:n2` labelled CI run (or local `pytest -m slow` on the touched suites) before merge. Never `git commit -a` — stage explicit paths.

**How the rtol=1e-12 pin is realized (read before Task 1).** A committed golden compared through a *fresh propagation* cannot hold 1e-12: the sparse-LU solves inside `make_pade_stepper` route through a multi-threaded BLAS whose summation order varies run-to-run — `libs/qscat/tests/test_core_td.py`'s module docstring records a *measured ~4e-9 relative cross-process drift* for exactly these paths, and its own regression pin is rtol=1e-6 for that reason. So the oracle splits into two layers:

- **Transform layer (the six kernels being consolidated), pinned at rtol=1e-12**: the extractors' `record`/`sigma` arithmetic is exercised on *deterministic, closed-form synthetic recorded series* — no propagation, no `SparseLU`, no threaded BLAS — so the golden comparison is stable at 1e-12 on the machine that generated it. Every task runs this layer with `QSCAT_KERNEL_ORACLE_RTOL=1e-12`.
- **Propagated/LU layer (end-to-end `td_*`, `ve/da/dr_cross_section`), pinned at rtol=1e-6**: same bound `test_core_td.py` already uses; it gates wiring, not last-digit arithmetic. Physics-level regressions are orders of magnitude larger than 1e-6.

The committed default tolerance is `1e-7` (transform) / `1e-6` (propagated) so the test also passes on the cross-architecture CI runner — the same portability bound the existing golden tests demonstrably clear (`test_td_extractors.py` pins ~1e-5-magnitude sigmas at atol=1e-12 ⇒ ~1e-7 relative, and CI is green). The env var `QSCAT_KERNEL_ORACLE_RTOL` tightens the transform layer to 1e-12 for the refactor gate on the developer machine.

## Findings addressed

From the 2026-08-25 release review. Each was re-verified against the current tree on 2026-08-25; verification notes in brackets.

- **lib-C4**: six near-identical S-matrix kernels — `time_dependent.py` `_s_vector_one_energy`, `td_extractors.py` `_tw_da_s_vector_one_energy`, `_dirac_s_vector_one_energy`, `_dirac_da_s_vector_one_energy`, `_flux_s_vector_one_energy`, `_flux_da_s_vector_one_energy` — sharing the skeleton (zeros-init → E<=0 early return → `_quadrature_weights` → `e_tot = E + eps[v_init]` → `k=sqrt(2E)` → `eta_in` → `phase = exp(1j*e_tot*t)` → per-channel loop with closed-channel continue → `s_raw = sum(weights*phase*…)*dt`), differing only in the outgoing factor (eta_out vs hankel_point_value vs Wronskian flux vs nuclear variants). Six matching σ kernels with the verbatim elastic-reference block (`ref = s_free[j] if ... else 1.0` / `sigma = pi*abs(s-ref)**2/(2E)`) repeated at three sites.
  [VERIFIED. The six S kernels are `time_dependent.py:312` and `td_extractors.py:188,456,532,827,915`. The quadrature helper's actual name is `quadrature_weights` (public, `time_dependent.py:226`), not `_quadrature_weights`. The elastic-reference block appears verbatim at `time_dependent.py:396-403`, `td_extractors.py:520-528` (`_dirac_sigma_one_energy`), `td_extractors.py:903-911` (`_flux_sigma_one_energy`); the three DA sigma kernels (`_tw_da_sigma_one_energy`, `_dirac_da_sigma_one_energy`, `_flux_da_sigma_one_energy`) repeat the same skeleton without the reference and with `_C_DA = np.pi` — numerically the *same* constant as the VE `np.pi`, so one kernel serves all six. Two extra facts that constrain the refactor: (1) `projects/n2_2d_td_cross_section/td_cross_section.py:20-21` imports `_s_vector_one_energy` and `_sigma_one_energy` from `qscat.core.time_dependent` by name, and `libs/qscat/tests/test_core_td.py:38` imports `_sigma_one_energy` — both names and positional signatures must survive; (2) nothing outside `td_extractors.py` imports any of its six private kernels, so five of them can be deleted outright.]
- **lib-M1**: `td_ve_cross_section` is three copy-pasted ~32-line method blocks (construct extractor → propagate → optional free-reference propagate → return sigma) where `td_da_cross_section` does the same job with one if/elif + one shared propagate.
  [VERIFIED: `time_dependent.py:544-640` vs `td_da_cross_section`'s single-dispatch shape at `time_dependent.py:791-855`.]
- **lib-M2**: `dr_cross_section` re-inlines `driven.ve_cross_section`'s lazy-LU+refactor sweep and `da_cross_section`'s exit-channel projection loop because the helper "cannot pass `charge` through to `channel_vector`" — but `channel_vector` already accepts `charge`. Forward the keyword and delete the inlined copies.
  [VERIFIED with one correction: the genuine duplication is the lazy-LU+refactor sweep (`dissociation.py:422-444` duplicates `driven.py:188-223`), and the docstring claim at `dissociation.py:383-387` is indeed stale — `channel_vector` (`channels.py:35-42`) has taken `charge: int = 0` since the H2+ work. The *exit-channel loop* (`dissociation.py:447-457`) is NOT a copy of `da_cross_section`'s: dr projects the volume T-matrix `<phi_ryd F_K | V_DR | Psi+>` while da reads a boundary flux value — only the closed-channel skeleton is shared, and that skeleton is 4 lines. The refactor therefore reuses `ve_cross_section(..., return_wavefunction=True)` for the Psi+ sweep (exactly as `da_cross_section` already does at `dissociation.py:255-265`) and keeps dr's own exit projection. `ResonanceModel` carries `charge` (`model/diatomic.py:70-74`), so `ve_cross_section` can forward `model.charge`; for every existing caller (`charge == 0`) the `channel_vector` branch is unchanged ⇒ bit-identical.]
- **lib-M17**: `riccati_bessel_en`/`riccati_bessel_en_mass` and `riccati_hankel_en`/`riccati_hankel_en_mass` are byte-identical except `2.0*k` vs `2.0*mu*k`; collapse to two functions with `mu: float = 1.0` (keep old names as thin deprecated aliases so nothing breaks).
  [VERIFIED: `special/radial.py:44-121`. `2.0*1.0 == 2.0` exactly, so the merged function at `mu=1.0` is bit-for-bit the old base function — `radial.py:99`'s own docstring already asserts this. The `_mass` names are imported by `correlation.py`, `dissociation.py`, `tuning/probes.py`, several test files, and re-exported by `special/__init__.py` — they stay as thin aliases; no call site changes.]
- **lib-m5**: `correlation.py` `_regular_coeffs`/`_outgoing_coeffs` differ only in the `f_vals` line.
  [VERIFIED: `correlation.py:122-155`; the shared 3 lines are `*sqrt(w)` → `.astype(complex128)` → zero past `R0`. Both names stay (a `test_correlation.py` docstring references `_outgoing_coeffs`); only the shared body moves into a helper.]
- **lib-m6**: `correlation.py` (~:279) re-derives the Riccati-Hankel VALUE inline from spherical_jn/yn where `riccati_hankel_en_mass` (imported in the same file) already provides it; only the derivative branch has a scipy reason.
  [VERIFIED: `outgoing_surface_wave`'s neutral branch, `correlation.py:276-283`. `phi = pref*r*h_l/2` is exactly `hankel_point_value(grid, r, k, l, 0, mass=mass)` (same multiplication order ⇒ bit-identical); `dphi` legitimately still needs `spherical_jn/yn(derivative=True)` because `qscat.special.radial` exposes no derivative primitive. The charged branch's finite-difference closure stays as-is (mpmath, out of scope for m6).]

**Non-goals (deliberate):** no `_map_energies` helper for the repeated scalar-vs-array `E` squeeze blocks (not in the findings; pure restructuring risk for zero reviewer value); no switching of internal `riccati_*_en_mass` call sites to the merged signature (the finding itself says keep the aliases so nothing breaks); no change to `outgoing_surface_wave`'s charged branch; no public API renames anywhere, so `docs/physics/td-extractors.md` and `td-da.md` need no edits (Task 8 verifies this by grep rather than asserting it).

---

## Task 1 — The differential oracle

**Files:**
- Create `libs/qscat/tests/test_kernel_consolidation_oracle.py`
- Generate + commit `libs/qscat/tests/kernel_consolidation_golden.npz` (note: `*.npz` is gitignored at `.gitignore:23` — stage it with `git add -f`)

**Steps:**

- [ ] Write `libs/qscat/tests/test_kernel_consolidation_oracle.py` with exactly this content:

```python
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
            TG, N2, EPS, CHI, V_INIT, [], NUCLEAR_WP_OUT,
            wp_in=WP_IN, dt=DT, axis="nuclear", n_channels=1,
        ),
        _synthetic_states(N_T, 0.6),
    )
    _check("da_tw_sigma", tw.sigma(E_DA), rtol=_RTOL)


def test_da_delta_transform_golden() -> None:
    di = _recorded(
        Dirac(
            TG, N2, EPS, CHI, V_INIT, [], NUCLEAR_SURFACE,
            wp_in=WP_IN, dt=DT, axis="nuclear", n_channels=1,
        ),
        _synthetic_states(N_T, 0.7),
    )
    _check("da_delta_sigma", di.sigma(E_DA), rtol=_RTOL)


def test_da_flow_transform_golden() -> None:
    fl = _recorded(
        Flux(
            TG, N2, EPS, CHI, V_INIT, [], NUCLEAR_SURFACE,
            wp_in=WP_IN, dt=DT, axis="nuclear", n_channels=1,
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
    free = PropagationResult(
        t=t, c=(0.7 * c + (0.1 + 0.05j)), norm=np.zeros(N_T), snapshots=[]
    )
    sigma = sigma_from_correlations(
        TG, N2, result, EPS, V_INIT, VPRIMES, E_VE,
        dt=DT, wp_in=WP_IN, wp_out=WP_OUT, free_result=free,
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
        TG, N2, EPS, CHI, V_INIT, VPRIMES, [0.10, 0.15],
        dt=DT, n_steps=N_STEPS, wp_in=WP_IN, wp_out=WP_OUT,
        method=method, position=POSITION, surface=POSITION,
    )
    _check(f"td_ve_{method}", sigma, rtol=_RTOL_PROPAGATED)


@pytest.mark.parametrize("method", ["flow", "delta", "tw"])
def test_td_da_end_to_end_golden(method: str) -> None:
    sigma = td_da_cross_section(
        TG, N2, EPS, CHI, V_INIT, E_DA,
        dt=DT, n_steps=N_STEPS, wp_in=WP_IN,
        method=method, surface=NUCLEAR_SURFACE, position=NUCLEAR_SURFACE,
        wp_out=NUCLEAR_WP_OUT, n_channels=1,
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
```

- [ ] Normalize the new file's layout before anything else (the plan's snippets are semantically exact, but the formatter has final say on line wrapping — this applies to every code snippet in this plan):
  ```bash
  uv run --no-sync ruff check --fix libs/qscat/tests/test_kernel_consolidation_oracle.py
  uv run --no-sync ruff format libs/qscat/tests/test_kernel_consolidation_oracle.py
  ```
- [ ] Generate the golden file (this MUST happen at the current, pre-refactor tree state — do not start Task 2 first):
  ```bash
  QSCAT_KERNEL_ORACLE_WRITE=1 uv run --no-sync pytest libs/qscat/tests/test_kernel_consolidation_oracle.py -q
  ```
  (No `-m` filter, no `-n` — the slow dr test must run so its key is written; expect a few minutes for the H2+ proxy.)
- [ ] Verify the pins hold from a *fresh process* at the default tolerances:
  ```bash
  uv run --no-sync pytest libs/qscat/tests/test_kernel_consolidation_oracle.py -q
  ```
- [ ] Verify the transform layer holds at the refactor gate tolerance, twice (two separate processes — this is the determinism check the whole plan leans on):
  ```bash
  QSCAT_KERNEL_ORACLE_RTOL=1e-12 uv run --no-sync pytest libs/qscat/tests/test_kernel_consolidation_oracle.py -q -m "not slow"
  QSCAT_KERNEL_ORACLE_RTOL=1e-12 uv run --no-sync pytest libs/qscat/tests/test_kernel_consolidation_oracle.py -q -m "not slow"
  ```
  If any transform-layer test fails here *before any refactoring*, the failing quantity is not process-stable at 1e-12 on this machine: measure the actual max relative deviation from the failure output, set the gate to one decade above it (e.g. `1e-11`), record the measured floor in the test's module docstring, and use that value as the gate in every later task. Do NOT weaken `_RTOL`'s committed default.
- [ ] Fast tier, types, lint (the standard gate — repeated verbatim so no task can skip it):
  ```bash
  uv run --no-sync pytest -m "not slow" -n auto --dist loadfile
  uv run --no-sync mypy libs/qscat/qscat apps/qscat-run/qscat_run
  uv run --no-sync ruff check . && uv run --no-sync ruff format --check .
  ```
- [ ] Commit (stage explicitly; the golden needs `-f` because `.gitignore` has `*.npz`):
  ```bash
  git add libs/qscat/tests/test_kernel_consolidation_oracle.py
  git add -f libs/qscat/tests/kernel_consolidation_golden.npz
  git commit -m "test(core): pin a differential oracle for the kernel consolidation"
  ```

---

## Task 2 — lib-C4 part 1: the shared kernels, wired into `time_dependent.py`

**File:** `libs/qscat/qscat/core/time_dependent.py`

**Steps:**

- [ ] Add `FemDvrEcsGrid` to the `qscat.dvr` import (line 71 becomes `from qscat.dvr import FemDvrEcsGrid, TensorGrid`) and add `from collections.abc import Callable` to the stdlib imports.
- [ ] Insert the shared kernel block immediately AFTER `quadrature_weights` (after line 246) and BEFORE `free_hamiltonian`. These four names are deliberately NOT added to `__all__`: they are the internal shared kernel `td_extractors.py` imports, not part of the advertised API.

```python
# --- Shared S-matrix / sigma transform kernels (kernel consolidation,
# 2026-08-25). Every TD energy-extraction route -- Tannor-Weeks, Dirac
# (delta), Flux (flow), on either the electronic (VE) or nuclear (DA) exit
# axis -- runs the same skeleton: zeros-init, E<=0 early return, Simpson/
# trapezoid quadrature weights, `e_tot = E + eps[v_init]`, the incident
# deconvolution `eta_in` (ALWAYS on the electronic incident axis, even for
# the nuclear DA extractors), `phase = exp(i*e_tot*t)`, then a per-channel
# loop that skips closed channels and forms one S-matrix element from the
# recorded series. Only the per-channel element differs, in two shapes
# (`correlation_channel_s` for TW/Dirac, `flux_channel_s` for Flux); the
# exit-axis mass enters ONLY through `kp = sqrt(2*exit_mass*excess)` and
# the strategy's own outgoing factor. `sigma_from_s` is the matching
# shared |S - ref|^2 step: the elastic (v'==v_init) reference applies on
# the electronic VE axis only (`elastic` mask); DA passes `elastic=None`
# and its prefactor is the SAME `pi` (the historical `_C_DA` -- see
# td_extractors.py's module docstring for the `S = 1 - 2*pi*i*T`
# reconciliation with the TI oracle's `4*pi^3`).
#
# Not in `__all__`: consumed by this module and `td_extractors.py`, kept
# out of the advertised public surface.

ChannelS = Callable[
    [int, float, npt.NDArray[np.float64], npt.NDArray[np.complex128], complex],
    complex,
]
"""Per-channel strategy for `s_vector_transform`, called as
`channel_s(j, kp, weights, phase, eta_in)` for each OPEN exit channel `j`
with its outgoing momentum `kp`; returns the S-matrix element `S_j`."""


def s_vector_transform(
    g_in: FemDvrEcsGrid,
    l_in: int,
    wp_in: _WpIn,
    t: npt.NDArray[np.float64],
    eps_init: float,
    thresholds: npt.NDArray[np.float64],
    E: float,
    exit_mass: float,
    channel_s: ChannelS,
) -> npt.NDArray[np.complex128]:
    """The shared raw-S-matrix skeleton (block comment above): one element
    per exit channel, `0` for closed channels (`E + eps_init -
    thresholds[j] <= 0`) and for `E <= 0`. `g_in`/`l_in`/`wp_in` are the
    incident electron's grid, partial wave and Gaussian parameters
    (`eta_incident` is always electronic); `thresholds` are the exit-channel
    threshold energies (`eps[vprimes]` for VE, `eps_e` for DA);
    `exit_mass` is the exit-axis reduced mass (1.0 electronic, `model.mu`
    nuclear), entering only `kp = sqrt(2*exit_mass*(e_tot - threshold))`."""
    S = np.zeros(len(thresholds), dtype=np.complex128)
    if E <= 0.0:
        return S
    weights = quadrature_weights(t.size)
    e_tot = E + eps_init
    k = float(np.sqrt(2.0 * E))
    eta_in = eta_incident(g_in, k, l_in, **wp_in)
    phase = np.exp(1j * e_tot * t)
    for j in range(len(thresholds)):
        excess = e_tot - thresholds[j]
        if excess <= 0.0:
            continue  # closed exit channel
        kp = float(np.sqrt(2.0 * exit_mass * excess))
        S[j] = channel_s(j, kp, weights, phase, eta_in)
    return S


def correlation_channel_s(
    outgoing_factor: Callable[[float], complex],
    c: npt.NDArray[np.complex128],
    dt: float,
) -> ChannelS:
    """The TW/Dirac-shaped per-channel element: `S_j = (sum_n w_n
    e^{i*e_tot*t_n} c_j(t_n)) * dt / (2*pi * conj(F_out(kp)) * eta_in)`.
    `outgoing_factor(kp)` is the method's outgoing deconvolution scalar --
    `eta_outgoing(...)` for Tannor-Weeks, `hankel_point_value(...)` for
    Dirac, on whichever axis/mass the caller closed over."""

    def s_element(
        j: int,
        kp: float,
        weights: npt.NDArray[np.float64],
        phase: npt.NDArray[np.complex128],
        eta_in: complex,
    ) -> complex:
        s_raw = np.sum(weights * phase * c[:, j]) * dt
        return complex(s_raw / (2.0 * np.pi * np.conj(outgoing_factor(kp)) * eta_in))

    return s_element


def flux_channel_s(
    outgoing_pair: Callable[[float], tuple[complex, complex]],
    b: npt.NDArray[np.complex128],
    d: npt.NDArray[np.complex128],
    dt: float,
    exit_mass: float,
) -> ChannelS:
    """The Flux-shaped per-channel element (the Wronskian transform --
    td_extractors.py's module docstring): `S_j = -i/(2*exit_mass*eta_in) *
    (sum_n w_n (conj(phi_out)*d_j(t_n) - b_j(t_n)*conj(dphi_out))
    e^{i*e_tot*t_n}) * dt`. `outgoing_pair(kp)` is `outgoing_surface_wave`'s
    `(phi_out, dphi_out)` on whichever axis/mass the caller closed over."""

    def s_element(
        j: int,
        kp: float,
        weights: npt.NDArray[np.float64],
        phase: npt.NDArray[np.complex128],
        eta_in: complex,
    ) -> complex:
        phi_out, dphi_out = outgoing_pair(kp)
        wronskian = np.conj(phi_out) * d[:, j] - b[:, j] * np.conj(dphi_out)
        s_raw = np.sum(weights * wronskian * phase) * dt
        return complex((-1j / (2.0 * exit_mass * eta_in)) * s_raw)

    return s_element


def sigma_from_s(
    s_full: npt.NDArray[np.complex128],
    s_free: npt.NDArray[np.complex128] | None,
    thresholds: npt.NDArray[np.float64],
    eps_init: float,
    E: float,
    elastic: npt.NDArray[np.bool_] | None,
) -> npt.NDArray[np.float64]:
    """The shared `sigma = pi*|S - ref|^2/(2E)` step, zeros for `E <= 0`
    and closed channels. `elastic` marks the diagonal (v'==v_init) VE
    channel(s), whose `ref` is `s_free[j]` when a free-reference S-vector
    is supplied, else the literal 1 (see `_sigma_one_energy`'s docstring
    for why callers should supply `s_free`); every other channel -- and
    every DA channel (`elastic=None`, DA has no diagonal) -- uses `ref=0`,
    where `pi*|S|^2/(2E)` is the historical `_C_DA = pi` convention."""
    sigma = np.zeros(len(thresholds), dtype=np.float64)
    if E <= 0.0:
        return sigma
    e_tot = E + eps_init
    for j in range(len(thresholds)):
        if e_tot - thresholds[j] <= 0.0:
            continue  # closed exit channel
        if elastic is not None and elastic[j]:
            ref = complex(s_free[j]) if s_free is not None else 1.0 + 0.0j
        else:
            ref = 0.0 + 0.0j
        sigma[j] = np.pi * abs(s_full[j] - ref) ** 2 / (2.0 * E)
    return sigma
```

- [ ] Replace the BODY of `_s_vector_one_energy` (keep the name, parameter list, and order exactly — `projects/n2_2d_td_cross_section/td_cross_section.py:20` imports it) with an assembly over the new kernel. Replace the whole function (lines 312-347) with:

```python
def _s_vector_one_energy(
    tgrid: TensorGrid,
    model: ResonanceModel,
    result: PropagationResult,
    eps: npt.NDArray[np.float64],
    v_init: int,
    vprimes: list[int],
    E: float,
    dt: float,
    wp_in: _WpIn,
    wp_out: _WpOut,
) -> npt.NDArray[np.complex128]:
    """The complex S-matrix `S_{v_init->v'}(E)` for each `v'`, shape `(len(vprimes),)`.

    `0` for closed channels (`E_tot - eps[v'] <= 0`) and for `E <= 0`. This
    is the raw Tannor-Weeks transform (module docstring) BEFORE the
    `|S - ref|^2` step -- now a thin assembly of the shared
    `s_vector_transform` skeleton with the TW outgoing factor
    (`correlation_channel_s` + `eta_outgoing` on the electronic axis).
    Kept under its original name/signature: the N2 project shim
    (`projects.n2_2d_td_cross_section.td_cross_section`) imports it.
    """
    g_elec = tgrid.grids[0]
    channel_s = correlation_channel_s(
        lambda kp: eta_outgoing(g_elec, kp, model.ell, **wp_out), result.c, dt
    )
    return s_vector_transform(
        g_elec,
        model.ell,
        wp_in,
        result.t,
        float(eps[v_init]),
        np.asarray([eps[vp] for vp in vprimes], dtype=np.float64),
        E,
        1.0,
        channel_s,
    )
```

- [ ] Replace the BODY of `_sigma_one_energy` (same constraint: name and positional signature stay — `test_core_td.py:38` and the N2 shim import it). Keep the existing docstring's ELASTIC-reference paragraphs verbatim (they are load-bearing physics documentation); replace only the code after the docstring (lines 386-404) with:

```python
    thresholds = np.asarray([eps[vp] for vp in vprimes], dtype=np.float64)
    elastic = np.asarray([vp == v_init for vp in vprimes], dtype=np.bool_)
    s_full = _s_vector_one_energy(tgrid, model, result, eps, v_init, vprimes, E, dt, wp_in, wp_out)
    s_free = None
    if free_result is not None:
        s_free = _s_vector_one_energy(
            tgrid, model, free_result, eps, v_init, vprimes, E, dt, wp_in, wp_out
        )
    return sigma_from_s(s_full, s_free, thresholds, float(eps[v_init]), E, elastic)
```

  (Behavior note, already accounted for: the old code early-returned zeros at `E <= 0` before building S-vectors; the new code calls `_s_vector_one_energy`, which itself early-returns zeros, then `sigma_from_s` returns zeros — identical output, negligible cost.)
- [ ] Verify (this exact block recurs in every task below as "**the standard gate**"):
  ```bash
  QSCAT_KERNEL_ORACLE_RTOL=1e-12 uv run --no-sync pytest libs/qscat/tests/test_kernel_consolidation_oracle.py -q -m "not slow"
  uv run --no-sync pytest -m "not slow" -n auto --dist loadfile
  uv run --no-sync mypy libs/qscat/qscat apps/qscat-run/qscat_run
  uv run --no-sync ruff check . && uv run --no-sync ruff format --check .
  ```
  The oracle's `ve_tw_*`, `sigma_from_correlations`, and `td_ve_tw` keys plus `test_td_extractors.py`'s pre-existing atol=1e-12 TW goldens are the direct gates for this task.
- [ ] Commit:
  ```bash
  git add libs/qscat/qscat/core/time_dependent.py
  git commit -m "refactor(core): one shared S-matrix/sigma transform kernel in time_dependent"
  ```

---

## Task 3 — lib-C4 part 2: drive all six `td_extractors.py` transforms through the shared kernel

**File:** `libs/qscat/qscat/core/td_extractors.py`

**Steps:**

- [ ] Update the imports (top of file):
  - From `.time_dependent`, import the kernels and drop what is no longer used. The import block (lines 152-157) becomes:
    ```python
    from .time_dependent import (
        Extractor,
        PropagationResult,
        correlation_channel_s,
        flux_channel_s,
        s_vector_transform,
        sigma_from_correlations,
        sigma_from_s,
    )
    ```
    (`quadrature_weights` is no longer used here after this task — ruff will confirm.)
  - From `.correlation` (lines 143-150), drop `eta_incident` (the shared kernel now computes it); keep `eta_outgoing`, `hankel_point_value`, `outgoing_channel`, `outgoing_channel_nuclear`, `outgoing_surface_wave`.
- [ ] DELETE these five functions entirely (verified: nothing outside this module references them): `_tw_da_s_vector_one_energy` (lines 188-228), `_dirac_s_vector_one_energy` (456-486), `_dirac_da_s_vector_one_energy` (532-570), `_flux_s_vector_one_energy` (827-869), `_flux_da_s_vector_one_energy` (915-954).
- [ ] Replace the five remaining per-energy sigma helpers with thin assemblies. Keep each function's name and parameter list unchanged (the extractor classes call them; no external importers). Replace `_tw_da_sigma_one_energy` with:

```python
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
    nuclear-axis Tannor-Weeks transform: the shared `s_vector_transform`
    skeleton with `eta_outgoing` moved to the NUCLEAR axis (mass `mu_r`,
    `l=0`) as the outgoing factor, then the shared `sigma_from_s` with
    `elastic=None` -- DA is a pure rearrangement channel with no
    `v'==v_init` diagonal, and its `pi*|S|^2/(2E)` prefactor is the same
    `_C_DA = pi` convention (see the module docstring)."""
    channel_s = correlation_channel_s(
        lambda k_r: eta_outgoing(g_nuc, k_r, 0, mass=mu_r, **wp_out), result.c, dt
    )
    s_full = s_vector_transform(
        g_elec, model.ell, wp_in, result.t, float(eps[v_init]), eps_e, E, mu_r, channel_s
    )
    return sigma_from_s(s_full, None, eps_e, float(eps[v_init]), E, None)
```

- [ ] Replace `_dirac_sigma_one_energy` with:

```python
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
    """`sigma_{v_init->v'}(E)` (bohr^2) via the delta transform, one energy:
    the shared `s_vector_transform` skeleton with `hankel_point_value` (the
    outgoing-Hankel-half VALUE at the fixed analysis point -- a delta test
    function's `F_out` is unintegrated) as the outgoing factor, then the
    shared `sigma_from_s`. Same elastic free-reference contract as
    `time_dependent._sigma_one_energy` (see that docstring): `S_free(E)`
    from a companion `V_int=0` propagation subtracts on the diagonal
    (`v'==v_init`) channel instead of a literal 1."""
    thresholds = np.asarray([eps[vp] for vp in vprimes], dtype=np.float64)
    elastic = np.asarray([vp == v_init for vp in vprimes], dtype=np.bool_)

    def s_vec(res: PropagationResult) -> npt.NDArray[np.complex128]:
        channel_s = correlation_channel_s(
            lambda kp: hankel_point_value(grid, z_position, kp, model.ell, model.charge),
            res.c,
            dt,
        )
        return s_vector_transform(
            grid, model.ell, wp_in, res.t, float(eps[v_init]), thresholds, E, 1.0, channel_s
        )

    s_full = s_vec(result)
    s_free = s_vec(free_result) if free_result is not None else None
    return sigma_from_s(s_full, s_free, thresholds, float(eps[v_init]), E, elastic)
```

- [ ] Replace `_dirac_da_sigma_one_energy` with:

```python
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
    nuclear-axis delta transform: `_dirac_sigma_one_energy`'s outgoing
    factor moved to the nuclear axis (`hankel_point_value` at mass `mu_r`,
    `l=0`); `elastic=None` (no free-reference -- DA has no `v'==v_init`
    diagonal), same `_C_DA = pi` convention."""
    channel_s = correlation_channel_s(
        lambda k_r: hankel_point_value(g_nuc, R_position, k_r, 0, model.charge, mass=mu_r),
        result.c,
        dt,
    )
    s_full = s_vector_transform(
        g_elec, model.ell, wp_in, result.t, float(eps[v_init]), eps_e, E, mu_r, channel_s
    )
    return sigma_from_s(s_full, None, eps_e, float(eps[v_init]), E, None)
```

- [ ] Replace `_flux_sigma_one_energy` with:

```python
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
    """`sigma_{v_init->v'}(E)` (bohr^2) via the flux transform, one energy:
    the shared `s_vector_transform` skeleton with the Wronskian element
    (`flux_channel_s` + `outgoing_surface_wave`'s `(phi_out, dphi_out)`
    pair, electronic mass 1.0), then the shared `sigma_from_s` -- same
    elastic free-reference pattern as `time_dependent._sigma_one_energy` /
    `_dirac_sigma_one_energy`."""
    thresholds = np.asarray([eps[vp] for vp in vprimes], dtype=np.float64)
    elastic = np.asarray([vp == v_init for vp in vprimes], dtype=np.bool_)

    def s_vec(
        t_a: npt.NDArray[np.float64],
        b_a: npt.NDArray[np.complex128],
        d_a: npt.NDArray[np.complex128],
    ) -> npt.NDArray[np.complex128]:
        channel_s = flux_channel_s(
            lambda kp: outgoing_surface_wave(grid, z_surface, kp, model.ell, model.charge),
            b_a,
            d_a,
            dt,
            1.0,
        )
        return s_vector_transform(
            grid, model.ell, wp_in, t_a, float(eps[v_init]), thresholds, E, 1.0, channel_s
        )

    s_full = s_vec(t, b, d)
    s_free = s_vec(*free) if free is not None else None
    return sigma_from_s(s_full, s_free, thresholds, float(eps[v_init]), E, elastic)
```

- [ ] Replace `_flux_da_sigma_one_energy` with:

```python
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
    nuclear-axis flux transform: `_flux_sigma_one_energy`'s Wronskian
    element with the outgoing wave moved to the nuclear axis (mass `mu_r`,
    `l=0`); `elastic=None` (DA has no `v'==v_init` diagonal)."""
    channel_s = flux_channel_s(
        lambda k_r: outgoing_surface_wave(g_nuc, R_surface, k_r, 0, model.charge, mass=mu_r),
        b,
        d,
        dt,
        mu_r,
    )
    s_full = s_vector_transform(
        g_elec, model.ell, wp_in, t, float(eps[v_init]), eps_e, E, mu_r, channel_s
    )
    return sigma_from_s(s_full, None, eps_e, float(eps[v_init]), E, None)
```

- [ ] Update the prose that references the deleted functions (the module's docstrings are load-bearing — do not leave dangling names):
  - Module docstring, `TannorWeeks(axis="nuclear")` section (lines 118-120): replace ``(this module's own `_tw_da_s_vector_one_energy`/`_tw_da_sigma_one_energy`, defined below, NOT routed through `sigma_from_correlations`: ...)`` with ``(this module's own `_tw_da_sigma_one_energy`, an assembly of `time_dependent`'s shared `s_vector_transform`/`sigma_from_s` kernels, NOT routed through `sigma_from_correlations`: ...)``.
  - Module docstring, `Dirac(axis="nuclear")` section (line 94): replace ``Its `sigma` is `_dirac_s_vector_one_energy`'s nuclear-axis twin`` with ``Its `sigma` is the electronic delta transform's nuclear-axis twin``.
  - `_C_DA` comment block (lines 169-173): unchanged — `_C_DA` itself STAYS defined (the docstrings and `docs/physics/td-da.md` reason about it by name), even though the shared `sigma_from_s` now carries the constant as `np.pi`; append one line to its comment: `# The shared sigma_from_s kernel realizes this constant as its np.pi prefactor.` If ruff flags `_C_DA` as unused after the rewrite, keep it and silence by referencing it where `Dirac`/`Flux`/`TannorWeeks` docstrings already do — it is module documentation; if an actual `F401`-style violation fires (it should not for module-level assignments), delete `_C_DA` and fold its comment into the `sigma_from_s` reference in this module's docstring instead.
- [ ] Verify — **the standard gate** (Task 2's four commands, verbatim). The oracle's `ve_delta_*`, `ve_flow_*`, `da_tw/delta/flow_sigma`, `td_ve_delta/flow`, and `td_da_*` keys, plus `test_td_extractors.py`'s whole fast set, are the direct gates.
- [ ] Commit:
  ```bash
  git add libs/qscat/qscat/core/td_extractors.py
  git commit -m "refactor(core): drive all six extractor transforms through the shared kernel"
  ```

---

## Task 4 — lib-m5 + lib-m6: `correlation.py` dedupe

**File:** `libs/qscat/qscat/core/correlation.py`

**Steps:**

- [ ] Add the shared coefficient helper immediately before `_regular_coeffs` (line 122):

```python
def _masked_coeffs(
    grid: FemDvrEcsGrid, f_vals: npt.NDArray[np.complex128] | npt.NDArray[np.float64]
) -> npt.NDArray[np.complex128]:
    """Function VALUES -> masked DVR coefficients: multiply by `sqrt(w)`
    (the bridge-summed complex weight -- the same conversion
    `channel_vector` applies; a coefficient vector like `chi_v` is never
    passed through here) and zero everything past the ECS pivot `R0`."""
    sqrt_w = np.sqrt(np.asarray(grid.weights, dtype=np.complex128))
    coeffs = (f_vals * sqrt_w).astype(np.complex128)
    coeffs[grid.real_points > grid.R0] = 0.0
    return coeffs
```

- [ ] Rewrite `_regular_coeffs` and `_outgoing_coeffs` as thin wrappers (both names stay — `test_correlation.py`'s docstrings reference `_outgoing_coeffs`, and the module docstring explains both):

```python
def _regular_coeffs(grid: FemDvrEcsGrid, k: float, l: int) -> npt.NDArray[np.complex128]:
    """`riccati_bessel_en(r, k, l) * sqrt(w_r)`, masked to the unscaled region
    (`_masked_coeffs`) -- the regular-function side of `eta_incident`."""
    return _masked_coeffs(grid, riccati_bessel_en(grid.real_points, k, l))


def _outgoing_coeffs(
    grid: FemDvrEcsGrid, k: float, l: int, *, mass: float = 1.0
) -> npt.NDArray[np.complex128]:
    """`h^{(1)}_{E,l}(r)/2 * sqrt(w_r)`, masked to the unscaled region
    (`_masked_coeffs`). `h^{(1)}` is the energy-normalized mass-`mass`
    outgoing Hankel half -- see the module docstring for why this (not the
    regular function) is `F_out`; `mass=1.0` reproduces the electronic path
    bit-for-bit (`riccati_hankel_en_mass`'s docstring), a nuclear (DA)
    caller passes `mass=model.mu`."""
    return _masked_coeffs(grid, riccati_hankel_en_mass(grid.real_points, k, l, mass) / 2.0)
```

  (Bit-identity note: the old `_outgoing_coeffs` also computed `f_vals = riccati/2.0` *before* the `* sqrt_w` — same operation order, identical bits.)
- [ ] In `outgoing_surface_wave`, replace the neutral branch's inline `phi` (lines 276-283) so the VALUE comes from `hankel_point_value` (defined above it in this file) and only the DERIVATIVE keeps the raw scipy pieces:

```python
    r = float(z_surface)
    if charge == 0:
        # The VALUE is exactly the outgoing-Hankel-half `hankel_point_value`
        # already provides (same multiplication order as the old inline
        # formula -- bit-identical); only the DERIVATIVE still needs
        # scipy's `derivative=True` pieces, because `qscat.special.radial`
        # exposes no derivative primitive (module docstring of `radial`).
        phi = hankel_point_value(grid, r, k, l, 0, mass=mass)
        x = k * r
        h_l = spherical_jn(l, x) + 1j * spherical_yn(l, x)
        h_l_prime = spherical_jn(l, x, derivative=True) + 1j * spherical_yn(l, x, derivative=True)
        dphi = np.sqrt(2.0 * mass * k / np.pi) * (h_l + x * h_l_prime) / 2.0
        return complex(phi), complex(dphi)
```

  Then update `outgoing_surface_wave`'s docstring paragraph that begins "Neutral (`charge == 0`): computed ANALYTICALLY..." to mention that the value half now comes from `hankel_point_value` and only the derivative is computed here. The `del grid` line at the top of the old function body (line 274) must be REMOVED — `grid` is now forwarded to `hankel_point_value` (which itself ignores it, keeping the call-site symmetry both docstrings describe; also update the "`grid` is accepted (unused)..." sentence accordingly). The charged branch (mpmath finite difference) stays byte-identical.
- [ ] Verify — **the standard gate**. The oracle's `correlation_factors` key plus `test_correlation.py` (which finite-difference-checks the analytic derivative) are the direct gates.
- [ ] Commit:
  ```bash
  git add libs/qscat/qscat/core/correlation.py
  git commit -m "refactor(core): dedupe correlation.py coefficient and Hankel-value helpers"
  ```

---

## Task 5 — lib-M2: `dr_cross_section` reuses `ve_cross_section`'s driven sweep

**Files:** `libs/qscat/qscat/core/driven.py`, `libs/qscat/qscat/core/dissociation.py`

**Steps:**

- [ ] In `driven.py`, thread the model's charge into the channel functions. `_sigma_at_one_energy` gains a keyword-only `charge: int = 0` (insert before `want_psi` in the signature) and forwards it to BOTH `channel_vector` calls:
  - line 94: `psi_i = channel_vector(tgrid, k, chi[v_init], l, charge=charge)`
  - line 103: `phi_f = channel_vector(tgrid, kp, chi[vp], l, charge=charge)`

  and `ve_cross_section`'s call site (line 210-221) passes `charge=model.charge`. For every existing caller `model.charge == 0`, and `channel_vector`'s `charge == 0` branch is the identical pre-existing `riccati_bessel_en` path — bit-for-bit unchanged. Add to `ve_cross_section`'s docstring, after the `model` paragraph (line 172-174):

  ```
  `model.charge` is forwarded to `channel_vector`, so an IONIC target's
  entrance/exit channels are the energy-normalized Coulomb functions
  (`coulomb_f_en`) rather than the free `riccati_bessel_en`; `charge == 0`
  (every neutral model) takes the identical pre-existing free-function
  branch. This is what lets `dissociation.dr_cross_section` reuse this
  solver's driven sweep for H2+ instead of re-inlining it.
  ```
- [ ] In `dissociation.py`, replace `dr_cross_section`'s inlined sweep with a `ve_cross_section` reuse. Replace the body from line 410 (`e_arr = ...`) through line 457 (end of the per-energy loop) with:

```python
    e_arr = np.atleast_1d(np.asarray(E, dtype=np.float64))
    mu = model.mu
    g_R = tgrid.grids[1]
    R_inf = g_R.R0

    eps_ryd, phi_ryd = anion_electronic_states(
        g_r=tgrid.grids[0], model=model, R_inf=R_inf, n_states=n_channels
    )
    v_dr = v_dr_diag(tgrid, model)
    mask = tgrid.real_mask()
    sqrt_w_R = tgrid.sqrt_weights()[1].ravel()

    # The driven Psi+ sweep is ve_cross_section's (analyze-once,
    # `SparseLU.refactor` per energy; `model.charge` forwarded to the
    # incident `channel_vector`, so H2+'s Coulomb entrance is built there)
    # -- exactly the reuse `da_cross_section` already performs. Its VE
    # sigma for the [v_init] channel is discarded; the extra cost is one
    # exit `channel_vector` + c-product per energy, marginal next to the
    # factorization.
    _, psis = ve_cross_section(
        tgrid,
        model,
        eps,
        chi,
        v_init,
        [v_init],
        e_arr,
        ordering=ordering,
        return_wavefunction=True,
    )
    psi_list = psis if isinstance(psis, list) else [psis]

    out = np.zeros((len(e_arr), n_channels), dtype=np.float64)
    amp = np.zeros((len(e_arr), n_channels), dtype=np.complex128)
    for ie, e in enumerate(e_arr):
        psi_plus = psi_list[ie]
        if psi_plus is None:  # E <= 0: below threshold, sigma == 0
            continue
        e_tot = float(e) + eps[v_init]
        v_psi = v_dr * psi_plus

        for n in range(n_channels):
            e_dr = e_tot - eps_ryd[n]
            if e_dr <= 0.0:
                continue  # closed Rydberg channel
            k_r = float(np.sqrt(2.0 * mu * e_dr))
            y_coeff = riccati_bessel_en_mass(g_R.real_points, k_r, 0, mu) * sqrt_w_R
            phi_f = tgrid.outer([phi_ryd[n], y_coeff])
            phi_f[~mask] = 0.0
            t = c_product(phi_f, v_psi)
            amp[ie, n] = t
            out[ie, n] = 4.0 * np.pi**3 * abs(t) ** 2 / (2.0 * float(e))
```

  The tail of the function (from `scalar = np.isscalar(E) ...`, line 459 onward) is unchanged — note it reads `psi_list`, which the new code still defines.
- [ ] Update `dr_cross_section`'s docstring: DELETE the stale paragraph at lines 383-387 ("The driven Lippmann-Schwinger solve for `Psi+` is replicated inline (rather than reusing `ve_cross_section`) because that helper cannot pass `charge` through to `channel_vector`; ...") and replace it with:

  ```
  The driven Lippmann-Schwinger solve for `Psi+` reuses
  `ve_cross_section(..., return_wavefunction=True)` -- the same
  analyze-once / `SparseLU.refactor`-per-energy sweep, with `model.charge`
  forwarded to `channel_vector` so the incident channel is Coulomb --
  exactly as `da_cross_section` does. Only the exit-channel read differs:
  DR projects the post-form volume T-matrix against `V_DR` (below), DA
  reads a boundary flux.
  ```
- [ ] Clean up `dissociation.py`'s now-unused imports: after the rewrite `sp` (`scipy.sparse`, only used for `sp.identity` in the deleted sweep), `SparseLU`, and `channel_vector` are unused — remove `import scipy.sparse as sp` (line 47), change line 50 to `from qscat.linalg import c_product`, and delete line 53 (`from .channels import channel_vector`). Keep `_Ordering` (still forwards) and `riccati_bessel_en_mass` (still used). Also delete the now-dangling local variables the old sweep needed (`H`, `v_diag`, `ident`, `lu`) — the replacement body above already omits them.
- [ ] Verify — **the standard gate**, PLUS the slow dr oracle and the slow dr suite (this task is the one that touches the H2+ path, so pay the minutes now):
  ```bash
  uv run --no-sync pytest libs/qscat/tests/test_kernel_consolidation_oracle.py::test_dr_golden -q
  uv run --no-sync pytest -m slow libs/qscat/tests/test_dissociation.py -q
  ```
  `test_dr_cproduct_matches_conjugated_dot_on_proxy`, `test_dr_amplitude_reproduces_the_returned_sigma`, and `test_dr_amplitude_matches_conjugated_oracle_value_and_phase` are independent mirrors of the old inlined arithmetic — they gate this rewrite at physics precision.
- [ ] Commit:
  ```bash
  git add libs/qscat/qscat/core/driven.py libs/qscat/qscat/core/dissociation.py
  git commit -m "refactor(core): dr_cross_section reuses ve_cross_section's driven sweep"
  ```

---

## Task 6 — lib-M1: collapse `td_ve_cross_section`'s three method blocks

**File:** `libs/qscat/qscat/core/time_dependent.py`

**Steps:**

- [ ] Replace `td_ve_cross_section`'s body after the docstring (lines 539-644: from `from .td_extractors import ...` to the closing `raise ValueError`) with the single-dispatch shape `td_da_cross_section` already uses:

```python
    from .td_extractors import Dirac, Flux, TannorWeeks  # deferred: avoids an import cycle

    psi0 = initial_state(tgrid, chi[v_init], **wp_in)
    hamiltonian = model.hamiltonian(tgrid)

    def build_extractor() -> Extractor:
        # One construction site for both the full run and the free-reference
        # run -- the two runs must use identically-configured extractors.
        if method == "tw":
            return TannorWeeks(tgrid, model, eps, chi, v_init, vprimes, wp_out, wp_in=wp_in, dt=dt)
        if method == "delta":
            if position is None:
                raise ValueError("td_ve_cross_section: method='delta' requires `position`")
            return Dirac(tgrid, model, eps, chi, v_init, vprimes, position, wp_in=wp_in, dt=dt)
        if method == "flow":
            if surface is None:
                raise ValueError("td_ve_cross_section: method='flow' requires `surface`")
            return Flux(tgrid, model, eps, chi, v_init, vprimes, surface, wp_in=wp_in, dt=dt)
        raise ValueError(
            f"td_ve_cross_section: unknown method {method!r} (must be one of 'tw', 'delta', 'flow')"
        )

    ext = build_extractor()
    propagate(
        tgrid,
        psi0,
        [],
        dt=dt,
        n_steps=n_steps,
        hamiltonian=hamiltonian,
        order=order,
        extractors=[ext],
    )

    ext_free: Extractor | None = None
    if subtract_free_reference and v_init in vprimes:
        ext_free = build_extractor()
        propagate(
            tgrid,
            psi0,
            [],
            dt=dt,
            n_steps=n_steps,
            hamiltonian=_free_hamiltonian(model, tgrid),
            order=order,
            extractors=[ext_free],
        )

    return ext.sigma(E, free=ext_free)
```

  Behavior notes (all verified against the old code): the unknown-`method` and missing-`position`/`surface` `ValueError` messages are byte-identical to the old ones (`test_td_extractors.py` matches on them); the only observable ordering change is that an unknown `method` now raises *before* the first propagation instead of after none-of-three-blocks matched — the old code also never propagated in that case, so nothing observable differs. `ext.sigma(E, free=ext_free)` through the `Extractor` protocol is the same call each block made concretely; `td_da_cross_section` already annotates `ext: Extractor` and calls `ext.sigma(E)`, so the protocol precedent is established.
- [ ] Trim the last paragraph of `td_ve_cross_section`'s docstring (lines 529-537, "`method="tw"` builds a `td_extractors.TannorWeeks` extractor ... mirror the SAME pattern ...") to describe the single dispatch:

  ```
  All three methods share one code path: the selected extractor is built
  (twice, identically, when the free-reference run applies), driven by
  `propagate(..., out_channels=[], extractors=[...])`, and asked for
  `sigma(E, free=...)`. `method="tw"` reproduces this function's
  pre-refactor implementation to machine precision (see
  `libs/qscat/tests/test_td_extractors.py`'s golden regression tests);
  any other `method` value raises `ValueError`.
  ```
- [ ] Verify — **the standard gate**. Direct gates: the oracle's `td_ve_tw/delta/flow` keys; `test_td_extractors.py`'s golden + `test_delta_method_matches_direct_dirac_construction` / `test_flow_method_matches_direct_flux_construction` / `test_unknown_method_raises`.
- [ ] Commit:
  ```bash
  git add libs/qscat/qscat/core/time_dependent.py
  git commit -m "refactor(core): collapse td_ve_cross_section's three method blocks"
  ```

---

## Task 7 — lib-M17: merge the four radial functions into two

**File:** `libs/qscat/qscat/special/radial.py`

**Steps:**

- [ ] Replace the four function definitions (lines 44-121) with two mass-general functions plus two thin alias wrappers. `__all__` (lines 36-41) is unchanged — all four names remain importable, and `special/__init__.py` needs no edit.

```python
def riccati_bessel_en(
    r: npt.NDArray[np.float64], k: float, l: int, mu: float = 1.0
) -> npt.NDArray[np.float64]:
    """`F_{E,l}(r) = sqrt(2 mu k / pi) r j_l(k r)`, energy-normalized at mass `mu`.

    The energy-normalized (`<F_E|F_E'> = delta(E-E')`) REGULAR radial
    solution for a particle of reduced mass `mu` and momentum
    `k = sqrt(2 mu E)`. `mu=1.0` (the default) is the electron case --
    bit-for-bit the historical mass-1 function, since `2.0*1.0 == 2.0`
    exactly and `mu` enters ONLY the normalization prefactor, never the
    momentum argument `k r`. `mu != 1` is the nuclear case, used for the
    OUTGOING NUCLEAR dissociation wave in the DA/DR exit channel (eMoScat
    `bessel::s_jEn(R, K, mu, l)`). `r` must be REAL (see module docstring);
    `k > 0`, `mu > 0`.
    """
    if k <= 0.0:
        raise ValueError(f"k must be positive, got {k}")
    if mu <= 0.0:
        raise ValueError(f"mu must be positive, got {mu}")
    rr = np.asarray(r, dtype=np.float64)
    out: npt.NDArray[np.float64] = np.sqrt(2.0 * mu * k / np.pi) * rr * spherical_jn(l, k * rr)
    return out


def riccati_hankel_en(
    r: npt.NDArray[np.float64], k: float, l: int, mu: float = 1.0
) -> npt.NDArray[np.complex128]:
    """`F^{(1)}_{E,l}(r) = sqrt(2 mu k / pi) r h_l^{(1)}(k r)`, energy-normalized
    at mass `mu`, `h_l^{(1)} = j_l + i y_l` the OUTGOING spherical Hankel
    function.

    The outgoing sibling of `riccati_bessel_en`, with the same mass
    convention: `mu` enters only the `sqrt(mu)` normalization prefactor
    (`2.0*1.0 == 2.0` exactly, so `mu=1.0` is bit-for-bit the historical
    mass-1 function), never the momentum argument `k r`. So
    `Re(F^{(1)}) == riccati_bessel_en` and `Im(F^{(1)}) = sqrt(2 mu k/pi)
    r y_l(k r)` at the same `(r, k, l, mu)`. `mu != 1` drives the OUTGOING
    NUCLEAR dissociation wave in the DA/DR flux (Wronskian) extractor
    (eMoScat `bessel::sphHankel1En(R, K, mu, l)` -- see
    `qscat.core.td_extractors.Flux`). `r` must be REAL (module docstring);
    `k > 0`, `mu > 0`.
    """
    if k <= 0.0:
        raise ValueError(f"k must be positive, got {k}")
    if mu <= 0.0:
        raise ValueError(f"mu must be positive, got {mu}")
    rr = np.asarray(r, dtype=np.float64)
    h1_l = spherical_jn(l, k * rr) + 1j * spherical_yn(l, k * rr)
    out: npt.NDArray[np.complex128] = (
        np.sqrt(2.0 * mu * k / np.pi) * rr.astype(np.complex128) * h1_l
    )
    return out


def riccati_bessel_en_mass(
    r: npt.NDArray[np.float64], k: float, l: int, mu: float
) -> npt.NDArray[np.float64]:
    """Deprecated alias for `riccati_bessel_en(r, k, l, mu)`: the mass
    generalization lives on the base name now. Kept so existing imports
    keep working; new code should call `riccati_bessel_en` directly."""
    return riccati_bessel_en(r, k, l, mu)


def riccati_hankel_en_mass(
    r: npt.NDArray[np.float64], k: float, l: int, mu: float
) -> npt.NDArray[np.complex128]:
    """Deprecated alias for `riccati_hankel_en(r, k, l, mu)`: the mass
    generalization lives on the base name now. Kept so existing imports
    keep working; new code should call `riccati_hankel_en` directly."""
    return riccati_hankel_en(r, k, l, mu)
```

  Two deliberate details: (1) `mu` in the aliases stays REQUIRED-positional, exactly matching the old `_mass` signatures; (2) the merged base functions gain the `mu <= 0` guard the `_mass` variants had — the default `mu=1.0` never trips it, so no existing caller can hit a new error.
- [ ] Update the module docstring's two formula lines (lines 6-16) to mention the `mu` parameter (e.g. append ", with an optional reduced mass `mu` entering only the normalization prefactor -- `mu=1.0` is the electron case" after the `riccati_bessel_en` description). Keep the eMoScat source anchors as they are.
- [ ] Verify — **the standard gate**. Direct gates: the oracle's `radial_functions` key and its two `np.array_equal` mu=1 identities; `test_special_radial.py`, `test_radial_mass.py`, `test_radial_hankel_mass.py` (all exercise the alias names, unchanged).
- [ ] Commit:
  ```bash
  git add libs/qscat/qscat/special/radial.py
  git commit -m "refactor(special): mass-generalize the Riccati radial functions in place"
  ```

---

## Task 8 — Closeout: leftovers audit, slow-tier gates, docs check

**Files:** none expected to change (this task verifies; it only edits if an audit step finds something).

**Steps:**

- [ ] Leftover audit — each of these greps must come back empty (excluding this plan file and `docs/_build`):
  ```bash
  grep -rn "_tw_da_s_vector_one_energy\|_dirac_s_vector_one_energy\|_dirac_da_s_vector_one_energy\|_flux_s_vector_one_energy\|_flux_da_s_vector_one_energy" \
      libs apps projects validation benchmarks docs/physics
  ```
  If a hit remains in prose, fix the prose (per-task docstring steps should already have covered every case).
- [ ] Docs check — no public name changed in this plan (`_s_vector_one_energy`/`_sigma_one_energy` kept, extractor classes and all `td_*`/`*_cross_section`/radial names kept), so `docs/physics/td-extractors.md` and `docs/physics/td-da.md` should need NO edits. Verify rather than assert:
  ```bash
  grep -rn "s_vector_one_energy\|sigma_one_energy\|riccati_bessel_en_mass\|riccati_hankel_en_mass" docs/physics/td-extractors.md docs/physics/td-da.md
  ```
  Any hit that names a *deleted* function gets its sentence updated to the surviving name (`s_vector_transform`/`sigma_from_s` or the kept sigma helper); hits naming kept functions/aliases stay. Also confirm the built API docs are not tracked (`git ls-files docs/_build | head` — expect empty; if anything IS tracked there, rebuild the docs per the Docs workflow instead of hand-editing).
- [ ] Full oracle, including the slow dr key, at the gate tolerance:
  ```bash
  QSCAT_KERNEL_ORACLE_RTOL=1e-12 uv run --no-sync pytest libs/qscat/tests/test_kernel_consolidation_oracle.py -q
  ```
- [ ] Slow-tier suites that exercise the touched paths (serial, no `-n` — these are the multi-GB decks; expect tens of minutes total):
  ```bash
  uv run --no-sync pytest -m slow libs/qscat/tests/test_td_extractors.py libs/qscat/tests/test_dissociation.py -q
  uv run --no-sync pytest -m slow projects/n2_2d_td_cross_section/test_td_cross_section.py -q
  uv run python -m validation.n2.experiment
  ```
  The named load-bearing gates inside those:
  - `libs/qscat/tests/test_td_extractors.py::test_delta_agrees_with_tw_same_trajectory`
  - `libs/qscat/tests/test_td_extractors.py::test_delta_agrees_with_ti_oracle_one_anchor`
  - `libs/qscat/tests/test_td_extractors.py::test_flux_agrees_with_tw_same_trajectory`
  - `libs/qscat/tests/test_td_extractors.py::test_flux_agrees_with_ti_oracle_one_anchor`
  - `libs/qscat/tests/test_td_extractors.py::test_nuclear_flux_da_converges_to_ti_oracle`
  - `libs/qscat/tests/test_td_extractors.py::test_nuclear_dirac_da_converges_to_ti_oracle`
  - `libs/qscat/tests/test_td_extractors.py::test_nuclear_tw_da_converges_to_ti_oracle`
  - `libs/qscat/tests/test_dissociation.py::test_dr_wellposed_and_threshold_ordered`
  - `libs/qscat/tests/test_dissociation.py::test_dr_cproduct_matches_conjugated_dot_on_proxy`
  - `libs/qscat/tests/test_dissociation.py::test_dr_amplitude_reproduces_the_returned_sigma`
  - `libs/qscat/tests/test_dissociation.py::test_dr_amplitude_composes_with_the_wavefunction_return`
  - `libs/qscat/tests/test_dissociation.py::test_dr_amplitude_matches_conjugated_oracle_value_and_phase`
  - `libs/qscat/tests/test_dissociation.py::test_da_shape_scalar_and_array` (and the other `@slow` da tests in that file)
  - `projects/n2_2d_td_cross_section/test_td_cross_section.py::test_v2a_td_matches_ti_at_e010` / `::test_v2a_td_matches_ti_at_e015`
  - the `validation.n2.experiment` harness's C5/D1/E1 group PASS lines
- [ ] **The standard gate**, one last time (fast tier + mypy + ruff + ruff format).
- [ ] On the PR, request the `validate:core` + `validate:n2` labelled CI run per the Global Constraints (the local slow runs above satisfy the same requirement if CI labels are unavailable).
- [ ] Commit only if an audit step changed a file:
  ```bash
  git add <exact paths changed by the audit, if any>
  git commit -m "chore(core): kernel-consolidation closeout fixes from the leftovers audit"
  ```

---

## Self-review notes (kept for the executor)

- Every task's code snippets are semantically exact but the formatter has final say on layout: after editing a file, run `uv run --no-sync ruff format <file>` before the standard gate so `ruff format --check .` passes. Formatting never changes the golden comparisons.
- Name/type consistency across tasks: Task 2 defines `ChannelS`, `s_vector_transform(g_in, l_in, wp_in, t, eps_init, thresholds, E, exit_mass, channel_s)`, `correlation_channel_s(outgoing_factor, c, dt)`, `flux_channel_s(outgoing_pair, b, d, dt, exit_mass)`, `sigma_from_s(s_full, s_free, thresholds, eps_init, E, elastic)`; Tasks 3 and 6 consume exactly these names/orders. `FemDvrEcsGrid` is imported into `time_dependent.py` in Task 2 and is already imported in `td_extractors.py` and `correlation.py`.
- Bit-identity of every substituted expression was checked during planning: `2.0*mass*x` at `mass=1.0` reduces to `2.0*x` exactly (left-associative `2.0*1.0 == 2.0`); `abs(s - (0+0j)) == abs(s)` for all IEEE inputs including signed zeros; `1.0 * v` (driven's `lam_scale`) is bit-identical; `hankel_point_value`'s value formula multiplies in the same order as the old inline `phi`; `float(np.float64)` conversions are exact.
- Survivors and their reasons: `_s_vector_one_energy` + `_sigma_one_energy` (imported by the N2 shim and `test_core_td.py`), `_free_hamiltonian` alias (imported by `test_td_extractors.py`), all four radial names (aliases; imported across `correlation.py`/`dissociation.py`/`tuning/probes.py`/tests), `_regular_coeffs`/`_outgoing_coeffs` (referenced by test docstrings), `_C_DA` (referenced by physics docs/docstrings).
- The oracle's two-layer tolerance design is a deliberate deviation from a naive single rtol=1e-12 pin; the measured ~4e-9 cross-process BLAS drift that forces it is documented in `libs/qscat/tests/test_core_td.py`'s module docstring and restated in Global Constraints above.
