# N₂ Time-Dependent VE Cross-Section Implementation Plan (sub-project #4)

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development. The TI cross section (#3) is the EXACT differential oracle — TD must converge to it. Internal checks (TD≈TI, convergence) are the gate.

**Goal:** Compute N₂ VE cross sections by Crank-Nicolson wavepacket propagation of the LCP model, cross-check against the TI result (#3), and flip harness Group D. Promote the CN propagator to `qscat.evolution`.

**Architecture:** `projects/n2_td_cross_section/`: `propagator.py` (Crank-Nicolson), `td_cross_section.py` (propagate + correlate + energy-transform). Reuses `projects.n2_ti_cross_section.{nuclear_grid,vibrational,vres,cross_section}` and `qscat.dvr`.

**Tech Stack:** Python 3.12, numpy; `qscat.dvr` (kinetic → T_nuc).

## Global Constraints

- Python `>=3.12`; `uv run pytest`. Atomic units (Hartree, bohr, bohr²). μ=12766.36.
- Packages are package-qualified now (post-refactor): `from projects.n2_ti_cross_section.vres import vres_on_grid`, etc. `projects/n2_td_cross_section/` HAS an `__init__.py` (it's a package). No sys.path hacks.
- **Formulas:** ψ(0)=d_{v_init}=√(Γ/2π)·χ_{v_init}. CN: `(I+iH·dt/2)ψ_{n+1}=(I−iH·dt/2)ψ_n`, `H=H_res=kinetic(grid,μ)+diag(V_d−iΓ/2)` (time-independent ⇒ LU-factor `A=I+iH·dt/2` once). Correlation `c_{v'}(t_n)=Σ_j d_{v'}[j]·ψ_n[j]` (c-product, NO conjugate). Energy transform `S_{v'}(E)=(1/i)·Σ_n w_n·e^{i(E+ε_{v_init})t_n}·c_{v'}(t_n)·dt` (Simpson weights `w_n`; `t_n=n·dt`). `σ=4π³|S|²/(2E)`, 0 if `E_tot−ε_{v'}≤0`.
- **The oracle:** `projects.n2_ti_cross_section.cross_section.ve_cross_section(...)` gives σ_TI. TD must match it (rtol ≤ 10%) at converged dt/n_steps.
- Commit trailer: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.

---

### Task 1: Crank-Nicolson propagator

**Files:** Create `projects/n2_td_cross_section/__init__.py` (empty), `projects/n2_td_cross_section/propagator.py`; Test `projects/n2_td_cross_section/test_propagator.py`.

**Interfaces:** `propagator.make_cn_stepper(H: ndarray, dt: float) -> Callable[[ndarray], ndarray]` — returns a `stepper(psi)` that advances one CN step, with `A=I+iH·dt/2` LU-factored once (use `scipy.linalg.lu_factor`/`lu_solve`, or `np.linalg.solve` with a cached factorization; numpy-only is fine via `np.linalg.solve(A, B@psi)` but prefer caching the factorization).

- [ ] **Step 1: Write failing tests — `test_propagator.py`**

```python
import numpy as np
from projects.n2_td_cross_section.propagator import make_cn_stepper


def test_cn_matches_exact_exp_for_hermitian():
    # small Hermitian H: CN step ≈ exp(-i H dt) to O(dt^3) per step
    rng = np.random.default_rng(0)
    A = rng.standard_normal((5, 5)) + 1j * rng.standard_normal((5, 5))
    H = A + A.conj().T                      # Hermitian
    dt = 1e-3
    step = make_cn_stepper(H, dt)
    psi0 = rng.standard_normal(5) + 1j * rng.standard_normal(5)
    exact = (np.linalg.matrix_power(np.eye(5), 0) @ psi0)  # placeholder
    from scipy.linalg import expm
    exact = expm(-1j * H * dt) @ psi0
    assert np.linalg.norm(step(psi0) - exact) < 1e-8       # O(dt^3) ~ 1e-9


def test_cn_unitary_for_hermitian():
    rng = np.random.default_rng(1)
    A = rng.standard_normal((6, 6)) + 1j * rng.standard_normal((6, 6))
    H = A + A.conj().T
    step = make_cn_stepper(H, 0.1)
    psi = rng.standard_normal(6) + 1j * rng.standard_normal(6)
    n0 = np.vdot(psi, psi).real
    for _ in range(50):
        psi = step(psi)
    assert abs(np.vdot(psi, psi).real - n0) < 1e-10        # CN preserves norm for Hermitian


def test_cn_decays_for_non_hermitian_decaying():
    # H with negative imaginary part (decaying) -> norm decreases
    H = np.diag([1.0 - 0.1j, 2.0 - 0.2j, 3.0 - 0.05j])
    step = make_cn_stepper(H, 0.05)
    psi = np.ones(3, dtype=complex)
    n0 = np.vdot(psi, psi).real
    for _ in range(100):
        psi = step(psi)
    assert np.vdot(psi, psi).real < n0                      # decays
```

- [ ] **Step 2: Run → fail.** `uv run pytest projects/n2_td_cross_section/test_propagator.py -q`.
- [ ] **Step 3: Implement `propagator.py`.**

```python
from __future__ import annotations
from typing import Callable
import numpy as np
from scipy.linalg import lu_factor, lu_solve


def make_cn_stepper(H: np.ndarray, dt: float) -> Callable[[np.ndarray], np.ndarray]:
    n = H.shape[0]
    ident = np.eye(n, dtype=complex)
    A = ident + 0.5j * dt * H          # (I + i H dt/2)
    B = ident - 0.5j * dt * H          # (I - i H dt/2)
    lu = lu_factor(A)

    def stepper(psi: np.ndarray) -> np.ndarray:
        return lu_solve(lu, B @ psi)
    return stepper
```
(scipy is already a dev dependency.)

- [ ] **Step 4: Run → pass** (3 tests). **Commit.**

---

### Task 2 (crux, low-risk): TD cross section + TD≈TI + convergence

**Files:** Create `projects/n2_td_cross_section/td_cross_section.py`; Test `test_td_cross_section.py`.

**Interfaces:** `td_cross_section.td_ve_cross_section(grid, mu, Vd, Gamma, eps, chi, v_init, vprimes, E, *, dt, n_steps) -> ndarray` — σ for the listed final channels at collision energy(ies) E, via CN propagation + energy transform.

- [ ] **Step 1: Write the TD≈TI + convergence tests — `test_td_cross_section.py`.** Build the shared setup ONCE at module scope (nuclear grid, `vibrational_states`, `vres_on_grid` — ~7s): 
  - **V1 (TD≈TI):** at anchor (E=0.1 Ha, v'=1) and one more (E=0.2, v'=2), compute σ_TD via `td_ve_cross_section` (converged dt/n_steps) and σ_TI via `projects.n2_ti_cross_section.cross_section.ve_cross_section`; assert `rtol ≤ 0.1` (10%). Also assert σ_TD real & ≥0.
  - **V2 (convergence):** σ_TD at (0.1, v'=1) with dt and dt/2 agree to `rtol ≤ 0.05`; and ‖ψ(T)‖ < 0.1·‖ψ(0)‖ (resonance depleted — expose or check via the propagation).
- [ ] **Step 2: Run → fail.**
- [ ] **Step 3: Implement `td_cross_section.py`.** doorways `d_v=√(Γ/2π)·chi[v]`; `H_res=kinetic(grid,μ)+diag(Vd−0.5j·Gamma)`; `step=make_cn_stepper(H_res, dt)`; `psi=d_{v_init}`; loop `n_steps`, at each `t_n=n·dt` record `c_{v'}(t_n)=chi/… ` → the c-product `Σ_j d_{v'}[j]·psi[j]` for each requested v'; after the loop, `S_{v'}(E)=(1/i)·Σ_n w_n·exp(1j·(E+eps[v_init])·t_n)·c_{v'}(t_n)·dt` with Simpson weights `w_n` (or trapezoidal if simpler and convergence holds); `sigma=4π³|S|²/(2E)`, 0 if closed. Accept scalar or array E (transform is cheap once c(t) is stored). Expose the final norm ratio for the V2 check.
- [ ] **Step 4: Run.** Tune dt (~0.5–5 a.u.) and n_steps (~1000–8000) until V1 (TD≈TI within 10%) and V2 (converged + norm decayed) pass. If TD is off from TI: check the (1/i) factor, the c-product (no conjugate), `E_tot=E+ε_{v_init}` in the exponential, the Simpson weights, and that the propagation runs long enough for ‖ψ‖ to decay (else the energy transform is truncated). If it won't converge to TI, STOP and report BLOCKED with the σ_TD vs σ_TI values and the norm-decay profile. **Commit** once V1/V2 pass.

---

### Task 3: Promote CN to `qscat.evolution` + wire Group D + docs

**Files:** Create `libs/qscat/qscat/evolution/crank_nicolson.py` + update `libs/qscat/qscat/evolution/__init__.py`; Test `libs/qscat/tests/test_crank_nicolson.py`; Modify `validation/n2/experiment.py` (+ maybe a small `validation/n2/td_check.py`), `validation/n2/reference.py` if needed; `docs/physics/`; `CLAUDE.md`.

- [ ] **Step 1: Promote the CN propagator to `qscat.evolution.make_cn_stepper`** (the generic primitive; mypy-strict clean). Have `projects/n2_td_cross_section/propagator.py` re-export or call it (keep the project copy thin, or import from qscat). Add `libs/qscat/tests/test_crank_nicolson.py` with the Task-1 propagator tests (exact-exp match, unitarity, decay), importing from `qscat.evolution`. `uv run mypy libs/qscat` → 0; `uv run pytest libs/qscat/tests/test_crank_nicolson.py -q` → pass.
- [ ] **Step 2: Wire harness Group D.** In `validation/n2/`, compute σ_TD at the gated anchors (reuse the vres already computed for C5 if practical, or accept the extra cost), and change Group D (`D1`) in `experiment.py` from PENDING to a real PASS/FAIL: PASS if σ_TD agrees with σ_TI within the V1 tolerance AND lands within the Houfek factor-3 bound. Guard with try/except → labeled FAIL on error (like B1/C5). Keep exit 0 when it passes; keep it CPU/docker-runnable (mind the added propagation cost — if the harness gets slow, use a modest dt/n_steps that still meets the 10% tolerance, documented).
- [ ] **Step 3: Docs.** `docs/physics/` — the TD method (CN, energy transform, TD=TI relation), the TD≈TI agreement, and that Group D is now green. `CLAUDE.md` — add `qscat.evolution` (Crank-Nicolson propagator) and the TD solver.
- [ ] **Step 4: Verify + commit.**
  - `uv run pytest projects/n2_td_cross_section libs/qscat validation/n2 -q` → all pass.
  - `uv run mypy libs/qscat` → 0.
  - `uv run python -m validation.n2.experiment` → Group D now PASS; ideally **0 PENDING** now, exit 0.
  - `docker/run-n2.sh` → completes, exit 0.

---

## Final verification

- [ ] TD σ matches TI σ at the gated anchors within 10%; converged in dt/n_steps; ‖ψ‖ decays.
- [ ] Harness Group D flips PENDING → PASS; `python -m validation.n2.experiment` exit 0 (0 FAIL); docker exit 0.
- [ ] CN propagator promoted to `qscat.evolution`, mypy-clean; `CLAUDE.md`/docs updated.
- [ ] No Rust; packages are properly package-qualified.
