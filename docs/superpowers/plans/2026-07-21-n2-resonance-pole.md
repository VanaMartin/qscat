# N₂ Electronic Resonance Pole Solver Implementation Plan (sub-project #2)

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development. Checkbox steps. Research-grade numerics — the pole must LAND in the physical window; tune resolution/window, never fake it.

**Goal:** Compute E_res(R), Γ(R), V_d(R) for the electron–N₂ ²Π_g resonance via two-angle ECS matching on the `qscat.dvr` grid, validated at R₀ (≈2.44 eV / 0.46 eV) and wired into the N₂ harness B1.

**Architecture:** `projects/n2_resonance/`: `potential.py` (N₂ electronic potentials, params from `validation/n2/config.json`), `grid_n2.py` (hand-built `FemDvrEcsGrid` factory), `pole.py` (Hamiltonian + two-angle pole finder + R-scan). Promote the general matcher to `qscat.ecs`.

**Tech Stack:** Python 3.12, numpy; `qscat.dvr` (FemDvrEcsGrid, kinetic, hamiltonian, eigen), `qscat.ecs`.

## Global Constraints

- Python `>=3.12`; `uv run pytest`. Atomic units; electron mass = 1; energy Hartree, length bohr.
- N₂ params (from `validation/n2/config.json`, verified in `validation/n2/model.py`):
  `l = impulsemomentum = 2`, `alpha_c = 0.4`, `lambda_inf=6.21066, lambda_1=1.05708, R_lambda=−27.9833, lambda_c=5.38022, R_c=2.405`, `R_0=2.01943, D_0=0.75102, alpha_0=1.1535`. `V_int(r,R)=−λ(R)exp(−α_c r²)`, `V_eff_el = l(l+1)/(2r²) + V_int`. λ(R) formula per `validation/n2/model.py`.
- **Expected pole at R₀** (port-scout FD prototype, sanity target — NOT the test oracle): `E_pole(R₀) ≈ 0.0898 − 0.00844i Ha` → E_res≈2.44 eV, Γ≈0.46 eV. `HARTREE_TO_EV = 27.211386245988`.
- `projects/n2_resonance/` has NO `__init__.py`; bare-name sibling imports. `qscat.dvr`/`qscat.ecs` import normally.
- Spec: `docs/superpowers/specs/2026-07-21-n2-resonance-pole-design.md`. Extraction: `.superpowers/sdd/n2-lcp-model-extraction.md`.
- Commit trailer: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.

---

### Task 1: N₂ electronic potential + hand-built ECS grid

**Files:** Create `projects/n2_resonance/potential.py`, `projects/n2_resonance/grid_n2.py`; Test `projects/n2_resonance/test_potential.py`, `test_grid_n2.py`.

**Interfaces:**
- Produces: `potential.v0(R)`, `potential.lam(R)`, `potential.v_int(r,R)`, `potential.v_eff_el(r,R)` (params loaded from `validation/n2/config.json` via an absolute path relative to repo root, or a copied constant block — match `validation/n2/model.py` exactly); `grid_n2.n2_electronic_grid(angle_deg, *, r_pivot=10.0, n_real=8, r_max=30.0, n_complex=8, quadrature=8) -> FemDvrEcsGrid`.

- [ ] **Step 1: Write failing potential tests — `test_potential.py`.** Assert `v_eff_el`/`v_int`/`lam`/`v0` match `validation/n2/model.py` to 1e-12 at sample (r,R) (import the reference: add the repo `validation/n2` dir to sys.path in the test, or recompute from the same formula and cross-check the CODATA constant). Concretely: `lam(R_c)==lambda_c`; `v0(R_0)==-D_0`; `v_int(1.0, R_0) < 0`; `v_eff_el(2.0,R_0) == v_int(2.0,R_0)+2*3/(2*2**2)`.

- [ ] **Step 2: Run → fail.** `uv run pytest projects/n2_resonance/test_potential.py -q`.

- [ ] **Step 3: Implement `potential.py`** — load params from `validation/n2/config.json` (path via `Path(__file__).resolve().parents[2] / "validation/n2/config.json"`), implement the four functions identical to `validation/n2/model.py` (v0 Morse, lam sigmoid, v_int Gaussian well, v_eff_el with l=2 centrifugal). Vectorized (numpy).

- [ ] **Step 4: Run → pass.**

- [ ] **Step 5: Write failing grid test — `test_grid_n2.py`.** `g = n2_electronic_grid(35.0)`; assert `g.R0 == r_pivot`; real-region points real, tail points rotated 35°; `g.n > 0`. Assert two grids at 35° and 44° have the SAME real_points (only the tail differs) — required for matching.

- [ ] **Step 6: Run → fail.**

- [ ] **Step 7: Implement `grid_n2.py`.** Build a `GridSpec`: `n_real` real elements tiling `[0, r_pivot]` (uniform length `r_pivot/n_real`), then `n_complex` complex elements tiling `[r_pivot, r_max]` (length `(r_max−r_pivot)/n_complex`, `angle_deg`), `quadrature=8`, `x_min=0`. Return `FemDvrEcsGrid(spec)`. (Uniform elements are fine for a converged toy; the well is at small r and the ECS tail captures the outgoing part.)

- [ ] **Step 8: Run → pass. Commit.**
```bash
git add projects/n2_resonance && git commit -m "feat(n2res): N2 electronic potential + hand-built ECS grid factory

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: Two-angle pole finder + V1/V2 (THE crux)

**Files:** Create `projects/n2_resonance/pole.py`; Test `projects/n2_resonance/test_pole.py`.

**Interfaces:**
- Consumes: Task 1; `qscat.dvr` (`kinetic`/`hamiltonian`/`eigen` or build H directly), `FemDvrEcsGrid`.
- Produces: `pole.electronic_hamiltonian(R, grid) -> ndarray`; `pole.find_pole(R, grid_a, grid_b, window) -> (E_pole: complex, residual: float)` where `window` is `(re_lo, re_hi, im_lo, im_hi)` bounding the search box; `resonance_curve` deferred to Task 3.

- [ ] **Step 1: Write the V1/V2 tests — `test_pole.py`**

```python
import numpy as np
from potential import v0
import pole
from grid_n2 import n2_electronic_grid

HARTREE_TO_EV = 27.211386245988
R0 = 2.01943
WINDOW = (0.04, 0.16, -0.05, 0.0)   # Re in [0.04,0.16] Ha, Im in [-0.05,0]


def test_V1_resonance_at_equilibrium():
    ga, gb = n2_electronic_grid(35.0), n2_electronic_grid(44.0)
    E, resid = pole.find_pole(R0, ga, gb, WINDOW)
    Eres_eV = E.real * HARTREE_TO_EV
    Gamma_eV = max(0.0, -2 * E.imag) * HARTREE_TO_EV
    assert 2.3 <= Eres_eV <= 2.5, Eres_eV
    assert 0.35 <= Gamma_eV <= 0.55, Gamma_eV


def test_V2_pole_is_stable():
    ga, gb = n2_electronic_grid(35.0), n2_electronic_grid(44.0)
    E, resid = pole.find_pole(R0, ga, gb, WINDOW)
    assert resid < 1e-3, resid                       # angle-stable (residual << Gamma)
    # resolution stability: coarser grid gives ~same pole (few %)
    ga2, gb2 = n2_electronic_grid(35.0, n_real=6, n_complex=6), n2_electronic_grid(44.0, n_real=6, n_complex=6)
    E2, _ = pole.find_pole(R0, ga2, gb2, WINDOW)
    assert abs(E - E2) / abs(E) < 0.05, (E, E2)
```

- [ ] **Step 2: Run → fail** (`ModuleNotFoundError: pole`).

- [ ] **Step 3: Implement `pole.py`.**
  - `electronic_hamiltonian(R, grid)`: `H = kinetic(grid, 1.0) + np.diag(v_eff_el(grid.points, R))` (use `qscat.dvr.kinetic`; `grid.points` are complex on the tail — `v_eff_el` must evaluate at complex r, which it does since it's analytic). Equivalent: `hamiltonian(grid, lambda z: v_eff_el(z, R), 1.0)`.
  - `find_pole(R, grid_a, grid_b, window)`: `Ea, _ = eigen(H_a)`, `Eb, _ = eigen(H_b)`; filter each to eigenvalues inside `window` (re/im box); for each `ea` in the filtered A-set, find nearest `eb`; the pole = the (ea, eb) pair with the SMALLEST `|ea − eb|` (the θ-stable one). Return `E_pole = 0.5*(ea+eb)`, `residual = |ea − eb|`. If the window is empty, raise a clear error (window too tight / grid too coarse).

- [ ] **Step 4: Run V1/V2.** `uv run pytest projects/n2_resonance/test_pole.py -q`.
  **If the pole is not in the window:** debug in this order — (a) print all eigenvalues near the window to see where the pole actually is; (b) widen `r_max`/increase `n_complex` (the outgoing tail may be under-resolved); (c) increase `quadrature`/`n_real` (the well region); (d) confirm the two angles give a matching pair (the continuum won't match, the pole will). Do NOT loosen the physical window [2.3,2.5]eV / [0.35,0.55]eV — that window IS the physics. If after genuine effort the pole won't land in the physical window, STOP and report BLOCKED with the eigenvalues you DO find near 0.09 Ha, so the controller can help.

- [ ] **Step 5: Commit** (once V1/V2 pass).
```bash
git add projects/n2_resonance && git commit -m "feat(n2res): two-angle ECS pole finder; N2 resonance at R0 ~2.4 eV validated

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: R-scan curves (V3) + promote matcher to qscat.ecs + wire harness B1 + docs

**Files:** Modify `projects/n2_resonance/pole.py` (add `resonance_curve`); Create `libs/qscat/qscat/ecs/` addition `find_resonance_pole`; Test `projects/n2_resonance/test_curve.py`; add `libs/qscat/tests/test_find_resonance_pole.py`; Modify `validation/n2/experiment.py` (+ maybe a small `validation/n2/resonance.py` helper); Modify `CLAUDE.md`, `docs/physics/`.

**Interfaces:**
- Produces: `pole.resonance_curve(R_grid, grid_a, grid_b) -> (E_res, Gamma, V_d)` (arrays); `qscat.ecs.find_resonance_pole(eigs_a, eigs_b, window) -> (E_pole, residual)` (the general matcher extracted from `find_pole`).

- [ ] **Step 1: Extract the general matcher to `qscat.ecs`.** Move the eigenvalue-matching core of `find_pole` into `qscat.ecs.find_resonance_pole(eigs_a, eigs_b, window)` (pure: takes two eigenvalue arrays + window, returns matched pole + residual). Have `pole.find_pole` call it. mypy-strict clean. Add `libs/qscat/tests/test_find_resonance_pole.py` with a synthetic case (two eigenvalue sets sharing one "pole" value + differing "continuum" values → returns the shared one). Run `uv run mypy libs/qscat` → 0 errors, `uv run pytest libs/qscat/tests/test_find_resonance_pole.py -q` → pass.

- [ ] **Step 2: Implement `resonance_curve` + V3 test — `test_curve.py`.** `resonance_curve(R_grid, grid_a, grid_b)`: for R in R_grid (e.g. `np.linspace(1.6, 3.0, 15)`), call `find_pole` seeded with a `window` recentered on the previous R's pole (continuation); collect E_res, Γ, V_d. V3 asserts: Γ ≥ 0 everywhere; E_res(R) and Γ(R) are smooth (no jump > a few× the local step — e.g. `max |ΔΔE_res| ` small, or monotone-ish through R₀). Run → pass. Debug mode-hops via the continuation window width, not by loosening physics.

- [ ] **Step 3: Wire the N₂ harness B1.** In `validation/n2/`, add a small helper that builds the N₂ electronic grids + calls `qscat.ecs.find_resonance_pole` (or imports the project solver) to compute `E_res(R₀)`, and change the B1 check in `experiment.py` from PENDING to a real PASS/FAIL against the literature window (`reference.LITERATURE["E_res_eV"]`). Keep it CPU/docker-runnable. Note: if importing from `projects/` into `validation/` is awkward, put the minimal grid+potential needed for B1 in `validation/n2/` (it already has the potential in `model.py`; add a tiny `resonance.py`). Run `uv run python validation/n2/experiment.py` → B1 now PASS, exit 0.

- [ ] **Step 4: Docs.** Add an N₂-resonance section to `docs/physics/` (method: two-angle matching, the R₀ result, the simplifications) and update `CLAUDE.md` for `qscat.ecs.find_resonance_pole`.

- [ ] **Step 5: Verify + commit.**
```bash
uv run pytest projects/n2_resonance libs/qscat -q     # all green
uv run python validation/n2/experiment.py             # B1 PASS, exit 0
git add projects/n2_resonance libs/qscat validation/n2 docs CLAUDE.md
git commit -m "feat(n2res): R-scan V_d(R)/Gamma(R); promote find_resonance_pole; N2 harness B1 green

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Final verification

- [ ] `uv run pytest projects/n2_resonance libs/qscat -q` green (V1–V3 + qscat suites).
- [ ] `E_res(R₀)` ∈ [2.3,2.5] eV, `Γ(R₀)` ∈ [0.35,0.55] eV; matching residual < 1e-3 Ha; resolution-stable; curves smooth.
- [ ] N₂ harness B1 flips PENDING → PASS; `validation/n2/experiment.py` exits 0; docker run still exits 0.
- [ ] `qscat.ecs.find_resonance_pole` promoted, mypy-strict clean; `CLAUDE.md`/docs updated.
- [ ] No Rust; no `__init__.py` under `projects/n2_resonance/`.
