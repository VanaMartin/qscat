# N₂ Time-Independent VE Cross-Section Implementation Plan (sub-project #3)

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development. Research-grade numerics — internal checks are the correctness gate; the Houfek anchors are a loose cross-MODEL comparison (LCP 1D vs Houfek 2D), reported not force-matched. If a computation won't converge to sane physics, report BLOCKED with what you see.

**Goal:** Compute N₂ VE cross sections σ_{0→v'}(E) via the time-independent LCP resolvent method on the nuclear FEM-DVR-ECS grid, validate internally (vib spacing, σ≥0, thresholds, resonance), compare the 6 C5 anchors to Houfek's data, and flip the harness C5 checks at a documented tolerance.

**Architecture:** `projects/n2_ti_cross_section/`: `nuclear_grid.py`, `vibrational.py`, `vres.py` (recompute V_d/Γ per nuclear-R via the #2 pole finder), `cross_section.py` (driven-equation solve). Reuses `qscat.dvr`, `qscat.ecs`, and `projects/n2_resonance/`.

**Tech Stack:** Python 3.12, numpy; `qscat.dvr` (FemDvrEcsGrid, kinetic, hamiltonian, eigen).

## Global Constraints

- Python `>=3.12`; `uv run pytest`. Atomic units (Hartree, bohr, bohr² cross sections). μ = 12766.36.
- **Formulas (port-scout-verified, `.superpowers/sdd/ti-cross-section-extraction.md` — READ IT):**
  - χ_v, ε_v = eigenpairs of `T_nuc(μ)+diag(V₀(R))` (bound states are real, θ-independent on the ECS grid).
  - doorway `d_v(R)=√(Γ(R)/2π)·χ_v(R)`; driven eqn `[E_tot−T_nuc(μ)−diag(V_d)+i·diag(Γ)/2]·ξ=d_v`, `E_tot=E+ε_v`.
  - `S_{v'←v}=⟨d_{v'}|ξ⟩` (DVR dot product; see the c-product note below); `σ_{v→v'}(E)=4π³·|S|²/(2E)`; open iff `E_tot−ε_{v'}>0` else 0.
  - `V_d(R)=V₀(R)+E_res(R)`, from the #2 pole finder (recompute per nuclear-R).
- **c-product note:** ECS matrices are complex-symmetric non-Hermitian; the correct inner product is the c-product (`Σ a_j b_j`, NO conjugate) for resolvent/doorway quantities, while the qscat `eigen` returns `v†v=1`-normalized vectors. Get this right — a wrong conjugate/normalization gives wrong (often complex or negative) σ. The internal checks (σ real & ≥0) will catch it.
- N₂ ω_e ≈ 2358 cm⁻¹ = **0.01074 Ha** (0.2924 eV) — the vibrational-spacing target.
  - **Note (superseded):** this real-N₂ target was superseded by a maintainer
    decision to accept eMoScat's model potential as-is (see
    `.superpowers/sdd/task1fix-report.md`). The shipped Part-1 check (`Step 5`
    below, `test_vibrational.py`) instead validates the FEM-DVR solver against
    the **analytic Morse spectrum of eMoScat's own potential** (matched to
    ~1e-14), not against this real-N₂ target. eMoScat's model gives
    `ε₁−ε₀ ≈ 0.0124` Ha — ~16% above real N₂'s `ω_e` — a documented property
    of the model potential, not a solver error; "ε₁−ε₀ ≈ N₂ ω_e" as written
    in this plan is not the criterion actually met.
- HARTREE_TO_EV = 27.211386245988. Houfek data: `validation/n2/data/CSVE.V00.J00` (energy Ha col 0; col j = v=0→(j−1); bohr²).
- No `__init__.py` under `projects/n2_ti_cross_section/`. Commit trailer: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.

---

### Task 1: Nuclear grid + neutral vibrational states

**Files:** Create `projects/n2_ti_cross_section/nuclear_grid.py`, `vibrational.py`; Test `test_vibrational.py`.

**Interfaces:** `nuclear_grid.n2_nuclear_grid(*, quadrature=14) -> FemDvrEcsGrid` (real 0→12 bohr per N2.json + 35° complex tail); `vibrational.vibrational_states(grid, mu, n) -> (eps: real ndarray (n,), chi: ndarray (n, nb))` — the n lowest-real-energy eigenpairs of `T_nuc(μ)+diag(V₀(grid.points))`, ordered by ε.

- [ ] **Step 1: Read** `.superpowers/sdd/ti-cross-section-extraction.md` §2,§7.
- [ ] **Step 2: Failing test — `test_vibrational.py`.** Build the grid; `eps, chi = vibrational_states(grid, 12766.36, 6)`. Assert: `eps` real (imag < 1e-6), ascending; **ε₁−ε₀ ≈ 0.01074 Ha within 5%** (N₂ vibrational quantum); ε roughly evenly spaced (anharmonicity small for low v). V₀ from `projects/n2_resonance/potential.v0` (import via sys.path, like earlier project cross-imports).
- [ ] **Step 3: Run → fail.**
- [ ] **Step 4: Implement.** `n2_nuclear_grid`: hand-built GridSpec, real elements tiling [0,12] (use N2.json's `lengths/points` [0.5,0.15,0.5,1.0]→[1.5,3.0,4.0,12.0] as a guide, or uniform enough to converge the low vib states), 35° tail to ~large R (e.g. r_max≈40, ~10 complex elements), quadrature 14. `vibrational_states`: `T = kinetic(grid, mu)`; `H0 = T + np.diag(v0(grid.points))`; `E, V = eigen(H0)`; select the n lowest eigenvalues with |Im| small (bound states), return (Re(eps), corresponding eigenvectors). Document the c-product/normalization convention chosen.
- [ ] **Step 5: Run → pass** (ε₁−ε₀ ≈ N₂ ω_e). If the spacing is wrong, the grid/mass/V₀ is wrong — debug there. **Commit.**

---

### Task 2: V_d(R)/Γ(R) recomputed per nuclear-R

**Files:** Create `projects/n2_ti_cross_section/vres.py`; Test `test_vres.py`.

**Interfaces:** `vres.vres_on_grid(grid) -> (Vd: ndarray (nb,), Gamma: ndarray (nb,))` — for each nuclear grid point R (`grid.points`), V_d(R)=V₀(R)+E_res(R), Γ(R), via the sub-project #2 two-angle pole finder.

- [ ] **Step 1: Failing test — `test_vres.py`.** `Vd, Gamma = vres_on_grid(grid)`; assert shapes = (nb,); at a REAL grid point near R₀=2.019, `E_res = Vd − v0(R) ≈ 0.0898 Ha` (the #2 result) and `Γ ≈ 0.0167 Ha` (2·0.00836) within a few %; `Γ ≥ 0` (clamp tiny negatives); Γ small/→0 at large R (beyond the ~2.4 bohr crossing).
- [ ] **Step 2: Run → fail.**
- [ ] **Step 3: Implement.** For each nuclear point R in `grid.points`: build the two electronic grids (`projects/n2_resonance/grid_n2.n2_electronic_grid(35)/(44)`), call `projects/n2_resonance/pole.find_pole(R, ga, gb, window)` with the window seeded near the last-found pole (continuation over the sorted real R points). `E_res=Re(E_pole)`, `Γ=max(0,−2 Im(E_pole))`, `Vd=v0(R)+E_res`. **Complex-tail R:** the electronic potential is analytic in R, so `find_pole` should work at complex R (Γ→0 there); if direct complex-R pole-finding is unstable, fall back to fitting V_d(R)/Γ(R) as smooth functions of real R and analytically continuing to the tail points (document which path you used). Reuse the electronic grids across R (only V_eff_el's R changes) for speed.
- [ ] **Step 4: Run → pass. Commit.** (If complex-tail pole-finding forces a design change, report DONE_WITH_CONCERNS describing it.)

---

### Task 3 (THE CRUX): VE cross section + internal checks + Houfek comparison

**Files:** Create `projects/n2_ti_cross_section/cross_section.py`; Test `test_cross_section.py`.

**Interfaces:** `cross_section.ve_cross_section(grid, mu, Vd, Gamma, eps, chi, v_init, vprimes, E) -> ndarray` — σ_{v_init→v'}(E) (bohr²) for the listed final channels `vprimes` at collision energy(ies) E. Builds `H_res=kinetic(grid,mu)+diag(Vd−i·Gamma/2)`, solves `(E_tot−H_res)ξ=d_{v_init}`, forms S and σ.

- [ ] **Step 1: Failing test — `test_cross_section.py`** implementing the internal Part-1 checks + the anchor comparison:
  - **Internal:** σ real & ≥ 0 at several (E, v'); a closed channel (E_tot−ε_{v'}<0) gives σ=0; σ_{0→1}(E) is larger in the ~2–3 eV region than near threshold (resonance enhancement).
  - **Anchors (report, loose gate):** compute σ_{0→v'}(E) at the 6 C5 coordinates (from `validation/n2/reference.ANCHOR_COORDS`), read Houfek's values from `CSVE.V00.J00`, and assert each is within a FACTOR of ~3 (documented loose bound; this is a cross-model check, tighten/loosen after seeing the numbers). PRINT the per-anchor ratio so the actual agreement is visible.
- [ ] **Step 2: Run → fail.**
- [ ] **Step 3: Implement `cross_section.py`** per the formulas. Doorways `d_v = sqrt(Gamma/(2*pi))*chi[v]`. For each E: `E_tot=E+eps[v_init]`; `M=(E_tot)*I − H_res`; `xi=np.linalg.solve(M, d_{v_init})`; `S=cprod(d_{v'}, xi)` using the correct DVR inner product (c-product, no conjugate, matching the vibrational-state normalization); `sigma=4*pi**3*abs(S)**2/(2*E)`; 0 if `E_tot−eps[v']<=0`. Reuse the LU/`M` across v' at fixed E.
- [ ] **Step 4: Run.** Debug order if σ is wrong: (a) confirm σ is REAL & positive — if complex/negative, the inner-product conjugation or eigenvector normalization is wrong (try c-product vs Hermitian, and vᵀv vs v†v normalization of χ); (b) check S has the right order of magnitude (print |S|); (c) confirm thresholds/units. The internal checks MUST pass (they're model-independent). For the anchors: if within a factor of ~3 and the trend/ordering across channels matches Houfek, that's success for a cross-model check — do NOT force exact agreement. If σ is qualitatively wrong (negative, or orders off, or no resonance), STOP and report BLOCKED with the |S|/σ values and what you tried.
- [ ] **Step 5:** Once internal checks pass and anchors are within the documented factor, **write the actual per-anchor agreement into the report** and **commit.**

---

### Task 4: Wire harness C5 + docs

**Files:** Modify `validation/n2/experiment.py` (+ maybe `validation/n2/cross_section.py` helper), `validation/n2/reference.py` (RTOL for C5 → documented cross-model bound); `docs/physics/`; `CLAUDE.md`. Optionally promote a generic `resolvent_cross_section` to `qscat`.

- [ ] **Step 1:** Set the C5 tolerance in `reference.py` (or the check) to the documented cross-model bound found in Task 3 (e.g. `C5_FACTOR = 3.0` with a comment: LCP 1D vs Houfek 2D). 
- [ ] **Step 2:** Add `validation/n2/cross_section.py` computing σ at the anchors (build grid + vibrational states + vres + ve_cross_section — import from `projects/n2_ti_cross_section` or replicate minimally per the existing lockstep pattern). Change C5 in `experiment.py` from PENDING to a real PASS/FAIL: each anchor passes if `1/factor ≤ σ_computed/σ_houfek ≤ factor`; detail string reports the ratio and notes "LCP vs Houfek 2D". Guard with try/except → labeled FAIL on solver error (matching the B1 pattern). Keep CPU/docker-runnable — NOTE: this adds real compute (many pole solves + linear solves) at harness time; if it's slow (>~30s), cache the anchor σ or reduce the vres recompute cost, but it must complete in the docker run.
- [ ] **Step 3:** `docs/physics/` section (TI-LCP method, the resolvent formula, the internal checks, the LCP-vs-2D caveat + actual agreement) and `CLAUDE.md` update.
- [ ] **Step 4: Verify + commit.**
  - `uv run pytest projects/n2_ti_cross_section libs/qscat validation/n2 -q` → green.
  - `uv run python validation/n2/experiment.py` → C5 anchors now PASS (or documented FAIL if a channel is off), exit reflects it.
  - `docker run --rm qmodeling:runtime python validation/n2/experiment.py` → completes, exit 0 if all pass.

---

## Final verification

- [ ] Internal checks pass (vib spacing ≈ N₂ ω_e; σ≥0; thresholds; resonance in 2–3 eV).
- [ ] The 6 C5 anchors compared to Houfek; per-anchor agreement documented; harness C5 flips PENDING→PASS at the documented cross-model tolerance.
- [ ] `validation/n2/experiment.py` runs (locally + docker); no `__init__.py` under the project; no Rust.
