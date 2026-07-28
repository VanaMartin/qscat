# Automatic discretisation tuner Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build `qscat.tuning` — deterministic primitives that analyze a potential and produce the minimal-cost FEM-DVR-ECS grid at a target precision (adaptive equidistribution mesh + h/p quadrature + a joint ECS-tail choice) — plus a `discretisation-tuner` skill that supervises the tuning loop, calibrated and gated against the eMoScat decks and the coarse-grid failures that bit us.

**Architecture:** `qscat.tuning` is a new library submodule of pure, unit-tested primitives: `analyze` (potential → local-wavenumber profile), `mesh` (equidistribution element layout + h/p sweep), `ecs` (the joint angle/exp-increment/tail-quadrature tail tuner, ≤~35° for double-ECS stability), `probes` (decoupled 1-D convergence at the energy extremes), `metrics`+`propose` (cost model + a-priori grid), `incident` (incident-state/test-function inputs + a TW auto-tune extension). The `discretisation-tuner` skill runs the analyze→propose→probe→refine→emit loop, making the judgment calls. See `docs/superpowers/specs/2026-07-28-discretisation-tuner-design.md`.

**Tech Stack:** Python ≥3.12, NumPy/SciPy, `qscat.dvr` (`ElementSpec`/`GridSpec`/`FemDvrEcsGrid`/`kinetic`/`eigen`/`TensorGrid`), `qscat.special` (`riccati_bessel_en`/`riccati_bessel_en_mass`/`coulomb_f_en`), `qscat.ecs.find_resonance_pole`, `qscat.core` (`vibrational_states`, `anion_electronic_states`, the cross sections for the 2-D spot-check), `qscat.model`. pytest, mypy --strict over `libs/qscat/qscat`, ruff.

## Global Constraints

- **Atomic units.** Every knob is derived from physics, never a magic number: the equidistribution phase-per-element constant `C` and the ECS safety fractions are **calibrated** (Task 8), not assumed; until calibrated they carry a documented provisional value.
- **Pure primitives.** `qscat.tuning.{analyze,mesh,ecs,metrics}` take plain callables/arrays (a potential `V(x)`, mass, energy range) — NOT models — so they are unit-testable on analytic potentials and reusable. Model-specific extraction (which V/mass/channels per coordinate) lives in `probes`/`propose` thin adapters.
- **`qscat.tuning` may import `qscat.model`** (it is tooling, not the model-independent core engine) — but keep the pure primitives model-free; only `probes`/`propose`/`incident` touch models. `qscat.core`'s no-model-import boundary is unaffected (tuning is a separate submodule).
- **Adaptive mesh = explicit `ElementSpec` length list.** The mesh generator emits a `list[ElementSpec]`; grids are built via `GridSpec(quadrature, elements, x_min)` → `FemDvrEcsGrid`. `GridSpec.R0 = x_min + Σ(real element lengths)`.
- **ECS tail** is the exp-decay regime, NOT equidistribution: exp-growth element lengths, angle `θ = min(per-coordinate stability limit, ~35° double-ECS cap)`, tail quadrature possibly ≠ real-region quadrature.
- **Convergence** = a quantity stable under ONE refinement step to `rtol` (default 1e-3), the repo convention. The tuner minimizes DVR points subject to that.
- mypy --strict over `libs/qscat/qscat` clean; ruff clean; the `l`-as-angular-momentum E741 per-file-ignore convention applies.

## File Structure

- `libs/qscat/qscat/tuning/__init__.py` (create) — public API re-exports.
- `libs/qscat/qscat/tuning/analyze.py`, `mesh.py`, `ecs.py`, `probes.py`, `metrics.py`, `propose.py`, `incident.py` (create) — the primitives.
- `libs/qscat/tests/test_tuning_*.py` (create) — per-module unit tests.
- `validation/tuning/` (create) — `test_emoscat_decks.py` (the reproduce-or-beat + flag-failures gate).
- `.claude/skills/discretisation-tuner/SKILL.md` (create) — the supervising skill.
- `docs/physics/discretisation-tuning.md` (create), `CLAUDE.md` (modify).

---

### Task 1: Potential analysis (`qscat.tuning.analyze`)

**Files:**
- Create: `libs/qscat/qscat/tuning/__init__.py`, `libs/qscat/qscat/tuning/analyze.py`
- Test: `libs/qscat/tests/test_tuning_analyze.py`

**Interfaces:**
- Produces:
  ```python
  @dataclass(frozen=True)
  class PotentialProfile:
      x: NDArray[float64]           # dense sample points (real)
      V: NDArray[float64]           # Re(V(x)) on x
      k: NDArray[float64]           # local wavenumber sqrt(2 m max(E_max - V, 0))
      kappa: NDArray[float64]       # forbidden-region decay sqrt(2 m max(V - E_max, 0))
      turning_points: NDArray[float64]
      singularities: NDArray[float64]   # x where |V| -> inf (e.g. the -1/r origin), from a growth test

  def analyze_potential(V: Callable[[NDArray], NDArray], x_min: float, x_max: float, m: float,
                        e_max: float, *, n_sample: int = 4000) -> PotentialProfile: ...
  ```
  `V` is evaluated on `n_sample` real points in `[x_min, x_max]` (its real part taken — potentials may return complex but the analysis uses `Re`). `k(x) = sqrt(2 m max(e_max - V, 0))`; `kappa` the complementary decay rate; turning points = sign changes of `e_max - V`; singularities = points where `|V|` exceeds a large threshold / grows without bound as `x -> x_min` (the `-1/r`).

**Design notes:** pure, no models. Sampling near a `1/r` singularity: start `x_min` a hair above 0 if `V(x_min)` is non-finite, and record the singularity. Turning points via `np.diff(np.sign(e_max - V))`. The profile is the sole input the mesh/ECS generators need.

- [ ] **Step 1: Write the failing test**

```python
# libs/qscat/tests/test_tuning_analyze.py
from __future__ import annotations

import numpy as np
from qscat.tuning import analyze_potential


def test_harmonic_turning_points_and_k():
    # V = 1/2 m w^2 x^2, m=1, w=1, E_max=2 -> turning points at x=+-2, k(0)=sqrt(2*2)=2
    m, w, E = 1.0, 1.0, 2.0
    V = lambda x: 0.5 * m * w**2 * np.asarray(x) ** 2
    p = analyze_potential(V, -5.0, 5.0, m, E)
    assert np.isclose(abs(p.turning_points).max(), 2.0, atol=1e-2)
    assert np.isclose(p.k[np.argmin(np.abs(p.x))], np.sqrt(2 * m * E), atol=1e-2)  # k at x=0
    # forbidden beyond the turning points: k=0, kappa>0
    assert p.k[np.argmax(p.x)] == 0.0 and p.kappa[np.argmax(p.x)] > 0.0


def test_detects_coulomb_singularity():
    V = lambda x: -1.0 / np.asarray(x)
    p = analyze_potential(V, 1e-3, 50.0, 1.0, 0.1)
    assert p.singularities.size >= 1 and p.singularities.min() < 0.05


def test_k_scales_with_mass_and_energy():
    V = lambda x: np.zeros_like(np.asarray(x, dtype=float))
    p = analyze_potential(V, 0.1, 10.0, 918.25, 0.05)   # heavy: large k
    assert np.allclose(p.k, np.sqrt(2 * 918.25 * 0.05), atol=1e-6)
```

- [ ] **Step 2: Run → fail** (`uv run pytest libs/qscat/tests/test_tuning_analyze.py -q` → ModuleNotFoundError).
- [ ] **Step 3: Implement `analyze.py`** per the Interfaces + Design notes; export `analyze_potential`, `PotentialProfile` from `tuning/__init__.py`.
- [ ] **Step 4: Run → pass.**
- [ ] **Step 5: mypy + ruff + commit** (`libs/qscat/qscat/tuning/`, the test) — `feat(tuning): potential analysis (local-wavenumber profile, turning points, singularities)`.

---

### Task 2: Equidistribution mesh + h/p sweep (`qscat.tuning.mesh`)

**Files:**
- Create: `libs/qscat/qscat/tuning/mesh.py`
- Modify: `libs/qscat/qscat/tuning/__init__.py`
- Test: `libs/qscat/tests/test_tuning_mesh.py`

**Interfaces:**
- Produces:
  ```python
  def equidistribution_elements(profile: PotentialProfile, order: int, *,
                                phase_per_element: float, min_len: float, max_len: float
                                ) -> list[float]: ...   # real-region element LENGTHS
  def optimal_real_mesh(profile: PotentialProfile, *, orders=(6, 8, 10, 14),
                        phase_coeff: float, min_len: float, max_len: float
                        ) -> tuple[list[float], int]: ...   # (element lengths, chosen order)
  ```
  `equidistribution_elements`: place element boundaries so `∫ₓₖ^{xₖ₊₁} k dx ≈ phase_per_element`; in classically-forbidden runs use a `kappa`-based length (resolve the decay), clamp to `[min_len, max_len]`, and refine near turning points/singularities (halve the local length). `optimal_real_mesh`: for each `order`, `phase_per_element = phase_coeff * (order - 1)`, build the mesh, count DVR points `≈ len(elements)*(order-1)`, return the `(mesh, order)` with the FEWEST points (the h/p optimum). `phase_coeff` is the calibrated `C` (Task 8; provisional value in a module constant `_PHASE_COEFF_PROVISIONAL`).

**Design notes:** the cumulative phase `Phi(x) = cumulative_trapezoid(k, x)`; boundaries at `Phi = j*phase_per_element` via `np.interp(targets, Phi, x)`. Element lengths = `np.diff(boundaries)`. Clamp + turning-point refinement as a post-pass. Return LENGTHS (not `ElementSpec`) — the caller wraps them + the ECS tail.

- [ ] **Step 1: Write the failing test**

```python
# libs/qscat/tests/test_tuning_mesh.py
from __future__ import annotations

import numpy as np
from qscat.tuning import analyze_potential, equidistribution_elements, optimal_real_mesh


def _flat_profile(k_const, x_max, m):
    E = k_const**2 / (2 * m)
    return analyze_potential(lambda x: np.zeros_like(np.asarray(x, float)), 0.0, x_max, m, E)


def test_uniform_k_gives_uniform_elements():
    # constant k -> equidistribution => (nearly) uniform element lengths
    p = _flat_profile(2.0, 20.0, 1.0)
    els = equidistribution_elements(p, 8, phase_per_element=5.0, min_len=1e-3, max_len=10.0)
    assert np.std(els) / np.mean(els) < 0.05                 # ~uniform
    assert abs(2.0 * np.mean(els) - 5.0) < 0.3               # k*h ~ phase_per_element


def test_finer_where_k_larger():
    # step: k large for x<5, small for x>5 -> smaller elements in the fast region
    def V(x):
        x = np.asarray(x, float); return np.where(x < 5.0, -2.0, 0.0)
    p = analyze_potential(V, 0.0, 20.0, 1.0, 0.5)
    els = equidistribution_elements(p, 8, phase_per_element=5.0, min_len=1e-3, max_len=10.0)
    # reconstruct boundaries; mean element length in [0,5] < mean in [5,20]
    b = np.concatenate([[0.0], np.cumsum(els)])
    lo = np.diff(b)[b[:-1] < 5.0].mean(); hi = np.diff(b)[b[:-1] >= 5.0].mean()
    assert lo < hi


def test_optimal_mesh_minimizes_points():
    p = _flat_profile(2.0, 20.0, 1.0)
    els, order = optimal_real_mesh(p, phase_coeff=1.0, min_len=1e-3, max_len=10.0)
    assert order in (6, 8, 10, 14)
    assert len(els) * (order - 1) <= 400                     # a sane point budget for this easy case
```

- [ ] **Step 2–5:** run→fail; implement `mesh.py` (export both); run→pass; mypy+ruff; commit `feat(tuning): equidistribution mesh + h/p quadrature sweep`.

---

### Task 3: ECS-tail tuner (`qscat.tuning.ecs`)

**Files:**
- Create: `libs/qscat/qscat/tuning/ecs.py`
- Modify: `libs/qscat/qscat/tuning/__init__.py`
- Test: `libs/qscat/tests/test_tuning_ecs.py`

**Interfaces:**
- Produces:
  ```python
  def max_stable_angle(V: Callable, R0: float, tail_extent: float, *,
                       angle_cap: float = 35.0, n_probe: int = 40) -> float: ...
  def tune_ecs_tail(K: float, R0: float, *, angle: float, order: int,
                    tail_alpha: float = 0.2, tail_skip: int = 2, decay_target: float = 1e-12
                    ) -> list[float]: ...   # tail element LENGTHS (exp-growth)
  ```
  `max_stable_angle`: increase θ from small to `angle_cap`; at each, evaluate `V` on the rotated contour `z = R0 + (x−R0)e^{iθ}` out to `R0+tail_extent`; reject θ where `|V|` GROWS along the tail (the analytic continuation diverges — e.g. a Gaussian past 45°). Return the largest non-diverging θ, capped at `angle_cap` (~35°, the double-ECS bound). `tune_ecs_tail`: exp-growth element lengths (`base·e^{α(i−skip+1)}`, `base` = the pivot real element length passed via `R0`-region context — accept `base` as a param) sized so the fastest wave `e^{−K·(x−R0)·sinθ}` reaches `decay_target` within the elements, using the existing `_ecs_tail`-style growth; return the length list.

**Design notes:** reuse `qscat.ecs.ecs_map` for the rotated contour. The tail length to reach `decay_target`: `L = −ln(decay_target)/(K sinθ)`; lay exp-growth elements spanning `[R0, R0+L]`. The number of elements follows from the growth + `L`. Tail quadrature: default = real-region order, but expose it so the skill can lower it for the decaying tail. Keep it a pure function of `(K, R0, angle, order, growth)`.

- [ ] **Step 1: Write the failing test**

```python
# libs/qscat/tests/test_tuning_ecs.py
from __future__ import annotations

import numpy as np
from qscat.tuning import max_stable_angle, tune_ecs_tail


def test_gaussian_angle_capped_below_45():
    # a Gaussian interaction exp(-alpha r^2) diverges on the rotated contour for theta>45
    V = lambda z: -np.exp(-0.4 * np.asarray(z) ** 2)
    ang = max_stable_angle(V, R0=12.0, tail_extent=40.0)
    assert ang <= 35.0 + 1e-9                    # never above the double-ECS cap
    # a bare -1/r (no Gaussian growth) is limited only by the cap
    ang2 = max_stable_angle(lambda z: -1.0 / np.asarray(z), R0=12.0, tail_extent=40.0)
    assert 30.0 <= ang2 <= 35.0 + 1e-9


def test_tail_absorbs_fast_wave():
    # K=58 (the F2 DA wave); at 35 deg the tail must reach ~1e-12 decay
    els = tune_ecs_tail(58.0, R0=10.7, angle=35.0, order=14, decay_target=1e-12)
    L = sum(els)
    assert np.exp(-58.0 * L * np.sin(np.deg2rad(35.0))) <= 1e-11   # decayed
    assert all(els[i] <= els[i + 1] + 1e-12 for i in range(len(els) - 1))  # exp-growth
```

- [ ] **Step 2–5:** run→fail; implement `ecs.py` (export both); run→pass; mypy+ruff; commit `feat(tuning): ECS-tail tuner (double-ECS-capped angle + exp-growth absorption)`.

---

### Task 4: 1-D convergence probes (`qscat.tuning.probes`)

**Files:**
- Create: `libs/qscat/qscat/tuning/probes.py`
- Modify: `libs/qscat/qscat/tuning/__init__.py`
- Test: `libs/qscat/tests/test_tuning_probes.py`

**Interfaces:**
- Produces (`ProbeResult = namedtuple("ProbeResult", "value converged cost detail")`):
  ```python
  def probe_nuclear(model, nuclear_grid, n_vib, *, rtol=1e-3) -> ProbeResult: ...
  def probe_electronic(model, elec_grid, R, *, window, rtol=1e-3) -> ProbeResult: ...
  def probe_channel_representation(grid, k, l, *, charge=0, mass=1.0, rtol=1e-3) -> ProbeResult: ...
  ```
  - `probe_nuclear`: `eps, _ = vibrational_states(nuclear_grid, model.mu, n_vib, model.v0)`; refine (bump quadrature by +2 OR halve the finest elements — a helper `refine(grid)`), recompute; `converged = max|Δeps|/|eps| < rtol`. `cost = nuclear_grid.n`.
  - `probe_electronic`: `E_pole = find_resonance_pole(...)` at fixed `R` on `elec_grid` and a refined grid; `converged` on the pole position+width. (Or `anion_electronic_states` for bound/Rydberg.)
  - `probe_channel_representation`: the channel function (`riccati_bessel_en`/`coulomb_f_en`) should be representable on the grid — project it onto the DVR basis and back, or check the quadrature integral of `|F|²` on the real region against the analytic energy-normalization; `converged = rel error < rtol`. This catches the "K≈58 unresolved" failure directly.

**Design notes:** provide a `refine(grid) -> grid` helper (one h- or p-refinement step). Probes are the empirical validators the skill reads. The channel-representation probe is the cheapest and most diagnostic (no eigensolve) — it is what would have caught the coarse-grid DA failure.

- [ ] **Steps (TDD):** test that on a KNOWN-GOOD grid the probes report `converged=True` (N₂ nuclear vibrational eps stable; the N₂ resonance pole stable; `riccati_bessel_en` at a modest k represented), and on a deliberately COARSE grid the channel-representation probe reports `converged=False` for a fast wave (K≈58 on a 1.0-bohr-element grid). Implement; mypy+ruff; commit `feat(tuning): decoupled 1-D convergence probes (nuclear/electronic/channel-representation)`.

---

### Task 5: Cost model + a-priori grid assembler (`qscat.tuning.metrics`, `propose`)

**Files:**
- Create: `libs/qscat/qscat/tuning/metrics.py`, `libs/qscat/qscat/tuning/propose.py`
- Modify: `libs/qscat/qscat/tuning/__init__.py`
- Test: `libs/qscat/tests/test_tuning_metrics.py`, `test_tuning_propose.py`

**Interfaces:**
- `metrics`: `grid_cost(grid) -> dict` (`n_points`); `tensor_cost(g_r, g_R) -> dict` (`n_unknowns = g_r.n*g_R.n`, `est_nnz`, `est_factor_gib`, `est_factor_seconds` — rough scalings anchored to the CLAUDE.md MUMPS/SuperLU measurements, e.g. the 143k-deck data points, clearly labelled ESTIMATES).
- `propose`: `propose_grid(model, coordinate, energy_range, *, rtol=1e-3, incident=None) -> FemDvrEcsGrid` — the ONE-SHOT a-priori grid: pick the per-coordinate `V`/mass/channel-`K` (a small model adapter: nuclear → `v0`, μ, `K=√(2μ E_DR_max)`; electronic → the effective electronic potential at the equilibrium R, mass 1, `k=√(2E_max)`, `charge=model.charge`), `analyze_potential` → `optimal_real_mesh` → `max_stable_angle`+`tune_ecs_tail` → assemble `ElementSpec` list → `FemDvrEcsGrid`. `incident` (Task 6) extends the real cutoff to contain the incident/test-function placement.

**Design notes:** `propose_grid` is the a-priori half of the hybrid; the skill runs the probe/refine loop on top. The model adapter is the only model-aware code; keep it small and explicit per coordinate. `tensor_cost`'s estimates are rough (anchored, labelled) — their job is relative ranking (is grid A cheaper than B), not absolute prediction.

- [ ] **Steps (TDD):** `grid_cost`/`tensor_cost` return the right point/unknown counts + monotone estimates; `propose_grid(N2, "nuclear", (0.04, 0.18))` builds a valid `FemDvrEcsGrid` with `R0` in a sane range and fewer-or-comparable points to the N₂ deck; `propose_grid(F2, "nuclear", (0.01, 0.05))` yields FINE outer elements (resolving the DA K≈58 wave — the fix the coarse grid missed). Implement; mypy+ruff; commit `feat(tuning): cost model + propose_grid (a-priori adaptive grid)`.

---

### Task 6: Incident/test-function inputs + TW analysis (`qscat.tuning.incident`)

**Files:**
- Create: `libs/qscat/qscat/tuning/incident.py`
- Modify: `libs/qscat/qscat/tuning/propose.py` (accept `incident`), `__init__.py`
- Test: `libs/qscat/tests/test_tuning_incident.py`

**Interfaces:**
- `IncidentSpec` (dataclass): for the TI route, the channel energies (from the energy range — a no-op beyond the channel-rep probe); for the TD route, the Gaussian wavepacket (`position`, `impulse`, `sigma`) + the observation boundary. `required_extent(incident) -> float` (how far the real region must reach to contain the wavepacket + observation).
- **BASELINE (this task):** `propose_grid(..., incident=IncidentSpec(...))` extends the real cutoff to `max(physics cutoff, required_extent(incident))` and ensures the resolution covers the incident's local wavenumber.
- **EXTENSION (deliver if tractable, else defer):** `tw_analysis(model, energy_range) -> IncidentSpec` — auto-place the TD Gaussian (position/impulse/width whose energy spectrum spans `[E_min,E_max]`, observation boundary) — PREPENDED to `propose_grid`. If the auto-tuning proves too complex, ship inputs-only and record the deferral.

**Design notes:** the wavepacket energy spectrum: a Gaussian `e^{−(r−r0)²/2σ²}e^{i p r}` has mean energy `p²/2` and an energy width set by `σ`; `tw_analysis` inverts `[E_min,E_max] → (p, σ, r0)` (r0 far enough that the packet starts outside the interaction and reaches the observation boundary within the propagation). Keep the baseline (inputs) solid; the auto-tune is best-effort.

- [ ] **Steps (TDD):** `required_extent` grows the cutoff for a far wavepacket; `propose_grid` with an `IncidentSpec` at position 45 extends the real region past 45; (extension) `tw_analysis(N2, (0.04,0.18))` returns a `p`/`sigma` whose spectrum brackets the range. Implement baseline (+ extension if clean); mypy+ruff; commit `feat(tuning): incident/test-function inputs (+ TW auto-tune extension)`.

---

### Task 7: The `discretisation-tuner` skill

**Files:**
- Create: `.claude/skills/discretisation-tuner/SKILL.md`
- Modify: `.claude/settings` skill listing if needed; `CLAUDE.md` skills table (Task 8).

**Interfaces:** the skill is a documented PROCEDURE (not code). It orchestrates: (1) `propose_grid` per coordinate for the energy range (+ incident spec); (2) run `probe_*` at `E_max` and near-threshold `E_min`, plus the channel-representation + ECS-absorption probes; (3) for any probe `converged=False`, refine that knob and re-probe; for any comfortably-over-converged knob, coarsen and re-probe; (4) once all probes hold `rtol` at minimal cost, run ONE 2-D spot-check (the actual cross section at the hardest energy) — on the proxy/Docker for huge decks; (5) emit the per-coordinate config (`ElementSpec` lists + quadrature + ECS) + a report (precision, cost vs the old grid, decisions).

**Design notes:** the skill encodes the loop + the judgment (which knob, when to stop, the cost/precision tradeoff, when to defer the 2-D check to Docker). It calls only `qscat.tuning` primitives. Include a worked example (tune the N₂ nuclear grid) and the STOP criteria. Follow the repo's SKILL.md conventions (frontmatter `name`/`description`, a checklist).

- [ ] **Steps:** write `SKILL.md` (procedure + worked example + stop criteria + output format); dry-run the procedure mentally against N₂; commit `feat(tuning): discretisation-tuner skill (supervised tuning loop)`.

---

### Task 8: Calibration + validation (the tuner's own gate)

**Files:**
- Create: `validation/tuning/__init__.py`, `validation/tuning/test_emoscat_decks.py`, `validation/tuning/calibrate.py`
- Modify: `libs/qscat/qscat/tuning/mesh.py` (set the calibrated `phase_coeff`), `docs/physics/discretisation-tuning.md`, `CLAUDE.md`.

**Design notes:** `calibrate.py` measures the phase-per-element constant `C` (and the ECS safety fractions) by requiring that `propose_grid` reproduce the KNOWN-GOOD eMoScat decks — sweep `C`, find the smallest that makes the probes converge on N₂/NO/F₂/H₂⁺, set `_PHASE_COEFF` in `mesh.py`. The gate `test_emoscat_decks.py`: (a) `propose_grid` for N₂/NO/F₂/(H₂⁺ proxy) reaches the same-or-better probe precision at same-or-fewer DVR points than the committed deck; (b) the tuner FLAGS the coarse shared N₂ grid as under-resolved for F₂ DA's K≈58 wave (the channel-representation probe returns `converged=False`) and the coarse H₂⁺ electronic grid for the Coulomb incident — the regression guards for the bugs that bit us.

- [ ] **Steps:** write `calibrate.py`, run it, set `_PHASE_COEFF`; write `test_emoscat_decks.py` (reproduce-or-beat + flag-failures, `@slow` for the 2-D confirmations); write `docs/physics/discretisation-tuning.md` + the `CLAUDE.md` skills-table + `qscat.tuning` entries; commit `feat(tuning): calibrate phase constant + gate vs eMoScat decks and coarse-grid failures`.

---

## Verification (whole sub-project)

- `uv run pytest -q -m "not slow"` pass; the `@slow` deck gates pass.
- `uv run mypy libs/qscat/qscat` 0; `uv run ruff check .` clean.
- The primitives are unit-tested on analytic potentials (turning points, uniform-k→uniform mesh, finer-where-k-larger, Gaussian angle cap, fast-wave absorption).
- `propose_grid` reproduces-or-beats the N₂/NO/F₂/(H₂⁺ proxy) eMoScat decks (same/better probe precision, same/fewer points) — with `C` calibrated to make this hold — and the channel-representation probe FLAGS the coarse-grid failures (F₂ DA K≈58, H₂⁺ Coulomb incident).
- The `discretisation-tuner` skill documents the end-to-end loop with a worked N₂ example; `docs/physics/discretisation-tuning.md` + `CLAUDE.md` updated.

## Out of scope (this plan)

- **Coupled-channel / non-adiabatic grid coupling** (the tensor-product decoupling suffices for current models).
- **Auto-running huge decks** (the H₂⁺ 2-D spot-check is proxy/Docker-deferred).
- **Rewiring N₂/NO/F₂/H₂⁺ production onto tuner grids** (the tuner EMITS configs; adoption is a separate opt-in follow-on).
- **The TW auto-tune** is best-effort (Task 6): inputs-only is the committed baseline; the auto-tuner ships if clean, else is a recorded follow-on.
