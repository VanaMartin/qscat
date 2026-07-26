# TI Energy-Sweep Symbolic Reuse + Cross-Section Display Implementation Plan (sub-project #9)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `SparseLU.refactor(A_new)` (reuse the symbolic analysis, re-run only the numeric factorization), wire `ve_cross_section_2d` to sweep energies analyze-once/refactor-per-energy, and deliver a dense exact-2-D σ(E) curve for N₂ displayed against Houfek's golden data via a reusable plotting utility.

**Architecture:** The reuse lives in `qscat.linalg.SparseLU` — MUMPS `factor(reuse_analysis=True)` (skips the SCOTCH ordering), a correct fresh-`splu` fallback on scipy. `ve_cross_section_2d` (#6) uses it; the physics is untouched. The display is a reusable `plot_cross_sections` utility + an N₂ figure driver. MUMPS-dependent tests skip on the Mac, run in Docker; SuperLU stays the differential oracle.

**Tech Stack:** Python 3.12, `scipy.sparse`, `qscat.linalg.SparseLU` (+MUMPS backend), matplotlib, `pytest-benchmark`.

**Design spec:** `docs/superpowers/specs/2026-07-26-ti-energy-sweep-reuse-design.md`

## Global Constraints

- Python `>=3.12`. Everything through `uv`; MUMPS work runs in the Docker `test` image (Mac has no MUMPS → MUMPS tests `@skipif(not mumps_available())`).
- **SuperLU stays the always-available fallback AND the differential oracle.** `refactor` on scipy re-runs `splu` (correct, no reuse); the reuse speedup materializes only on MUMPS.
- **The reuse is bit-identical:** `SparseLU(A0); lu.refactor(A1); lu.solve(b)` must equal `SparseLU(A1).solve(b)` to round-off (V1). Physics through the reuse sweep is unchanged (V2).
- **`A(E) = E_tot·I − H` has a constant sparsity pattern across E** (H fixed, identity keeps the diagonal). `refactor` requires an identical pattern to the analyzed matrix — a guard raises on a mismatch (MUMPS `reuse_analysis` is only valid for the same structure).
- `refactor` keeps the SAME backend and `symmetric` decision as construction (the analysis assumed them).
- Core `qscat` stays numpy/scipy-only; `uv run mypy libs/qscat` at **0** (the MUMPS module stays import-guarded); `uv run ruff check .` clean (line length 100).
- The N₂ harness must not regress: **23 PASS / 0 PENDING / 6 NOTE / 0 FAIL**, exit 0. `projects/` must not import `validation/` (the Houfek-reading display driver lives on the `validation/` side or reads a path).
- Commit trailer: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`

---

## File Structure

| File | Responsibility |
|---|---|
| `libs/qscat/qscat/linalg/sparse_lu.py` | **Modify.** `SparseLU.refactor()`; `_ScipyBackend.refactor()` |
| `libs/qscat/qscat/linalg/_mumps_backend.py` | **Modify.** `_MumpsBackend.refactor()` (reuse_analysis) + pattern guard |
| `libs/qscat/tests/test_sparse_lu.py` | **Modify.** V1 scipy refactor differential + pattern-guard tests |
| `libs/qscat/tests/test_mumps_backend.py` | **Modify.** V1 MUMPS reuse differential (`@skipif`, Docker) |
| `projects/n2_2d_cross_section/cross_section_2d.py` | **Modify.** sweep analyze-once/refactor-per-energy |
| `projects/n2_2d_cross_section/test_cross_section_2d.py` | **Modify.** V2 sweep == per-energy (physics unchanged) |
| `projects/n2_2d_cross_section/cross_section_plot.py` | **Create.** `plot_cross_sections(...)` reusable utility |
| `validation/n2/ti_curve.py` | **Create.** dense N₂ exact-2D σ(E) sweep + Houfek overlay driver |
| `validation/n2/test_ti_curve.py` | **Create.** V4 curve-vs-Houfek |
| `benchmarks/sweep_reuse.py` | **Create.** V3 sweep reuse-vs-no-reuse speedup (Docker) |
| `docs/physics/figures/n2-2d-ti-cross-section.png` | **Create.** the committed σ(E)-vs-Houfek figure |
| `docs/physics/ti-energy-sweep-reuse.md`, `CLAUDE.md` | **Create/Modify.** method + measured results |

---

### Task 1: `SparseLU.refactor()` — reuse the symbolic analysis

**Files:** Modify `libs/qscat/qscat/linalg/sparse_lu.py`, `_mumps_backend.py`, `test_sparse_lu.py`, `test_mumps_backend.py`.

**Interfaces:**
- Produces: `SparseLU.refactor(A_new: sp.spmatrix) -> None`; `_ScipyBackend.refactor(csc)`, `_MumpsBackend.refactor(csc)`.

**Background you need.** `_MumpsBackend.__init__` already does `set_matrix(a); analyze(); factor()` on a persistent `self._ctx`. `refactor` re-supplies a same-pattern matrix and calls `factor(reuse_analysis=True)` — the Context keeps the analysis, so the SCOTCH ordering is skipped. `_ScipyBackend` has no clean symbolic reuse, so its `refactor` just re-runs `splu` (correct, no speedup). **The pattern guard is a correctness requirement:** `reuse_analysis=True` is only valid for the identical structure, so each backend stores the analyzed matrix's pattern (of the SUPPLIED matrix — `triu` for symmetric, full for not) and `refactor` raises `ValueError` on a mismatch. Verified in a Docker spike: MUMPS analyze-once + `factor(reuse_analysis=True)` per diagonal shift matches fresh SuperLU to ~1e-15.

- [ ] **Step 1: Write the failing scipy-path tests (Mac, no MUMPS) — append to `test_sparse_lu.py`.**

```python
def test_refactor_scipy_matches_fresh_factorization() -> None:
    """refactor(A1) then solve == a fresh SparseLU(A1).solve, on the scipy path."""
    n = 200
    A0 = _complex_symmetric(n, seed=40)
    # A1 = A0 with a different diagonal shift -> SAME sparsity pattern
    A1 = (A0 + (2.0 + 1.0j) * sp.identity(n, dtype=complex)).tocsc()
    rng = np.random.default_rng(41)
    b = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    lu = SparseLU(A0, backend="scipy")
    lu.refactor(A1)
    x = lu.solve(b)
    x_fresh = SparseLU(A1, backend="scipy").solve(b)
    assert np.linalg.norm(x - x_fresh) / np.linalg.norm(x_fresh) < 1e-10
    assert np.linalg.norm(A1 @ x - b) / np.linalg.norm(b) < 1e-10


def test_refactor_rejects_pattern_mismatch() -> None:
    """reuse_analysis is only valid for an identical pattern -> guard raises."""
    A0 = _complex_symmetric(80, seed=42)
    B = _complex_symmetric(80, seed=43)  # different random pattern, same shape
    lu = SparseLU(A0, backend="scipy")
    with pytest.raises(ValueError, match="pattern"):
        lu.refactor(B)
    with pytest.raises(ValueError, match="shape|pattern"):
        lu.refactor(_complex_symmetric(70, seed=44))  # different shape too


def test_refactor_reuses_backend_and_symmetry() -> None:
    A0 = _complex_symmetric(100, seed=45)
    A1 = (A0 + (1.0 + 1.0j) * sp.identity(100, dtype=complex)).tocsc()
    lu = SparseLU(A0, backend="scipy")
    assert lu.symmetric is True
    lu.refactor(A1)
    assert lu.backend_used == "scipy" and lu.symmetric is True
```

- [ ] **Step 2: Run → fail** (`AttributeError: refactor`).

- [ ] **Step 3: Implement.** In `_ScipyBackend`: store the analyzed pattern and add `refactor`:

```python
def __init__(self, csc: sp.csc_matrix[np.complex128], ordering: _Ordering) -> None:
    self._ordering = ordering
    self._pattern = (csc.indices.copy(), csc.indptr.copy())   # for the refactor guard
    self._lu: spla.SuperLU[np.complex128] = spla.splu(csc, permc_spec=ordering)

def refactor(self, csc: sp.csc_matrix[np.complex128]) -> None:
    _check_pattern(self._pattern, csc)   # raises ValueError on mismatch
    self._lu = spla.splu(csc, permc_spec=self._ordering)   # scipy: fresh factor, no reuse
```

In `_MumpsBackend`: store the supplied-matrix pattern (`triu` or full) and add `refactor`:

```python
def refactor(self, csc: sp.csc_matrix[np.complex128]) -> None:
    a = sp.triu(csc).tocsc() if self._symmetric else csc
    _check_pattern(self._pattern, a)     # a's pattern must match the analyzed one
    self._ctx.set_matrix(a, symmetric=self._symmetric)
    self._ctx.factor(reuse_analysis=True)   # <-- reuse the SCOTCH ordering
```

(store `self._pattern = (a.indices.copy(), a.indptr.copy())` in `__init__` after building `a`.)

A shared helper `_check_pattern(pattern, csc)` raising `ValueError("...shape..." / "...pattern...")` on shape or `indices/indptr` mismatch (canonicalize `csc` with `.sorted_indices()`/`.sum_duplicates()` first so the comparison is well-defined).

In `SparseLU`: add

```python
def refactor(self, A_new: sp.spmatrix) -> None:
    """Re-factorize `A_new` reusing this object's symbolic analysis.

    `A_new` MUST share the original matrix's sparsity pattern (e.g. a diagonal
    shift `E*I - H` across energies). On the MUMPS backend this reuses the
    analysis (skips re-ordering); on scipy it re-runs `splu` (correct, no
    reuse). Keeps the original backend and symmetry decision. Raises
    `ValueError` on a shape/pattern mismatch.
    """
    if A_new.shape != self._shape:
        raise ValueError(f"refactor shape {A_new.shape} != {self._shape}")
    csc = sp.csc_matrix(A_new, dtype=np.complex128)
    self._impl.refactor(csc)
    self._nnz = int(csc.nnz)
```

- [ ] **Step 4: Run scipy tests → pass** (Mac). `uv run mypy libs/qscat` → 0; `uv run ruff check .` → clean; full core suite bit-identical.

- [ ] **Step 5: Write the MUMPS reuse differential test (`@skipif`, Docker) — `test_mumps_backend.py`.**

```python
def test_mumps_refactor_reuses_analysis_matches_fresh() -> None:
    """analyze once, refactor(A_shift) per shift == fresh SuperLU each time."""
    n = 400
    A0 = _complex_symmetric(n, seed=300)
    lu = SparseLU(A0, backend="mumps")
    rng = np.random.default_rng(9)
    b = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    for shift in (2.0 + 1.0j, -3.0 + 0.5j, 5.0 - 2.0j):
        A = (A0 + shift * sp.identity(n, dtype=complex)).tocsc()
        lu.refactor(A)
        x = lu.solve(b)
        x_ref = SparseLU(A, backend="scipy").solve(b)
        assert np.linalg.norm(x - x_ref) / np.linalg.norm(x_ref) < 1e-9
    assert lu.backend_used == "mumps"
```

- [ ] **Step 6: Run in Docker** (`docker run ... uv sync --extra mumps ... pytest test_mumps_backend.py`) → pass; skips on the Mac. **Commit** with the in-Docker reuse rel-err.

---

### Task 2: sweep `ve_cross_section_2d` with reuse + benchmark

**Files:** Modify `projects/n2_2d_cross_section/cross_section_2d.py`, `test_cross_section_2d.py`; create `benchmarks/sweep_reuse.py`.

**Background you need.** The current loop builds a fresh `SparseLU((e_tot*ident - H))` per energy (`cross_section_2d.py:194-202`). Rewrite it to construct the solver ONCE (analyze+factor at the first open energy) and `refactor((e_tot*ident - H))` for each subsequent energy — the returned σ and its scalar/array shape contract are unchanged. Note `_sigma_at_one_energy` takes the `lu` object; keep that, just feed it the reused-then-refactored solver. Skip closed/`E<=0` energies as today (no factorization needed there — but the pattern is E-independent, so the solver can still be built once at the first energy that needs it).

- [ ] **Step 1: V2 test — sweep result == per-energy result.** In `test_cross_section_2d.py`, assert `ve_cross_section_2d(..., E=[array of energies])` (now reuse-swept) equals the same call computed one energy at a time (scalar calls) to round-off — the reuse must not change the physics. (This runs on the Mac via the scipy path AND in Docker via MUMPS; both must hold.)
- [ ] **Step 2: Run → (may already pass if the refactor is transparent; if you change the loop, confirm it stays green).**
- [ ] **Step 3: Rewrite the sweep loop** in `ve_cross_section_2d` to analyze-once/refactor-per-energy. Preserve the `return_wavefunction`/scalar/array contract exactly.
- [ ] **Step 4: Run V2 → pass** (Mac scipy + the existing anchor tests). mypy 0, ruff clean.
- [ ] **Step 5: V3 benchmark — `benchmarks/sweep_reuse.py`.** On the MUMPS backend, sweep `M` energies (e.g. 50–100) over a real N₂ working-grid matrix `(E·I − H)` two ways: (a) reuse — one `SparseLU` + `refactor` per energy; (b) no-reuse — a fresh `SparseLU` per energy. Report total wall time each and the per-energy analysis saving. **Measure — do not assert a speedup;** the analysis is a bounded fraction of the factorization, so the saving is real but modest and grows with `M`. Also record the scipy path (no reuse → no speedup, expected). Run it in Docker; save the numbers.
- [ ] **Step 6: Commit** with the measured sweep reuse timing (reuse vs no-reuse total, and the fraction saved).

---

### Task 3: dense σ(E) + `plot_cross_sections` + N₂ figure vs Houfek

**Files:** Create `projects/n2_2d_cross_section/cross_section_plot.py`, `validation/n2/ti_curve.py`, `validation/n2/test_ti_curve.py`; commit `docs/physics/figures/n2-2d-ti-cross-section.png`.

**Interfaces:**
- `cross_section_plot.plot_cross_sections(E_grid, sigma, *, channels=None, reference=None, thresholds=None, usable=None, title=None, path) -> None` — experiment-agnostic: plots `sigma[E, channel]` (bohr²) vs `E_grid` per VE channel, overlaying `reference` golden data (an `(E_ref, sigma_ref[E, channel])` bundle) if given, marking thresholds; `matplotlib.use("Agg")`.
- `validation/n2/ti_curve.py`: `compute_ti_curve(E_grid, vprimes) -> sigma[E, channel]` (the dense exact-2D sweep via `ve_cross_section_2d`, which now reuses the analysis), `houfek_reference()` (load `CSVE.V00.J00` via `validation.n2.loader`), and a `main()` that computes the curve, plots it vs Houfek, saves the PNG + a `.npz`.

**Background you need.** `plot_cross_sections` lives in `projects/` (generic, no Houfek dependency — it takes reference data as an argument). The N₂ *driver* that reads Houfek lives in `validation/n2/` (which MAY import `projects/` and `validation.n2.loader`; `projects/` must NOT import `validation/`). The dense curve is now affordable because of Task 1/2's reuse. The exact-2D solver matched Houfek to ~1e-6 at the #6 gated anchors, so the dense curve should track Houfek tightly through the resonance, with the known elastic-background / near-threshold caveats (which the exact solver largely closes — show honestly).

- [ ] **Step 1: Tests.** `test_ti_curve.py`: `compute_ti_curve` returns real, ≥0 σ of the right shape; V4 — at the gated anchor energies the dense-curve σ matches Houfek within the #6 documented bound (reuse `reference.ANCHOR_FACTOR`/the exact2d tolerance); `plot_cross_sections` produces a non-empty PNG (tmp path). Keep the test grid/energy-count modest (or mark `@slow`) so it is not minutes-long; the committed figure is generated separately at full density.
- [ ] **Step 2: Implement** `plot_cross_sections` (generic) and `ti_curve.py` (the N₂ driver).
- [ ] **Step 3: Generate the real figure at full density** — a dense `E_grid` over [0, 0.2] Ha, the VE channels, exact-2D σ(E) via the reuse sweep, overlaid on Houfek. Commit `docs/physics/figures/n2-2d-ti-cross-section.png`; note the `.npz` location.
- [ ] **Step 4: Commit** with the curve-vs-Houfek agreement summary.

---

### Task 4: docs + CLAUDE.md + final verification

**Files:** Create `docs/physics/ti-energy-sweep-reuse.md`; modify `CLAUDE.md`.

- [ ] **Step 1: Write `docs/physics/ti-energy-sweep-reuse.md`:** the symbolic/numeric split (analysis pattern-only, reused; numeric per energy); `SparseLU.refactor` (MUMPS `reuse_analysis`; scipy fresh-splu fallback; the pattern guard); the wired sweep; the **measured** sweep-reuse speedup (honestly — bounded by the analysis fraction); the dense exact-2D σ(E) curve and its Houfek agreement, embedding the figure; that SuperLU stays the differential oracle; the deferred levers (complex64, PARDISO, packaging, Rust). Cross-reference `docs/physics/mumps-sparse-backend.md`.
- [ ] **Step 2: Update `CLAUDE.md`:** `SparseLU.refactor` (symbolic reuse); the reuse-enabled energy sweep in `ve_cross_section_2d`; the `plot_cross_sections` utility + `validation/n2/ti_curve`; the new doc + figure.
- [ ] **Step 3: Full verification:**
```
uv run pytest -q -m "not slow"        # Mac: MUMPS tests skip, all else pass, bit-identical
uv run mypy libs/qscat                # 0
uv run ruff check .                   # clean
uv run python -m validation.n2.experiment   # 23/0/6/0, exit 0, unchanged
docker/build.sh test                  # MUMPS reuse tests RUN + pass in-container
```
Report the in-container reuse test result + the sweep benchmark + the curve-vs-Houfek numbers.
- [ ] **Step 4: Commit.**

---

## Final verification

- [ ] `SparseLU.refactor()` reuses the symbolic analysis (MUMPS `reuse_analysis`; scipy fresh-splu fallback), with a pattern guard; solve after refactor is bit-identical to a fresh factorization (V1), both backends.
- [ ] `ve_cross_section_2d` sweeps analyze-once/refactor-per-energy; σ is unchanged vs per-energy (V2); the N₂ harness stays 23/0/6/0.
- [ ] The sweep reuse speedup is MEASURED (V3, honestly — bounded by the analysis fraction), reuse vs no-reuse, on MUMPS.
- [ ] The dense exact-2D σ(E) curve matches Houfek across the resonance (V4); the figure + `.npz` committed; `plot_cross_sections` is experiment-agnostic (N₂ now, F₂/NO later).
- [ ] SuperLU stays the differential oracle; core `qscat` numpy/scipy-only; no `projects/`→`validation/` import; mypy 0; ruff clean; docker green.
- [ ] Deferred (not done here): complex64, MKL PARDISO, packaging, Rust non-LU kernel, interactive display.
