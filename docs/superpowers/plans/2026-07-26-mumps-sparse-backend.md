# MUMPS Complex-Symmetric Sparse Backend Implementation Plan (sub-project #8)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a MUMPS complex-symmetric (`SYM=2`) sparse-solver backend behind `qscat.linalg.SparseLU`, so the driven solve (#6) and the CN propagation (#7) factorize a `H = Hᵀ` matrix in ~half the work/memory — with SuperLU kept as the always-available fallback and differential oracle, and the speedup **measured**, not assumed.

**Architecture:** MUMPS is a *vendor library* provisioned in the Docker base (apt `libmumps-seq-dev`); `python-mumps` is an **optional** `qscat[mumps]` extra; `SparseLU` gains a backend dispatch (`auto` → MUMPS if present, else SuperLU). MUMPS-dependent tests skip gracefully when MUMPS is absent (green on the Mac dev box) and run in the Docker `test` image (MUMPS present).

**Tech Stack:** Python 3.12, `scipy.sparse` (SuperLU baseline/oracle), `python-mumps` (Kwant binding) + system MUMPS, `pytest-benchmark`, Docker.

**Design spec:** `docs/superpowers/specs/2026-07-26-mumps-sparse-backend-design.md`
**Optimization-target memory:** the sparse LU hot path; MKL PARDISO did all N₂/NO/F₂ in <1 hr.

## Global Constraints

- Python `>=3.12`. Everything through `uv` — never bare `python`/`pip`/conda. MUMPS-dependent work runs **in the Docker `test` image** (MUMPS present); the Mac dev box has no MUMPS, so MUMPS tests must `skipif` there.
- **Core `qscat` stays numpy/scipy-only.** MUMPS is an OPTIONAL extra (`qscat[mumps]`); the core wheel and the core test suite must pass with MUMPS absent.
- **SuperLU remains the always-available fallback AND the differential oracle.** With MUMPS absent, `backend="auto"` behaves bit-identically to today; every existing result is unchanged.
- **The matrices are complex-symmetric** (`A = Aᵀ ≠ A†`, ~1e-13). MUMPS `SYM=2` for symmetric `A`, `SYM=0` otherwise.
- **MUMPS `SYM=2` correctness trap: it takes only the UPPER TRIANGLE of `A`.** Supplying the full matrix double-counts off-diagonals and silently corrupts the factorization. This is the V1 differential-test target.
- `uv run mypy libs/qscat` at **0 errors** (strict; type-checks tests). `uv run ruff check .` clean (line length 100; rules `E, F, I, UP, B, NPY`).
- The N₂ harness must not regress: **23 PASS / 0 PENDING / 6 NOTE / 0 FAIL**, exit 0.
- Package-absolute imports; `projects/` must not import `validation/`.
- Commit trailer: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`

---

## File Structure

| File | Responsibility |
|---|---|
| `docker/base.Dockerfile` | **Modify.** Add MUMPS (`libmumps-seq-dev` + ordering libs) + pkg-config gate |
| `libs/qscat/pyproject.toml` | **Modify.** Add the `[project.optional-dependencies] mumps` extra |
| `libs/qscat/qscat/linalg/sparse_lu.py` | **Modify.** Backend dispatch; keep SuperLU path intact |
| `libs/qscat/qscat/linalg/_mumps_backend.py` | **Create.** The MUMPS `SYM=2` backend (import-guarded) |
| `libs/qscat/tests/test_sparse_lu.py` | **Modify.** Dispatch + fallback tests (no MUMPS needed) |
| `libs/qscat/tests/test_mumps_backend.py` | **Create.** V1 differential, V3 benchmark (`skipif` no MUMPS) |
| `validation/n2/test_backend_equivalence.py` | **Create.** V2: #6/#7 physics unchanged through MUMPS (`skipif`) |
| `docs/physics/mumps-sparse-backend.md` | **Create.** The measured benchmark + method note |
| `CLAUDE.md` | **Modify.** MUMPS backend, the `qscat[mumps]` extra, the base-image addition |

---

### Task 1 (SPIKE): Provision MUMPS in the Docker base + verify `python-mumps` + discover the complex-symmetric API

**This is a de-risking spike in the target environment.** `python-mumps` builds against a *system* MUMPS via pkg-config/cmake; it is NOT installable on the Mac dev box. The exact complex-symmetric (`SYM=2`, upper-triangle) API and whether Debian's `libmumps-seq-dev` cleanly satisfies the binding's build are UNKNOWN and must be established here. **Deliverable: a base image with a working MUMPS + python-mumps, and a documented, verified minimal complex-symmetric solve recipe that Tasks 3–4 build on.**

**Files:** Modify `docker/base.Dockerfile`, `libs/qscat/pyproject.toml`. Scratch spike script under the scratchpad (not committed).

- [ ] **Step 1: Add MUMPS to `docker/base.Dockerfile`.** Extend the apt install with `libmumps-seq-dev` and its ordering libraries (`libmetis-dev`, `libscotch-dev` if needed for fill-reducing ordering), and extend the pkg-config sanity gate. Example (adjust package names to what Debian bookworm actually provides):

```dockerfile
RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential curl ca-certificates pkg-config \
      libopenblas-dev liblapacke-dev libfftw3-dev \
      libmumps-seq-dev libmetis-dev libscotch-dev \
    && rm -rf /var/lib/apt/lists/*
```
Extend the gate to also report MUMPS. **Note the real risk:** Debian's MUMPS packaging may not ship a `.pc` pkg-config file named as `python-mumps` expects (`dmumps_seq`). If the pkg-config gate or the later `uv pip install python-mumps` build fails to find MUMPS, resolve it in this task — options in order: (a) install the package that provides the `.pc`/headers (`libmumps-dev`, `coinor-libmumps-dev`, or similar); (b) write a small `dmumps_seq.pc` into `PKG_CONFIG_PATH`; (c) pass CMake hints to the build; (d) if Debian packaging is intractable, document it and fall back (build MUMPS from source in the base, OR switch the base to a conda-forge `mumps-seq`). Record what actually worked.

- [ ] **Step 2: Build the base image and confirm the gate.** `DOCKER_BUILDKIT=0 docker/build.sh test` (or build just the base). Confirm the pkg-config gate prints a MUMPS version — i.e. the system library is discoverable.

- [ ] **Step 3: Declare the optional extra** in `libs/qscat/pyproject.toml`:

```toml
[project.optional-dependencies]
mumps = ["python-mumps>=0.0.6"]
```
(Core `dependencies` stay `["numpy>=2", "scipy>=1.14"]` — unchanged.)

- [ ] **Step 4: Get `python-mumps` to build/import in the container.** In the `test` image (or a throwaway container off the base), install the extra and import it. Resolve any build failure per Step 1's fallback list. Success = `python -c "import mumps"` works in the container.

- [ ] **Step 5: Discover + verify the complex-symmetric API (the spike's core).** Write a scratch script (run IN the container) that, using the real `mumps` module:
  1. builds a small complex-symmetric matrix `A = Aᵀ` (not Hermitian) and a RHS `b`;
  2. factors and solves it **via MUMPS in `SYM=2` mode, supplying only the upper triangle** (discover exactly how the binding wants symmetry declared and the matrix supplied — `mumps.Context(...)` analyze/factor/solve, or `mumps.spsolve`, and the symmetry/SYM parameter);
  3. compares to `scipy.sparse.linalg.splu` on the SAME full matrix: `‖x_mumps − x_scipy‖/‖x‖ < 1e-10`;
  4. reads MUMPS's stats (INFOG: L+U nnz / factorization memory / the ordering it chose) so Tasks 3–4 know how to surface `fill_factor`/`memory_bytes`/`ordering_used`.

  **Record the verified recipe** (exact calls, how the upper triangle is supplied, how SYM is set, how stats are read) in the Task-1 report — Tasks 3–4 consume it verbatim. If a correct `SYM=2` solve cannot be achieved, report **BLOCKED** with what failed (this is the point of a spike).

- [ ] **Step 6: Commit** the Dockerfile + pyproject changes (the scratch script is not committed; the recipe goes in the report). Message records what provisioning path worked and the verified round-off agreement.

---

### Task 2: `SparseLU` backend dispatch (SuperLU behind it; MUMPS not yet)

**Files:** Modify `libs/qscat/qscat/linalg/sparse_lu.py`, `libs/qscat/tests/test_sparse_lu.py`.

**Interfaces:**
- Consumes: nothing new (pure refactor).
- Produces: `SparseLU(A, *, ordering="COLAMD", backend="auto", symmetric=None)`, with new read-only `backend_used: str` and `ordering_used: str`; existing `shape`/`ordering`/`fill_factor`/`memory_bytes()`/`solve()` unchanged.

**Background you need.** This task is a **backward-compatible refactor with NO MUMPS yet** — it is fully developable and testable on the Mac. Extract the current SuperLU logic into a private `_ScipyBackend`, add the dispatch that (for now) always selects SuperLU, and add the `backend`/`symmetric` parameters and the `backend_used`/`ordering_used` properties. `backend="mumps"` raises `NotImplementedError` until Task 3 (or, cleaner: raise the "MUMPS not available" error path Task 3 will fill in). Every existing test and caller (`ve_cross_section_2d`, `make_sparse_cn_stepper`, `test_sparse_lu.py`) must keep passing **bit-identically**.

- [ ] **Step 1: Write the failing tests** (append to `libs/qscat/tests/test_sparse_lu.py`):

```python
def test_default_backend_is_scipy_and_bit_identical() -> None:
    """backend='auto' with MUMPS absent == the old SuperLU behaviour exactly."""
    A = _complex_symmetric(200, seed=20)
    rng = np.random.default_rng(21)
    b = rng.standard_normal(200) + 1j * rng.standard_normal(200)
    lu = SparseLU(A)  # unchanged call site
    assert lu.backend_used == "scipy"
    x = lu.solve(b)
    assert np.linalg.norm(A @ x - b) / np.linalg.norm(b) < 1e-12


def test_force_scipy_backend() -> None:
    A = _complex_symmetric(120, seed=22)
    lu = SparseLU(A, backend="scipy")
    assert lu.backend_used == "scipy"
    assert lu.ordering_used in {"COLAMD", "NATURAL", "MMD_ATA", "MMD_AT_PLUS_A"}


def test_ordering_still_applies_on_scipy_path() -> None:
    A = _complex_symmetric(300, seed=23)
    lu = SparseLU(A, ordering="MMD_AT_PLUS_A", backend="scipy")
    assert lu.ordering == "MMD_AT_PLUS_A"
    assert lu.backend_used == "scipy"


def test_symmetric_autodetect_flag_is_recorded() -> None:
    """A == A.T is detected (used by the MUMPS path later); scipy ignores it."""
    A = _complex_symmetric(80, seed=24)  # symmetric fixture
    lu = SparseLU(A)  # symmetric=None => auto-detect
    # exposed for the MUMPS path; on scipy it is informational only
    assert lu.backend_used == "scipy"
```

(These need no MUMPS. Reuse the existing `_complex_symmetric` fixture in that file.)

- [ ] **Step 2: Run → fail** (`AttributeError: backend_used` / unexpected kwarg).

- [ ] **Step 3: Refactor `sparse_lu.py`.** Keep the module docstring's SuperLU memory guidance. Structure:
  - `_Backend = Literal["auto", "scipy", "mumps"]`.
  - A private `_factor_scipy(csc, ordering)` (or a `_ScipyBackend` class) holding the exact current `splu` logic and its `fill_factor`/`memory_bytes` semantics.
  - `SparseLU.__init__(self, A, *, ordering="COLAMD", backend="auto", symmetric=None)`: convert to CSC/complex128 as today; detect `symmetric` (if `None`, `symmetric = (abs(A - A.T)).max() == 0` — but do it on the sparse structure cheaply; document the O(nnz) cost); **select the backend**: `auto` → try MUMPS (Task 3 fills this in; for now falls straight through to scipy), else scipy; `scipy` → scipy; `mumps` → the MUMPS path (Task 3) or a clear `RuntimeError("MUMPS backend requested but qscat[mumps]/system MUMPS not available")`.
  - `backend_used` records which ran; `ordering_used` returns the scipy `permc_spec` on the scipy path (MUMPS's chosen ordering on the MUMPS path, Task 3).
  - `solve`, `fill_factor`, `memory_bytes`, `shape`, `ordering` unchanged for the scipy path.

- [ ] **Step 4: Run → pass** (new + all existing `test_sparse_lu.py` tests). `uv run mypy libs/qscat` → 0; `uv run ruff check .` → clean.

- [ ] **Step 5: Regression — the whole repo is bit-identical.** `uv run pytest -q -m "not slow"` passes; `uv run python -m validation.n2.experiment` → 23/0/6/0 exit 0 (unchanged — the dispatch selects scipy everywhere, no MUMPS present locally).

- [ ] **Step 6: Commit.**

---

### Task 3 (CRUX): The MUMPS complex-symmetric backend

**Files:** Create `libs/qscat/qscat/linalg/_mumps_backend.py`; modify `sparse_lu.py` (wire the MUMPS path into the dispatch); create `libs/qscat/tests/test_mumps_backend.py`.

**Interfaces:**
- Consumes: **Task 1's verified `python-mumps` recipe** (exact `SYM=2` upper-triangle calls + stats reads); `qscat.linalg.SparseLU`'s dispatch (Task 2).
- Produces: `_mumps_backend.mumps_available() -> bool` (import guard); a `_MumpsBackend` with `factor(csc, symmetric) → handle`, `solve(handle, b)`, `fill_factor`, `memory_bytes`, `ordering_used`.

**Background you need.** This is the crux and it is **MUMPS-dependent — developed and tested in the Docker `test` image**, not on the Mac. Use Task 1's verified recipe for the exact binding calls. The essential requirements:
- **Import-guarded:** `mumps_available()` returns False (never raises) when `import mumps` fails, so `backend="auto"` cleanly falls to SuperLU and the core suite is unaffected.
- **`SYM=2` for symmetric `A`, supplying ONLY the upper triangle** (`scipy.sparse.triu(A)`), `SYM=0` for non-symmetric. This is THE correctness trap — the V1 differential test exists to catch it.
- **Factor once, solve many** (matches the `SparseLU` contract — MUMPS analyze+factor at construction, solve on demand).
- **Diagnostics from MUMPS's INFOG** (Task 1 found which entries): `fill_factor` (L+U nnz / A.nnz), `memory_bytes` (MUMPS's reported factor memory), `ordering_used` (the ordering MUMPS selected, e.g. "metis"/"amd"/"pord").

- [ ] **Step 1: Write the V1 differential test — `libs/qscat/tests/test_mumps_backend.py`.** All tests `@pytest.mark.skipif(not mumps_available(), reason="system MUMPS / qscat[mumps] not installed")` so they SKIP on the Mac and RUN in Docker.

```python
import numpy as np
import pytest
import scipy.sparse as sp
from qscat.linalg import SparseLU
from qscat.linalg._mumps_backend import mumps_available

pytestmark = pytest.mark.skipif(
    not mumps_available(), reason="system MUMPS / qscat[mumps] not installed"
)


def _complex_symmetric(n, seed):
    rng = np.random.default_rng(seed)
    nnz = 5 * n
    r = rng.integers(0, n, nnz); c = rng.integers(0, n, nnz)
    v = rng.standard_normal(nnz) + 1j * rng.standard_normal(nnz)
    m = sp.coo_matrix((v, (r, c)), shape=(n, n), dtype=complex).tocsr()
    m = m + m.T                                    # complex SYMMETRIC
    m = m + sp.identity(n, format="csr", dtype=complex) * (10.0 + 3.0j)
    return sp.csc_matrix(m)


@pytest.mark.parametrize("n", [50, 400])
def test_mumps_solve_matches_scipy_to_roundoff(n: int) -> None:
    """THE gate: MUMPS SYM=2 (upper triangle) == SuperLU on the full matrix."""
    A = _complex_symmetric(n, seed=100 + n)
    rng = np.random.default_rng(7)
    b = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    x_mumps = SparseLU(A, backend="mumps").solve(b)
    x_scipy = SparseLU(A, backend="scipy").solve(b)
    assert np.linalg.norm(x_mumps - x_scipy) / np.linalg.norm(x_scipy) < 1e-10
    assert np.linalg.norm(A @ x_mumps - b) / np.linalg.norm(b) < 1e-10


def test_mumps_used_and_reports_sym2() -> None:
    A = _complex_symmetric(120, seed=200)
    lu = SparseLU(A, backend="mumps")
    assert lu.backend_used == "mumps"
    assert lu.fill_factor >= 1.0
    assert lu.ordering_used  # non-empty; the ordering MUMPS chose


def test_mumps_upper_triangle_trap_would_be_caught() -> None:
    """A deliberately asymmetric matrix must NOT be silently treated as symmetric:
    forcing symmetric=True on a non-symmetric A must differ from the truth,
    proving the symmetric path really uses only the upper triangle."""
    rng = np.random.default_rng(3)
    n = 60
    m = sp.random(n, n, density=0.1, format="csc", random_state=rng, dtype=float) + \
        1j * sp.random(n, n, density=0.1, format="csc", random_state=rng, dtype=float)
    A = (m + sp.identity(n) * (10 + 3j)).tocsc()   # NOT symmetric
    b = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    x_true = SparseLU(A, backend="scipy").solve(b)
    x_sym_wrong = SparseLU(A, backend="mumps", symmetric=True).solve(b)  # wrongly claims sym
    # forcing SYM=2 on a non-symmetric A (upper triangle only) gives a DIFFERENT,
    # wrong answer -- confirming the symmetric path genuinely drops the lower triangle
    assert np.linalg.norm(x_sym_wrong - x_true) / np.linalg.norm(x_true) > 1e-3
    # and the correct (auto/unsymmetric) MUMPS path matches scipy
    x_unsym = SparseLU(A, backend="mumps", symmetric=False).solve(b)
    assert np.linalg.norm(x_unsym - x_true) / np.linalg.norm(x_true) < 1e-10
```

- [ ] **Step 2: Run in Docker → fail** (`_mumps_backend` missing). On the Mac these SKIP — that is expected and correct.

- [ ] **Step 3: Implement `_mumps_backend.py`** using Task 1's verified recipe: the import guard, the `SYM=2` upper-triangle factorization, `SYM=0` fallback, factor-once/solve-many, and the INFOG-based diagnostics. Wire the MUMPS path into `SparseLU`'s dispatch (`auto` selects it when `mumps_available()`; `mumps` forces it).

- [ ] **Step 4: Run the tests IN DOCKER** (`docker run --rm qmodeling:test uv run --no-sync pytest libs/qscat/tests/test_mumps_backend.py -q`) → pass. Confirm they still SKIP cleanly on the Mac (`uv run pytest libs/qscat/tests/test_mumps_backend.py -q` → skipped, not failed).

- [ ] **Step 5: mypy/ruff** (on the Mac; the import-guarded MUMPS module must type-check without MUMPS installed — guard the `import mumps` so mypy/ruff pass with it absent). Full core suite still green.

- [ ] **Step 6: Commit** with the differential-agreement numbers from the in-Docker run.

---

### Task 4: V2 physics-unchanged + V4 fallback + V3 benchmark

**Files:** Create `validation/n2/test_backend_equivalence.py` (V2); create the benchmark in `libs/qscat/tests/test_mumps_backend.py` or a `benchmarks/` module (V3); extend `test_sparse_lu.py` for V4.

- [ ] **Step 1 (V4 — fallback/absence, no MUMPS needed, runs on Mac):** add tests: `backend="auto"` with MUMPS absent → `backend_used == "scipy"` and results bit-identical; `backend="mumps"` with MUMPS absent → a clear `RuntimeError` naming the missing extra. (On the Mac these exercise the absent-MUMPS path directly.)

- [ ] **Step 2 (V2 — physics unchanged, `skipif` no MUMPS, runs in Docker):** `validation/n2/test_backend_equivalence.py` — recompute a #6 exact anchor (via `ve_cross_section_2d`, which uses `SparseLU`) and a #7 σ_TD (short propagation via `make_sparse_cn_stepper`) with the factorization forced through MUMPS vs SuperLU, and assert the physics agrees to tight rtol (the backend must not change the answer). Keep it modest-grid so it is not a 250s test; mark `@pytest.mark.slow` if needed. **Note:** `SparseLU`/`make_sparse_cn_stepper` may need a `backend=` pass-through, OR use a context/monkeypatch to force the backend — pick the cleaner: a thread-local/default-backend override is tidiest, but a simple `backend` kwarg threaded through `ve_cross_section_2d`/`make_sparse_cn_stepper` (defaulting to `"auto"`) is acceptable and explicit. Decide and document.

- [ ] **Step 3 (V3 — the benchmark, the deliverable, in Docker):** a `pytest-benchmark` (or a runnable script) that factors + solves a real N₂ 2-D matrix (build via `build_h2d` on a representative grid — the working grid N≈27k first; the production N≈143k if the container has the RAM) with `backend="scipy"` vs `backend="mumps"`, reporting **factor time, solve time, peak memory (RSS delta), fill_factor, and the ordering MUMPS used**. Print a table. **Measure the win — do not assert a specific speedup** (if MUMPS-seq does not beat SuperLU, that is a real finding; the fallback means nothing regresses). Save the numbers for the docs.

- [ ] **Step 4: Run V4 on the Mac (pass), V2+V3 in Docker (pass/measured).** Record the benchmark table.

- [ ] **Step 5: Commit** with the measured factor/solve/memory numbers (MUMPS vs SuperLU, vs the <1 hr bar) in the message.

---

### Task 5: Docs + CLAUDE.md + final verification

**Files:** Create `docs/physics/mumps-sparse-backend.md`; modify `CLAUDE.md`.

- [ ] **Step 1: Write `docs/physics/mumps-sparse-backend.md`:** why the sparse LU is the hot path; why the matrices are complex-symmetric and SuperLU can't exploit it; the MUMPS `SYM=2` backend and the upper-triangle trap; the dispatch (`auto`/`scipy`/`mumps`) with SuperLU as fallback+oracle; the **measured benchmark table** (MUMPS vs SuperLU: factor/solve/memory/fill/ordering, on the real N₂ matrix, vs the <1 hr bar — honestly, including if MUMPS did NOT win); how to provision MUMPS (Docker apt = clean; the Mac convenience path); and that this is stage-4 optimization with the Python SuperLU kept as the differential oracle. Note the deferred levers (complex64, symbolic reuse, MKL PARDISO, packaging, Rust non-LU kernels).

- [ ] **Step 2: Update `CLAUDE.md`:** the MUMPS backend + `backend=` dispatch in `qscat.linalg.SparseLU`; the `qscat[mumps]` optional extra; the `docker/base.Dockerfile` MUMPS vendor addition and the "MUMPS tests run in the container, skip on a MUMPS-less box" convention; the new doc.

- [ ] **Step 3: Full verification:**
```
uv run pytest -q -m "not slow"                 # Mac: MUMPS tests SKIP, all else pass, bit-identical
uv run mypy libs/qscat                         # 0
uv run ruff check .                            # clean
uv run python -m validation.n2.experiment      # 23/0/6/0 exit 0, unchanged
docker/build.sh test                           # base has MUMPS; in-container suite (incl MUMPS tests) passes
docker run --rm qmodeling:test uv run --no-sync pytest libs/qscat/tests/test_mumps_backend.py -q   # MUMPS tests RUN + pass
```
Report the in-container MUMPS test results and the benchmark table.

- [ ] **Step 4: Commit.**

---

## Final verification

- [ ] MUMPS provisioned in the Docker base (pkg-config gated); `python-mumps` builds/imports there; the verified `SYM=2` recipe is documented.
- [ ] `SparseLU` dispatch: `auto` → MUMPS-if-present else SuperLU; `scipy`/`mumps` forced; **SuperLU is the fallback AND the differential oracle**; existing call sites bit-identical.
- [ ] V1: MUMPS `SYM=2` (upper triangle) solve matches SuperLU to round-off; the upper-triangle trap is actively tested.
- [ ] V2: #6 anchors and #7 σ_TD unchanged through the MUMPS backend (physics is backend-independent).
- [ ] V3: the benchmark MEASURES factor/solve/memory/fill/ordering (MUMPS vs SuperLU, vs the <1 hr bar) — win reported honestly, not assumed.
- [ ] V4: MUMPS-absent → silent SuperLU fallback, bit-identical; `backend="mumps"` absent → clear error. Core suite (numpy/scipy-only) passes with no MUMPS.
- [ ] Core `qscat` stays numpy/scipy-only; MUMPS is an opt-in extra. mypy 0; ruff clean; harness 23/0/6/0; docker green.
- [ ] Deferred (not done here): complex64+refinement, symbolic/numeric reuse, MKL PARDISO, the publishing pipeline, a Rust non-LU kernel.
