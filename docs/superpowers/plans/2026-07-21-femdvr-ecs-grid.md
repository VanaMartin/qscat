# FEM-DVR-ECS Grid + Kinetic Operator Implementation Plan (sub-project #1)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development. Steps use checkbox (`- [ ]`) syntax. This is correctness-critical numerical code — follow the extraction doc exactly.

**Goal:** A validated FEM-DVR-ECS radial grid + kinetic operator in `projects/femdvr_ecs/`, passing analytic benchmarks (particle-in-box, harmonic oscillator, ECS continuum rotation, bound-state θ-independence), then promoted to `qscat.dvr`/`qscat.ecs`.

**Architecture:** Pure Python/numpy. `spec.py` (grid spec dataclasses) → `grid.py` (`FemDvrEcsGrid`: GLL nodes, differentiation matrix, complex element boundaries, bridge-summed weights, Dirichlet drop) → `kinetic.py` (FEM-DVR −½d²/dz² assembly) → `operators.py` (`hamiltonian`, `eigen`). Correctness over speed (dense).

**Tech Stack:** Python 3.12, numpy, scipy; pytest.

## Global Constraints

- Python `>=3.12`; `uv run pytest` for tests. Atomic units (Hartree, bohr); electron mass = 1.
- **Authoritative construction reference:** `.superpowers/sdd/femdvr-ecs-extraction.md` (port-scout's extraction of `reference/eMoScat`). Implementers MUST read it before coding — it has the exact bridge-weight, normalization, and ECS-Jacobian recipe.
- Design spec: `docs/superpowers/specs/2026-07-21-femdvr-ecs-grid-design.md`.
- `projects/femdvr_ecs/` has **no `__init__.py`** — modules import by bare name (pytest/py prepend the file's dir to sys.path). Tests live in `projects/femdvr_ecs/` beside the code.
- Gauss-Lobatto nodes/weights and the differentiation matrix are built with numpy/scipy — NOT the reference's hand-rolled QL solver.
- Key correctness traps (from extraction): bridge weight = sum of both elements' local Lobatto weight at the shared node; kinetic basis normalization uses that **global bridge-summed** weight; drop the two outermost global points (Dirichlet); H is complex-symmetric non-Hermitian → `np.linalg.eig`, sort by `Re(E)`.
- Commit trailer: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.

---

### Task 1: Grid spec + FEM-DVR-ECS grid geometry

**Files:**
- Create: `projects/femdvr_ecs/spec.py`
- Create: `projects/femdvr_ecs/gll.py` (Gauss-Lobatto-Legendre nodes/weights + differentiation matrix)
- Create: `projects/femdvr_ecs/grid.py`
- Test: `projects/femdvr_ecs/test_gll.py`, `projects/femdvr_ecs/test_grid.py`

**Interfaces:**
- Produces: `spec.ElementSpec(length: float, angle_deg: float = 0.0)`, `spec.GridSpec(quadrature: int, elements: list[ElementSpec], x_min: float = 0.0)` (validates complex elements are contiguous at the end). `gll.gll_nodes_weights(n) -> (x, w)` on (−1,1); `gll.diff_matrix(x) -> D` with `D[j,i] = L_i'(x_j)`. `grid.FemDvrEcsGrid(spec)` with `.points` (complex (nb,)), `.weights` (complex (nb,)), `.real_points` (real (nb,)), `.R0: float`, `.n: int`.

- [ ] **Step 1: Read the extraction doc.** Read `.superpowers/sdd/femdvr-ecs-extraction.md` sections 1–3 fully before writing code.

- [ ] **Step 2: Write failing GLL tests — `test_gll.py`**

```python
import numpy as np
import gll


def test_gll_nodes_endpoints_and_count():
    x, w = gll.gll_nodes_weights(6)
    assert x.shape == (6,) and w.shape == (6,)
    assert np.isclose(x[0], -1.0) and np.isclose(x[-1], 1.0)
    assert np.all(np.diff(x) > 0)               # ascending
    assert np.isclose(w.sum(), 2.0)             # weights integrate 1 over (-1,1)


def test_gll_integrates_polynomials_exactly():
    # n-point GLL is exact for degree <= 2n-3
    x, w = gll.gll_nodes_weights(6)             # exact to degree 9
    for p in range(0, 10):
        approx = np.sum(w * x**p)
        exact = 0.0 if p % 2 else 2.0 / (p + 1)
        assert abs(approx - exact) < 1e-12, (p, approx, exact)


def test_diff_matrix_differentiates_exactly():
    x, _ = gll.gll_nodes_weights(7)
    D = gll.diff_matrix(x)                       # D[j,i] = L_i'(x_j); (D @ f)[j] = f'(x_j)
    for c in ([0, 1, 0, 0], [0, 0, 1, 0], [1, -2, 3, -1]):   # polynomials up to degree 3
        f = sum(ci * x**i for i, ci in enumerate(c))
        df = sum(i * ci * x**(i - 1) for i, ci in enumerate(c) if i >= 1)
        assert np.allclose(D @ f, df, atol=1e-10)
```

- [ ] **Step 3: Run to verify fail** — `uv run pytest projects/femdvr_ecs/test_gll.py -q` → `ModuleNotFoundError: No module named 'gll'`.

- [ ] **Step 4: Implement `gll.py`**

Gauss-Lobatto-Legendre: interior nodes are the roots of `P'_{n-1}`; endpoints are ±1. Weights `w_i = 2 / (n(n-1) [P_{n-1}(x_i)]²)`. Use `numpy.polynomial.legendre`:
```python
import numpy as np
from numpy.polynomial import legendre as L


def gll_nodes_weights(n: int):
    """Gauss-Lobatto-Legendre nodes and weights on (-1, 1), n points incl. endpoints."""
    if n < 2:
        raise ValueError("n >= 2")
    # interior nodes = roots of P'_{n-1}
    coeff = np.zeros(n)          # P_{n-1}
    coeff[n - 1] = 1.0
    dcoeff = L.legder(coeff)     # P'_{n-1}
    interior = np.sort(L.legroots(dcoeff)) if n > 2 else np.array([])
    x = np.concatenate(([-1.0], interior, [1.0]))
    Pn1 = L.legval(x, coeff)     # P_{n-1}(x_i)
    w = 2.0 / (n * (n - 1) * Pn1**2)
    return x, w


def diff_matrix(x: np.ndarray) -> np.ndarray:
    """Collocation differentiation matrix D with (D @ f)[j] = f'(x_j) for f sampled at x.
    D[j, i] = L_i'(x_j), via barycentric weights (robust for any node set)."""
    n = x.size
    # barycentric weights
    bw = np.ones(n)
    for i in range(n):
        bw[i] = 1.0 / np.prod([x[i] - x[k] for k in range(n) if k != i])
    D = np.zeros((n, n))
    for j in range(n):
        for i in range(n):
            if i != j:
                D[j, i] = (bw[i] / bw[j]) / (x[j] - x[i])
        D[j, j] = -np.sum(D[j, :])   # negative sum trick
    return D
```

- [ ] **Step 5: Run GLL tests** — `uv run pytest projects/femdvr_ecs/test_gll.py -q` → 3 passed.

- [ ] **Step 6: Write failing grid tests — `test_grid.py`**

```python
import numpy as np
import gll
from spec import ElementSpec, GridSpec
from grid import FemDvrEcsGrid


def _real_grid(nq=6, lengths=(1.0, 1.0, 1.0)):
    return GridSpec(quadrature=nq, elements=[ElementSpec(L) for L in lengths], x_min=0.0)


def test_point_count_and_dirichlet_drop():
    nq, nel = 6, 3
    g = FemDvrEcsGrid(_real_grid(nq, (1.0,) * nel))
    assert g.n == nel * (nq - 1) + 1 - 2          # bridge sharing + 2 endpoints dropped
    # outermost points (x_min=0 and x_max=3) are NOT in .points (Dirichlet)
    assert g.real_points.min() > 0.0 and g.real_points.max() < 3.0


def test_real_region_points_are_real():
    g = FemDvrEcsGrid(_real_grid())
    assert np.allclose(g.points.imag, 0.0)         # no ECS -> purely real
    assert np.allclose(g.points, g.real_points)


def test_ecs_tail_points_are_rotated():
    # 2 real elements then 1 complex element at 30 deg; pivot R0 = 2.0
    spec = GridSpec(quadrature=6, elements=[ElementSpec(1.0), ElementSpec(1.0), ElementSpec(1.0, 30.0)])
    g = FemDvrEcsGrid(spec)
    assert np.isclose(g.R0, 2.0)
    tail = g.points[g.real_points > g.R0 + 1e-9]
    # z = R0 + (x-R0) e^{i theta}; arg of (z-R0) ~ 30 deg
    ang = np.degrees(np.angle(tail - g.R0))
    assert np.allclose(ang, 30.0, atol=1e-6)


def test_weights_bridge_summed_at_shared_nodes():
    # interior element-boundary points carry a weight ~ sum of two half-element contributions;
    # they should be (roughly) larger than the small end-weights within an element.
    g = FemDvrEcsGrid(_real_grid(nq=6, lengths=(1.0, 1.0)))
    assert np.all(g.weights.real > 0)


def test_spec_rejects_noncontiguous_complex():
    import pytest
    with pytest.raises(ValueError):
        GridSpec(quadrature=6, elements=[ElementSpec(1.0, 30.0), ElementSpec(1.0, 0.0)])  # complex before real
```

- [ ] **Step 7: Run to verify fail** — `uv run pytest projects/femdvr_ecs/test_grid.py -q` → import error.

- [ ] **Step 8: Implement `spec.py` and `grid.py`** following extraction §1–§3. Key points:
  - `GridSpec.__post_init__`: validate that once an element has `angle_deg != 0`, all later elements are complex too (ECS tail contiguous); `R0 = x_min + sum(real element lengths)`.
  - `grid.py`: per element, complex length `Lk = length * exp(i·deg2rad(angle))`; complex cumulative boundaries `az`; place local GLL nodes `z = hz·ξ + hz + az[k]` with `hz = 0.5·Lk`; accumulate global points (shared boundary node written once) and **bridge-summed weights** `w[shared] = w_left_end + w_right_start`; keep `real_points` from the unscaled lengths; drop the first and last global point.
  - Store `self.points, self.weights, self.real_points, self.R0, self.n`, and also `self._dLp` (from `gll.diff_matrix`) and per-element bookkeeping (`hz` per element, global index ranges) that `kinetic.py` will need — expose via attributes `self.nq`, `self.element_ranges`, `self.hz` (list of complex half-lengths), `self.dLp`.

- [ ] **Step 9: Run grid tests** — `uv run pytest projects/femdvr_ecs/test_grid.py -q` → 5 passed.

- [ ] **Step 10: Commit**

```bash
git add projects/femdvr_ecs && git commit -m "feat(dvr): FEM-DVR-ECS grid geometry + GLL nodes/diff-matrix

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: Kinetic operator + Hamiltonian/eigen + B1 (particle-in-box) + B2 (harmonic oscillator)

**Files:**
- Create: `projects/femdvr_ecs/kinetic.py`
- Create: `projects/femdvr_ecs/operators.py`
- Test: `projects/femdvr_ecs/test_kinetic_benchmarks.py`

**Interfaces:**
- Consumes: `FemDvrEcsGrid` (Task 1), its `.dLp`, `.hz`, `.element_ranges`, `.weights`, `.points`, `.nq`.
- Produces: `kinetic.kinetic(grid, mass) -> (n,n) complex ndarray`; `operators.hamiltonian(grid, V, mass) -> ndarray` (V callable `V(z)` or ndarray at `grid.points`); `operators.eigen(H) -> (E, vecs)` with `E` sorted by ascending `Re`.

- [ ] **Step 1: Write the failing benchmark tests — `test_kinetic_benchmarks.py`**

```python
import numpy as np
from spec import ElementSpec, GridSpec
from grid import FemDvrEcsGrid
from operators import hamiltonian, eigen


def _box(L=1.0, nel=4, nq=10):
    return FemDvrEcsGrid(GridSpec(quadrature=nq, elements=[ElementSpec(L / nel)] * nel))


def test_B1_particle_in_box():
    L, m = 1.0, 1.0
    g = _box(L=L, nel=4, nq=12)
    H = hamiltonian(g, lambda z: 0.0 * z, mass=m)   # V = 0
    E, _ = eigen(H)
    exact = np.array([n**2 * np.pi**2 / (2 * m * L**2) for n in range(1, 6)])
    assert np.allclose(E[:5].real, exact, rtol=1e-6), (E[:5].real, exact)
    assert np.allclose(E[:5].imag, 0.0, atol=1e-9)


def test_B1_spectral_convergence():
    L, m = 1.0, 1.0
    err = []
    for nq in (6, 9, 12):
        g = _box(L=L, nel=3, nq=nq)
        E, _ = eigen(hamiltonian(g, lambda z: 0.0 * z, mass=m))
        err.append(abs(E[0].real - np.pi**2 / (2 * m * L**2)))
    assert err[0] > err[1] > err[2]                 # error falls as order rises


def test_B2_harmonic_oscillator():
    m, omega, L = 1.0, 1.0, 20.0
    xc = L / 2
    g = FemDvrEcsGrid(GridSpec(quadrature=10, elements=[ElementSpec(L / 10)] * 10))
    H = hamiltonian(g, lambda z: 0.5 * m * omega**2 * (z - xc) ** 2, mass=m)
    E, _ = eigen(H)
    exact = np.array([omega * (n + 0.5) for n in range(5)])
    assert np.allclose(E[:5].real, exact, rtol=1e-6), (E[:5].real, exact)
```

- [ ] **Step 2: Run to verify fail** — `uv run pytest projects/femdvr_ecs/test_kinetic_benchmarks.py -q` → import error.

- [ ] **Step 3: Implement `kinetic.py`** following extraction §2. Per element k (global index range `[a, b)` in the retained basis, local GLL weights `wl`, diff matrix `dLp`, complex half-length `hz`):
  - local scaled weights `wze[l] = hz * wl[l]`; scaled derivatives `dBF[i, l] = dLp[l, i] / hz` (note `dLp[j,i]=L_i'(x_j)`, so `dBF` indexed [basis i, node l]);
  - normalize each basis row by `1/sqrt(w_global[i])` (the bridge-summed global weight from the grid);
  - `T_local[i, j] = (1/(2*mass)) * sum_l wze[l] * dBF[i, l] * dBF[j, l]`;
  - scatter-add `T_local` into the global matrix over the element's retained indices, so shared boundary rows/cols accumulate across adjacent elements (the bridge coupling). Respect the Dirichlet trim (first element drops local index 0, last element drops local index nq−1).

  The exact index bookkeeping mirrors `KineticEnergy.cpp:15-87` as summarized in the extraction — implement to match, and let the B1 test (which is extremely sensitive to any assembly error) be the arbiter.

- [ ] **Step 4: Implement `operators.py`**

```python
import numpy as np


def hamiltonian(grid, V, mass):
    from kinetic import kinetic
    T = kinetic(grid, mass)
    Vals = V(grid.points) if callable(V) else np.asarray(V)
    Vals = np.broadcast_to(Vals, (grid.n,)).astype(complex)
    return T + np.diag(Vals)


def eigen(H):
    E, vecs = np.linalg.eig(H)           # complex, non-Hermitian
    order = np.argsort(E.real)
    return E[order], vecs[:, order]
```

- [ ] **Step 5: Run the benchmark tests** — `uv run pytest projects/femdvr_ecs/test_kinetic_benchmarks.py -q` → 3 passed. If B1 fails, the bug is in the kinetic assembly (bridge weight or normalization) — debug there, do NOT loosen the tolerance.

- [ ] **Step 6: Commit**

```bash
git add projects/femdvr_ecs && git commit -m "feat(dvr): FEM-DVR kinetic operator; particle-in-box & HO benchmarks pass

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: ECS validation — B3 (continuum rotation) + B4 (bound-state θ-independence)

**Files:**
- Test: `projects/femdvr_ecs/test_ecs_benchmarks.py`

(No new production code — ECS is already built into `grid.py`/`kinetic.py`; this task is the numerical-validation gate for the ECS machinery.)

**Interfaces:** consumes the Task 1–2 modules.

- [ ] **Step 1: Write the ECS benchmark tests — `test_ecs_benchmarks.py`**

```python
import numpy as np
from spec import ElementSpec, GridSpec
from grid import FemDvrEcsGrid
from operators import hamiltonian, eigen


def _ecs_grid(theta_deg, real_len=10.0, nreal=8, ncomplex=8, nq=8):
    els = [ElementSpec(real_len / nreal)] * nreal + [ElementSpec(real_len / nreal, theta_deg)] * ncomplex
    return FemDvrEcsGrid(GridSpec(quadrature=nq, elements=els))


def test_B3_continuum_rotation():
    # Free particle on an ECS grid: continuum eigenvalues lie on arg(E) ~ -2*theta.
    theta = 30.0
    g = _ecs_grid(theta)
    E, _ = eigen(hamiltonian(g, lambda z: 0.0 * z, mass=1.0))
    # pick mid-spectrum eigenvalues with sizeable |E| (avoid ~0 and top-of-grid edge)
    mag = np.abs(E)
    sel = E[(mag > 0.2) & (mag < 5.0)]
    ang = np.degrees(np.angle(sel))
    # most such eigenvalues cluster near -2*theta
    near = np.abs(ang - (-2 * theta)) < 5.0
    assert near.mean() > 0.5, (np.median(ang), -2 * theta)


def test_B4_bound_state_theta_independence():
    # Square well V=-V0 on [0,a], deep enough for a bound state; energy invariant under theta.
    a, V0 = 3.0, 5.0
    def Vwell(z):
        return np.where(np.real(z) <= a, -V0, 0.0).astype(complex)
    Eb = []
    for theta in (20.0, 35.0):
        g = _ecs_grid(theta, real_len=12.0, nreal=6, ncomplex=6, nq=10)
        E, _ = eigen(hamiltonian(g, Vwell, mass=1.0))
        bound = E[E.real < 0].real
        assert bound.size >= 1, "expected a bound state"
        Eb.append(bound.min())
    assert abs(Eb[0] - Eb[1]) < 1e-4, Eb        # theta-independent
```

Note for the implementer: if `test_B4` finds the well's real cutoff at `a=3.0` must sit on an element boundary for a clean result, adjust `nreal`/`real_len` so an element edge lands at `a` (document the choice). The physics target — a θ-independent bound level — is fixed; the grid layout may be tuned to express it cleanly.

- [ ] **Step 2: Run** — `uv run pytest projects/femdvr_ecs/test_ecs_benchmarks.py -q` → 2 passed. If B3/B4 fail, the bug is in the ECS Jacobian handling in the kinetic assembly (the `e^{iθ}` in `hz`) — fix in `kinetic.py`/`grid.py`, do not weaken the test.

- [ ] **Step 3: Commit**

```bash
git add projects/femdvr_ecs && git commit -m "test(ecs): continuum-rotation & bound-state theta-independence benchmarks pass

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 4: Promote to `qscat.dvr` / `qscat.ecs` + docs

**Files:**
- Create: `libs/qscat/qscat/dvr/femdvr_ecs.py` (grid + kinetic + operators), `libs/qscat/qscat/dvr/__init__.py` (re-exports)
- Create: `libs/qscat/qscat/ecs/__init__.py` (ECS contour/transform helper + re-export)
- Create: `libs/qscat/tests/test_femdvr_ecs.py` (the benchmark suite, importing from `qscat`)
- Create: `docs/physics/femdvr-ecs.md`
- Modify: `CLAUDE.md` (note the new capability)

**Interfaces:**
- Produces: `from qscat.dvr import FemDvrEcsGrid, kinetic, hamiltonian, eigen`, `from qscat.dvr import ElementSpec, GridSpec`; `qscat.ecs.ecs_map(x, R0, theta_deg)` helper (the `z(x)` transform, extracted for reuse).

- [ ] **Step 1: Move the validated modules into `qscat.dvr`.** Consolidate `spec.py`, `gll.py`, `grid.py`, `kinetic.py`, `operators.py` into `libs/qscat/qscat/dvr/` (keep the module split or merge into `femdvr_ecs.py` — prefer keeping `gll.py` and `spec.py` separate, they are reusable). Fix imports to package-relative (`from .gll import ...`). `qscat/dvr/__init__.py` re-exports the public names.

- [ ] **Step 2: Extract the ECS map into `qscat.ecs`**

```python
# libs/qscat/qscat/ecs/__init__.py
"""Exterior Complex Scaling utilities."""
import numpy as np


def ecs_map(x, R0: float, theta_deg: float):
    """z(x) = x for x<=R0, else R0 + (x-R0) e^{i theta}. x may be array; theta in degrees."""
    x = np.asarray(x, dtype=float)
    eit = np.exp(1j * np.deg2rad(theta_deg))
    return np.where(x <= R0, x.astype(complex), R0 + (x - R0) * eit)
```
Have `grid.py` use `qscat.ecs.ecs_map` for its tail coordinates (single source of the ECS transform).

- [ ] **Step 3: Port the benchmark tests to `libs/qscat/tests/test_femdvr_ecs.py`** — same B1–B4 assertions, importing `from qscat.dvr import ...`. Run: `uv run pytest libs/qscat/tests/test_femdvr_ecs.py -q` → all pass (B1, B1-convergence, B2, B3, B4).

- [ ] **Step 4: Leave `projects/femdvr_ecs/` in place** as the development toy (the spec keeps it as the origin of the promoted code). Add a one-line note at the top of `projects/femdvr_ecs/` (a `README.md`) pointing to the promoted `qscat.dvr` location.

- [ ] **Step 5: Write `docs/physics/femdvr-ecs.md`** — the method (FEM-DVR-ECS), the kinetic assembly, the ECS map, the diagonal-potential assumption, and the four validation benchmarks with their analytic references. Cross-reference the extraction doc.

- [ ] **Step 6: Update `CLAUDE.md`** — add `qscat.dvr` (FEM-DVR-ECS grid + kinetic) and `qscat.ecs` (complex-scaling map) to the standard-library capability description.

- [ ] **Step 7: Verify + commit**

Run: `uv run pytest libs/qscat -q` → all pass (units + new FEM-DVR-ECS suite). Then:
```bash
git add libs/qscat docs/physics projects/femdvr_ecs CLAUDE.md
git commit -m "feat(qscat): promote FEM-DVR-ECS grid + kinetic to qscat.dvr/qscat.ecs

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Final verification

- [ ] `uv run pytest projects/femdvr_ecs libs/qscat/tests/test_femdvr_ecs.py -q` green (B1–B4 in both locations).
- [ ] B1/B2 hit `rtol ≤ 1e-6` and B1 shows spectral convergence; B3 clusters on `arg(E)≈−2θ`; B4 bound level θ-invariant `< 1e-4`.
- [ ] `from qscat.dvr import FemDvrEcsGrid, kinetic, hamiltonian, eigen` works; `qscat.ecs.ecs_map` works.
- [ ] `docs/physics/femdvr-ecs.md` and the `CLAUDE.md` update describe the capability accurately.
- [ ] No Rust, no N₂ physics, no `__init__.py` under `projects/femdvr_ecs/`.
