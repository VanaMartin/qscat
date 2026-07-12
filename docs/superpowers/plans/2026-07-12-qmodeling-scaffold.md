# qModeling Scaffold Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Bootstrap the `qModeling` monorepo — Python-first standard library, an example Rust kernel, Docker runtime, read-only reference sources, and the Claude Code skills/agents/`CLAUDE.md` that encode the development lifecycle — with nothing physics-specific ported yet.

**Architecture:** A `uv` workspace at the repo root. Validated Python lives in `libs/qscat` (the standard library); compiled hot paths live as PyO3/maturin crates in `native/`; research toy models live in `projects/`; validation harnesses in `validation/`; old codebases in read-only `reference/`. Claude Code assets in `.claude/` teach the design → toy → validate → optimize → promote lifecycle.

**Tech Stack:** Python 3.12 + uv; Rust (stable) + PyO3 + maturin; pytest + hypothesis + pytest-benchmark; ruff + mypy + clippy + rustfmt; Docker (multi-stage, CPU); mpmath/numpy/scipy.

## Global Constraints

- Python floor: **`requires-python = ">=3.12"`** everywhere.
- Package manager: **`uv`** for all Python env/build/lock operations. No pip/conda in instructions.
- Compiled kernel default: **Rust** via **maturin** build backend + **PyO3 >= 0.22**.
- Standard-library import name: **`qscat`** (domain-reserved; keep verbatim).
- Everything must **run locally on CPU** and be **buildable in Docker**. No GPU/CUDA code this phase.
- Units convention: **atomic units** (Hartree energy, Bohr length) throughout physics/numerics.
- `reference/` is **read-only** — never edited or imported as a live dependency.
- Every commit message ends with the trailer: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.

---

### Task 1: Repo foundation & workspace tooling

**Files:**
- Create: `pyproject.toml` (workspace root)
- Create: `.gitignore`
- Create: `README.md`
- Create: `.pre-commit-config.yaml`

**Interfaces:**
- Produces: a resolvable `uv` workspace with members `libs/*` and `native/*`; a dev dependency group named `dev`.

- [ ] **Step 1: Create `pyproject.toml`**

```toml
[project]
name = "qmodeling"
version = "0.0.0"
description = "Quantum-mechanics research monorepo (QSCAT)"
requires-python = ">=3.12"
readme = "README.md"

[tool.uv]
package = false

[tool.uv.workspace]
members = ["libs/*", "native/*"]

[dependency-groups]
dev = [
    "pytest>=8",
    "pytest-benchmark>=4",
    "hypothesis>=6",
    "ruff>=0.6",
    "mypy>=1.11",
    "maturin>=1.7",
    "numpy>=2",
    "scipy>=1.14",
    "mpmath>=1.3",
    "matplotlib>=3.9",
]

[tool.ruff]
line-length = 100
target-version = "py312"

[tool.ruff.lint]
select = ["E", "F", "I", "UP", "B", "NPY"]

[tool.mypy]
python_version = "3.12"
strict = true
```

- [ ] **Step 2: Create `.gitignore`**

```gitignore
# Python
__pycache__/
*.py[cod]
.venv/
.uv/
*.egg-info/
.pytest_cache/
.mypy_cache/
.ruff_cache/
dist/
build/
# Rust
target/
Cargo.lock
# maturin
*.so
*.pyd
# OS / editors
.DS_Store
.idea/
.vscode/
# Data / outputs
*.npz
*.h5
outputs/
```

- [ ] **Step 3: Create `.pre-commit-config.yaml`**

```yaml
repos:
  - repo: https://github.com/astral-sh/ruff-pre-commit
    rev: v0.6.9
    hooks:
      - id: ruff
        args: [--fix]
      - id: ruff-format
  - repo: local
    hooks:
      - id: cargo-fmt
        name: cargo fmt
        entry: bash -c 'cd native/qscat-kernels && cargo fmt --check'
        language: system
        files: \.rs$
        pass_filenames: false
```

- [ ] **Step 4: Create `README.md`** (concise top-level; full manual lives in CLAUDE.md, Task 8)

```markdown
# qModeling

Quantum-mechanics research monorepo (QSCAT). Python-first, HPC-capable.

Development lifecycle: **design → Python toy model → validate → optimize hot paths
in Rust → promote to the `qscat` standard library.** Everything runs on CPU locally
and is containerizable.

## Layout
- `libs/qscat/` — the standard library (validated, reusable)
- `native/` — Rust kernels (PyO3/maturin)
- `projects/` — research toy models
- `validation/` — analytic benchmarks & golden data
- `reference/` — read-only old codebases (eMoScat, libXcuda)
- `docs/` — specs, physics notes, ADRs
- `docker/` — CPU test / runtime images

## Quickstart
```bash
uv sync
uv run pytest
```

See `CLAUDE.md` for the full operating manual.
```

- [ ] **Step 5: Verify the workspace resolves**

Run: `uv sync`
Expected: completes without error and creates `.venv/` (members don't exist yet, so this only installs the dev group + root — acceptable; re-run after Tasks 2–3).

Note: if `uv sync` errors because workspace members are empty, proceed to Task 2 and re-run there. Do not add stub members here.

- [ ] **Step 6: Commit**

```bash
git add pyproject.toml .gitignore README.md .pre-commit-config.yaml
git commit -m "chore: workspace root, tooling, gitignore

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: `qscat` standard-library skeleton with a validated units module

**Files:**
- Create: `libs/qscat/pyproject.toml`
- Create: `libs/qscat/qscat/__init__.py`
- Create: `libs/qscat/qscat/units.py`
- Create: `libs/qscat/qscat/{special,dvr,ecs,evolution,linalg}/__init__.py` (empty package stubs)
- Create: `libs/qscat/tests/test_units.py`

**Interfaces:**
- Produces: `qscat.__version__: str`; `qscat.units.HARTREE_TO_EV: float`; `qscat.units.hartree_to_ev(x: float | np.ndarray) -> float | np.ndarray`; `qscat.units.ev_to_hartree(x) -> ...`. Later tasks import `qscat`.

- [ ] **Step 1: Create `libs/qscat/pyproject.toml`**

```toml
[project]
name = "qscat"
version = "0.0.0"
description = "QSCAT standard library — validated quantum-scattering numerics"
requires-python = ">=3.12"
dependencies = ["numpy>=2"]

[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[tool.hatch.build.targets.wheel]
packages = ["qscat"]
```

- [ ] **Step 2: Write the failing test — `libs/qscat/tests/test_units.py`**

```python
import numpy as np
import pytest

import qscat
from qscat import units


def test_version_is_string():
    assert isinstance(qscat.__version__, str)


def test_hartree_to_ev_known_value():
    # CODATA 2018: 1 Hartree = 27.211386245988 eV
    assert units.hartree_to_ev(1.0) == pytest.approx(27.211386245988, rel=1e-12)


def test_roundtrip_is_identity():
    x = np.array([0.0, 0.5, 2.5, -1.3])
    np.testing.assert_allclose(units.ev_to_hartree(units.hartree_to_ev(x)), x, rtol=1e-14)
```

- [ ] **Step 3: Run the test to verify it fails**

Run: `uv run pytest libs/qscat/tests/test_units.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'qscat'`.

- [ ] **Step 4: Create the package — `libs/qscat/qscat/__init__.py`**

```python
"""QSCAT standard library."""

__version__ = "0.0.0"
```

- [ ] **Step 5: Create `libs/qscat/qscat/units.py`**

```python
"""Atomic-unit conversions. All physics in qModeling uses atomic units."""

from __future__ import annotations

# CODATA 2018
HARTREE_TO_EV: float = 27.211386245988
EV_TO_HARTREE: float = 1.0 / HARTREE_TO_EV


def hartree_to_ev(x):
    """Convert energy in Hartree to electron-volts."""
    return x * HARTREE_TO_EV


def ev_to_hartree(x):
    """Convert energy in electron-volts to Hartree."""
    return x * EV_TO_HARTREE
```

- [ ] **Step 6: Create empty subpackage stubs**

Create each of these files containing a single docstring line, e.g. `libs/qscat/qscat/special/__init__.py`:

```python
"""Special functions (Coulomb, Bessel, ...). Populated during method porting."""
```

Repeat for `dvr` ("FEM-DVR basis and operators."), `ecs` ("Exterior Complex Scaling."), `evolution` ("Time propagation (Crank–Nicolson, ...)."), `linalg` ("Linear-algebra helpers.").

- [ ] **Step 7: Run the tests to verify they pass**

Run: `uv sync && uv run pytest libs/qscat/tests/ -v`
Expected: 3 passed.

- [ ] **Step 8: Commit**

```bash
git add libs/qscat pyproject.toml
git commit -m "feat(qscat): standard-library skeleton with validated units module

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: Example Rust kernel (`qscat-kernels`) with differential test

**Files:**
- Create: `native/qscat-kernels/pyproject.toml`
- Create: `native/qscat-kernels/Cargo.toml`
- Create: `native/qscat-kernels/src/lib.rs`
- Create: `native/qscat-kernels/tests/test_l2_norm.py`

**Interfaces:**
- Consumes: nothing.
- Produces: Python module `qscat_kernels` exposing `l2_norm(v: Sequence[float]) -> float`. This is the reference pattern for all future kernels: a Rust implementation differential-tested against a Python/numpy oracle.

- [ ] **Step 1: Create `native/qscat-kernels/Cargo.toml`**

```toml
[package]
name = "qscat-kernels"
version = "0.0.0"
edition = "2021"

[lib]
name = "qscat_kernels"
crate-type = ["cdylib"]

[dependencies]
pyo3 = { version = "0.22", features = ["extension-module"] }
```

- [ ] **Step 2: Create `native/qscat-kernels/pyproject.toml`**

```toml
[project]
name = "qscat-kernels"
version = "0.0.0"
description = "QSCAT compiled kernels (Rust/PyO3)"
requires-python = ">=3.12"

[build-system]
requires = ["maturin>=1.7"]
build-backend = "maturin"

[tool.maturin]
module-name = "qscat_kernels"
```

- [ ] **Step 3: Write the failing differential test — `native/qscat-kernels/tests/test_l2_norm.py`**

```python
import numpy as np

import qscat_kernels


def test_l2_norm_matches_numpy():
    v = [3.0, 4.0]
    assert qscat_kernels.l2_norm(v) == 5.0


def test_l2_norm_differential_vs_numpy():
    rng = np.random.default_rng(0)
    for _ in range(100):
        v = rng.standard_normal(rng.integers(1, 50)).tolist()
        assert abs(qscat_kernels.l2_norm(v) - float(np.linalg.norm(v))) < 1e-12
```

- [ ] **Step 4: Run the test to verify it fails**

Run: `uv run pytest native/qscat-kernels/tests/ -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'qscat_kernels'`.

- [ ] **Step 5: Implement the kernel — `native/qscat-kernels/src/lib.rs`**

```rust
use pyo3::prelude::*;

/// Euclidean (L2) norm of a vector.
#[pyfunction]
fn l2_norm(v: Vec<f64>) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

#[pymodule]
fn qscat_kernels(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(l2_norm, m)?)?;
    Ok(())
}
```

- [ ] **Step 6: Build the wheel into the environment**

Run: `uv run maturin develop --manifest-path native/qscat-kernels/Cargo.toml`
Expected: builds and installs `qscat_kernels` into `.venv`.

- [ ] **Step 7: Run the tests to verify they pass**

Run: `uv run pytest native/qscat-kernels/tests/ -v`
Expected: 2 passed.

- [ ] **Step 8: Commit**

```bash
git add native/qscat-kernels
git commit -m "feat(kernels): example PyO3 l2_norm kernel with differential test

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 4: Read-only reference sources

**Files:**
- Create: `reference/README.md`
- Add: `reference/eMoScat/` (copy of tracked source, no `.git`)
- Add: `reference/libXcuda/` (git submodule from CodeCommit)
- Modify: `.gitmodules` (created by submodule add)

**Interfaces:**
- Produces: read-only oracle trees under `reference/`. Nothing imports these.

- [ ] **Step 1: Create `reference/README.md`**

```markdown
# reference/ — read-only porting oracles

These trees are the prior QSCAT codebases, kept ONLY as sources of algorithms and
math for clean reimplementation. **Do not edit, build, or import them as dependencies.**

- `eMoScat/` — C++/CUDA electron–molecule scattering (FEM-DVR-ECS, Coulomb/Bessel,
  Crank–Nicolson, LCP/NRM). Snapshot copy.
- `libXcuda/` — recovered CUDA layer (submodule). Reference for eventual GPU kernels.

Use the `port-scout` agent to extract a method's algorithm before reimplementing it.
```

- [ ] **Step 2: Copy eMoScat tracked source (excluding its git history)**

```bash
mkdir -p reference/eMoScat
git -C /Users/martin/src/eMoScat archive HEAD | tar -x -C reference/eMoScat
```
Expected: `reference/eMoScat/README.md`, `source/`, `include/`, etc. present; no `.git`.

- [ ] **Step 3: Add libXcuda as a submodule**

```bash
git submodule add ssh://git-codecommit.eu-central-1.amazonaws.com/v1/repos/libXcuda reference/libXcuda
```
Expected: `.gitmodules` created, `reference/libXcuda/` populated.
If SSH auth fails, STOP and report — do not fall back silently; the user must fix credentials.

- [ ] **Step 4: Verify**

Run: `ls reference/eMoScat/source | head && ls reference/libXcuda | head && cat .gitmodules`
Expected: eMoScat source files listed; libXcuda contents listed; submodule entry present.

- [ ] **Step 5: Commit**

```bash
git add reference .gitmodules
git commit -m "chore(reference): add read-only eMoScat snapshot and libXcuda submodule

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 5: Docker CPU test/runtime image

**Files:**
- Create: `docker/Dockerfile`
- Create: `.dockerignore`
- Create: `docker/README.md`

**Interfaces:**
- Produces: a multi-stage image that builds the Rust kernel + installs the workspace and runs `pytest` on CPU.

- [ ] **Step 1: Create `.dockerignore`**

```dockerignore
.git
.venv
target
reference/libXcuda
**/__pycache__
.pytest_cache
.mypy_cache
.ruff_cache
dist
build
```

- [ ] **Step 2: Create `docker/Dockerfile`**

```dockerfile
# syntax=docker/dockerfile:1
FROM ghcr.io/astral-sh/uv:python3.12-bookworm-slim AS build
RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential curl ca-certificates && rm -rf /var/lib/apt/lists/*
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"
WORKDIR /app
COPY . .
RUN uv sync --frozen || uv sync
RUN uv run maturin develop --release --manifest-path native/qscat-kernels/Cargo.toml

FROM build AS test
RUN uv run pytest -q

FROM ghcr.io/astral-sh/uv:python3.12-bookworm-slim AS runtime
WORKDIR /app
COPY --from=build /app /app
ENV PATH="/app/.venv/bin:${PATH}"
CMD ["python", "-c", "import qscat, qscat_kernels; print('qscat', qscat.__version__, 'ready')"]
```

- [ ] **Step 3: Create `docker/README.md`**

```markdown
# docker/

CPU-only images. AWS deployment (later) builds on the `runtime` stage.

```bash
# Run the test stage (fails the build if tests fail):
docker build --target test -f docker/Dockerfile .

# Build and run the runtime image:
docker build --target runtime -t qmodeling:latest -f docker/Dockerfile .
docker run --rm qmodeling:latest
```
```

- [ ] **Step 4: Verify the image builds and tests pass inside it**

Run: `docker build --target test -f docker/Dockerfile .`
Expected: build completes; `pytest` stage passes. (If Docker is unavailable in the exec environment, note it and defer to the user.)

- [ ] **Step 5: Commit**

```bash
git add docker .dockerignore
git commit -m "build: multi-stage CPU Docker image (build/test/runtime)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 6: Custom skills

**REQUIRED SUB-SKILL:** Use `superpowers:writing-skills` to author each skill (correct frontmatter, structure, and a smoke check).

**Files:** create `SKILL.md` under each of:
- `.claude/skills/qm-method-lifecycle/`
- `.claude/skills/numerical-validation/`
- `.claude/skills/python-to-rust-kernel/`
- `.claude/skills/containerize-and-run/`
- `.claude/skills/qscat-conventions/`

**Interfaces:**
- Produces: five discoverable skills. Each `SKILL.md` has YAML frontmatter with `name` and `description`, then a body. Use the exact `name`/`description` below verbatim (they drive discovery). Bodies are authored via writing-skills; required content is itemized per skill.

- [ ] **Step 1: `qm-method-lifecycle/SKILL.md`**

Frontmatter:
```yaml
---
name: qm-method-lifecycle
description: Use when adding or porting any quantum-mechanics method or numerical capability to qModeling — enforces the design → Python toy model → validate → optimize-in-Rust → promote-to-qscat lifecycle.
---
```
Body must cover, as a numbered checklist the worker turns into todos:
1. **Design** — write/point to a spec (`docs/superpowers/specs/`), state physics, units (atomic), interfaces, and success criteria.
2. **Toy model** — pure-Python implementation in `projects/<name>/`; correctness over speed.
3. **Validate** — invoke `numerical-validation`; do not proceed until it passes.
4. **Optimize** — only if profiling shows a hot path; invoke `python-to-rust-kernel`.
5. **Promote** — move validated, reusable code into `libs/qscat/qscat/<subpackage>/`, keep tests, update `CLAUDE.md`/docs.
Include a rule: the Python version is always retained as the differential-test oracle after optimization.

- [ ] **Step 2: `numerical-validation/SKILL.md`**

Frontmatter:
```yaml
---
name: numerical-validation
description: Use when validating quantum/numerical code where exact equality does not apply — analytic benchmarks, convergence studies, conservation checks, and differential testing against references.
---
```
Body must cover: (a) **analytic benchmarks** (harmonic oscillator eigenvalues, hydrogen energies, Coulomb phase shifts) with `pytest.approx`/`np.testing.assert_allclose` and explicit tolerances; (b) **convergence studies** — refine grid/basis, assert monotone error decrease toward a known rate; (c) **conservation checks** — norm preservation and unitarity of time evolution; (d) **differential testing** vs `reference/` outputs and vs `mpmath` high-precision oracles; (e) tolerance conventions (state `rtol`/`atol` explicitly, never bare `==` on floats). Cross-reference `superpowers:test-driven-development`.

- [ ] **Step 3: `python-to-rust-kernel/SKILL.md`**

Frontmatter:
```yaml
---
name: python-to-rust-kernel
description: Use when a validated Python method has a proven hot path worth moving to a compiled kernel — profile, scaffold a PyO3/maturin crate in native/, mirror the Python API, benchmark, and keep Python as the differential oracle.
---
```
Body must cover: (1) **profile first** (`pytest-benchmark`/`cProfile`) — no kernel without evidence; (2) scaffold crate in `native/<name>/` mirroring `native/qscat-kernels` (Cargo.toml, pyproject with maturin backend, `src/lib.rs`); (3) API parity — Rust function signature mirrors the Python one; (4) `uv run maturin develop`; (5) differential test vs the retained Python implementation (pattern in `native/qscat-kernels/tests/`); (6) `criterion` bench in Rust; (7) record the measured speedup.

- [ ] **Step 4: `containerize-and-run/SKILL.md`**

Frontmatter:
```yaml
---
name: containerize-and-run
description: Use when packaging a qModeling capability to run locally in Docker or prepare it for AWS — reproducible multi-stage builds via uv + maturin, CPU-only.
---
```
Body must cover: reuse `docker/Dockerfile` stages (`build`/`test`/`runtime`); `uv sync --frozen` for reproducibility; building kernels with `maturin develop --release`; the `docker build --target test` gate; note that AWS deployment extends the `runtime` stage (deferred). Reference `docker/README.md`.

- [ ] **Step 5: `qscat-conventions/SKILL.md`**

Frontmatter:
```yaml
---
name: qscat-conventions
description: Reference for qModeling/QSCAT shared conventions — atomic units, FEM-DVR-ECS notation, tolerance defaults, and standard-library layout. Consult when unsure how the project names or measures things.
---
```
Body must cover: atomic units (Hartree, Bohr) as the default; `qscat` subpackage map (`special`, `dvr`, `ecs`, `evolution`, `linalg`, `units`); default tolerance conventions; naming (snake_case Python, kernel module names). Keep concise — this is a lookup, not a process.

- [ ] **Step 6: Smoke-check all skills**

Run: `for f in .claude/skills/*/SKILL.md; do echo "== $f =="; head -5 "$f"; done`
Expected: each prints valid `---` frontmatter with `name:` and `description:`. Confirm every `name:` matches its directory name.

- [ ] **Step 7: Commit**

```bash
git add .claude/skills
git commit -m "feat(skills): add qModeling lifecycle, validation, kernel, docker, conventions skills

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 7: Custom agents

**Files:** create one Markdown file per agent under `.claude/agents/`:
- `.claude/agents/port-scout.md`
- `.claude/agents/physics-reviewer.md`
- `.claude/agents/rust-kernel-engineer.md`

**Interfaces:**
- Produces: three subagent definitions with frontmatter (`name`, `description`, `tools`) and a system-prompt body. Use the exact frontmatter below.

- [ ] **Step 1: `.claude/agents/port-scout.md`**

```markdown
---
name: port-scout
description: Read-only archaeologist over reference/. Given a method or module, extracts the underlying math and algorithm (not the C++) so it can be reimplemented cleanly in Python. Use before porting anything from eMoScat/libXcuda.
tools: Read, Grep, Glob, Bash
---

You explore ONLY `reference/` (eMoScat, libXcuda). You never edit those files and never
edit code outside your report. Given a target method, produce: (1) the mathematical
formulation and key equations; (2) the algorithm/control flow; (3) inputs, outputs,
units, and boundary/edge conditions; (4) numerical pitfalls the original code handled
(e.g. ECS contour, singularities); (5) a proposed clean Python interface. Return a
concise, structured report — the caller reimplements from your report, not the C++.
```

- [ ] **Step 2: `.claude/agents/physics-reviewer.md`**

```markdown
---
name: physics-reviewer
description: Reviews quantum/numerical code for physical and numerical correctness — units, conservation laws, boundary conditions, ECS contour handling, and convergence — the issues a generic code reviewer misses. Use before promoting a method into qscat.
tools: Read, Grep, Glob, Bash
---

You review for physics/numerics correctness, not style. Check: atomic-unit consistency;
conservation (norm, probability, energy where applicable) and unitarity of time
evolution; boundary conditions and asymptotics; ECS contour correctness; basis/grid
convergence and stated tolerances; differential agreement with references/oracles.
Report findings ranked by severity with a concrete failing scenario for each. Do not
rubber-stamp; if something is unverified, say so.
```

- [ ] **Step 3: `.claude/agents/rust-kernel-engineer.md`**

```markdown
---
name: rust-kernel-engineer
description: Specialist for building PyO3/Rust kernels in native/ — mirrors a validated Python API, adds criterion benchmarks, and differential-tests against the Python oracle. Use during the optimize stage of the lifecycle.
tools: Read, Edit, Write, Bash, Grep, Glob
---

You implement compiled kernels following native/qscat-kernels as the reference pattern.
Requirements: Rust signature mirrors the Python function; build with
`uv run maturin develop`; a pytest differential test vs the retained Python
implementation must pass (tight tolerance); add a criterion benchmark and report the
measured speedup. Keep the Python version intact as the oracle. Never optimize without
a profile showing the path is hot.
```

- [ ] **Step 4: Smoke-check agents**

Run: `for f in .claude/agents/*.md; do echo "== $f =="; head -5 "$f"; done`
Expected: each prints valid frontmatter with `name`, `description`, `tools`. Confirm `name:` matches filename.

- [ ] **Step 5: Commit**

```bash
git add .claude/agents
git commit -m "feat(agents): add port-scout, physics-reviewer, rust-kernel-engineer

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 8: CLAUDE.md, docs scaffolding, and settings

**Files:**
- Create: `CLAUDE.md`
- Create: `.claude/settings.json`
- Create: `docs/physics/README.md`
- Create: `docs/adr/0001-record-architecture-decisions.md`
- Create: `docs/adr/0002-python-first-rust-kernels.md`
- Create: `validation/README.md`
- Create: `projects/README.md`

**Interfaces:**
- Consumes: everything built in Tasks 1–7 (paths, skills, agents).
- Produces: the operating manual and empty-dir READMEs that keep `projects/`, `validation/`, `docs/physics`, `docs/adr` tracked.

- [ ] **Step 1: Create `CLAUDE.md`**

Content must include, accurately reflecting the built repo:
- **Purpose** — Python-first QM research monorepo (QSCAT); resuscitate/extend prior work.
- **The lifecycle** — design → toy → validate → optimize (Rust) → promote to `qscat`; invariants: CPU-runnable + containerizable at every stage.
- **Repo map** — the `libs/ native/ projects/ validation/ reference/ docs/ docker/` tree with one line each.
- **Tech decisions** — Python 3.12 + uv; Rust + PyO3/maturin; pytest/hypothesis; ruff/mypy/clippy; atomic units; Docker CPU; GPU/AWS deferred.
- **Skills & agents** — table mapping each custom skill and agent (Tasks 6–7) to when to use it, plus a pointer to Superpowers for brainstorming/TDD/debugging/plans.
- **Reference oracles** — `reference/` is read-only; use `port-scout` before porting.
- **Common commands** — `uv sync`, `uv run pytest`, `uv run maturin develop --manifest-path native/<crate>/Cargo.toml`, `docker build --target test -f docker/Dockerfile .`.

- [ ] **Step 2: Create `.claude/settings.json`**

```json
{
  "permissions": {
    "allow": [
      "Bash(uv sync)",
      "Bash(uv run pytest:*)",
      "Bash(uv run ruff:*)",
      "Bash(uv run mypy:*)",
      "Bash(uv run maturin develop:*)",
      "Bash(cargo test:*)",
      "Bash(cargo clippy:*)",
      "Bash(git status:*)",
      "Bash(git diff:*)"
    ]
  }
}
```

- [ ] **Step 3: Create `docs/adr/0001-record-architecture-decisions.md`**

```markdown
# 1. Record architecture decisions

Date: 2026-07-12
Status: accepted

## Context
We need a durable record of significant technical choices as qModeling grows.

## Decision
Use lightweight ADRs in `docs/adr/`, numbered sequentially.

## Consequences
Each significant decision gets a short, immutable record. Supersede rather than edit.
```

- [ ] **Step 4: Create `docs/adr/0002-python-first-rust-kernels.md`**

```markdown
# 2. Python-first development with Rust kernels

Date: 2026-07-12
Status: accepted

## Context
Prior work (eMoScat) is C++/CUDA, judged increasingly costly to maintain. Research
velocity favors Python; heavy numerics need a compiled language. GPU (libXcuda) is
recovered but deferred.

## Decision
Python (uv) is the primary language. Validated methods live in `qscat`. Proven hot
paths move to Rust (PyO3/maturin) with the Python version kept as the differential
oracle. Everything runs on CPU locally and is containerizable. GPU/AWS deferred.

## Consequences
Fast iteration; safe optimization; reproducible builds. Two toolchains (Python + Rust)
must be maintained. See the `qm-method-lifecycle` skill.
```

- [ ] **Step 5: Create the tracked-dir READMEs**

`docs/physics/README.md`:
```markdown
# docs/physics
Method derivations, equations, unit conventions (atomic units), and literature
references. One file per method; link from the method's spec.
```

`validation/README.md`:
```markdown
# validation
Analytic benchmarks, golden datasets, and convergence studies shared across methods.
See the `numerical-validation` skill for how these are used.
```

`projects/README.md`:
```markdown
# projects
Per-problem research and toy models — lifecycle stages 1–2 (design + Python toy).
Validated, reusable code is promoted from here into `libs/qscat`.
```

- [ ] **Step 6: Verify links and paths in CLAUDE.md resolve**

Run: `uv run pytest -q && ls CLAUDE.md .claude/settings.json docs/physics/README.md docs/adr/*.md validation/README.md projects/README.md`
Expected: tests pass; all files listed. Manually confirm every path named in `CLAUDE.md` exists.

- [ ] **Step 7: Commit**

```bash
git add CLAUDE.md .claude/settings.json docs projects validation
git commit -m "docs: add CLAUDE.md operating manual, ADRs, and dir scaffolding

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Final verification (whole scaffold)

- [ ] `uv sync` resolves the full workspace.
- [ ] `uv run maturin develop --manifest-path native/qscat-kernels/Cargo.toml` builds the kernel.
- [ ] `uv run pytest` is green (units + kernel differential tests).
- [ ] `docker build --target test -f docker/Dockerfile .` passes (or is explicitly deferred if Docker unavailable).
- [ ] `import qscat` and `import qscat_kernels` both succeed.
- [ ] All 5 skills and 3 agents present with valid frontmatter; every `name:` matches its file/dir.
- [ ] `reference/eMoScat` populated (no `.git`); `reference/libXcuda` submodule present.
- [ ] `CLAUDE.md` accurately describes the built repo; every path it names exists.
