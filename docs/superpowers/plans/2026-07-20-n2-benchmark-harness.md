# N₂ LCP Benchmark Harness Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax.

**Goal:** A CPU-Docker-runnable N₂ ²Π_g benchmark harness under `validation/n2/` — closed-form model checks and Houfek golden-data integrity green now, solver-dependent checks reported PENDING.

**Architecture:** A self-contained `validation/n2/` package: `loader.py` (parse the golden data), `model.py` (closed-form potentials), `reference.py` (literature + anchor coordinates + tolerances), `experiment.py` (PASS/PENDING/FAIL table + exit code), `test_n2.py` (pytest for the green checks). Modules are siblings imported by name (pytest/`python` prepend the file's directory to `sys.path`; `validation/n2/` has **no** `__init__.py`).

**Tech Stack:** Python 3.12, numpy, `qscat.units`; pytest; the existing layered Docker image.

## Global Constraints

- Python `>=3.12`; `uv` for all Python ops; run tests with `uv run pytest`.
- Atomic units throughout: energy in **Hartree**, cross sections in **bohr²**, lengths in bohr.
- Golden data `CSVE.V00.J00`: 400 rows × 32 cols, Fortran `E` notation, whitespace-separated; col 1 energy (Ha), col 2 elastic (v=0→0), cols 3–32 v=0→1…30. Source: Karel Houfek, time-independent calc (external to eMoScat).
- N₂ model params (a.u.): `D_0=0.75102, alpha_0=1.1535, R_0=2.01943, lambda_inf=6.21066, lambda_1=1.05708, R_lambda=-27.9833, lambda_c=5.38022, R_c=2.405, alpha_c=0.4, l=2, reduced_mass=12766.36`.
- Reference/anchor cross-section values are **looked up from `CSVE.V00.J00`**, never hardcoded.
- `experiment.py` exits **0** when no FAIL (PENDING ≠ FAIL), non-zero on any integrity/closed-form FAIL.
- No `__init__.py` in `validation/n2/`. Design spec: `docs/superpowers/specs/2026-07-20-n2-benchmark-harness-design.md`.
- Commit trailer: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.

---

### Task 1: Golden data + loader + integrity checks

**Files:**
- Move: `CSVE.V00.J00` (repo root) → `validation/n2/data/CSVE.V00.J00`
- Create: `validation/n2/data/MANIFEST.md`
- Create: `validation/n2/loader.py`
- Test: `validation/n2/test_loader.py`

**Interfaces:**
- Produces: `loader.DATA_PATH: Path`; `loader.load() -> CrossSectionData` where `CrossSectionData` is a dataclass with `energy: np.ndarray (400,)` (Ha), `sigma: np.ndarray (400, 31)` (bohr²; column j = v=0→j, j=0 elastic), and `n_energies`, `n_channels` properties. `loader.integrity_checks() -> list[tuple[str, bool, str]]` returning `(name, ok, detail)` for C1–C4.

- [ ] **Step 1: Move the data file and write the manifest**

```bash
mkdir -p validation/n2/data
git mv CSVE.V00.J00 validation/n2/data/CSVE.V00.J00 2>/dev/null || { mkdir -p validation/n2/data && mv CSVE.V00.J00 validation/n2/data/CSVE.V00.J00; }
```
`validation/n2/data/MANIFEST.md`:
```markdown
# CSVE.V00.J00 — N₂ vibrational-excitation cross sections (golden data)

- **Source:** Karel Houfek, time-independent calculation. External to this repo (not from eMoScat).
- **System:** electron–N₂, ²Π_g resonance (LCP model). Initial state v=0, J=0.
- **Format:** 400 rows × 32 whitespace-separated columns, Fortran `E` notation.
  - Column 1: collision energy, **Hartree** (5e-4 … 0.2, step 5e-4).
  - Column 2: elastic / vibrationally-elastic (v=0→0).
  - Columns 3–32: v=0→1, v=0→2, …, v=0→30.
- **Units:** cross sections in **atomic units (bohr²)**.
- Higher-v channels are exactly 0 below their energetic threshold.
- Used as regression anchors for a future time-independent solver (see reference.py).
```

- [ ] **Step 2: Write the failing loader test — `validation/n2/test_loader.py`**

```python
import numpy as np
import loader


def test_shape_and_grid():
    d = loader.load()
    assert d.energy.shape == (400,)
    assert d.sigma.shape == (400, 31)
    assert d.energy[0] == 5e-4
    assert abs(d.energy[-1] - 0.2) < 1e-12
    # strictly increasing
    assert np.all(np.diff(d.energy) > 0)


def test_nonnegative_and_elastic_column():
    d = loader.load()
    assert np.all(d.sigma >= 0.0)
    # elastic (v=0->0) is column index 0, grows into the resonance region
    assert d.sigma[-1, 0] > d.sigma[0, 0]


def test_threshold_ordering():
    # channel v=0->(j+1) opens at an energy >= where v=0->j opens (higher channels open later)
    d = loader.load()
    def first_open(j):
        nz = np.nonzero(d.sigma[:, j] > 0)[0]
        return d.energy[nz[0]] if nz.size else np.inf
    opens = [first_open(j) for j in range(1, 31)]  # skip elastic
    finite = [o for o in opens if np.isfinite(o)]
    assert finite == sorted(finite)


def test_integrity_checks_all_pass():
    results = loader.integrity_checks()
    assert results, "expected integrity checks"
    assert all(ok for _name, ok, _detail in results), results
```

- [ ] **Step 3: Run the test to verify it fails**

Run: `uv run pytest validation/n2/test_loader.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'loader'`.

- [ ] **Step 4: Implement `validation/n2/loader.py`**

```python
"""Load and integrity-check the Houfek time-independent N₂ VE cross-section data."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

DATA_PATH = Path(__file__).parent / "data" / "CSVE.V00.J00"


@dataclass(frozen=True)
class CrossSectionData:
    energy: np.ndarray  # (N,) Hartree
    sigma: np.ndarray   # (N, 31) bohr²; column j = v=0->j (j=0 elastic)

    @property
    def n_energies(self) -> int:
        return int(self.energy.shape[0])

    @property
    def n_channels(self) -> int:
        return int(self.sigma.shape[1])


def load(path: Path = DATA_PATH) -> CrossSectionData:
    raw = np.loadtxt(path)
    return CrossSectionData(energy=raw[:, 0].copy(), sigma=raw[:, 1:].copy())


def integrity_checks(path: Path = DATA_PATH) -> list[tuple[str, bool, str]]:
    d = load(path)
    out: list[tuple[str, bool, str]] = []
    out.append(("C1 shape 400x32", d.n_energies == 400 and d.n_channels == 31,
                f"{d.n_energies} energies, {d.n_channels} channels"))
    inc = bool(np.all(np.diff(d.energy) > 0))
    out.append(("C2 energy strictly increasing (Ha)", inc,
                f"[{d.energy[0]:.4g}, {d.energy[-1]:.4g}]"))
    out.append(("C3 cross sections non-negative (bohr^2)", bool(np.all(d.sigma >= 0.0)),
                f"min={d.sigma.min():.3e}"))

    def first_open(j: int) -> float:
        nz = np.nonzero(d.sigma[:, j] > 0)[0]
        return float(d.energy[nz[0]]) if nz.size else float("inf")

    opens = [first_open(j) for j in range(1, d.n_channels)]
    finite = [o for o in opens if np.isfinite(o)]
    out.append(("C4 channel thresholds ordered", finite == sorted(finite),
                f"{sum(np.isfinite(opens))} channels open in range"))
    return out
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `uv run pytest validation/n2/test_loader.py -q`
Expected: 4 passed.

- [ ] **Step 6: Commit**

```bash
git add validation/n2 && git commit -m "feat(n2): golden-data loader + integrity checks; relocate CSVE.V00.J00

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: Closed-form model + checks

**Files:**
- Create: `validation/n2/model.py`
- Create: `validation/n2/config.json`
- Test: `validation/n2/test_model.py`

**Interfaces:**
- Consumes: nothing from Task 1.
- Produces: `model.PARAMS: dict` (loaded from config.json); `model.v0(R)`, `model.lam(R)`, `model.v_int(r, R)`, `model.v_eff_el(r, R)` (all accept float or np.ndarray, return same); `model.model_checks() -> list[tuple[str, bool, str]]` for A1–A5.

- [ ] **Step 1: Create `validation/n2/config.json`**

```json
{
  "provenance": "reference/eMoScat/input/experimental/N2-model.json (model.potential + model.reduced_mass)",
  "units": "atomic units (Hartree, bohr)",
  "reduced_mass": 12766.36,
  "impulsemomentum": 2,
  "potential": {
    "D_0": 0.75102, "alpha_0": 1.1535, "R_0": 2.01943,
    "lambda_inf": 6.21066, "lambda_1": 1.05708, "R_lambda": -27.9833,
    "lambda_c": 5.38022, "R_c": 2.405, "alpha_c": 0.4
  }
}
```

- [ ] **Step 2: Write the failing test — `validation/n2/test_model.py`**

```python
import numpy as np
import model


def test_morse_minimum_and_depth():
    p = model.PARAMS["potential"]
    R0, D0 = p["R_0"], p["D_0"]
    # V0(R0) == -D0 exactly
    assert model.v0(R0) == np.float64(-D0) or abs(model.v0(R0) + D0) < 1e-12
    # argmin over a fine grid is at R0
    R = np.linspace(1.0, 6.0, 200001)
    assert abs(R[np.argmin(model.v0(R))] - R0) < 1e-3
    # asymptote -> 0
    assert abs(model.v0(20.0)) < 1e-6


def test_lambda_at_Rc():
    p = model.PARAMS["potential"]
    assert abs(model.lam(p["R_c"]) - p["lambda_c"]) < 1e-12


def test_v_int_is_negative_decaying_well():
    R0 = model.PARAMS["potential"]["R_0"]
    assert model.v_int(1.0, R0) < 0.0
    assert abs(model.v_int(10.0, R0)) < abs(model.v_int(1.0, R0))  # decays in r


def test_v_eff_has_centrifugal_term():
    R0 = model.PARAMS["potential"]["R_0"]
    r = 2.0
    l = model.PARAMS["impulsemomentum"]
    assert abs(model.v_eff_el(r, R0) - (model.v_int(r, R0) + l * (l + 1) / (2 * r**2))) < 1e-14


def test_model_checks_all_pass():
    results = model.model_checks()
    assert results and all(ok for _n, ok, _d in results), results
```

- [ ] **Step 3: Run to verify it fails**

Run: `uv run pytest validation/n2/test_model.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'model'`.

- [ ] **Step 4: Implement `validation/n2/model.py`**

```python
"""Closed-form N₂ LCP potentials (extracted from reference/eMoScat, verified).

E_res(R)/Γ(R) are NOT closed form (ECS eigenvalue pole) and are out of scope here.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

PARAMS: dict = json.loads((Path(__file__).parent / "config.json").read_text())


def v0(R):
    """Neutral N₂ Morse potential (Hartree). Minimum -D_0 at R_0."""
    p = PARAMS["potential"]
    a, R0, D0 = p["alpha_0"], p["R_0"], p["D_0"]
    return D0 * (np.exp(-2 * a * (R - R0)) - 2 * np.exp(-a * (R - R0)))


def lam(R):
    """Interaction strength λ(R); λ(R_c) == λ_c."""
    p = PARAMS["potential"]
    li, l1, Rl, lc, Rc = (p["lambda_inf"], p["lambda_1"], p["R_lambda"],
                          p["lambda_c"], p["R_c"])
    lam0 = (lc - li) * (1 + np.exp(l1 * (Rc - Rl)))
    return li + lam0 / (1 + np.exp(l1 * (R - Rl)))


def v_int(r, R):
    """Electron–molecule interaction potential (Hartree)."""
    return -lam(R) * np.exp(-PARAMS["potential"]["alpha_c"] * np.asarray(r) ** 2)


def v_eff_el(r, R):
    """Fixed-R electronic effective potential incl. l(l+1)/2r² centrifugal term."""
    l = PARAMS["impulsemomentum"]
    r = np.asarray(r, dtype=float)
    return v_int(r, R) + l * (l + 1) / (2 * r**2)


def model_checks() -> list[tuple[str, bool, str]]:
    p = PARAMS["potential"]
    R0, D0, Rc, lc = p["R_0"], p["D_0"], p["R_c"], p["lambda_c"]
    out: list[tuple[str, bool, str]] = []
    out.append(("A1 V0(R0) == -D_0", abs(float(v0(R0)) + D0) < 1e-12, f"{float(v0(R0)):.6f} Ha"))
    Rg = np.linspace(1.0, 6.0, 200001)
    out.append(("A2 Morse minimum at R0", abs(Rg[np.argmin(v0(Rg))] - R0) < 1e-3, f"R0={R0}"))
    out.append(("A3 V0(inf) -> 0", abs(float(v0(20.0))) < 1e-6, f"{float(v0(20.0)):.2e}"))
    out.append(("A4 lambda(Rc) == lambda_c", abs(float(lam(Rc)) - lc) < 1e-12, f"{float(lam(Rc)):.6f}"))
    out.append(("A5 V_int negative well", float(v_int(1.0, R0)) < 0.0, f"{float(v_int(1.0, R0)):.4f} Ha"))
    return out
```

- [ ] **Step 5: Run to verify pass**

Run: `uv run pytest validation/n2/test_model.py -q`
Expected: 5 passed.

- [ ] **Step 6: Commit**

```bash
git add validation/n2 && git commit -m "feat(n2): closed-form LCP model (V0, lambda, V_int) + model checks

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: Reference anchors + experiment harness + docker wrapper + README

**Files:**
- Create: `validation/n2/reference.py`
- Create: `validation/n2/experiment.py`
- Create: `validation/n2/test_n2.py`
- Create: `validation/n2/README.md`
- Create: `docker/run-n2.sh`

**Interfaces:**
- Consumes: `loader` (Task 1), `model` (Task 2).
- Produces: `experiment.run_checks() -> list[Check]` where `Check` is `(group, name, status, detail)` with status in `{"PASS","PENDING","FAIL"}`; `experiment.main() -> int` (exit code). `reference.LITERATURE`, `reference.ANCHORS`, `reference.RTOL`.

- [ ] **Step 1: Create `validation/n2/reference.py`**

```python
"""Reference values: literature resonance + Houfek golden-data anchor coordinates.

Anchor *values* are looked up from CSVE.V00.J00, never hardcoded. RTOL applies when a
future time-independent solver's output is compared at the anchor coordinates.
"""

from __future__ import annotations

import numpy as np

import loader

# Literature electron–N₂ ²Π_g shape resonance (eV) — Schulz; Berman/Domcke.
LITERATURE = {"E_res_eV": (2.3, 2.4), "Gamma_eV": (0.35, 0.55)}

# Anchor coordinates: (energy_Ha, channel_index). channel 0 = elastic, j = v=0->j.
# Chosen near E=0.2 Ha (resonance region), one mid-range, one near-threshold.
ANCHOR_COORDS = [
    (0.2, 0), (0.2, 1), (0.2, 2), (0.2, 3),
    (0.1, 1), (0.02, 1),
]
RTOL = 0.05  # 5% — tune when the TI solver lands


def anchors() -> list[tuple[float, int, float]]:
    """Resolve (energy, channel) -> reference sigma (bohr²) from the golden data."""
    d = loader.load()
    out: list[tuple[float, int, float]] = []
    for e, ch in ANCHOR_COORDS:
        i = int(np.argmin(np.abs(d.energy - e)))
        out.append((float(d.energy[i]), ch, float(d.sigma[i, ch])))
    return out
```

- [ ] **Step 2: Write the failing harness test — `validation/n2/test_n2.py`**

```python
import experiment


def test_green_groups_pass_pending_never_fail():
    checks = experiment.run_checks()
    statuses = {c[2] for c in checks}
    assert statuses <= {"PASS", "PENDING", "FAIL"}
    # Every Group A and Group C-integrity check must PASS; nothing may FAIL.
    assert not any(c[2] == "FAIL" for c in checks), [c for c in checks if c[2] == "FAIL"]
    assert any(c[2] == "PASS" for c in checks)
    assert any(c[2] == "PENDING" for c in checks)  # resonance / cross-section anchors


def test_main_exits_zero():
    assert experiment.main() == 0
```

- [ ] **Step 3: Run to verify it fails**

Run: `uv run pytest validation/n2/test_n2.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'experiment'`.

- [ ] **Step 4: Implement `validation/n2/experiment.py`**

```python
"""N₂ LCP benchmark harness. Prints a PASS/PENDING/FAIL table; exits non-zero on FAIL.

Run: python validation/n2/experiment.py  (or in Docker: docker run --rm qmodeling:runtime
python validation/n2/experiment.py)
"""

from __future__ import annotations

import sys

import loader
import model
import reference

Check = tuple[str, str, str, str]  # (group, name, status, detail)


def run_checks() -> list[Check]:
    checks: list[Check] = []

    # Group A — closed-form model (green now)
    for name, ok, detail in model.model_checks():
        checks.append(("A model", name, "PASS" if ok else "FAIL", detail))

    # Group C1–C4 — golden-data integrity (green now)
    for name, ok, detail in loader.integrity_checks():
        checks.append(("C data", name, "PASS" if ok else "FAIL", detail))

    # Group B — resonance position (needs ECS eigensolver): PENDING
    lo, hi = reference.LITERATURE["E_res_eV"]
    checks.append(("B resonance", "B1 E_res(R0) in literature window", "PENDING",
                   f"expect {lo}-{hi} eV; needs ECS eigensolver"))

    # Group C5 — cross-section value anchors vs Houfek data (needs TI solver): PENDING
    for e, ch, ref in reference.anchors():
        lbl = "elastic" if ch == 0 else f"v=0->{ch}"
        checks.append(("C anchors", f"C5 sigma({e:.4g} Ha, {lbl})", "PENDING",
                       f"ref={ref:.4e} bohr^2, rtol={reference.RTOL:.0%}; needs TI solver"))

    # Group D — time-dependent model: PENDING (later)
    checks.append(("D time-dependent", "D1 TD cross sections", "PENDING",
                   "needs time-dependent LCP propagation"))
    return checks


def main() -> int:
    checks = run_checks()
    width = max(len(f"{g}: {n}") for g, n, _s, _d in checks)
    print("N2 LCP benchmark harness")
    print("=" * (width + 30))
    for group, name, status, detail in checks:
        print(f"[{status:7}] {group}: {name:<{width - len(group) - 2}}  {detail}")
    n_pass = sum(c[2] == "PASS" for c in checks)
    n_pend = sum(c[2] == "PENDING" for c in checks)
    n_fail = sum(c[2] == "FAIL" for c in checks)
    print("=" * (width + 30))
    print(f"{n_pass} PASS, {n_pend} PENDING, {n_fail} FAIL")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 5: Run tests to verify pass**

Run: `uv run pytest validation/n2 -q`
Expected: all tests pass (loader + model + n2).

- [ ] **Step 6: Create `docker/run-n2.sh` and `validation/n2/README.md`**

`docker/run-n2.sh`:
```bash
#!/usr/bin/env bash
# Run the N2 LCP benchmark harness inside the CPU runtime image.
set -euo pipefail
export DOCKER_BUILDKIT=0
docker build -t qmodeling-base:latest -f docker/base.Dockerfile . >/dev/null
docker build --target runtime -t qmodeling:runtime -f docker/Dockerfile . >/dev/null
docker run --rm qmodeling:runtime python validation/n2/experiment.py
```
Make executable: `chmod +x docker/run-n2.sh`.

`validation/n2/README.md`: describe the physics (electron–N₂ ²Π_g LCP model), the check
groups (A/C green now, B/C5/D pending which solver), the Houfek data provenance and units,
and how to run: `uv run python validation/n2/experiment.py`, `uv run pytest validation/n2`,
and `docker/run-n2.sh`.

- [ ] **Step 7: Verify locally and in Docker**

Run: `uv run python validation/n2/experiment.py; echo "exit=$?"`
Expected: table printed, Group A + C integrity PASS, B/C5/D PENDING, `exit=0`.

Run: `docker/run-n2.sh`
Expected: same table inside the runtime image (validation/ is copied into the image), exit 0.
(If Docker is unavailable, note it — do not fake.)

- [ ] **Step 8: Commit**

```bash
git add validation/n2 docker/run-n2.sh && git commit -m "feat(n2): reference anchors + benchmark harness + docker/run-n2.sh

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Final verification

- [ ] `uv run pytest validation/n2` green.
- [ ] `uv run python validation/n2/experiment.py` prints the table, exit 0, with B/C5/D PENDING.
- [ ] `docker run --rm qmodeling:runtime python validation/n2/experiment.py` produces the same table, exit 0.
- [ ] `CSVE.V00.J00` is under `validation/n2/data/` (not repo root); MANIFEST credits Houfek.
- [ ] Anchor reference values come from the data file, not hardcoded.
