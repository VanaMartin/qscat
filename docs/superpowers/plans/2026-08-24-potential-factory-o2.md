# Potential Factory — O₂ (image match → benchmark on the tables → promotion) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** The first real-molecule model out of the factory: an O₂-like `FlexibleDiatomicModel` fitted to Alt & Houfek's published nonlocal-model data (PRA **103**, 032829 (2021)), first from digitised figures ("image match", to learn what the 2-D form can do), then benchmarked against the authors' tabulated curves; promoted into `qscat.model`/`qscat.factory` with a tuner-built `MoleculePreset` so `qscat-run` can run it.

**Architecture:** Three phases in one plan. (1) Data: a `mastering-references` note for the paper, digitised `V_0`/`V_ion`/`Γ` tables with a stated precision, and a loader that checks them. (2) Fit: an O₂ `Target` (T0 from the curve + spectroscopic constants, T1 from the digitised curves + EA(O), T3 from Table II), the factory run, a report and an overlay figure; the acceptance test gates the report against the measured budget. (3) Promotion: the ansatz moves to `qscat.model.flexible`, the factory to `qscat.factory`, an `O2` preset with `propose_grid`-built grids and the 2-D spot-check, and the docs.

**Tech Stack:** as the core plan, plus `qscat.tuning.propose_grid` / `refine_to_2d_convergence`, `apps/qscat-run` presets, matplotlib (only in the figure driver, `pytest.importorskip`), WebPlotDigitizer or equivalent for the figures (manual step, recorded).

**Spec:** `docs/superpowers/specs/2026-08-24-potential-factory-design.md` — "Data acquisition for O₂", "Validation (b)–(d)", "Placement and lifecycle". Depends on the core plan and the budget plan being merged.

## Global Constraints

- Every published number used carries a page-anchored locator in `reference/literature/alt-houfek-2021-pra103-032829.md` (the `mastering-references` skill; page-offset check before any locator is written; parity checks *run*, not asserted).
- Digitised data is labelled as such everywhere (file header, loader docstring, report provenance) with its precision; it is never presented as the paper's numbers.
- Phase 2 (tables from the authors) is gated on the data being present: tests `pytest.skip` when `validation/factory/data/o2/tables/` is absent — never fake it.
- Atomic units in every file; the paper's eV/Å values are converted in the loader with `qscat.units` (state the constants used).
- The factory never fits an experimental observable; the paper's VE cross sections (Figs. 5–9) are **not** read by any code here.
- `qscat.core` must not import `qscat.model`/`qscat.factory` (`test_core_no_model_import.py`).
- Long runs are `@pytest.mark.slow`; the O₂ 2-D spot-check runs in Docker/MUMPS if the laptop cannot hold it.

---

### Task 1: Reference note for Alt & Houfek (2021)

**Files:**
- Create: `reference/literature/alt-houfek-2021-pra103-032829.md`
- Modify: `reference/literature/README.md` (table row)
- (PDF: `reference/literature/alt-houfek-2021-pra103-032829.pdf`, gitignored — obtain from <https://utf.mff.cuni.cz/librtfy/papers/0000842/physreva.103.032829_alt.houfek.pdf>)

**Interfaces:**
- Produces: locators the O₂ target cites — Table I (EA(O) 1.4611 eV, EA(O₂) 0.450 eV expt / 0.474 calc+shift, D₀ 5.165 eV; p. 032829-3), Table II (`a0 = 13.836690`, `a1 = 0.892095`, `a2 = −0.935987`, `b0 = 3.015014`, `b1 = 0.718160`, `c = −14.260279`; p. 032829-4), Eqs. (20)–(31) (p. 032829-4), Fig. 2 (curves; p. 032829-3), Fig. 3 (eigenphase sums at `R = 1.80 … 2.25 a₀`; p. 032829-5), Fig. 4 (`Γ̃`, `Δ̃`; p. 032829-5), the `l = 2`, `α = 2.5` statement (p. 032829-4), and the statistical factors `g = 2/3` / `1/3` (p. 032829-4).

- [ ] **Step 1: Invoke the skill and write the note**

Run `Skill: mastering-references` with the PDF path; follow it: page-offset check (this paper is journal-paginated `032829-N`, extractor page `N` = printed `032829-N` — verify on page 1), then the note in the house form (Role, Locator table, Equations transcribed, Parameters and numeric values with the parity check against the repo — none yet; state "no code consumes these values until `validation/factory/targets/o2.py`", Findings and limits, Pagination).

- [ ] **Step 2: Add the README row**

```markdown
| [`alt-houfek-2021-pra103-032829.md`](alt-houfek-2021-pra103-032829.md) | V. Alt, K. Houfek, *Resonant collisions of electrons with O₂ via the lowest-lying ²Π_g state of O₂⁻*, Phys. Rev. A **103**, 032829 (2021). [DOI](https://doi.org/10.1103/PhysRevA.103.032829) | **The O₂ target of the potential factory**: the published nonlocal-model construction (MRCI curves + R-matrix eigenphase sums → `A(R)`, `B(R)`, `α`, `c`, Table II). |
```

- [ ] **Step 3: Commit**

```bash
git add reference/literature/alt-houfek-2021-pra103-032829.md reference/literature/README.md
git commit -m "docs(refs): Alt & Houfek 2021 (O2 nonlocal model) reference note"
```

---

### Task 2: Digitised curves and their loader

**Files:**
- Create: `validation/factory/data/o2/README.md`, `validation/factory/data/o2/v0.csv`, `validation/factory/data/o2/v_ion.csv`, `validation/factory/data/o2/gamma.csv`
- Create: `validation/factory/targets/__init__.py` (empty), `validation/factory/targets/o2_data.py`
- Test: `validation/factory/targets/test_o2_data.py`

**Interfaces:**
- Produces: `O2Curves(R: ndarray, v0: ndarray, v_ion: ndarray, gamma: ndarray, precision_Ha: float, source: str)` and `load_o2(kind: Literal["digitised", "tables"] = "digitised") -> O2Curves` (`"tables"` reads `validation/factory/data/o2/tables/{v0,v_ion,gamma}.csv` with the same columns and raises `FileNotFoundError` when absent). CSV columns: `R_bohr, value_eV` (Fig. 2's axes are bohr and eV; `Γ` is plotted ×2 — the loader halves it and says so). Energies are converted with `qscat.units` (`EV_TO_HARTREE`, or whatever name `units.py` exposes — read it). `v0` is shifted so `v0(R → ∞) = 0` (the paper plots the O(³P)+O(³P) limit at 0 eV already; assert the last point is within precision of 0); `v_ion` is on the same scale.

- [ ] **Step 1: Digitise (manual, recorded)**

Digitise Fig. 2 (p. 032829-3) with WebPlotDigitizer: the ³Σ_g⁻ curve (`v0`), the ²Π_g full+dashed curve (`v_ion`: bound region full line, `Re` of the resonance region dashed), and `Γ(R)×2` (dash-dotted), 40–60 points each over `R ∈ [1.6, 6.0]` bohr. Record in `README.md`: the tool and version, the axis calibration points, the estimated precision (state it as a number, e.g. `±0.03 eV` → `precision_Ha`), and the date. Write the three CSVs with a two-line header (`# digitised from Alt & Houfek 2021 Fig. 2, see README.md` / `R_bohr,value_eV`).

- [ ] **Step 2: Write the failing tests**

```python
# validation/factory/targets/test_o2_data.py
from __future__ import annotations

import numpy as np
import pytest

from validation.factory.targets.o2_data import load_o2

EV = 27.211386
EA_O_EV = 1.4611  # Alt & Houfek 2021 Table I (expt), p. 032829-3


def test_digitised_curves_are_monotone_where_they_must_be():
    c = load_o2("digitised")
    assert np.all(np.diff(c.R) > 0)
    i_min = int(np.argmin(c.v0))
    assert np.all(np.diff(c.v0[:i_min]) < 0) and np.all(np.diff(c.v0[i_min:]) > 0)
    assert np.all(c.gamma >= 0.0)


def test_asymptote_matches_ea_of_oxygen_within_precision():
    c = load_o2("digitised")
    assert abs(c.v0[-1]) < 3 * c.precision_Ha
    assert abs((c.v_ion[-1] - c.v0[-1]) + EA_O_EV / EV) < 3 * c.precision_Ha


def test_gamma_is_zero_beyond_the_crossing():
    c = load_o2("digitised")
    cross = c.R[np.flatnonzero(np.sign(c.v_ion - c.v0)[:-1] * np.sign(c.v_ion - c.v0)[1:] < 0)[0]]
    assert np.all(c.gamma[c.R > cross + 0.1] < 3 * c.precision_Ha)


def test_tables_kind_skips_cleanly_when_absent():
    try:
        load_o2("tables")
    except FileNotFoundError:
        pytest.skip("authors' tables not present")
```

- [ ] **Step 3: Run tests to verify they fail**

Run: `uv run --no-sync pytest validation/factory/targets/test_o2_data.py -q`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 4: Implement the loader**

```python
# validation/factory/targets/o2_data.py
"""Alt & Houfek (2021) O2 curves. `digitised` = read off Fig. 2 (precision stated in
data/o2/README.md); `tables` = the authors' tabulated curves when present."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np
import numpy.typing as npt
from qscat import units  # use the eV->Hartree constant units.py exposes

_DATA = Path(__file__).resolve().parents[1] / "data" / "o2"
_PRECISION_EV = {"digitised": 0.03, "tables": 1e-4}  # digitised value from README.md


@dataclass(frozen=True)
class O2Curves:
    R: npt.NDArray[np.float64]
    v0: npt.NDArray[np.float64]
    v_ion: npt.NDArray[np.float64]
    gamma: npt.NDArray[np.float64]
    precision_Ha: float
    source: str


def _read(path: Path) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    arr = np.loadtxt(path, delimiter=",", comments="#", skiprows=0)
    order = np.argsort(arr[:, 0])
    return arr[order, 0], arr[order, 1] * units.EV_TO_HARTREE


def load_o2(kind: Literal["digitised", "tables"] = "digitised") -> O2Curves:
    base = _DATA if kind == "digitised" else _DATA / "tables"
    R0, v0 = _read(base / "v0.csv")
    R1, v_ion = _read(base / "v_ion.csv")
    R2, g2 = _read(base / "gamma.csv")
    R = R0
    v_ion = np.interp(R, R1, v_ion)
    gamma = np.clip(np.interp(R, R2, g2, left=0.0, right=0.0) / 2.0, 0.0, None)  # Fig. 2 plots Gamma x2
    return O2Curves(R, v0, v_ion, gamma, _PRECISION_EV[kind] * units.EV_TO_HARTREE, f"Alt & Houfek 2021 Fig. 2 ({kind})")
```

- [ ] **Step 5: Run tests, commit**

Run: `uv run --no-sync pytest validation/factory/targets/test_o2_data.py -q`
Expected: 3 passed, 1 skipped

```bash
uv run ruff check validation/factory && uv run ruff format validation/factory
git add validation/factory/data/o2/README.md validation/factory/data/o2/v0.csv validation/factory/data/o2/v_ion.csv validation/factory/data/o2/gamma.csv validation/factory/targets/__init__.py validation/factory/targets/o2_data.py validation/factory/targets/test_o2_data.py
git commit -m "data(o2): digitised Alt & Houfek Fig. 2 curves with stated precision + loader"
```

---

### Task 3: The O₂ `Target` and seed

**Files:**
- Create: `validation/factory/targets/o2.py`
- Test: `validation/factory/targets/test_o2_target.py`

**Interfaces:**
- Consumes: `projects.potential_factory.target.{Target, NeutralTarget, ResonanceTarget, CouplingTarget, Curve, Provenance}`, `projects.potential_factory.ansatz.{FlexibleDiatomicModel, SmoothR}`, `load_o2`.
- Produces: `O2_MU = 14582.6` (¹⁶O₂ reduced mass, `m(¹⁶O)/2` in electron masses — compute from `15.99491461956 u × 1822.888486` and state it), `o2_target(kind="digitised") -> Target` with `ell=2`, `charge=0`, `coordinates=("R",)`, `neutral` = the `v0` curve over `[1.8, 5.5]` plus `constants={"R_e": 2.282, "D_e": 0.1898}` (R_e from spectroscopy — verify the value against NIST/Huber–Herzberg before committing and record the locator; `D_e = D_0 + ω_e/2` from Table I's 5.165 eV), `resonance` = `(v_ion, gamma, ea=1.4611 eV)`, `coupling = CouplingTarget.from_alt_houfek(a0=13.836690, a1=0.892095, a2=-0.935987, b0=3.015014, b1=0.718160, alpha=2.5, R_range=(1.8, 2.25), eps_window=(0.002, 0.22))`, provenance per slot. Also `o2_seed() -> FlexibleDiatomicModel`: `mu=O2_MU, ell=2, D_e=0.19, R_e=2.28, betas=(1.4,), p=3, lam=SmoothR(6.0, 1.0, 5.0, 2.3, R_e=2.28), alpha=SmoothR(0.5, 0.0, 1.0, 0.0, R_e=2.28), shell=None, alpha_b=2.0, r_b=3.0` (a d-wave well of N₂-like depth; the fit moves it).

- [ ] **Step 1: Write the failing tests**

```python
# validation/factory/targets/test_o2_target.py
from __future__ import annotations

import numpy as np
from qscat.model import ResonanceModel

from validation.factory.targets.o2 import O2_MU, o2_seed, o2_target


def test_o2_target_is_complete_and_has_provenance():
    t = o2_target()
    assert t.ell == 2 and t.mu == O2_MU and t.coordinates == ("R",)
    assert t.neutral is not None and t.resonance is not None and t.coupling is not None
    assert {"neutral", "resonance", "coupling"} <= set(t.provenance)
    assert t.coupling.alpha_exponent == 2.5
    assert t.resonance.ea > 0.05


def test_o2_target_coupling_is_positive_on_the_window():
    t = o2_target()
    eps = np.geomspace(*t.coupling.eps_window, 5)
    R = np.linspace(*t.coupling.R_range, 4)
    assert np.all(t.coupling.gamma_tilde(eps[:, None], R[None, :]) > 0)


def test_o2_seed_is_a_resonance_model():
    assert isinstance(o2_seed(), ResonanceModel)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run --no-sync pytest validation/factory/targets/test_o2_target.py -q`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Implement**

```python
# validation/factory/targets/o2.py
"""The O2 Target for the potential factory (Alt & Houfek, PRA 103, 032829 (2021))."""

from __future__ import annotations

from qscat import units

from projects.potential_factory.ansatz import FlexibleDiatomicModel, SmoothR
from projects.potential_factory.target import CouplingTarget, Curve, NeutralTarget, Provenance, ResonanceTarget, Target
from validation.factory.targets.o2_data import load_o2

__all__ = ["O2_MU", "o2_target", "o2_seed"]

_NOTE = "reference/literature/alt-houfek-2021-pra103-032829.md"
O2_MU = 15.99491461956 * 1822.888486 / 2.0  # m(16O)/2 in m_e
_EA_O = 1.4611 * units.EV_TO_HARTREE  # Table I (expt), p. 032829-3
_D0 = 5.165 * units.EV_TO_HARTREE  # Table I (expt), p. 032829-3
_OMEGA_E = 1580.19 / 219474.63  # cm^-1 -> Ha; Huber-Herzberg (verify + cite before commit)
_R_E = 2.282  # bohr; Huber-Herzberg (verify + cite before commit)


def o2_target(kind: str = "digitised") -> Target:
    c = load_o2(kind)  # type: ignore[arg-type]
    rng = (1.8, min(5.5, float(c.R.max())))
    return Target(
        name=f"O2 ({kind})",
        mu=O2_MU,
        ell=2,
        charge=0,
        coordinates=("R",),
        neutral=NeutralTarget(Curve.from_table(c.R, c.v0), {"R_e": _R_E, "D_e": _D0 + 0.5 * _OMEGA_E, "omega_e": _OMEGA_E}, rng),
        resonance=ResonanceTarget(Curve.from_table(c.R, c.v_ion), Curve.from_table(c.R, c.gamma), _EA_O, rng),
        coupling=CouplingTarget.from_alt_houfek(
            a0=13.836690, a1=0.892095, a2=-0.935987, b0=3.015014, b1=0.718160, alpha=2.5,
            R_range=(1.8, 2.25), eps_window=(0.002, 0.22),
        ),
        provenance={
            "neutral": Provenance(_NOTE, f"Fig. 2 p. 032829-3 ({c.source}); Table I p. 032829-3"),
            "resonance": Provenance(_NOTE, f"Fig. 2 p. 032829-3 ({c.source}); EA(O) Table I"),
            "coupling": Provenance(_NOTE, "Table II p. 032829-4; Eqs. (24)-(27) p. 032829-4"),
        },
    )


def o2_seed() -> FlexibleDiatomicModel:
    return FlexibleDiatomicModel(
        mu=O2_MU, ell=2, D_e=0.19, R_e=2.28, betas=(1.4,), p=3,
        lam=SmoothR(f_inf=6.0, f_0=1.0, f_1=5.0, R_f=2.3, R_e=2.28),
        alpha=SmoothR(f_inf=0.5, f_0=0.0, f_1=1.0, R_f=0.0, R_e=2.28),
        shell=None, alpha_b=2.0, r_b=3.0,
    )
```

- [ ] **Step 4: Run tests, verify the spectroscopic constants, commit**

Run: `uv run --no-sync pytest validation/factory/targets/test_o2_target.py -q` → 3 passed.
Before committing, check `_R_E` and `_OMEGA_E` against the NIST WebBook (Huber–Herzberg) and write the locator into the reference note's "Parameters" section under a "spectroscopic constants used alongside" heading.

```bash
uv run ruff check validation/factory && uv run ruff format validation/factory
git add validation/factory/targets/o2.py validation/factory/targets/test_o2_target.py reference/literature/alt-houfek-2021-pra103-032829.md
git commit -m "feat(validation/factory): the O2 Target (Table II coupling, Fig. 2 curves, EA(O)) and seed"
```

---

### Task 4: Fit O₂ (phase 1: image match), report, figure, and the acceptance gate

**Files:**
- Create: `validation/factory/fit_o2.py` (`__main__`), `validation/factory/results/o2_digitised.json` (committed output), `docs/physics/figures/o2-factory-fit.png`
- Test: `validation/factory/test_o2_fit.py`

**Interfaces:**
- Produces: `run_o2(kind: str, *, n_nodes=24, out: Path) -> tuple[FlexibleDiatomicModel, FitReport]` — `fit(o2_target(kind), o2_seed(), pair=ElectronicPair(), n_beta=3, n_nodes=n_nodes, lam_coeffs=2, alpha_coeffs=2)`, writes `out / f"o2_{kind}.json"` and, when matplotlib is importable, the overlay figure (three panels: `v0`, `v_ion` and `Γ` vs `R` — target points, fitted curves, digitisation precision as an error band; `Γ̃(ε)` at `R = 2.0` — Table II form vs model); CLI `python -m validation.factory.fit_o2 --kind digitised|tables`.

- [ ] **Step 1: Write the failing test**

```python
# validation/factory/test_o2_fit.py
from __future__ import annotations

from pathlib import Path

import pytest

from projects.potential_factory.report import FitReport, Tolerances

RESULT = Path("validation/factory/results/o2_digitised.json")


def test_committed_o2_report_meets_the_budget_or_says_which_tier_does_not():
    rep = FitReport.from_json(RESULT)
    names = [t.name for t in rep.tiers]
    assert names == ["T0", "T1", "T3"]
    assert rep.tiers[0].status == "met", rep.tiers[0].detail
    assert rep.tiers[1].status in ("met", "not met")
    assert rep.da_threshold_sign == -1  # O2 DA to O- + O is EXOTHERMIC for... (verify from Fig. 2: O(3P)+O-(2P) lies below v=0? if not, +1) -- set from the figure and cite
    assert rep.ecs_bounds_deg["tail_growth"] < 10


@pytest.mark.slow
def test_o2_fit_reproduces_the_committed_report(tmp_path):
    from validation.factory.fit_o2 import run_o2

    _, rep = run_o2("digitised", n_nodes=12, out=tmp_path)
    ref = FitReport.from_json(RESULT)
    assert [t.status for t in rep.tiers] == [t.status for t in ref.tiers]
    tol = Tolerances()
    assert rep.tiers[1].rms <= max(ref.tiers[1].rms * 1.5, tol.e_res_rms)


@pytest.mark.slow
def test_o2_fit_on_authors_tables_if_present(tmp_path):
    from validation.factory.fit_o2 import run_o2

    try:
        _, rep = run_o2("tables", n_nodes=24, out=tmp_path)
    except FileNotFoundError:
        pytest.skip("authors' tables not present")
    assert all(t.status == "met" for t in rep.tiers), [t.detail for t in rep.tiers]
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run --no-sync pytest validation/factory/test_o2_fit.py -q -m "not slow"`
Expected: FAIL (`FileNotFoundError` for the result JSON)

- [ ] **Step 3: Implement the driver**

```python
# validation/factory/fit_o2.py
"""python -m validation.factory.fit_o2 --kind digitised|tables -- fit the O2 target."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from projects.potential_factory.ansatz import FlexibleDiatomicModel
from projects.potential_factory.fit import fit, model_gamma_tilde
from projects.potential_factory.report import FitReport
from projects.potential_factory.tracker import ElectronicPair
from validation.factory.targets.o2 import o2_seed, o2_target
from validation.factory.targets.o2_data import load_o2

EV = 27.211386


def run_o2(kind: str, *, n_nodes: int = 24, out: Path) -> tuple[FlexibleDiatomicModel, FitReport]:
    target = o2_target(kind)
    pair = ElectronicPair()
    model, report = fit(target, o2_seed(), pair=pair, n_beta=3, n_nodes=n_nodes, lam_coeffs=2, alpha_coeffs=2)
    out.mkdir(parents=True, exist_ok=True)
    report.to_json(out / f"o2_{kind}.json")
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return model, report
    c = load_o2(kind)  # type: ignore[arg-type]
    R = np.linspace(c.R.min(), c.R.max(), 300)
    fig, ax = plt.subplots(2, 2, figsize=(10, 7))
    ax[0, 0].plot(c.R, c.v0 * EV, "k.", label=f"target ({c.source})")
    ax[0, 0].plot(R, model.v0(R).real * EV, "-", label="factory v0")
    ax[0, 0].fill_between(c.R, (c.v0 - c.precision_Ha) * EV, (c.v0 + c.precision_Ha) * EV, alpha=0.2)
    ax[0, 0].set(xlabel="R (bohr)", ylabel="V0 (eV)")
    ax[0, 1].plot(c.R, c.v_ion * EV, "k.", label="target V_ion")
    from projects.potential_factory.extract import extract_target

    t_fit = extract_target(model, pair=pair, R_desc=np.linspace(c.R.max(), c.R.min(), 24), n_eps=3)
    Rd = np.linspace(c.R.max(), c.R.min(), 24)
    ax[0, 1].plot(Rd, t_fit.resonance.v_ion(Rd) * EV, "-", label="factory V_ion")
    ax[0, 1].set(xlabel="R (bohr)", ylabel="V_ion (eV)")
    ax[1, 0].plot(c.R, c.gamma * EV, "k.", label="target Gamma")
    ax[1, 0].plot(Rd, t_fit.resonance.gamma(Rd) * EV, "-", label="factory Gamma")
    ax[1, 0].set(xlabel="R (bohr)", ylabel="Gamma (eV)")
    eps = np.geomspace(*target.coupling.eps_window, 30)
    ax[1, 1].loglog(eps * EV, target.coupling.gamma_tilde(eps, 2.0) * EV, "k-", label="Table II, R=2.0")
    ax[1, 1].loglog(eps * EV, model_gamma_tilde(model, pair, eps, np.array([2.0]))[:, 0] * EV, "--", label="factory, R=2.0")
    ax[1, 1].set(xlabel="eps (eV)", ylabel="Gamma~ (eV)")
    for a in ax.ravel():
        a.legend(fontsize=8)
    fig.suptitle(f"O2 factory fit ({kind}): " + ", ".join(f"{t.name} {t.status}" for t in report.tiers))
    fig.tight_layout()
    fig.savefig(Path("docs/physics/figures") / f"o2-factory-fit{'' if kind == 'digitised' else '-tables'}.png", dpi=130)
    return model, report


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kind", choices=["digitised", "tables"], default="digitised")
    ap.add_argument("--n-nodes", type=int, default=24)
    ap.add_argument("--out", type=Path, default=Path("validation/factory/results"))
    a = ap.parse_args()
    _, rep = run_o2(a.kind, n_nodes=a.n_nodes, out=a.out)
    for t in rep.tiers:
        print(f"{t.name}: {t.status}  rms={t.rms:.3e} max={t.max:.3e}  {t.detail}")
    print("crossing R =", rep.crossing_R, " DA sign =", rep.da_threshold_sign)


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run the fit, inspect, commit the report + figure**

Run: `uv run --no-sync python -m validation.factory.fit_o2 --kind digitised` (tens of minutes). Read the printed tier lines. If T1 is `not met`, that is a **result** — record the detail in Task 6's note; do not loosen the budget. Set the `da_threshold_sign` assertion in the test from Fig. 2 (and cite the page in the note), then:

Run: `uv run --no-sync pytest validation/factory/test_o2_fit.py -q -m "not slow"` → 1 passed.

```bash
uv run ruff check validation/factory && uv run ruff format validation/factory
git add validation/factory/fit_o2.py validation/factory/test_o2_fit.py validation/factory/results/o2_digitised.json docs/physics/figures/o2-factory-fit.png
git commit -m "feat(validation/factory): O2 image-match fit — report, overlay figure, acceptance gate"
```

- [ ] **Step 5: Phase 2 when the tables arrive**

Put the authors' tables under `validation/factory/data/o2/tables/{v0,v_ion,gamma}.csv` (same columns, header stating the provenance and the precision), set `_PRECISION_EV["tables"]` from what they state, run `--kind tables`, commit `results/o2_tables.json` + `docs/physics/figures/o2-factory-fit-tables.png`, and un-skip `test_o2_fit_on_authors_tables_if_present` by presence. The acceptance verdict of the spec's "Validation (b)" is *this* run.

---

### Task 5: Promotion — `qscat.model.flexible`, `qscat.factory`, and the O₂ preset

**Files:**
- Create: `libs/qscat/qscat/model/flexible.py` (move `ansatz.py`'s content), `libs/qscat/qscat/factory/__init__.py`, `libs/qscat/qscat/factory/{target,tracker,extract,fit,report}.py` (moved), `libs/qscat/tests/test_model_flexible.py`, `libs/qscat/tests/test_factory_*.py` (moved tests, with the slow round trip kept slow)
- Modify: `libs/qscat/qscat/model/__init__.py` (export `FlexibleDiatomicModel`, `SmoothR`), `libs/qscat/qscat/model/library.py` (add `O2 = FlexibleDiatomicModel(**<the fitted parameters from results/o2_<best>.json>)` with a comment naming the report and the paper), `projects/potential_factory/*.py` (thin re-export shims: `from qscat.factory.fit import *  # noqa`), `apps/qscat-run/qscat_run/presets.py` (an `O2` `MoleculePreset`), `apps/qscat-run/README.md` (registry row), `apps/qscat-run/examples/o2-ve.yaml`
- Test: `apps/qscat-run/tests/test_presets_o2.py` (or wherever the preset tests live — read the app's `tests/`)

**Interfaces:**
- Produces: `qscat.model.O2` (the fitted model — from the `tables` fit if it exists, else `digitised`, and the docstring says which), `qscat.factory.fit(...)` with the core plan's signature, and `presets` entry `O2` whose `ti_grid()` builds `TensorGrid([propose_grid(O2, "electronic", (E_min, E_max)), propose_grid(O2, "nuclear", (E_min, E_max), channel="ve")])` for the energy window `(0.002, 0.10)` Ha (Fig. 5's 0–2.5 eV), `valid_observables=frozenset({"ve", "eigenstates", "lcp"})` (add `"da"` only after the DA channel is checked open at these energies — Fig. 2 places O⁻+O relative to v=0; state it), `n_vib=12`.

- [ ] **Step 1: Move the code (git mv), fix imports, run everything**

```bash
git mv projects/potential_factory/ansatz.py libs/qscat/qscat/model/flexible.py
mkdir -p libs/qscat/qscat/factory
for f in target tracker extract fit report; do git mv projects/potential_factory/$f.py libs/qscat/qscat/factory/$f.py; done
for f in ansatz tracker target fit roundtrip report; do git mv projects/potential_factory/test_$f.py libs/qscat/tests/test_factory_$f.py; done
```

Replace every `projects.potential_factory.ansatz` import with `qscat.model.flexible`, every `projects.potential_factory.X` with `qscat.factory.X`; recreate `projects/potential_factory/{ansatz,target,tracker,extract,fit,report}.py` as one-line re-exports so `validation/factory` keeps working; add the `qscat.factory` package docstring (what it is, that it imports `qscat.model` and `qscat.core`, and that `qscat.core` must never import it).

Run: `uv run --no-sync pytest libs/qscat/tests/test_core_no_model_import.py libs/qscat/tests/test_factory_ansatz.py libs/qscat/tests/test_factory_tracker.py -q -m "not slow"` → pass; `uv run mypy libs/qscat/qscat` → clean (add annotations where strict mode complains; the moved code was written typed).

- [ ] **Step 2: Add `O2` to the library and the preset; write the failing preset test**

```python
# apps/qscat-run/tests/test_presets_o2.py
from __future__ import annotations

import pytest
from qscat.model import O2, ResonanceModel
from qscat_run import presets


def test_o2_is_registered_and_is_a_resonance_model():
    assert isinstance(O2, ResonanceModel) and O2.ell == 2
    p = presets.get("O2")  # the real accessor name from presets.py
    assert "ve" in p.valid_observables


@pytest.mark.slow
def test_o2_ti_grid_builds_and_2d_spot_check_converges():
    from qscat.tuning.refine2d import refine_to_2d_convergence  # read its signature first

    p = presets.get("O2")
    tg = p.ti_grid()
    assert tg.grids[0].n * tg.grids[1].n < 400_000
    # one nuclear h-refinement must change sigma_VE(0->1) at E=0.05 Ha by < 2 %
    rel = refine_to_2d_convergence(O2, tg, energy=0.05, vprime=1)  # adapt to the real API
    assert rel < 0.02
```

Then in `library.py`:

```python
O2 = FlexibleDiatomicModel(
    mu=O2_MU_VALUE, ell=2, D_e=..., R_e=..., betas=(...,), p=3,
    lam=SmoothR(...), alpha=SmoothR(...), shell=SmoothR(...) or None, alpha_b=..., r_b=...,
)
"""O2-like model FITTED by the potential factory to Alt & Houfek, PRA 103, 032829
(2021) -- report validation/factory/results/o2_<kind>.json; NOT a hand-tuned deck.
Tiers met: <copy from the report>."""
```

with every `...` replaced by the numbers from the committed report (`parameters` dict → constructor; a small test asserts `params(O2) == FitReport.from_json(...).parameters`).

- [ ] **Step 3: Build the preset with the tuner and run the spot-check**

Follow the `discretisation-tuner` skill: `propose_grid(O2, "nuclear", (0.002, 0.10), channel="ve")`, `propose_grid(O2, "electronic", (0.002, 0.10))`, the 1-D probes at both energy extremes, then the 2-D spot-check; if the spot-check fails, `refine_to_2d_convergence` and record the refinement in the preset's docstring. Register the preset (`molecule="O2"`, `variant="default"`, TD grid = TI grid for now, `default_incident` copied from N₂'s with the packet centred inside the electronic box, `da_test_function=None`).

Run: `uv run --no-sync pytest apps/qscat-run/tests/test_presets_o2.py -q -m "not slow"` → pass; `-m slow` → pass (Docker if needed).

- [ ] **Step 4: Example config and README row**

```yaml
# apps/qscat-run/examples/o2-ve.yaml
molecule: O2
methods: [ti]
observables:
  ve:
    v_init: 0
    vprimes: [0, 1, 2, 3]
energies:
  min: 0.002
  max: 0.10
  n: 120
artifacts:
  csv: true
  png: true
```

Run it: `uv run --no-sync qscat-run apps/qscat-run/examples/o2-ve.yaml` (or the app's entry point); commit the produced PNG as `docs/physics/figures/o2-2d-ti-cross-section.png`. Add the README registry row (`O2 | 0 | ve, eigenstates, lcp | default`).

- [ ] **Step 5: Commit**

```bash
uv run ruff check . && uv run ruff format . && uv run mypy libs/qscat/qscat
git add -A libs/qscat/qscat/model libs/qscat/qscat/factory libs/qscat/tests projects/potential_factory apps/qscat-run docs/physics/figures/o2-2d-ti-cross-section.png
git commit -m "feat(qscat): promote the potential factory (qscat.factory, qscat.model.flexible); O2 fitted model + tuner-built preset"
```

(`git add -A <paths>` is scoped to those paths; check `git status` first that nothing else is staged.)

---

### Task 6: Physics review and documentation

**Files:**
- Modify: `docs/physics/potential-factory.md` (final: promoted paths, the O₂ results, the phase-1/phase-2 distinction, the limits found), `docs/molecules/` (add an O₂ page in the same form as the existing molecule pages — read one first), `CLAUDE.md` (move the `projects/potential_factory` bullet to a `qscat.factory` entry under `libs/`, add `O2` to the `qscat.model` sentence), `docs/physics/potential-factory-options.md` (status line → "implemented as `qscat.factory`; this note remains the survey")

- [ ] **Step 1: Dispatch the physics reviewer**

`Agent: physics-reviewer` over `libs/qscat/qscat/factory/`, `libs/qscat/qscat/model/flexible.py`, `validation/factory/targets/o2.py`, with the spec's checklist: units; the `Γ`-support condition on the fitted O₂ (`Γ ≠ 0` only where `V_0 < V_ion` — check on the report's crossing); the DA threshold sign vs Fig. 2; ECS bounds; the c-product gradient's sign convention; that nothing reads an experimental cross section. Fix what it finds; record anything deliberately left in the note's Limitations.

- [ ] **Step 2: Write the docs**

`docs/physics/potential-factory.md` gains: the O₂ section with the tier table from `o2_digitised.json` (and `o2_tables.json` when present), the figure, the sentence "phase 1 is an image match at ±0.03 eV; the benchmark is phase 2", the measured wall time, and the limits found (which tier could not be met and why, in the 2-D form's terms). The molecule page lists the observables run and links the cross-section figure. CLAUDE.md: the `libs/` entry (6–10 lines in the house style) and the `qscat.model` sentence.

- [ ] **Step 3: Commit; then prepare the branch for review**

```bash
git add docs/physics/potential-factory.md docs/physics/potential-factory-options.md docs/molecules CLAUDE.md
git commit -m "docs(factory): O2 results, molecule page, repo map; physics review applied"
```

Then run the `review-ready` skill (main must stand alone: no plan/spec references from shipped files except the `docs/superpowers` cross-links the house style allows; `tidy-history` if the branch is fix-on-fix).

---

## Self-review against the spec

- **Data acquisition, two phases, digitisation precision as an uncertainty floor, tables from the authors, NIST fallback, locators via `mastering-references`** → Tasks 1–3, Task 4 Step 5.
- **Validation (b)** O₂ acceptance against the budget → Task 4 (the gate is the budget-set `Tolerances`; a `not met` tier is committed as a result, never loosened). **(c)** tuner grids + 2-D spot-check → Task 5. **(d)** physics-reviewer → Task 6.
- **Placement**: `qscat.factory` sibling of `qscat.model`, `FlexibleDiatomicModel` in `qscat.model`, `validation/factory/`, `projects/` shim → Task 5.
- **`MoleculePreset` / `qscat-run` by registry name** → Task 5.
- **Never experiment** → Global Constraints; Figs. 5–9 unread by code.
- Unknown-name flags left deliberately for the executor to read the source: the presets accessor, `refine_to_2d_convergence`'s signature, the app's test directory, `units.py`'s eV constant, and `da_threshold_sign` for O₂ (set from Fig. 2 with a citation — the plan does not assert physics it has not checked).
