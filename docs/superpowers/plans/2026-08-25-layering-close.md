# Layering Close (Phase 1.4) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close the last seven `projects/` → `validation/` layering violations (two direct imports, four `config.json` path traversals, one anchor-coordinate dependency), then lock the boundary with an enforcement test so a new violation fails CI instead of accumulating.

**Architecture:** No new module is needed as the layer-neutral home for the N₂ parameters — it already exists: `qscat.model.N2` (`libs/qscat/qscat/model/library.py:28-40`) carries exactly the contents of `validation/n2/config.json` (`reduced_mass` → `N2.mu`, `impulsemomentum` → `N2.ell`, the ten `potential` values → the ten dataclass fields), and `libs/qscat/tests/test_model.py` already gates those literals against an independent transcription of the same deck constants. Every `projects/` read of `config.json` therefore becomes a read of `qscat.model.N2` (the allowed direction: projects → qscat). The Houfek anchor coordinates stay on the **validation** side (see Task 4 for the justification), and the projects-side test that consumed them is retired in favor of a validation-side gate over the already-general `validation/n2/cross_section.py` classification.

**Tech Stack:** Python 3.12 / uv, pytest (fast tier), `ast` for the import-graph enforcement test, numpy.

**Spec:** the "Findings addressed" section below (self-contained; from the 2026-08-25 release review)

**Sequencing:** This plan runs AFTER the kernel-consolidation and api-surface plans (same tree). File-overlap check: this plan touches only `projects/n2_resonance/{potential.py,test_potential.py}`, `projects/n2_ti_cross_section/{test_vibrational.py,test_cross_section.py}`, `projects/n2_td_cross_section/test_td_cross_section.py`, `validation/n2/{cross_section.py,test_anchor_gate.py}`, and `tests/test_layering.py` — the kernel-consolidation plan works in `native/` + the `libs/qscat` kernel mirrors, and the api-surface plan works on the `libs/qscat` public API; none of those files appear here, so the three plans do not overlap (re-verify against their final file lists when they land, since all three are being written concurrently). The experiment-lifecycle plan (2026-08-25) runs after THIS plan and deliberately builds on its end state (it later deletes `projects/n2_ti_cross_section/{cross_section.py,test_cross_section.py}`, which this plan edits first).

## Global Constraints

PyPI release DEFERRED until the peer-reviewed article publishes — repo-only distribution, no publishing tasks. After every task: `uv run --no-sync pytest -m "not slow" -n auto --dist loadfile` green; `uv run --no-sync mypy libs/qscat/qscat apps/qscat-run/qscat_run` clean; `uv run --no-sync ruff check .` + `ruff format --check .` clean. Physics-bearing moves are IDENTITY-PRESERVING (differential test pinning old vs new output before deleting the old). Tasks touching validation/n2 or validation/diatomic solver paths need a `validate:n2` / `validate:diatomic` labelled run before merge. Never `git commit -a`. Layering rule: **validation may import projects and qscat; projects may import qscat only; qscat imports neither.**

## Findings addressed

> **exp-C1 (remainder).** The code still violates the one-directional rule in seven places: (a) `projects/n2_resonance/test_potential.py` imports `validation.n2.model`; (b) `projects/n2_ti_cross_section/test_cross_section.py` imports `validation.n2.loader` and `validation.n2.reference` (`ANCHOR_COORDS`); (c) four modules read `validation/n2/config.json` by relative path traversal (`projects/n2_resonance/potential.py`, `projects/n2_ti_cross_section/test_vibrational.py` + `test_cross_section.py`, `projects/n2_td_cross_section/test_td_cross_section.py`). Design the layer-neutral home and plan the moves; final task: an enforcement test that fails on any projects→validation import, including config-path traversal.

Verified locations (2026-08-25 tree):

| # | Violation | Where |
|---|---|---|
| 1 | `from validation.n2 import model as ref_model` | `projects/n2_resonance/test_potential.py:11` |
| 2 | `from validation.n2 import loader` | `projects/n2_ti_cross_section/test_cross_section.py:37` |
| 3 | `from validation.n2.reference import ANCHOR_COORDS` | `projects/n2_ti_cross_section/test_cross_section.py:38` |
| 4 | `Path(__file__).resolve().parents[2] / "validation" / "n2" / "config.json"` | `projects/n2_resonance/potential.py:39-41` |
| 5 | same traversal | `projects/n2_ti_cross_section/test_vibrational.py:30-32` |
| 6 | same traversal | `projects/n2_ti_cross_section/test_cross_section.py:40-42` |
| 7 | same traversal | `projects/n2_td_cross_section/test_td_cross_section.py:46-48` |

**Design decision — the layer-neutral parameter home is `qscat.model.N2`.** `config.json` holds `reduced_mass: 12766.36`, `impulsemomentum: 2`, and a ten-value `potential` block; `qscat.model.N2` holds the identical values as dataclass fields (`mu`, `ell`, `D0`, `alpha0`, `R0`, `lambda_inf`, `lambda_1`, `R_lambda`, `lambda_c`, `R_c`, `alpha_c`), and `libs/qscat/tests/test_model.py` locks them to an independent transcription of the deck. `config.json` itself stays where it is, as validation-side provenance (its `provenance`/`note` fields document the eMoScat origin); nothing in `projects/` reads it after this plan.

**Design decision — `ANCHOR_COORDS` stay in `validation` (golden-data-side); the projects-side test is retired, not relocated.** Read of the two imports: `test_cross_section.py::test_houfek_anchor_agreement` uses `loader.load()` to fetch Houfek's `CSVE.V00.J00` table and `ANCHOR_COORDS` to index into it, then gates 4 of 6 anchors at `ANCHOR_FACTOR = 3` with the other two hardcoded in `_KNOWN_MODEL_LIMITATION_ANCHORS`. The anchor coordinates are meaningless without the golden dataset they index (which lives under `validation/n2/data/`), so a "neutral home" would strand coordinates away from their data — they belong with it. And the projects test itself is the *less general* duplicate of machinery validation already owns: `validation/n2/cross_section.py::classify` (lines 114-133) derives the same GATED/DOCUMENTED-LIMITED split generally from `(energy, channel)` — its own docstring (lines 106-151 of the projects test, lines 19-41 of the validation module) says the hardcoded two-coordinate set is the thing `classify` supersedes. What is missing on the validation side is only a *pytest* gate over `compute_anchor_results()` (today it is exercised by the manual `python -m validation.n2.experiment` harness and the `@slow` exact-2D `test_anchors.py`). So: add that fast-tier gate in validation (Task 4a), then delete the projects duplicate (Task 4b). Coverage is preserved, generality improves, and violations 2, 3 and 6 disappear together.

---

## Task 1 — Re-anchor `projects/n2_resonance/test_potential.py` on `qscat.model.N2`

Removes violation 1. The cross-check the test performs (two independent N₂ potential implementations agree to 1e-12) is preserved; only the oracle changes from `validation.n2.model` to the layer-legal `qscat.model.N2` — which is itself locked to the same constants by `libs/qscat/tests/test_model.py`, so the check is not weakened. Note `potential.v_eff_el` EXCLUDES `v0` while `N2.surface` INCLUDES it, so the v_eff_el comparison is against `N2.surface(r, R) - N2.v0(R)`.

**Files:**
- Modify: `projects/n2_resonance/test_potential.py`
- Test: `projects/n2_resonance/test_potential.py` (self)

**Interfaces:**
- Consumes: `qscat.model.N2` (`DiatomicResonanceModel` instance: `.v0(R) -> NDArray[complex128]`, `.lam(R) -> NDArray[complex128]`, `.v_int(r, R) -> NDArray[complex128]`, `.surface(r, R) -> NDArray[complex128]`); `projects.n2_resonance.potential` (`v0`, `lam`, `v_int`, `v_eff_el`, `PARAMS`).
- Produces: no API change — a test-only edit.

**Steps:**

- [ ] Edit `projects/n2_resonance/test_potential.py`: delete line 11 (`from validation.n2 import model as ref_model`), add `from qscat.model import N2`, and replace `test_matches_reference_model_to_1e_12` with:

  ```python
  def test_matches_library_model_to_1e_12():
      """Cross-check against qscat.model.N2 -- the same eMoScat constants,
      independently implemented (DiatomicResonanceModel) and independently
      gated (libs/qscat/tests/test_model.py). potential.v_eff_el excludes
      v0(R) while N2.surface includes it, hence the subtraction."""
      rs = np.array([0.5, 1.0, 2.0, 3.0, 5.0, 10.0])
      Rs = np.array([1.5, 2.01943, 2.405, 3.0, 4.0])
      for R in Rs:
          assert abs(potential.v0(R) - complex(N2.v0(R))) < 1e-12
          assert abs(potential.lam(R) - complex(N2.lam(R))) < 1e-12
          for r in rs:
              assert abs(potential.v_int(r, R) - complex(N2.v_int(r, R))) < 1e-12
              assert abs(
                  potential.v_eff_el(r, R) - (complex(N2.surface(r, R)) - complex(N2.v0(R)))
              ) < 1e-12
  ```

  Update the module docstring (lines 1-6): the cross-check target is now `qscat.model.N2`, not `validation/n2/model.py`.
- [ ] Run: `uv run --no-sync pytest projects/n2_resonance/test_potential.py -q` — expect all 6 tests PASS (values are equal; only the oracle import moved).
- [ ] Run: `grep -rn "validation" projects/n2_resonance/test_potential.py` — expect no matches.
- [ ] Run the full gate (fast tier + mypy + ruff, per Global Constraints).
- [ ] Commit: `git add projects/n2_resonance/test_potential.py && git commit -m "test(n2-resonance): cross-check the potential against qscat.model.N2, not validation"`

---

## Task 2 — `projects/n2_resonance/potential.py`: parameters from `qscat.model.N2`, not `config.json`

Removes violation 4. `potential.py`'s `PARAMS` dict shape is kept (its only in-repo consumer is `test_potential.py`, which reads `PARAMS["potential"][...]` and `PARAMS["impulsemomentum"]`), but it is now BUILT from the `N2` dataclass fields instead of parsed from `validation/n2/config.json`. The formulas themselves are untouched here — full delegation of the function bodies to `N2` is the experiment-lifecycle plan's exp-M2 task, which runs after this plan.

**Files:**
- Modify: `projects/n2_resonance/potential.py`
- Test: `projects/n2_resonance/test_potential.py` (unchanged from Task 1), `projects/n2_resonance/test_pole.py`, `projects/n2_ti_cross_section/test_vres.py` (existing consumers exercise the values)

**Interfaces:**
- Consumes: `qscat.model.N2` fields (`mu: float`, `ell: int`, `D0`, `alpha0`, `R0`, `lambda_inf`, `lambda_1`, `R_lambda`, `lambda_c`, `R_c`, `alpha_c`: all `float`).
- Produces: `projects.n2_resonance.potential.PARAMS: dict` (same keys/values as before: `"reduced_mass"`, `"impulsemomentum"`, `"potential": {"D_0", "alpha_0", "R_0", "lambda_inf", "lambda_1", "R_lambda", "lambda_c", "R_c", "alpha_c"}`); `v0`/`lam`/`v_int`/`v_eff_el` signatures unchanged.

**Steps:**

- [ ] In `projects/n2_resonance/potential.py`, delete the `import json` / `from pathlib import Path` lines and the `PARAMS = json.loads((Path(__file__).resolve().parents[2] / "validation" / "n2" / "config.json").read_text())` block (lines 34-41); replace with:

  ```python
  from qscat.model import N2

  # The N2 deck constants, read from the layer-neutral single source
  # `qscat.model.N2` (locked to the eMoScat deck by
  # libs/qscat/tests/test_model.py). The dict shape mirrors the historical
  # validation/n2/config.json layout so existing PARAMS consumers are
  # unaffected.
  PARAMS: dict = {
      "reduced_mass": N2.mu,
      "impulsemomentum": N2.ell,
      "potential": {
          "D_0": N2.D0,
          "alpha_0": N2.alpha0,
          "R_0": N2.R0,
          "lambda_inf": N2.lambda_inf,
          "lambda_1": N2.lambda_1,
          "R_lambda": N2.R_lambda,
          "lambda_c": N2.lambda_c,
          "R_c": N2.R_c,
          "alpha_c": N2.alpha_c,
      },
  }
  ```

  Update the module docstring paragraph (lines 3-11): parameters now come from `qscat.model.N2`; the 1e-12 cross-check named there is `test_matches_library_model_to_1e_12` against `qscat.model.N2`.
- [ ] Run: `uv run --no-sync pytest projects/n2_resonance/ projects/n2_ti_cross_section/test_vres.py -q -m "not slow"` — expect PASS (identical values, different source).
- [ ] Run the full gate.
- [ ] Commit: `git add projects/n2_resonance/potential.py && git commit -m "refactor(n2-resonance): read the N2 deck constants from qscat.model.N2, not validation/n2/config.json"`

---

## Task 3 — Replace the `config.json` traversal in the three project test modules

Removes violations 5 and 7, and the traversal half of 6 (the import half of 6 falls in Task 4). Each module reads only scalars that exist on `N2`: `MU = _CONFIG["reduced_mass"]` → `N2.mu`; `test_vibrational.py` additionally reads `D_0` → `N2.D0` and `alpha_0` → `N2.alpha0`.

**Files:**
- Modify: `projects/n2_ti_cross_section/test_vibrational.py`, `projects/n2_ti_cross_section/test_cross_section.py`, `projects/n2_td_cross_section/test_td_cross_section.py`
- Test: the three files themselves

**Interfaces:**
- Consumes: `qscat.model.N2.mu -> float` (12766.36), `N2.D0 -> float` (0.75102), `N2.alpha0 -> float` (1.1535).
- Produces: no API change — test-only edits; module constants `MU`, `D0`, `ALPHA0` keep their names and values.

**Steps:**

- [ ] In `projects/n2_ti_cross_section/test_vibrational.py`: delete the `import json` / `from pathlib import Path` imports and lines 30-35 (`_CONFIG = ...`, `MU = ...`, `D0 = ...`, `ALPHA0 = ...`); add `from qscat.model import N2` and

  ```python
  MU = N2.mu  # N2 nuclear reduced mass (a.u.), 12766.36
  D0 = N2.D0  # 0.75102 Ha
  ALPHA0 = N2.alpha0  # 1.1535 bohr^-1
  ```

- [ ] In `projects/n2_ti_cross_section/test_cross_section.py`: delete `import json` / `from pathlib import Path` and lines 40-43 (`_CONFIG` + `MU`); add `from qscat.model import N2` and `MU = N2.mu`. (Leave lines 37-38's validation imports in place — Task 4 removes them together with the test that uses them.)
- [ ] In `projects/n2_td_cross_section/test_td_cross_section.py`: same edit — delete `import json` / `from pathlib import Path` and lines 46-49, add `from qscat.model import N2` and `MU = N2.mu`.
- [ ] Run: `uv run --no-sync pytest projects/n2_ti_cross_section/test_vibrational.py projects/n2_ti_cross_section/test_cross_section.py projects/n2_td_cross_section/test_td_cross_section.py -q -m "not slow"` — expect PASS.
- [ ] Run: `grep -rn '"validation"' projects/` — expect the only remaining hit to be `projects/n2_ti_cross_section/test_cross_section.py`'s import lines (gone after Task 4) and docstring prose.
- [ ] Run the full gate.
- [ ] Commit: `git add projects/n2_ti_cross_section/test_vibrational.py projects/n2_ti_cross_section/test_cross_section.py projects/n2_td_cross_section/test_td_cross_section.py && git commit -m "test: read the N2 constants from qscat.model.N2 instead of traversing into validation/"`

---

## Task 4 — Move the Houfek anchor gate to validation; retire the projects duplicate

Removes violations 2, 3, and the last of 6. Order matters for coverage: the validation-side gate lands and passes FIRST, then the projects-side duplicate is deleted.

### Task 4a — validation-side anchor gate over `compute_anchor_results()`

**Files:**
- Create: `validation/n2/test_anchor_gate.py`
- Test: itself

**Interfaces:**
- Consumes: `validation.n2.cross_section.compute_anchor_results() -> list[AnchorResult]` (fields: `energy_ha: float`, `channel: int`, `sigma_computed: float`, `sigma_houfek: float`, `ratio: float`, `gated: bool`, `mechanism: str`); `validation.n2.reference.ANCHOR_FACTOR: float` (3.0).
- Produces: the fast-tier pytest gate that previously lived (in less general form) in `projects/n2_ti_cross_section/test_cross_section.py::test_houfek_anchor_agreement`.

**Steps:**

- [ ] Create `validation/n2/test_anchor_gate.py`:

  ```python
  """The C5 Houfek anchor gate, as a pytest test.

  Replaces `projects/n2_ti_cross_section/test_cross_section.py`'s
  `test_houfek_anchor_agreement`, which hardcoded the two DOCUMENTED-LIMITED
  coordinates; here the GATED / DOCUMENTED-LIMITED split comes from
  `validation.n2.cross_section.classify`'s general `(energy, channel)` rule,
  so the gate and the `python -m validation.n2.experiment` harness cannot
  disagree. `_build_system` is lru_cached, so the ~7 s vres_on_grid walk is
  paid once per process however many tests read it.
  """

  from __future__ import annotations

  from validation.n2 import reference
  from validation.n2.cross_section import compute_anchor_results


  def test_gated_anchors_within_factor_band():
      results = compute_anchor_results()
      assert len(results) == 6
      gated = [r for r in results if r.gated]
      assert len(gated) == 4  # the general rule must still gate exactly 4 of 6
      for r in gated:
          assert 1.0 / reference.ANCHOR_FACTOR <= r.ratio <= reference.ANCHOR_FACTOR, (
              f"anchor (E={r.energy_ha}, v'={r.channel}) ratio {r.ratio:.3f} "
              f"outside factor-of-{reference.ANCHOR_FACTOR} band"
          )


  def test_documented_limited_anchors_carry_a_mechanism():
      for r in compute_anchor_results():
          if not r.gated:
              assert r.mechanism  # never silently excluded
          print(
              f"E={r.energy_ha:.4f} Ha, v'={r.channel}: computed={r.sigma_computed:.4e} "
              f"houfek={r.sigma_houfek:.4e} ratio={r.ratio:.3f} "
              f"{'GATED' if r.gated else '[' + r.mechanism + ']'}"
          )
  ```

- [ ] Run: `uv run --no-sync pytest validation/n2/test_anchor_gate.py -q` — expect 2 PASS (the 4 gated anchors sit at ratios ~0.8-1.2 today).
- [ ] Commit: `git add validation/n2/test_anchor_gate.py && git commit -m "test(validation/n2): gate the C5 Houfek anchors via the general classification"`

### Task 4b — delete the projects-side duplicate and its validation imports

**Files:**
- Modify: `projects/n2_ti_cross_section/test_cross_section.py`, `validation/n2/cross_section.py`
- Test: `projects/n2_ti_cross_section/test_cross_section.py`, `validation/n2/test_anchor_gate.py`

**Interfaces:**
- Consumes: nothing new.
- Produces: `projects/n2_ti_cross_section/test_cross_section.py` keeps exactly its three INTERNAL tests (`test_sigma_real_and_nonnegative`, `test_closed_channel_is_exactly_zero`, `test_v0_to_v1_resonance_enhancement`) plus the module fixture; no validation imports remain anywhere in `projects/`.

**Steps:**

- [ ] In `projects/n2_ti_cross_section/test_cross_section.py`, delete: lines 37-38 (the two `validation.n2` imports), the `ANCHOR_FACTOR` constant (lines 45-49), the `_KNOWN_MODEL_LIMITATION_ANCHORS` block with its long comment (lines 106-151), and `test_houfek_anchor_agreement` (lines 154-187). Trim the module docstring's "HOUFEK anchors" bullet to a pointer: "the Houfek anchor comparison lives in `validation/n2/test_anchor_gate.py` (validation may import projects; not the reverse)."
- [ ] In `validation/n2/cross_section.py`, rewrite the docstring's "KNOWN EXCEPTIONS" sentence (lines 10-17): the exceptions no longer exist — state that the boundary is now clean and enforced by `tests/test_layering.py` (landed in Task 5; forward reference is fine within one branch).
- [ ] Run: `uv run --no-sync pytest projects/n2_ti_cross_section/test_cross_section.py validation/n2/test_anchor_gate.py -q -m "not slow"` — expect PASS (3 + 2 tests).
- [ ] Run: `grep -rn "validation" projects/ --include="*.py"` — read the hits: every remaining one must be docstring/comment prose; no `import` lines, no path fragments in code.
- [ ] Run the full gate. This task touches `validation/n2` — run the `validate:n2` labelled slow tier before merging the branch.
- [ ] Commit: `git add projects/n2_ti_cross_section/test_cross_section.py validation/n2/cross_section.py && git commit -m "test: retire the projects-side Houfek anchor duplicate in favor of the validation gate"`

---

## Task 5 — Enforcement test: the boundary can never silently regress

The import-graph test that fails on any `projects/` → `validation/` import AND on any path-string traversal into `validation/` from non-docstring code. It lands last (the tree must be green after every task, and the violations it polices were just removed), so its detection power is demonstrated with a temporary canary file rather than against real violations.

**Files:**
- Create: `tests/test_layering.py` (the existing top-level `tests/` package — home of `test_docs_portability.py` etc.)
- Test: itself

**Interfaces:**
- Consumes: the filesystem under `projects/` (every `*.py`), Python `ast`.
- Produces: `tests/test_layering.py::test_projects_never_import_validation` and `::test_projects_never_reference_validation_paths`.

**Steps:**

- [ ] Create `tests/test_layering.py`:

  ```python
  """Layering enforcement: projects/ must not depend on validation/.

  The rule (CLAUDE.md, docs/adr): validation may import projects and qscat;
  projects may import qscat only; qscat imports neither. The qscat and
  qscat_run sides are already enforced (test_core_no_model_import.py,
  apps/qscat-run/tests/test_no_validation_import.py); this test closes the
  projects side, in BOTH forms the 2026-08-25 release review found in the
  wild: a literal `import validation...`, and a filesystem traversal into
  `validation/` via a path string (the config.json pattern). Docstrings may
  mention validation/ in prose; string literals in CODE may not.
  """

  from __future__ import annotations

  import ast
  from pathlib import Path

  REPO = Path(__file__).resolve().parents[1]
  PROJECTS = REPO / "projects"


  def _project_files() -> list[Path]:
      files = sorted(PROJECTS.rglob("*.py"))
      assert files, f"no Python files found under {PROJECTS}"
      return files


  def _docstring_nodes(tree: ast.Module) -> set[int]:
      """id()s of every Constant node that is a module/class/function docstring."""
      out: set[int] = set()
      for node in ast.walk(tree):
          if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
              body = node.body
              if (
                  body
                  and isinstance(body[0], ast.Expr)
                  and isinstance(body[0].value, ast.Constant)
                  and isinstance(body[0].value.value, str)
              ):
                  out.add(id(body[0].value))
      return out


  def test_projects_never_import_validation() -> None:
      offenders: list[str] = []
      for path in _project_files():
          tree = ast.parse(path.read_text(), filename=str(path))
          for node in ast.walk(tree):
              if isinstance(node, ast.Import):
                  if any(
                      a.name == "validation" or a.name.startswith("validation.")
                      for a in node.names
                  ):
                      offenders.append(f"{path.relative_to(REPO)}:{node.lineno}")
              elif isinstance(node, ast.ImportFrom):
                  mod = node.module or ""
                  if node.level == 0 and (mod == "validation" or mod.startswith("validation.")):
                      offenders.append(f"{path.relative_to(REPO)}:{node.lineno}")
      assert not offenders, (
          "projects/ imports validation/ (forbidden direction; move the "
          f"dependency into validation/ or qscat): {offenders}"
      )


  def test_projects_never_reference_validation_paths() -> None:
      offenders: list[str] = []
      for path in _project_files():
          tree = ast.parse(path.read_text(), filename=str(path))
          doc_ids = _docstring_nodes(tree)
          for node in ast.walk(tree):
              if (
                  isinstance(node, ast.Constant)
                  and isinstance(node.value, str)
                  and id(node) not in doc_ids
                  and ("validation/" in node.value or node.value == "validation")
              ):
                  offenders.append(f"{path.relative_to(REPO)}:{node.lineno} {node.value!r}")
      assert not offenders, (
          "projects/ code contains a path string reaching into validation/ "
          "(the config.json-traversal pattern; read qscat.model instead): "
          f"{offenders}"
      )
  ```

- [ ] Run: `uv run --no-sync pytest tests/test_layering.py -q` — expect 2 PASS (Tasks 1-4 removed every violation).
- [ ] Canary check (proves the test detects both violation forms; leaves no trace):
  - `printf 'import validation.n2\nP = "validation/n2/config.json"\n' > projects/_layering_canary.py`
  - `uv run --no-sync pytest tests/test_layering.py -q` — expect **2 FAILED**, each assertion message naming `projects/_layering_canary.py:1` / `:2`.
  - `rm projects/_layering_canary.py`
  - `uv run --no-sync pytest tests/test_layering.py -q` — expect 2 PASS again.
- [ ] Run the full gate.
- [ ] Commit: `git add tests/test_layering.py && git commit -m "test: enforce the projects->validation layering boundary (imports and path traversal)"`

---

## Completion checklist

- [ ] `grep -rn "from validation\|import validation" projects/ --include="*.py"` → no matches.
- [ ] `grep -rn '"validation"' projects/ --include="*.py"` → no matches outside docstrings.
- [ ] Fast tier, mypy, ruff, ruff-format all clean (Global Constraints).
- [ ] `validate:n2` labelled slow-tier run green (Task 4 touched `validation/n2`).
- [ ] All seven findings-table rows resolved; `tests/test_layering.py` green with the canary check demonstrated in the task log.
