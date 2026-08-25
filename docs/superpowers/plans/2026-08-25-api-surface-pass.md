# API Surface Pass Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Freeze the public surface the eventual release ships: give `ScatteringProblem` real typed signatures and full observable coverage, resolve the two intra-repo name collisions, replace free-form `str` parameters with `Literal` types, type the tuner's `incident` parameter, make `td_ve_cross_section`'s `wp_out` optional the way `td_da_cross_section`'s already is, write the argument-order convention down as an ADR, and apply ADR 0004's provisional marker to every wide functional solver signature.

**Architecture:** No numerics change anywhere in this plan — every task is typing, naming, wiring, or documentation. The facade (`qscat.core.problem.ScatteringProblem`) mirrors the functional solvers' parameter names, defaults, and `@overload` sets exactly, and delegates positionally; new observables reach the facade the same way (`lcp_da_cross_section`, `resonance_levels`, `exact_resonance_states`, `nrm_ve_cross_section`, `nrm_da_cross_section`), with `qscat.core.nrm` imported *deferred inside the method bodies* so `import qscat.core` still never pulls `nrm` in (the hard boundary `qscat.core.__init__` documents). Renamed public names keep a one-cycle deprecated alias via module-level `__getattr__`. Literals are the mypy/IDE layer; every existing runtime `raise ValueError` validation path stays.

**Tech Stack:** Python 3.12, numpy/scipy, pytest (`uv run --no-sync pytest`), mypy --strict, ruff. No new dependencies, no Rust changes, no Docker changes.

**Spec:** the "Findings addressed" section below (self-contained; from the 2026-08-25 release review)

## Global Constraints

PyPI release DEFERRED until the peer-reviewed article publishes — repo-only distribution, no publishing tasks; this pass is still what freezes the surface the eventual release ships. Executes AFTER the kernel-consolidation plan; every task on a shared file begins with a re-verify step. After every task: `uv run --no-sync pytest -m "not slow" -n auto --dist loadfile` green; `uv run --no-sync mypy libs/qscat/qscat apps/qscat-run/qscat_run` clean; `uv run --no-sync ruff check .` + `ruff format --check .` clean; `uv run --no-sync pytest tests/test_api_docs_coverage.py -q` green (it gates docs/api against `__all__` — renames must update docs/api/ in the same task). Renamed public names keep a deprecated re-export for one cycle via module-level `__getattr__` emitting DeprecationWarning (show the shim code fully once, in the first task that uses it). Never `git commit -a`.

**The standard gate** (referenced by every task below, run verbatim):

```bash
uv run --no-sync pytest -m "not slow" -n auto --dist loadfile
uv run --no-sync mypy libs/qscat/qscat apps/qscat-run/qscat_run
uv run --no-sync ruff check . && uv run --no-sync ruff format --check .
uv run --no-sync pytest tests/test_api_docs_coverage.py -q
```

Two standing notes for the executor:

- **Formatter has final say on layout.** Every code snippet in this plan is semantically exact; after editing a file run `uv run --no-sync ruff format <file>` before the gate so `ruff format --check .` passes.
- **Sequencing with sibling 2026-08-25 plans.** This plan runs after `2026-08-25-kernel-consolidation.md` and its tasks reference the *post-consolidation* shapes of `core/time_dependent.py`, `core/td_extractors.py`, and `core/dissociation.py` (re-verify steps below confirm them). Three things this plan deliberately does NOT own: the public `Ordering` Literal export (`2026-08-25-library-structure.md`, lib-M12 — Task 6 there deletes the private copies; our facade adds its own private copy with a pointer comment, and its Task 6 says to adapt if lib-M12 landed first), the `core/lcp.py` package split (same plan), and the `lcp_ve_cross_section` graduation into `qscat.core.lcp` (`2026-08-25-experiment-lifecycle.md` — our facade leaves a marked extension point, Task 8).

## Findings addressed

From the 2026-08-25 release review. Each was re-verified against the current tree on 2026-08-25; verification notes in brackets.

- **lib-C2**: `ScatteringProblem` (core/problem.py) is the documented "recommended, stable" entry point but all methods are `(..., **kwargs: Any) -> Any`, and it covers ~6 of ~15 observables. Give every existing method a real typed signature (mirror the functional solver's parameters — read each one; where the functional signature is wide, the facade should expose the same names/types, not `**kwargs`) and real return types; then ADD the missing observables as typed methods: `lcp_da_cross_section`, `resonance_levels`, `lcp_resonance_levels`, `exact_resonance_states`, `nrm_ve_cross_section`, `nrm_da_cross_section` (decide which belong on the facade vs are legitimately expert-only — justify each exclusion in the plan; the docstring must stop calling the facade complete if you exclude any). Tests: a signature-coverage test asserting no facade method has `**kwargs`, plus per-method delegation tests on a tiny grid (the pattern test_scattering_problem.py already uses post-PR#39).
  [VERIFIED: `core/problem.py:88-130` — all seven methods are `**kwargs: Any) -> Any`. Functional signatures read in full: `driven.ve_cross_section` (2 overloads + impl, `ordering`/`lam_scale`/`return_wavefunction`), `dissociation.da_cross_section` (2 overloads, `n_channels=1`), `dissociation.dr_cross_section` (4 overloads, `n_channels=3`, `return_amplitude`), the four `td_*` (keyword-only `dt`/`n_steps`/`wp_in`/…). Facade scope decision and exclusion justifications are in Task 8's design notes.]
- **lib-C3**: rename one of each colliding pair. Decide the names after reading both (suggested: `qscat.core.nrm.scattering.free_hamiltonian` → `electronic_free_hamiltonian`; `qscat.tuning.resonance.resonance_curve` → `resonance_curve_arrays` — but choose what reads best in each module's vocabulary and say why). Deprecated aliases one cycle; update all callers + docs/api + any docs/physics mentions (grep).
  [VERIFIED: both collisions carry explicit `NAME COLLISION` docstring notes — `time_dependent.py:261` ↔ `nrm/scattering.py:39`, and `bo.py:233` ↔ `tuning/resonance.py:130`. Callers mapped: nrm's `free_hamiltonian` is used only inside `nrm/scattering.py:98` and `libs/qscat/tests/test_nrm_scattering.py` (never re-exported by `nrm/__init__`, never documented in docs/api); tuning's `resonance_curve` is used by `tuning/__init__.py:84,105`, `tuning/propose.py:36,312` (+ 4 docstring mentions), `libs/qscat/tests/test_tuning_propose.py:18`, `test_tuning_resonance.py:41,45`, `docs/api/tuning.md:32`, `docs/physics/discretisation-tuning.md:170`. Name decisions and rationale in Tasks 4 and 5.]
- **lib-M14**: `Literal` the free-form str params — `method` ("tw"|"delta"|"flow") and `axis` ("electronic"|"nuclear") in core/time_dependent.py + core/td_extractors.py; `channel` in tuning/propose.py (match the existing Coordinate Literal style); `verdict` in core/assignment.py (the seven documented values — read :36-46 and make the Literal the single source, with the docstring pointing at it). Keep runtime validation (the raise paths) — the Literal is for mypy/IDE, the raise for runtime strings.
  [VERIFIED: `method: str = "tw"` / `"flow"` in `td_ve_cross_section`/`td_da_cross_section`; `axis: str = "electronic"` in all three extractor `__init__`s, `_AXES = ("electronic", "nuclear")` at `td_extractors.py:167`; `channel: str = "ve"` at `propose.py:360` with the runtime check at `:447-450` and `Coordinate = Literal[...]` at `:43` as the style precedent; the seven verdicts (`ok`, `spurious`, `basis-limited`, `box-limited`, `weak`, `mixed`, `distant`) in the table at `assignment.py:36-46` and the assignments at `:354-371`, field `verdict: str` at `:154`. One extra consumer found: `apps/qscat-run/qscat_run/runner.py`'s `_build_extractor(..., axis: str, ...)` (`:474`) passes `axis` into the extractor constructors, so it must adopt the Literal in the same task or mypy fails.]
- **lib-M15**: tuning/propose.py `incident: object | None` → `IncidentSpec | None` (TYPE_CHECKING import per the file's existing ResonanceModel pattern); replace the silent `getattr(..., lambda: 0.0)` reads with direct attribute access; check nothing but IncidentSpec is ever passed (grep call sites).
  [VERIFIED: `propose.py:358` and the two getattr reads at `:466-469`. Repo-wide grep of `propose_grid(` call sites (libs, apps, projects, validation, benchmarks): the only callers passing `incident=` are in `libs/qscat/tests/test_tuning_incident.py`, and every one passes a real `IncidentSpec`. `IncidentSpec.required_extent`'s docstring (`incident.py:93-101`) explicitly documents the getattr protocol and must be rewritten with it.]
- **lib-M18**: make `wp_out` optional in `td_ve_cross_section` (`_WpOut | None = None`), raising only when method=="tw" — exactly the td_da_cross_section pattern; update callers that pass a dummy.
  [VERIFIED: `td_ve_cross_section` takes `wp_out: _WpOut` required; `td_da_cross_section` already has `wp_out: _WpOut | None = None` + `raise ValueError("td_da_cross_section: method='tw' requires `wp_out`")`. Dummy-passers: four tests in `libs/qscat/tests/test_td_extractors.py` pass `wp_out=WP_OUT` with `method="delta"`/`"flow"` (`test_delta_method_requires_position`, `test_flow_method_requires_surface`, `test_delta_method_matches_direct_dirac_construction`, `test_flow_method_matches_direct_flux_construction`). `td_ve_cross_sections_all` correctly keeps `wp_out` REQUIRED (it always builds the TW extractor). The kernel-consolidation plan's Task 6 rewrote `td_ve_cross_section`'s body into a `build_extractor` closure — the edit in Task 2 below is written against that shape.]
- **lib-M3**: settle the argument-order convention. Recommendation to evaluate: cross-section functions keep grids-first (the majority + the qscat-run runner's call pattern), constructors/analysis keep model-first; `lcp_da_cross_section(nuclear_grid, mu, Vd, Gamma, ...)` is the real outlier (takes bare mu instead of model) — decide whether to add a model-accepting signature now or document the exception; either way the RULE gets written into docs (a short ADR — read docs/adr/ for the next number, likely 0006 or 0007 depending on the infrastructure plan's ADR; check docs/superpowers/plans/2026-08-25-phase2-infrastructure.md which reserves ADR 0006 for the Rust decision, so use the next free number) and into the qscat-conventions skill.
  [VERIFIED: the split is real and consistent — every σ/observable solver is grids-first (`ve/da/dr_cross_section`, all four `td_*`, `nrm_ve/nrm_da`, `lcp_da`), every model-derived builder is model-first (`resonance_levels`, `local_complex_potential`, `exact_resonance_states`, `propose_grid`, `tw_analysis`, `bo.resonance_curve`). `docs/adr/` holds 0001-0005; `2026-08-25-phase2-infrastructure.md` reserves 0006 (Rust stub decision); no sibling plan claims 0007 (grepped all five 2026-08-25 plans). Decision, taken in Task 6: DOCUMENT the `lcp_da_cross_section` exception, do not add a model-accepting signature — rationale in the task.]
- **docs-N15**: apply ADR 0004's provisional-signature marker to every wide-signature solver it names (read the ADR; grep "provisional" — currently 2 sites) — one docstring line each.
  [VERIFIED: ADR 0004 §3 defines the marker for "the wide solver signatures targeted for the context-object refactor" — the `(grids, model, eps, chi, v_init, …)` group. `grep -rn provisional` over code finds exactly ONE site (`problem.py:23`, which *claims* the functional solvers "are marked provisional" — currently false) plus the docs sites (`adr/0004`, `api/index.md:12`). Ten solver docstrings get the one-line marker in Task 9, which also makes `problem.py:23`'s claim true.]

**Ordering rationale** (per the review): Literals + `wp_out` + `IncidentSpec` first (small, independent), renames next, the ADR while its decisions are fresh, `ScatteringProblem` last (biggest, and its docstrings/tests reference the settled names), provisional markers after the facade exists so each marker can truthfully point at its stable facade route.

---

## Task 1 — lib-M14: `Method`, `Axis`, `Channel`, `Verdict` Literals

**Files:**
- `libs/qscat/qscat/core/time_dependent.py`
- `libs/qscat/qscat/core/td_extractors.py`
- `libs/qscat/qscat/core/assignment.py`
- `libs/qscat/qscat/core/__init__.py`
- `libs/qscat/qscat/tuning/propose.py`
- `apps/qscat-run/qscat_run/runner.py`
- `docs/api/core.md`
- Create `libs/qscat/tests/test_api_literals.py`

**Interfaces:**
- Produces: `qscat.core.time_dependent.Method = Literal["tw", "delta", "flow"]` (module-level, NOT in `__all__`); `qscat.core.td_extractors.Axis = Literal["electronic", "nuclear"]` (NOT in `__all__`; `_AXES = get_args(Axis)`); `qscat.tuning.propose.Channel = Literal["ve", "dissociation"]` (NOT in `__all__`, matching `Coordinate`'s precedent); `qscat.core.assignment.Verdict = Literal["ok", "spurious", "basis-limited", "box-limited", "weak", "mixed", "distant"]` (PUBLIC: added to `assignment.__all__` and re-exported from `qscat.core`).
- Consumes (annotation changes only, no behavior): `td_ve_cross_section(..., method: Method = "tw", ...)`, `td_da_cross_section(..., method: Method = "flow", ...)`; `TannorWeeks/Dirac/Flux.__init__(..., axis: Axis = "electronic", ...)`; `propose_grid(..., channel: Channel = "ve", ...)`; `OverlapPair.verdict: Verdict`; runner's `_build_extractor(..., axis: Literal["electronic", "nuclear"], ...)`.

**Steps:**

- [ ] Re-verify the post-consolidation shapes this task edits (the kernel plan moved/rewrote code in both core files):
  ```bash
  grep -n "method: str" libs/qscat/qscat/core/time_dependent.py         # expect 2 hits: td_ve (= "tw"), td_da (= "flow")
  grep -n "_AXES = \|axis: str" libs/qscat/qscat/core/td_extractors.py   # expect _AXES tuple + 3 extractor __init__ hits + _check_axis/_axis_grid_index
  grep -n "channel: str" libs/qscat/qscat/tuning/propose.py              # expect 1 hit in propose_grid
  grep -n "verdict: str" libs/qscat/qscat/core/assignment.py             # expect 1 hit (OverlapPair)
  ```
  If a count differs, read the surrounding code and apply the same edits to every hit of the same role.
- [ ] Write the failing test, `libs/qscat/tests/test_api_literals.py`:

```python
"""lib-M14 (2026-08-25 API surface pass): the free-form str parameters carry
Literal types, and each Literal is the single source of its legal values --
the runtime `raise` paths (which stay, for callers holding a plain str)
validate against the SAME tuple the type checker sees."""

from __future__ import annotations

from typing import get_args


def test_method_literal_names_the_three_extraction_methods() -> None:
    from qscat.core.time_dependent import Method

    assert get_args(Method) == ("tw", "delta", "flow")


def test_axis_literal_is_the_single_source_of_the_axes_tuple() -> None:
    from qscat.core.td_extractors import _AXES, Axis

    assert get_args(Axis) == ("electronic", "nuclear")
    assert _AXES == get_args(Axis)


def test_channel_literal_names_the_two_mesh_channels() -> None:
    from qscat.tuning.propose import Channel

    assert get_args(Channel) == ("ve", "dissociation")


def test_verdict_literal_is_public_and_names_the_seven_verdicts() -> None:
    # Public through both the module and the qscat.core re-export (it is the
    # vocabulary `OverlapPair.verdict` speaks; typed user code needs it).
    from qscat.core import Verdict as core_verdict
    from qscat.core.assignment import OverlapPair, Verdict

    assert core_verdict is Verdict
    assert get_args(Verdict) == (
        "ok",
        "spurious",
        "basis-limited",
        "box-limited",
        "weak",
        "mixed",
        "distant",
    )
    assert OverlapPair.__dataclass_fields__["verdict"].type == "Verdict"
```

- [ ] Run it, expect ImportError on every test (none of the Literals exist yet):
  ```bash
  uv run --no-sync pytest libs/qscat/tests/test_api_literals.py -q
  ```
- [ ] `core/time_dependent.py`: extend the typing import to `from typing import TYPE_CHECKING, Literal, Protocol`, and add immediately after the `_WpIn = dict[str, float]` / `_WpOut = dict[str, float]` pair:

```python
# The three TD energy-extraction methods (see td_extractors.py's module
# docstring for the transforms). The Literal is the mypy/IDE layer; the
# dispatch `raise ValueError` paths in td_ve_cross_section /
# td_da_cross_section stay as the runtime layer for callers holding a plain
# str. Not in `__all__`: an annotation vocabulary, not an advertised name.
Method = Literal["tw", "delta", "flow"]
```

  Then change the two signatures: `method: str = "tw"` → `method: Method = "tw"` (`td_ve_cross_section`) and `method: str = "flow"` → `method: Method = "flow"` (`td_da_cross_section`). Nothing else in either function changes — the `raise ValueError(...unknown method...)` paths stay verbatim (mypy strict does not enable `warn_unreachable`, so the post-narrowing raise is fine).
- [ ] `core/td_extractors.py`: extend the typing import (currently `from typing import TYPE_CHECKING`) to `from typing import TYPE_CHECKING, Literal, get_args`, and replace the `_AXES = ("electronic", "nuclear")` line with:

```python
# The two exit axes every extractor implements (electronic = VE, nuclear =
# DA; module docstring). The Literal is the single source: `_AXES` is derived
# from it, so `_check_axis`'s runtime raise and the type checker can never
# disagree. Not in `__all__` (same reasoning as `time_dependent.Method`).
Axis = Literal["electronic", "nuclear"]
_AXES = get_args(Axis)
```

  Then change `axis: str = "electronic"` → `axis: Axis = "electronic"` in all three extractor `__init__` signatures (`TannorWeeks`, `Dirac`, `Flux`). `_check_axis(axis: str, cls_name: str)` and `_axis_grid_index(axis: str)` keep `str` deliberately — they are the runtime validation layer and are called before/after narrowing.
- [ ] `apps/qscat-run/qscat_run/runner.py`: `_build_extractor`'s `axis: str` parameter becomes `axis: Literal["electronic", "nuclear"]` (`Literal` is already imported at `runner.py:71`; both call sites pass the literal strings `"electronic"` / `"nuclear"`, so nothing else changes).
- [ ] `libs/qscat/qscat/tuning/propose.py`: add directly under the existing `Coordinate = Literal["nuclear", "electronic"]`:

```python
# The two physical channels the mesh can target (propose_grid's docstring).
# Same style/status as `Coordinate` above: an annotation vocabulary, not in
# `__all__`. The runtime check in propose_grid validates against get_args of
# this Literal so the two layers cannot drift.
Channel = Literal["ve", "dissociation"]
```

  Change `channel: str = "ve"` → `channel: Channel = "ve"`, add `get_args` to the file's typing import, and rewrite the runtime check to read the Literal while keeping the error message byte-identical (a test matches on it):

```python
    if channel not in get_args(Channel):
        raise ValueError(f"channel must be 've' or 'dissociation', got {channel!r}")
```

- [ ] `core/assignment.py`: add `from typing import Literal` to the imports, and insert immediately before the `NO_PARTNER` constant block:

```python
# The seven verdicts `pair_by_overlap` can return, in priority order -- THE
# single source of the legal values. The module docstring's table ("What the
# verdicts mean") explains each one; `OverlapPair.verdict` is typed with this
# and the assignment chain in `pair_by_overlap` is checked against it.
Verdict = Literal[
    "ok",
    "spurious",
    "basis-limited",
    "box-limited",
    "weak",
    "mixed",
    "distant",
]
```

  Then: (1) `OverlapPair`'s field `verdict: str` → `verdict: Verdict`; (2) in `pair_by_overlap`, immediately before the verdict if/elif chain, add the two annotations `verdict: Verdict` and `level: int | None` (the chain's tuple-unpack assignments then type-check against the Literal — a typo'd verdict string becomes a mypy error, which is the point); (3) in the module docstring's "What the verdicts mean" section, after the sentence introducing the seven values, add: ``The seven strings are the `Verdict` Literal (this module's single source of the legal values); the table below explains them.``; (4) add `"Verdict"` to `assignment.__all__` (sorted position: after `"PeakAlignment"`, before `"overlap"` — the list is case-sensitive sorted).
- [ ] `core/__init__.py`: add `Verdict` to the `from .assignment import (...)` block and to `__all__` (sorted: after `"VibrationalBasis"`).
- [ ] `docs/api/core.md`: in the "Assignment and verification" section's `eval-rst` block, add one line ABOVE `.. autofunction:: qscat.core.pair_by_overlap`:

```
.. autodata:: qscat.core.Verdict
```

  (Required in the same task: `tests/test_api_docs_coverage.py` gates `qscat.core.__all__` against docs/api pages, and `Verdict` just became public.)
- [ ] Run the new test, expect pass; then the standard gate. The pre-existing runtime-validation tests are the direct behavior gates and must stay green untouched: `test_td_extractors.py::test_unknown_method_raises`, `test_tuning_propose.py::test_propose_grid_rejects_unknown_channel`, and `test_assignment.py`'s verdict tests.
- [ ] Commit:
  ```bash
  git add libs/qscat/qscat/core/time_dependent.py libs/qscat/qscat/core/td_extractors.py \
      libs/qscat/qscat/core/assignment.py libs/qscat/qscat/core/__init__.py \
      libs/qscat/qscat/tuning/propose.py apps/qscat-run/qscat_run/runner.py \
      docs/api/core.md libs/qscat/tests/test_api_literals.py
  git commit -m "feat(core,tuning): Literal types for method/axis/channel/verdict parameters"
  ```

---

## Task 2 — lib-M18: `wp_out` optional in `td_ve_cross_section`

**File:** `libs/qscat/qscat/core/time_dependent.py`, `libs/qscat/tests/test_td_extractors.py`

**Interfaces:**
- Produces: `td_ve_cross_section(..., *, dt: float, n_steps: int, wp_in: _WpIn, wp_out: _WpOut | None = None, order: int = 3, subtract_free_reference: bool = True, method: Method = "tw", position: int | None = None, surface: int | None = None) -> npt.NDArray[np.float64]` — raising `ValueError("td_ve_cross_section: method='tw' requires `wp_out`")` only when `method == "tw"` and `wp_out is None`, exactly `td_da_cross_section`'s pattern.
- Consumes: the post-kernel-consolidation `build_extractor` closure inside `td_ve_cross_section` (kernel plan Task 6). `td_ve_cross_sections_all` keeps `wp_out: _WpOut` REQUIRED (it always builds all three extractors including TW).

**Steps:**

- [ ] Re-verify the post-consolidation body: `grep -n "def build_extractor" libs/qscat/qscat/core/time_dependent.py` must hit inside `td_ve_cross_section` (kernel plan Task 6 landed). If it does not (the consolidation was skipped or reshaped), apply the same `wp_out is None` guard to whatever construction site builds the `TannorWeeks` extractor instead — the guard must run before the constructor, and only on the tw path.
- [ ] Add the failing test to `libs/qscat/tests/test_td_extractors.py` (next to `test_delta_method_requires_position`):

```python
def test_tw_method_requires_wp_out() -> None:
    with pytest.raises(ValueError, match="requires `wp_out`"):
        td_ve_cross_section(
            TG,
            N2,
            EPS,
            CHI,
            V_INIT,
            VPRIMES,
            0.10,
            dt=DT,
            n_steps=N_STEPS,
            wp_in=WP_IN,
            method="tw",
        )
```

- [ ] Run it, expect FAIL — today the call fails at the *signature* level (`TypeError: missing ... 'wp_out'`), not with the ValueError the test demands:
  ```bash
  uv run --no-sync pytest libs/qscat/tests/test_td_extractors.py::test_tw_method_requires_wp_out -q
  ```
- [ ] Implement in `td_ve_cross_section`:
  1. Signature: `wp_out: _WpOut,` → `wp_out: _WpOut | None = None,` (it is already keyword-only; position in the keyword block unchanged).
  2. In `build_extractor`'s tw branch, guard before construction:

```python
        if method == "tw":
            if wp_out is None:
                raise ValueError("td_ve_cross_section: method='tw' requires `wp_out`")
            return TannorWeeks(tgrid, model, eps, chi, v_init, vprimes, wp_out, wp_in=wp_in, dt=dt)
```

  (mypy narrows `wp_out` to `_WpOut` after the raise; the `"delta"`/`"flow"` branches never touch it. Note the guard sits inside `build_extractor`, so the free-reference second construction re-checks it identically — it cannot have become None in between, but the symmetry keeps the closure self-contained.)
  3. Docstring: replace the sentence fragment ``are the outgoing test function's parameters -- used only by `method="tw"` (`outgoing_channel`/`eta_outgoing`); ignored by `"delta"`/`"flow"`.`` with ``are the outgoing test function's parameters -- required by `method="tw"` (omitting them raises `ValueError`, mirroring `td_da_cross_section`'s contract) and unused by `"delta"`/`"flow"`, which need no propagated test packet.``
- [ ] Update the four dummy-passers in `libs/qscat/tests/test_td_extractors.py`: delete the `wp_out=WP_OUT,` line from `test_delta_method_requires_position`, `test_flow_method_requires_surface`, `test_delta_method_matches_direct_dirac_construction`, and `test_flow_method_matches_direct_flux_construction` — these calls now exercise the new optionality (delta/flow with NO `wp_out` at all), which is exactly the behavior lib-M18 buys. Leave `test_unknown_method_raises` as is (its `wp_out` is not a dummy for any particular method and the test targets the dispatch error).
- [ ] Run the touched tests, expect pass; then the standard gate:
  ```bash
  uv run --no-sync pytest libs/qscat/tests/test_td_extractors.py -q -m "not slow"
  ```
- [ ] Commit:
  ```bash
  git add libs/qscat/qscat/core/time_dependent.py libs/qscat/tests/test_td_extractors.py
  git commit -m "feat(core): td_ve_cross_section only requires wp_out when method='tw'"
  ```

---

## Task 3 — lib-M15: `incident: IncidentSpec | None`, direct attribute access

**Files:** `libs/qscat/qscat/tuning/propose.py`, `libs/qscat/qscat/tuning/incident.py`, `libs/qscat/tests/test_tuning_incident.py`

**Interfaces:**
- Produces: `propose_grid(..., incident: IncidentSpec | None = None, ...)` with the body reading `incident.required_extent()` and `incident.incident_energy()` directly.
- Consumes: `qscat.tuning.incident.IncidentSpec` (frozen dataclass; both methods exist — `incident.py:94` and `:107`). All existing call sites pass `IncidentSpec` or omit the argument (verified repo-wide; only `test_tuning_incident.py` passes it).

**Steps:**

- [ ] Add the failing test to `libs/qscat/tests/test_tuning_incident.py` (add `import inspect` to its imports):

```python
def test_propose_grid_incident_parameter_is_typed_incident_spec() -> None:
    """lib-M15: the parameter is `IncidentSpec | None`, not a duck-typed
    `object | None` read through silent getattr defaults."""
    ann = inspect.signature(propose_grid).parameters["incident"].annotation
    assert ann == "IncidentSpec | None"
```

- [ ] Run it, expect FAIL (annotation is currently `"object | None"`).
- [ ] Implement in `tuning/propose.py`:
  1. In the existing `if TYPE_CHECKING:` block (which already imports `ResonanceModel` — the pattern the finding names), add `from .incident import IncidentSpec`.
  2. Signature: `incident: object | None = None,` → `incident: IncidentSpec | None = None,`.
  3. Replace the getattr block in the body:

```python
    if incident is not None:
        x_max = max(x_max, float(incident.required_extent()))
        e_max_mesh = max(e_max_mesh, float(incident.incident_energy()))
```

  4. Docstring: the `incident` section currently documents the getattr protocol at length. Replace its two bullet leads — ``EXTENT: `getattr(incident, "required_extent", lambda: 0.0)()` extends`` → ``EXTENT: `incident.required_extent()` extends`` and ``RESOLUTION: `getattr(incident, "incident_energy", lambda: 0.0)()` raises`` → ``RESOLUTION: `incident.incident_energy()` raises`` — and DELETE the paragraph beginning ``Both getattrs default to `0.0` (a no-op) for any duck-typed `incident` that does not define the corresponding method; `IncidentSpec` defines both, precisely so these duck-typed calls work against it unchanged -- see `qscat.tuning.incident`'s docstring for the reconciliation.`` replacing it with one sentence: ``The parameter is a real `IncidentSpec` (no duck typing -- a wrong object now fails loudly at the call instead of silently contributing `0.0`); the placement logic itself (impulse/width/observation boundary, `tw_analysis`) lives in `qscat.tuning.incident` -- this is only the extent/resolution floor.``
- [ ] Implement in `tuning/incident.py`: `IncidentSpec.required_extent`'s docstring paragraph ``A METHOD (not a free function) so `propose_grid`'s existing `getattr(incident, "required_extent", lambda: 0.0)()` call works against an `IncidentSpec` unchanged.`` becomes ``A METHOD on the spec (`propose_grid` calls it directly on its typed `incident: IncidentSpec | None` parameter); the module-level `required_extent(spec)` below remains for callers who prefer a free function.``
- [ ] Grep for stragglers — must come back empty:
  ```bash
  grep -rn "getattr(incident" libs apps validation projects
  ```
- [ ] Run `libs/qscat/tests/test_tuning_incident.py` (all of it — the behavioral tests `test_propose_grid_nuclear_extends_past_a_far_incident_position`, `test_propose_grid_incident_energy_above_range_refines_resolution`, etc. are the proof the direct calls behave identically), expect pass; then the standard gate.
- [ ] Commit:
  ```bash
  git add libs/qscat/qscat/tuning/propose.py libs/qscat/qscat/tuning/incident.py \
      libs/qscat/tests/test_tuning_incident.py
  git commit -m "refactor(tuning): propose_grid takes a typed IncidentSpec, no duck-typed getattr"
  ```

---

## Task 4 — lib-C3 (pair 1): `qscat.core.nrm.scattering.free_hamiltonian` → `electronic_free_hamiltonian`

**Files:** `libs/qscat/qscat/core/nrm/scattering.py`, `libs/qscat/qscat/core/time_dependent.py`, `libs/qscat/tests/test_nrm_scattering.py`

**Name decision (and why this side of the pair renames):** `qscat.core.time_dependent.free_hamiltonian` keeps its name — it is in `time_dependent.__all__`, imported by `apps/qscat-run/qscat_run/runner.py:101`, and aliased as `_free_hamiltonian` for `test_td_extractors.py`; renaming it churns three consumers plus docs. The nrm one is consumed only by its own module and one test file, and its collision note already has to explain in prose what the new name says directly: it is the 1-D **electronic** free Hamiltonian (`T_r + centrifugal`, fixed-R, no molecular potential), in a package whose whole vocabulary revolves around the electronic/nuclear split. So: `electronic_free_hamiltonian`. (`bare_electronic_hamiltonian` was considered and rejected — "free" is the physics term both docstrings and PRA 77 use for the V=0 reference.)

**Interfaces:**
- Produces: `qscat.core.nrm.scattering.electronic_free_hamiltonian(grid: FemDvrEcsGrid, ell: int) -> npt.NDArray[np.complex128]` (in `scattering.__all__`); deprecated module-`__getattr__` alias `free_hamiltonian` (NOT in `__all__` — ADR 0004 pre-1.0 allows the surface change; the alias keeps imports working one cycle).
- Consumes: `scattering.py`'s own `scattering_state` (the one internal caller) and `test_nrm_scattering.py` (import + 9 call sites).

**Steps:**

- [ ] Add the failing test to `libs/qscat/tests/test_nrm_scattering.py`:

```python
def test_free_hamiltonian_is_a_deprecated_alias_for_electronic_free_hamiltonian():
    """lib-C3 rename (2026-08-25): the old name warns and resolves to the new
    function for one release cycle (ADR 0004)."""
    from qscat.core.nrm import scattering

    with pytest.warns(DeprecationWarning, match="electronic_free_hamiltonian"):
        old = scattering.free_hamiltonian
    assert old is scattering.electronic_free_hamiltonian
```

- [ ] Run it, expect FAIL (`electronic_free_hamiltonian` does not exist; no warning fires).
- [ ] Implement in `core/nrm/scattering.py`:
  1. Rename the function: `def free_hamiltonian(` → `def electronic_free_hamiltonian(`; in `__all__`, `"free_hamiltonian"` → `"electronic_free_hamiltonian"` (keep the list sorted).
  2. Update the one internal call site (`h_free = free_hamiltonian(grid, ell)` inside `scattering_state`) to the new name; grep the file for any second use (`scattering_state_minus` documents but may not call it).
  3. Rewrite the function's NAME COLLISION paragraph — the collision is now resolved, so the note becomes a disambiguation: ``Renamed from `free_hamiltonian` (2026-08-25 API surface pass) to end the collision with `qscat.core.time_dependent.free_hamiltonian`, which is a different function: that one is the FULL 2-D `model.hamiltonian` with only the electron-molecule interaction removed (the elastic free-reference propagation); this one is the bare 1-D electronic kinetic-plus-centrifugal operator, carrying no molecular potential.``
  4. Add `import warnings` to the module imports and the deprecation shim at the END of the module — **this is the plan's one full showing of the shim; Task 5 reuses it with only the names/module changed**:

```python
# --- Deprecated aliases (2026-08-25 API surface pass) ------------------------
# One release cycle per ADR 0004, then delete this block. Not in `__all__`:
# the public surface is the new name; the alias only keeps old imports alive.

_DEPRECATED = {"free_hamiltonian": "electronic_free_hamiltonian"}


def __getattr__(name: str) -> object:
    if name in _DEPRECATED:
        new = _DEPRECATED[name]
        warnings.warn(
            f"{__name__}.{name} was renamed to {new} in the 2026-08-25 API "
            "surface pass; the old name is a deprecated alias for one release "
            "cycle (docs/adr/0004-public-api-stability-policy.md)",
            DeprecationWarning,
            stacklevel=2,
        )
        return globals()[new]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
```

  5. In `core/time_dependent.py`, update `free_hamiltonian`'s NAME COLLISION paragraph to the resolved form: ``Not to be confused with `qscat.core.nrm.scattering.electronic_free_hamiltonian` (renamed from `free_hamiltonian` in the 2026-08-25 API surface pass precisely to end this collision): that one is the bare 1-D electronic kinetic-plus-centrifugal operator with NO molecular potential at all; this one is the FULL 2-D `model.hamiltonian` with only the electron-molecule interaction removed.``
  6. Update `test_nrm_scattering.py`: the import at the top (`free_hamiltonian` → `electronic_free_hamiltonian`) and every call site (grep shows 9 `h_free = free_hamiltonian(grid, ell=1)`-shaped lines plus one docstring mention at `:50` — update the docstring sentence too).
- [ ] Verify no stragglers (the only remaining hits must be the shim itself and the two rewritten docstring notes):
  ```bash
  grep -rn "free_hamiltonian" libs/qscat/qscat/core/nrm apps validation projects benchmarks docs/physics docs/api
  ```
  (Expected survivors outside nrm: `time_dependent.free_hamiltonian` + its `_free_hamiltonian` alias and the runner import — those keep the name by decision.)
- [ ] Run the test, expect pass; then the standard gate (no docs/api change needed: the nrm scattering submodule is not documented in docs/api and not gated by the coverage test — verified).
- [ ] Commit:
  ```bash
  git add libs/qscat/qscat/core/nrm/scattering.py libs/qscat/qscat/core/time_dependent.py \
      libs/qscat/tests/test_nrm_scattering.py
  git commit -m "refactor(nrm): rename free_hamiltonian to electronic_free_hamiltonian"
  ```

---

## Task 5 — lib-C3 (pair 2): `qscat.tuning.resonance.resonance_curve` → `resonance_curve_arrays`

**Files:** `libs/qscat/qscat/tuning/resonance.py`, `libs/qscat/qscat/tuning/__init__.py`, `libs/qscat/qscat/tuning/propose.py`, `libs/qscat/qscat/core/bo.py`, `libs/qscat/tests/test_tuning_resonance.py`, `libs/qscat/tests/test_tuning_propose.py`, `docs/api/tuning.md`, `docs/physics/discretisation-tuning.md`

**Name decision (and why this side renames):** `qscat.core.bo.resonance_curve` keeps its name — it is re-exported in `qscat.core.__all__`, documented at `docs/api/core.md:132`, consumed by `validation/n2/pole_verification.py`, and it is the physics-facing object (an `ElectronicCurves` carrying phase-aligned eigenvectors for building BO basis states). The tuning one exists only to size a grid, and its own NAME COLLISION note (`tuning/resonance.py:130`) states the distinguishing trait in one clause: "This one returns plain `(R, V_d, Gamma)` arrays". The rename encodes exactly that clause: `resonance_curve_arrays`. (`sampled_resonance_curve` was considered — the sparse sampling is an implementation economy, not the contract; the arrays-not-states return IS the contract.)

**Interfaces:**
- Produces: `qscat.tuning.resonance.resonance_curve_arrays(model, elec_grid_a, elec_grid_b, *, R_max: float = 22.0, n_dense: int = 25, region: tuple[float, float] | None = None) -> tuple[FloatArray, FloatArray, FloatArray]` (signature unchanged, name only); exported as `qscat.tuning.resonance_curve_arrays` (in both `__all__`s); deprecated `__getattr__` alias `resonance_curve` at BOTH module levels (`qscat.tuning.resonance` and `qscat.tuning` — both import paths were public).
- Consumes: `tuning/propose.py:312`'s call `R, Vd, Gamma = resonance_curve(model, ga, gb, R_max=x_max, n_dense=resonance_n_dense)`.

**Steps:**

- [ ] Add the failing tests to `libs/qscat/tests/test_tuning_resonance.py` (it already imports `pytest`; add if not):

```python
def test_resonance_curve_is_a_deprecated_alias_at_both_import_paths():
    """lib-C3 rename (2026-08-25): both public paths to the old name warn and
    resolve to `resonance_curve_arrays` for one release cycle (ADR 0004)."""
    import qscat.tuning
    from qscat.tuning import resonance

    with pytest.warns(DeprecationWarning, match="resonance_curve_arrays"):
        old_mod = resonance.resonance_curve
    with pytest.warns(DeprecationWarning, match="resonance_curve_arrays"):
        old_pkg = qscat.tuning.resonance_curve
    assert old_mod is resonance.resonance_curve_arrays
    assert old_pkg is resonance.resonance_curve_arrays
```

- [ ] Run it, expect FAIL.
- [ ] Implement:
  1. `tuning/resonance.py`: rename `def resonance_curve(` → `def resonance_curve_arrays(`; `__all__ = ["interaction_region", "resonance_curve"]` → `["interaction_region", "resonance_curve_arrays"]`; add `import warnings`; append the Task-4 shim with `_DEPRECATED = {"resonance_curve": "resonance_curve_arrays"}`. Rewrite the docstring's NAME COLLISION paragraph to the resolved form: ``Renamed from `resonance_curve` (2026-08-25 API surface pass) to end the collision with `qscat.core.bo.resonance_curve`, which shares the same underlying pole walk but returns an `ElectronicCurves` carrying the eigenVECTORS, for building Born-Oppenheimer basis states. This one returns plain `(R, V_d, Gamma)` arrays -- the name says so -- and exists only to size a grid, so it discards the states and samples as sparsely as it can.``
  2. `tuning/__init__.py`: `from .resonance import interaction_region, resonance_curve` → `from .resonance import interaction_region, resonance_curve_arrays`; in `__all__`, `"resonance_curve"` → `"resonance_curve_arrays"`; update the docstring bullet at `:51` to the new name; add `import warnings` and the same shim with `_DEPRECATED = {"resonance_curve": "resonance_curve_arrays"}` (the package `__init__` needs its own copy because `from qscat.tuning import resonance_curve` resolved against the package namespace).
  3. `tuning/propose.py`: update the import (`from .resonance import interaction_region, resonance_curve_arrays`), the call at the `_resonant_nuclear_mesh` site (`R, Vd, Gamma = resonance_curve_arrays(model, ga, gb, R_max=x_max, n_dense=resonance_n_dense)`), and the four prose mentions (grep `resonance_curve` in the file: the comment at ~`:128`, the docstrings at ~`:242` and ~`:276`, and `propose_grid`'s channel docstring at ~`:423` — each becomes `resonance_curve_arrays`, keeping the `qscat.tuning.resonance.` qualification where present).
  4. `core/bo.py:233`: rewrite its NAME COLLISION paragraph to the resolved form: ``The grid-sizing sibling `qscat.tuning.resonance.resonance_curve_arrays` (renamed from `resonance_curve` in the 2026-08-25 API surface pass) runs the same underlying pole walk but returns plain `(R, V_d, Gamma)` arrays with no states -- use that one to size a grid, this one to build BO basis states.``
  5. Tests: `test_tuning_propose.py:18` import → `from qscat.tuning.resonance import resonance_curve_arrays` (and its use sites in that file, grep); `test_tuning_resonance.py:41,45` → `from qscat.tuning import interaction_region, resonance_curve_arrays` / `R, Vd, G = resonance_curve_arrays(F2, ga, gb, R_max=22.0, n_dense=20)`.
  6. `docs/api/tuning.md:32`: autosummary entry `resonance_curve` → `resonance_curve_arrays` (this is the edit the api-docs coverage gate enforces in the same task).
  7. `docs/physics/discretisation-tuning.md:170`: `qscat.tuning.resonance.resonance_curve` → `qscat.tuning.resonance.resonance_curve_arrays`.
- [ ] Verify no stragglers — after this, every remaining `resonance_curve` hit must be either `qscat.core.bo.resonance_curve` (kept), a `resonance_curve_arrays` substring match, one of the two shims, or `projects/n2_resonance/pole.resonance_curve` (a project-local toy predating promotion, not part of the qscat surface — left alone):
  ```bash
  grep -rn "resonance_curve" libs apps validation projects docs/api docs/physics --include="*.py" --include="*.md" | grep -v "resonance_curve_arrays" | grep -v "bo\." | grep -v "n2_resonance"
  ```
  Read every remaining hit and resolve it to one of the categories above (validation/n2/pole_verification.py's hits are the KEPT `qscat.core` import — untouched).
- [ ] Run the test, expect pass; then the standard gate.
- [ ] Commit:
  ```bash
  git add libs/qscat/qscat/tuning/resonance.py libs/qscat/qscat/tuning/__init__.py \
      libs/qscat/qscat/tuning/propose.py libs/qscat/qscat/core/bo.py \
      libs/qscat/tests/test_tuning_resonance.py libs/qscat/tests/test_tuning_propose.py \
      docs/api/tuning.md docs/physics/discretisation-tuning.md
  git commit -m "refactor(tuning): rename resonance_curve to resonance_curve_arrays"
  ```

---

## Task 6 — lib-M3: the argument-order convention — ADR 0007 + qscat-conventions skill

**Files:** create `docs/adr/0007-solver-argument-order.md`; `.claude/skills/qscat-conventions/SKILL.md`; `libs/qscat/qscat/core/lcp.py` (one docstring paragraph)

**Decision (evaluated as the finding asks):** the recommendation holds and matches what the code already does with one exception. RULE: (1) **Observable solvers** — functions that compute a physical observable on a prepared discretisation (`*_cross_section`, level/state solvers taking a basis) — take the discretisation first: `(grid(s), model, eps, chi, v_init, <per-call physics: vprimes/E/phi_d>, *, options)`. This is the majority (`ve/da/dr`, all four `td_*`, `nrm_ve/nrm_da`, `lcp_da`) and the qscat-run runner's call pattern. (2) **Model-derived builders** — functions whose output is derived from the model itself (curves, level pipelines, grids, tuning) — take the model first: `resonance_levels`, `local_complex_potential`, `exact_resonance_states`, `bo.resonance_curve`, `propose_grid`, `tw_analysis`. (3) **The `lcp_da_cross_section` exception is documented, NOT papered over with a model-accepting signature**: the LCP solve is defined by the *curve* `(V_d, Γ)` plus `mu` — the model appears nowhere in its 1-D equation, and curves legitimately arrive from elsewhere (a `resonance_levels(return_curve=True)` run, a fit, a file). A model-first overload would either hide an expensive electronic pole walk inside a cross-section call or take a redundant argument it ignores. The ergonomic cost is paid down instead by the facade (Task 8): `ScatteringProblem.lcp_da_cross_section` supplies `mu`/`eps`/`chi`/`v_init` from the bundle, leaving only the curve per-call.

**Steps:**

- [ ] Re-verify the ADR number is still free: `ls docs/adr/` — expect 0001-0005 plus (if the phase2-infrastructure plan has executed) `0006-rust-kernels-*.md`. Use 0007 either way: 0006 is reserved by that plan and a numbering gap is harmless if it never lands.
- [ ] Write `docs/adr/0007-solver-argument-order.md`:

```markdown
# 7. Solver argument order

Date: 2026-08-25

## Status

Accepted

## Context

The public solvers grew across eight sub-projects and the 2026-08-25 release
review asked whether their argument orders follow a rule or an accident.
Reading them all: they follow a rule with one deliberate exception, but the
rule was never written down, so every new solver re-derived it.

## Decision

1. **Observable solvers** — functions computing a physical observable on a
   prepared discretisation (`*_cross_section`, level/state solvers that take a
   vibrational basis) — take the discretisation first:

       f(grid(s), model, eps, chi, v_init, <per-call physics>, *, options)

   `<per-call physics>` is what varies between calls on one problem
   (`vprimes`, `E`, `phi_d`); everything tunable-but-stable is keyword-only.
   Examples: `ve_cross_section`, `da_cross_section`, `dr_cross_section`, the
   four `td_*` solvers, `nrm_ve_cross_section`, `nrm_da_cross_section`.

2. **Model-derived builders** — functions whose output is derived from the
   model itself (potential curves, resonance-level pipelines, proposed grids,
   incident placement) — take the model first, then the grids they need:
   `resonance_levels`, `local_complex_potential`, `exact_resonance_states`,
   `bo.resonance_curve`, `propose_grid`, `tw_analysis`.

3. **`lcp_da_cross_section` is a documented exception to (1)**: it takes
   `(nuclear_grid, mu, Vd, Gamma, eps, chi, v_init, E, ...)` — a bare reduced
   mass and a curve, not a model — because the LCP equation contains no model:
   its physics input IS the curve, which legitimately arrives from
   `resonance_levels(return_curve=True)`, from a fit, or from a file. A
   model-accepting signature would either hide an expensive electronic pole
   walk inside a cross-section call or carry a redundant argument. The
   `ScatteringProblem.lcp_da_cross_section` method supplies
   `mu`/`eps`/`chi`/`v_init` from its bundle, so the exception costs facade
   users nothing.

4. New solvers follow (1) or (2); a new exception needs its own recorded
   reason, in its docstring and here.

## Consequences

- A reader can predict any solver's leading arguments from its kind, and the
  `ScatteringProblem` facade can bundle the shared `(grid, model, eps, chi,
  v_init)` group mechanically.
- The rule is also recorded in the `qscat-conventions` skill for lookup.
```

- [ ] Add to `.claude/skills/qscat-conventions/SKILL.md`, in the `## Naming` section (as a new bullet after the `snake_case` bullet):

```markdown
- **Solver argument order** (docs/adr/0007): observable solvers take the
  discretisation first — `f(grid(s), model, eps, chi, v_init, <per-call
  physics>, *, options)` — while model-derived builders (curves, level
  pipelines, tuning) take the model first. `lcp_da_cross_section` is the one
  documented exception (curve-first, no model — see the ADR).
```

- [ ] Add the exception note to `lcp_da_cross_section`'s docstring in `core/lcp.py` (re-verify first: if the library-structure plan's lcp split already landed, the function lives in the `lcp` package's dissociation module — edit it there; the docstring is the same). After the first paragraph, insert: ``Argument-order note (docs/adr/0007): this solver deliberately takes `(nuclear_grid, mu, Vd, Gamma, ...)` rather than a `model` -- the LCP equation contains no model; its physics input IS the curve, which may come from `resonance_levels(return_curve=True)`, a fit, or a file. `ScatteringProblem.lcp_da_cross_section` supplies `mu`/`eps`/`chi`/`v_init` from its bundle.`` (Written here, before Task 8 adds that method, but committed docs may reference it — Task 8 lands in the same PR; if executing strictly incrementally, this sentence's last clause is forward-looking for exactly two tasks.)
- [ ] No failing test for a documentation task; the verification is: `uv run --no-sync pytest tests/test_docs_portability.py -q` (the docs link checker — the new ADR must not break it) plus the standard gate.
- [ ] Commit:
  ```bash
  git add docs/adr/0007-solver-argument-order.md .claude/skills/qscat-conventions/SKILL.md \
      libs/qscat/qscat/core/lcp.py
  git commit -m "docs(adr): record the solver argument-order convention (ADR 0007)"
  ```

---

## Task 7 — lib-C2 (part 1): real typed signatures for `ScatteringProblem`'s seven existing methods

**Files:** `libs/qscat/qscat/core/problem.py`, `libs/qscat/tests/test_scattering_problem.py`

**Interfaces:**
- Consumes (verified against the code, all post-consolidation): `driven.ve_cross_section(tgrid, model, eps, chi, v_init, vprimes, E, *, ordering="COLAMD", lam_scale=1.0, return_wavefunction=False)` with `Literal[True/False]` overloads; `dissociation.da_cross_section(..., v_init, E, *, n_channels=1, ordering="COLAMD", return_wavefunction=False)` (2 overloads); `dissociation.dr_cross_section(..., v_init, E, *, n_channels=3, ordering="COLAMD", return_wavefunction=False, return_amplitude=False)` (4 overloads); `td_ve_cross_section(..., *, dt, n_steps, wp_in, wp_out=None, order=3, subtract_free_reference=True, method="tw", position=None, surface=None)` (post-Task-2); `td_ve_cross_sections_all(..., *, dt, n_steps, wp_in, wp_out, position, surface, order=3, subtract_free_reference=True)`; `td_da_cross_section(..., *, dt, n_steps, wp_in, method="flow", surface=None, position=None, wp_out=None, n_channels=1, order=3)`; `td_da_cross_sections_all(..., *, dt, n_steps, wp_in, surface, position, wp_out, n_channels=1, order=3)`.
- Produces: the same seven methods on `ScatteringProblem`, minus the bundled `(grid, model, eps, chi, v_init)` group, with identical parameter names, defaults, overload sets, and return types. No `**kwargs` anywhere.

**Steps:**

- [ ] Re-verify the functional signatures this task mirrors (the kernel-consolidation plan touched `time_dependent.py` and `dissociation.py`; it kept all public signatures, but confirm):
  ```bash
  grep -n "def ve_cross_section\|def da_cross_section\|def dr_cross_section\|def td_ve_cross_section\|def td_da_cross_section\|def td_ve_cross_sections_all\|def td_da_cross_sections_all" \
      libs/qscat/qscat/core/driven.py libs/qscat/qscat/core/dissociation.py libs/qscat/qscat/core/time_dependent.py
  ```
  Read each def's full parameter list and compare against the facade signatures below; if any drifted, the facade mirrors the CODE, not this plan.
- [ ] Add the failing signature-coverage test to `libs/qscat/tests/test_scattering_problem.py` (add `import inspect` to its imports):

```python
def test_facade_methods_have_real_signatures() -> None:
    """lib-C2: every public facade method exposes the functional solver's
    real parameters -- no `**kwargs`, no `*args`, and a real return type."""
    for name, fn in inspect.getmembers(ScatteringProblem, predicate=inspect.isfunction):
        if name.startswith("_"):
            continue
        sig = inspect.signature(fn)
        kinds = {p.kind for p in sig.parameters.values()}
        assert inspect.Parameter.VAR_KEYWORD not in kinds, f"{name} takes **kwargs"
        assert inspect.Parameter.VAR_POSITIONAL not in kinds, f"{name} takes *args"
        assert sig.return_annotation is not inspect.Signature.empty, f"{name} has no return type"
```

- [ ] Run it, expect FAIL (all seven current methods take `**kwargs: Any`).
- [ ] Rewrite `core/problem.py`'s import block and type aliases (module docstring untouched in this task; Task 8 rewrites it):

```python
from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Literal, overload

import numpy as np
import numpy.typing as npt

from .dissociation import da_cross_section, dr_cross_section
from .driven import ve_cross_section
from .time_dependent import (
    Method,
    td_da_cross_section,
    td_da_cross_sections_all,
    td_ve_cross_section,
    td_ve_cross_sections_all,
)
from .vibrational import VibrationalBasis, vibrational_states

if TYPE_CHECKING:
    from qscat.dvr import TensorGrid
    from qscat.model import ResonanceModel

__all__ = ["ScatteringProblem"]

# Mirrors the private copies in driven.py / dissociation.py / lcp.py (scipy
# splu's permc_spec). The library-structure pass (lib-M12) consolidates all of
# them into a public `qscat.linalg.Ordering`; if that has already landed,
# import that name here instead of redefining.
_Ordering = Literal["NATURAL", "MMD_ATA", "MMD_AT_PLUS_A", "COLAMD"]

# Return/parameter conventions, identical to the functional solvers mirrored
# below (see driven.py / dissociation.py / time_dependent.py).
_Sigma = npt.NDArray[np.float64]
_Psi = npt.NDArray[np.complex128] | None
_PsiOut = _Psi | list[_Psi]
_Amp = npt.NDArray[np.complex128]
_WpIn = dict[str, float]
_WpOut = dict[str, float]
```

- [ ] Replace the seven method bodies. The class header, `__post_init__`, and the `eps`/`chi` properties stay exactly as they are. The methods (full replacement — note the overloaded functional solvers are re-dispatched through explicit `Literal` branches, because a call with a plain `bool` flag matches none of their overloads under mypy):

```python
    # --- time-independent observables ---------------------------------------

    @overload
    def ve_cross_section(
        self,
        vprimes: list[int],
        E: float | npt.ArrayLike,
        *,
        ordering: _Ordering = ...,
        lam_scale: float = ...,
        return_wavefunction: Literal[False] = ...,
    ) -> _Sigma: ...

    @overload
    def ve_cross_section(
        self,
        vprimes: list[int],
        E: float | npt.ArrayLike,
        *,
        ordering: _Ordering = ...,
        lam_scale: float = ...,
        return_wavefunction: Literal[True],
    ) -> tuple[_Sigma, _PsiOut]: ...

    def ve_cross_section(
        self,
        vprimes: list[int],
        E: float | npt.ArrayLike,
        *,
        ordering: _Ordering = "COLAMD",
        lam_scale: float = 1.0,
        return_wavefunction: bool = False,
    ) -> _Sigma | tuple[_Sigma, _PsiOut]:
        """Vibrational-excitation cross section; same parameters, defaults and
        return convention as `qscat.core.ve_cross_section` (which this
        delegates to with the bundled grid/model/basis)."""
        if return_wavefunction:
            return ve_cross_section(
                self.grid, self.model, self.eps, self.chi, self.v_init, vprimes, E,
                ordering=ordering, lam_scale=lam_scale, return_wavefunction=True,
            )
        return ve_cross_section(
            self.grid, self.model, self.eps, self.chi, self.v_init, vprimes, E,
            ordering=ordering, lam_scale=lam_scale,
        )

    @overload
    def da_cross_section(
        self,
        E: float | npt.ArrayLike,
        *,
        n_channels: int = ...,
        ordering: _Ordering = ...,
        return_wavefunction: Literal[False] = ...,
    ) -> _Sigma: ...

    @overload
    def da_cross_section(
        self,
        E: float | npt.ArrayLike,
        *,
        n_channels: int = ...,
        ordering: _Ordering = ...,
        return_wavefunction: Literal[True],
    ) -> tuple[_Sigma, _PsiOut]: ...

    def da_cross_section(
        self,
        E: float | npt.ArrayLike,
        *,
        n_channels: int = 1,
        ordering: _Ordering = "COLAMD",
        return_wavefunction: bool = False,
    ) -> _Sigma | tuple[_Sigma, _PsiOut]:
        """Dissociative-attachment cross section; see `qscat.core.da_cross_section`."""
        if return_wavefunction:
            return da_cross_section(
                self.grid, self.model, self.eps, self.chi, self.v_init, E,
                n_channels=n_channels, ordering=ordering, return_wavefunction=True,
            )
        return da_cross_section(
            self.grid, self.model, self.eps, self.chi, self.v_init, E,
            n_channels=n_channels, ordering=ordering,
        )

    @overload
    def dr_cross_section(
        self,
        E: float | npt.ArrayLike,
        *,
        n_channels: int = ...,
        ordering: _Ordering = ...,
        return_wavefunction: Literal[False] = ...,
        return_amplitude: Literal[False] = ...,
    ) -> _Sigma: ...

    @overload
    def dr_cross_section(
        self,
        E: float | npt.ArrayLike,
        *,
        n_channels: int = ...,
        ordering: _Ordering = ...,
        return_wavefunction: Literal[True],
        return_amplitude: Literal[False] = ...,
    ) -> tuple[_Sigma, _PsiOut]: ...

    @overload
    def dr_cross_section(
        self,
        E: float | npt.ArrayLike,
        *,
        n_channels: int = ...,
        ordering: _Ordering = ...,
        return_wavefunction: Literal[False] = ...,
        return_amplitude: Literal[True],
    ) -> tuple[_Sigma, _Amp]: ...

    @overload
    def dr_cross_section(
        self,
        E: float | npt.ArrayLike,
        *,
        n_channels: int = ...,
        ordering: _Ordering = ...,
        return_wavefunction: Literal[True],
        return_amplitude: Literal[True],
    ) -> tuple[_Sigma, _PsiOut, _Amp]: ...

    def dr_cross_section(
        self,
        E: float | npt.ArrayLike,
        *,
        n_channels: int = 3,
        ordering: _Ordering = "COLAMD",
        return_wavefunction: bool = False,
        return_amplitude: bool = False,
    ) -> _Sigma | tuple[_Sigma, _PsiOut] | tuple[_Sigma, _Amp] | tuple[_Sigma, _PsiOut, _Amp]:
        """Dissociative-recombination cross section (ionic target); see
        `qscat.core.dr_cross_section`."""
        if return_wavefunction and return_amplitude:
            return dr_cross_section(
                self.grid, self.model, self.eps, self.chi, self.v_init, E,
                n_channels=n_channels, ordering=ordering,
                return_wavefunction=True, return_amplitude=True,
            )
        if return_wavefunction:
            return dr_cross_section(
                self.grid, self.model, self.eps, self.chi, self.v_init, E,
                n_channels=n_channels, ordering=ordering, return_wavefunction=True,
            )
        if return_amplitude:
            return dr_cross_section(
                self.grid, self.model, self.eps, self.chi, self.v_init, E,
                n_channels=n_channels, ordering=ordering, return_amplitude=True,
            )
        return dr_cross_section(
            self.grid, self.model, self.eps, self.chi, self.v_init, E,
            n_channels=n_channels, ordering=ordering,
        )

    # --- time-dependent observables -----------------------------------------

    def td_ve_cross_section(
        self,
        vprimes: list[int],
        E: float | npt.ArrayLike,
        *,
        dt: float,
        n_steps: int,
        wp_in: _WpIn,
        wp_out: _WpOut | None = None,
        order: int = 3,
        subtract_free_reference: bool = True,
        method: Method = "tw",
        position: int | None = None,
        surface: int | None = None,
    ) -> npt.NDArray[np.float64]:
        """Time-dependent VE cross section; see `qscat.core.td_ve_cross_section`
        (same method/`wp_out`/`position`/`surface` contract)."""
        return td_ve_cross_section(
            self.grid, self.model, self.eps, self.chi, self.v_init, vprimes, E,
            dt=dt, n_steps=n_steps, wp_in=wp_in, wp_out=wp_out, order=order,
            subtract_free_reference=subtract_free_reference, method=method,
            position=position, surface=surface,
        )

    def td_ve_cross_sections_all(
        self,
        vprimes: list[int],
        E: float | npt.ArrayLike,
        *,
        dt: float,
        n_steps: int,
        wp_in: _WpIn,
        wp_out: _WpOut,
        position: int,
        surface: int,
        order: int = 3,
        subtract_free_reference: bool = True,
    ) -> dict[str, npt.NDArray[np.float64]]:
        """All three TD-VE extractors from ONE propagation; see
        `qscat.core.td_ve_cross_sections_all`."""
        return td_ve_cross_sections_all(
            self.grid, self.model, self.eps, self.chi, self.v_init, vprimes, E,
            dt=dt, n_steps=n_steps, wp_in=wp_in, wp_out=wp_out,
            position=position, surface=surface, order=order,
            subtract_free_reference=subtract_free_reference,
        )

    def td_da_cross_section(
        self,
        E: float | npt.ArrayLike,
        *,
        dt: float,
        n_steps: int,
        wp_in: _WpIn,
        method: Method = "flow",
        surface: int | None = None,
        position: int | None = None,
        wp_out: _WpOut | None = None,
        n_channels: int = 1,
        order: int = 3,
    ) -> npt.NDArray[np.float64]:
        """Time-dependent DA cross section; see `qscat.core.td_da_cross_section`."""
        return td_da_cross_section(
            self.grid, self.model, self.eps, self.chi, self.v_init, E,
            dt=dt, n_steps=n_steps, wp_in=wp_in, method=method, surface=surface,
            position=position, wp_out=wp_out, n_channels=n_channels, order=order,
        )

    def td_da_cross_sections_all(
        self,
        E: float | npt.ArrayLike,
        *,
        dt: float,
        n_steps: int,
        wp_in: _WpIn,
        surface: int,
        position: int,
        wp_out: _WpOut,
        n_channels: int = 1,
        order: int = 3,
    ) -> dict[str, npt.NDArray[np.float64]]:
        """All three TD-DA extractors from ONE propagation; see
        `qscat.core.td_da_cross_sections_all`."""
        return td_da_cross_sections_all(
            self.grid, self.model, self.eps, self.chi, self.v_init, E,
            dt=dt, n_steps=n_steps, wp_in=wp_in, surface=surface,
            position=position, wp_out=wp_out, n_channels=n_channels, order=order,
        )
```

- [ ] Add the per-method delegation tests to `libs/qscat/tests/test_scattering_problem.py`. The existing `test_problem_ve_matches_functional_api` already covers `ve`; add the rest. Real tiny solves for everything cheap and in-process-deterministic (the file's own `_grid()` is byte-identical to the kernel oracle's deck, so the TD constants below are its proven-fast settings — `POSITION = 37` is a real-region electronic DVR index, `NUCLEAR_SURFACE = 90` a real-region nuclear index on this exact grid); argument-capture for `dr` (its real solve is slow-tier — mpmath Coulomb channel functions dominate):

```python
from qscat.core import (
    da_cross_section,
    td_da_cross_section,
    td_da_cross_sections_all,
    td_ve_cross_section,
    td_ve_cross_sections_all,
)

WP_IN = {"r0": 4.0, "p0": -0.5, "sigma": 1.2}
WP_OUT = {"r0_out": 6.0, "p0_out": 0.5, "sigma_out": 1.0}
NUCLEAR_WP_OUT = {"r0_out": 7.0, "p0_out": 5.0, "sigma_out": 1.0}
POSITION = 37
NUCLEAR_SURFACE = 90
DT = 0.2
N_STEPS = 3


def _problem_and_basis() -> tuple[ScatteringProblem, np.ndarray, np.ndarray, TensorGrid]:
    tg = _grid()
    prob = ScatteringProblem(grid=tg, model=N2, n_vib=4, v_init=0)
    eps, chi = vibrational_states(tg.grids[1], N2.mu, 4, N2.v0)
    return prob, eps, chi, tg


def test_problem_da_matches_functional_api() -> None:
    prob, eps, chi, tg = _problem_and_basis()
    E = np.array([0.10, 0.60])  # one closed, one open DA channel on N2
    expected = da_cross_section(tg, N2, eps, chi, 0, E, n_channels=1)
    assert np.array_equal(prob.da_cross_section(E, n_channels=1), expected)


def test_problem_dr_delegates_exact_arguments(monkeypatch) -> None:
    """dr's real solve is slow-tier (mpmath Coulomb channels), so delegation
    is checked by argument capture instead of a re-solve."""
    prob, _, _, tg = _problem_and_basis()
    sentinel = np.zeros((2, 2))
    seen: dict[str, object] = {}

    def fake_dr(tgrid, model, eps, chi, v_init, E, **kw):
        seen.update(tgrid=tgrid, model=model, eps=eps, chi=chi, v_init=v_init, E=E, **kw)
        return sentinel

    monkeypatch.setattr("qscat.core.problem.dr_cross_section", fake_dr)
    got = prob.dr_cross_section([0.01, 0.03], n_channels=2)
    assert got is sentinel
    assert seen["tgrid"] is tg
    assert seen["model"] is N2
    assert seen["v_init"] == 0
    assert np.array_equal(seen["eps"], prob.eps)
    assert np.array_equal(seen["chi"], prob.chi)
    assert seen["n_channels"] == 2
    assert seen["ordering"] == "COLAMD"


def test_problem_td_ve_matches_functional_api() -> None:
    prob, eps, chi, tg = _problem_and_basis()
    E = [0.10, 0.15]
    expected = td_ve_cross_section(
        tg, N2, eps, chi, 0, [0, 1], E, dt=DT, n_steps=N_STEPS, wp_in=WP_IN, wp_out=WP_OUT
    )
    got = prob.td_ve_cross_section([0, 1], E, dt=DT, n_steps=N_STEPS, wp_in=WP_IN, wp_out=WP_OUT)
    assert np.array_equal(got, expected)


def test_problem_td_ve_all_matches_functional_api() -> None:
    prob, eps, chi, tg = _problem_and_basis()
    expected = td_ve_cross_sections_all(
        tg, N2, eps, chi, 0, [0, 1], 0.10, dt=DT, n_steps=N_STEPS,
        wp_in=WP_IN, wp_out=WP_OUT, position=POSITION, surface=POSITION,
    )
    got = prob.td_ve_cross_sections_all(
        [0, 1], 0.10, dt=DT, n_steps=N_STEPS,
        wp_in=WP_IN, wp_out=WP_OUT, position=POSITION, surface=POSITION,
    )
    assert set(got) == {"tw", "delta", "flow"}
    for key in expected:
        assert np.array_equal(got[key], expected[key])


def test_problem_td_da_matches_functional_api() -> None:
    prob, eps, chi, tg = _problem_and_basis()
    expected = td_da_cross_section(
        tg, N2, eps, chi, 0, 0.60, dt=DT, n_steps=N_STEPS, wp_in=WP_IN,
        method="flow", surface=NUCLEAR_SURFACE, n_channels=1,
    )
    got = prob.td_da_cross_section(
        0.60, dt=DT, n_steps=N_STEPS, wp_in=WP_IN,
        method="flow", surface=NUCLEAR_SURFACE, n_channels=1,
    )
    assert np.array_equal(got, expected)


def test_problem_td_da_all_matches_functional_api() -> None:
    prob, eps, chi, tg = _problem_and_basis()
    kw = dict(
        dt=DT, n_steps=N_STEPS, wp_in=WP_IN, surface=NUCLEAR_SURFACE,
        position=NUCLEAR_SURFACE, wp_out=NUCLEAR_WP_OUT, n_channels=1,
    )
    expected = td_da_cross_sections_all(tg, N2, eps, chi, 0, 0.60, **kw)
    got = prob.td_da_cross_sections_all(0.60, **kw)
    assert set(got) == {"flow", "delta", "tw"}
    for key in expected:
        assert np.array_equal(got[key], expected[key])
```

  (`TensorGrid` is already imported in this test file. The in-process `array_equal` pattern is proven by the file's existing ve test and by `test_td_extractors.py`'s `atol=0` method-vs-direct tests — the sparse solves and propagations are deterministic within one process.)
- [ ] Run the whole file, expect pass; then the standard gate:
  ```bash
  uv run --no-sync pytest libs/qscat/tests/test_scattering_problem.py -q
  ```
- [ ] Commit:
  ```bash
  git add libs/qscat/qscat/core/problem.py libs/qscat/tests/test_scattering_problem.py
  git commit -m "feat(core): real typed signatures on all seven ScatteringProblem methods"
  ```

---

## Task 8 — lib-C2 (part 2): the missing observables on the facade

**Files:** `libs/qscat/qscat/core/problem.py`, `libs/qscat/tests/test_scattering_problem.py`, `docs/api/core.md`

**Facade scope decision (as the finding demands, with each exclusion justified):**

*Added* (5): `lcp_da_cross_section`, `resonance_levels`, `exact_resonance_states`, `nrm_ve_cross_section`, `nrm_da_cross_section` — each one's required `(grid(s), model, eps, chi, v_init)`-shaped inputs come from the bundle; only genuinely per-call physics remains in the signature.

*Excluded, and why:*
- `lcp_resonance_levels` — legitimately expert-only. Of its seven required parameters the bundle supplies exactly one (`mu`); the other six (two angle-paired nuclear grids, the curve laid onto each, `Gamma`) are expert inputs the facade cannot derive from its single `TensorGrid`. The model-first route that computes those inputs internally is `resonance_levels`, which IS on the facade — a user who hand-builds curves is already operating below the facade's abstraction.
- `td_nrm_ve_cross_section` / `td_nrm_da_cross_section` — not in the finding's list; the newest capability, whose knob surface (`markovian`, `Vd`/`Gamma` overrides, `rank_tol`, `unabsorbed_tol`) is still settling. They remain functional-only until the surface stabilizes.
- `local_complex_potential` and the other curve/state builders (`resonance_pole_walk`, `bo.*`) — ingredient builders, not observables; the facade is an observable surface (and `resonance_levels(return_curve=True)`, on the facade, is the documented route to the curve).
- `lcp_ve_cross_section` — does not exist in `qscat.core` at this plan's execution time (the experiment-lifecycle plan graduates it); a marked extension point is left in the code (below).

**Interfaces:**
- Consumes (verified): `lcp.lcp_da_cross_section(nuclear_grid, mu, Vd, Gamma, eps, chi, v_init, E, *, ordering="COLAMD", return_wavefunction=False)` (2 overloads); `lcp.resonance_levels(model, nuclear_grid_a, nuclear_grid_b, elec_grid_a, elec_grid_b, *, re_half_width=0.05, im_half_width=0.05, resid_tol=1e-3, window=None, n_levels=None, rel_tol=1e-4, atol=1e-8, golden_rule=True, return_curve=False)` (2 overloads on `return_curve`); `resonance.exact_resonance_states(model, grid_base, grid_electronic, grid_nuclear, *, shifts, window, k=8, rel_tol=1e-4, atol=1e-8) -> ExactResonanceStates`; `nrm.nrm_ve_cross_section(nuclear_grid, elec_grid, model, phi_d, eps, chi, v_init, vprimes, E, *, ingredients=None, n_states=None, include_background=True)`; `nrm.nrm_da_cross_section(nuclear_grid, elec_grid, model, phi_d, eps, chi, v_init, E, *, ingredients=None, n_states=None)` — note nrm takes NUCLEAR grid first, electronic second.
- Produces: five new typed methods (signatures below), a rewritten module/class docstring that states the coverage honestly, and a marked `lcp_ve_cross_section` extension point.

**Steps:**

- [ ] Re-verify the five functional signatures (`lcp.py` may have become a package under the library-structure plan — the import path `from .lcp import ...` is preserved by that plan's own constraint, so only confirm):
  ```bash
  grep -rn "def lcp_da_cross_section\|def resonance_levels\|def exact_resonance_states" libs/qscat/qscat/core/
  grep -n "def nrm_ve_cross_section\|def nrm_da_cross_section" libs/qscat/qscat/core/nrm/vibrational_excitation.py libs/qscat/qscat/core/nrm/dissociation.py
  ```
- [ ] Add the failing delegation tests to `libs/qscat/tests/test_scattering_problem.py`:

```python
def test_problem_lcp_da_matches_functional_api() -> None:
    """Synthetic curve arrays: delegation equality needs a well-posed solve,
    not converged physics (docs/adr/0005 point 7)."""
    from qscat.core import lcp_da_cross_section

    prob, eps, chi, tg = _problem_and_basis()
    g_R = tg.grids[1]
    Vd = (0.2 * (1.0 - np.exp(-(g_R.points - 2.0))) ** 2 + 0.05).astype(np.complex128)
    Gamma = 0.01 * np.exp(-np.abs(g_R.points - 2.4) ** 2).astype(np.float64)
    E = np.array([0.02, 0.05])
    expected = lcp_da_cross_section(g_R, N2.mu, Vd, Gamma, eps, chi, 0, E)
    got = prob.lcp_da_cross_section(E, Vd=Vd, Gamma=Gamma)
    assert np.array_equal(got, expected)


def test_problem_resonance_levels_delegates_exact_arguments(monkeypatch) -> None:
    """The electronic pole walk is minutes-scale; delegation is checked by
    argument capture (the walk's own gates live in test_lcp_resonance_levels
    and the validation harness)."""
    prob, _, _, tg = _problem_and_basis()
    nuc_b = nuclear_grid(r_max=14.0, quadrature=6, n_complex=3, angle_deg=25.0)
    elec_b = electronic_grid(r_max=12.0, order=5, n_complex=3, angle_deg=25.0)
    sentinel = object()
    seen: dict[str, object] = {}

    def fake_levels(model, nuclear_grid_a, nuclear_grid_b, elec_grid_a, elec_grid_b, **kw):
        seen.update(
            model=model, nuc_a=nuclear_grid_a, nuc_b=nuclear_grid_b,
            elec_a=elec_grid_a, elec_b=elec_grid_b, **kw,
        )
        return sentinel

    monkeypatch.setattr("qscat.core.problem.resonance_levels", fake_levels)
    got = prob.resonance_levels(nuc_b, elec_b, n_levels=3)
    assert got is sentinel
    assert seen["model"] is N2
    assert seen["nuc_a"] is tg.grids[1] and seen["elec_a"] is tg.grids[0]
    assert seen["nuc_b"] is nuc_b and seen["elec_b"] is elec_b
    assert seen["n_levels"] == 3 and seen["return_curve"] is False


def test_problem_exact_resonance_states_delegates_exact_arguments(monkeypatch) -> None:
    """2-D pole searches are minutes of sparse factorizations; argument
    capture checks the wiring (the solver's own gates live in
    test_exact_resonance_states.py)."""
    prob, _, _, tg = _problem_and_basis()
    g_elec = TensorGrid(
        [electronic_grid(r_max=12.0, order=5, n_complex=3, angle_deg=40.0), tg.grids[1]]
    )
    g_nuc = TensorGrid(
        [tg.grids[0], nuclear_grid(r_max=14.0, quadrature=6, n_complex=3, angle_deg=25.0)]
    )
    sentinel = object()
    seen: dict[str, object] = {}

    def fake_exact(model, grid_base, grid_electronic, grid_nuclear, **kw):
        seen.update(model=model, base=grid_base, ge=grid_electronic, gn=grid_nuclear, **kw)
        return sentinel

    monkeypatch.setattr("qscat.core.problem.exact_resonance_states", fake_exact)
    got = prob.exact_resonance_states(
        g_elec, g_nuc, shifts=[-0.66 - 0.004j], window=(-0.75, -0.55, -0.05, 0.0)
    )
    assert got is sentinel
    assert seen["model"] is N2 and seen["base"] is tg
    assert seen["ge"] is g_elec and seen["gn"] is g_nuc
    assert seen["shifts"] == [-0.66 - 0.004j]
    assert seen["k"] == 8


def test_problem_nrm_methods_delegate_exact_arguments(monkeypatch) -> None:
    """The NRM ingredient build (fixed-R electronic eigenbases) is the
    expensive part and its physics gates are validation/diatomic's; argument
    capture checks the facade wiring, including the NUCLEAR-grid-first
    argument order nrm uses."""
    import qscat.core.nrm as nrm_pkg

    prob, _, _, tg = _problem_and_basis()
    phi_d = object()  # any DiscreteState; never touched by the fake
    sentinel = object()
    seen: dict[str, object] = {}

    def fake_ve(nuclear_grid, elec_grid, model, phi_d_got, eps, chi, v_init, vprimes, E, **kw):
        seen.update(
            nuc=nuclear_grid, elec=elec_grid, model=model, phi_d=phi_d_got,
            v_init=v_init, vprimes=vprimes, E=E, **kw,
        )
        return sentinel

    monkeypatch.setattr(nrm_pkg, "nrm_ve_cross_section", fake_ve)
    got = prob.nrm_ve_cross_section(phi_d, [0, 1], 0.05, n_states=20)
    assert got is sentinel
    assert seen["nuc"] is tg.grids[1] and seen["elec"] is tg.grids[0]
    assert seen["model"] is N2 and seen["phi_d"] is phi_d
    assert seen["vprimes"] == [0, 1] and seen["n_states"] == 20
    assert seen["include_background"] is True

    def fake_da(nuclear_grid, elec_grid, model, phi_d_got, eps, chi, v_init, E, **kw):
        seen.clear()
        seen.update(nuc=nuclear_grid, elec=elec_grid, phi_d=phi_d_got, E=E, **kw)
        return sentinel

    monkeypatch.setattr(nrm_pkg, "nrm_da_cross_section", fake_da)
    got = prob.nrm_da_cross_section(phi_d, 0.05)
    assert got is sentinel
    assert seen["nuc"] is tg.grids[1] and seen["elec"] is tg.grids[0]
    assert seen["phi_d"] is phi_d and seen["ingredients"] is None
```

  (Add `from qscat.core.grids import electronic_grid, nuclear_grid` to the file's imports if not present — it already imports both.)
- [ ] Run them, expect FAIL (`AttributeError`: the methods do not exist).
- [ ] Implement in `core/problem.py`:
  1. Extend the runtime imports: `from .lcp import lcp_da_cross_section, resonance_levels` and `from .resonance import exact_resonance_states`. Extend the `TYPE_CHECKING` block: `from qscat.dvr import FemDvrEcsGrid, TensorGrid`, `from .lcp import ResonanceLevels`, `from .nrm import DiscreteState, NrmIngredients`, `from .resonance import ExactResonanceStates`. Add the alias `_Window = tuple[float, float, float, float]` next to the others. (Importing `.lcp`/`.resonance` at runtime is already `qscat.core.__init__`'s behavior — no new import edge. `.nrm` is TYPE_CHECKING-only here; the runtime imports are deferred inside the two methods, preserving the hard boundary `core/__init__.py`'s docstring documents: `import qscat.core` must not pull `nrm` in.)
  2. Append the new methods to the class:

```python
    # --- LCP / resonance observables ----------------------------------------

    @overload
    def lcp_da_cross_section(
        self,
        E: float | npt.ArrayLike,
        *,
        Vd: npt.NDArray[np.complex128],
        Gamma: npt.NDArray[np.float64],
        ordering: _Ordering = ...,
        return_wavefunction: Literal[False] = ...,
    ) -> _Sigma: ...

    @overload
    def lcp_da_cross_section(
        self,
        E: float | npt.ArrayLike,
        *,
        Vd: npt.NDArray[np.complex128],
        Gamma: npt.NDArray[np.float64],
        ordering: _Ordering = ...,
        return_wavefunction: Literal[True],
    ) -> tuple[_Sigma, _PsiOut]: ...

    def lcp_da_cross_section(
        self,
        E: float | npt.ArrayLike,
        *,
        Vd: npt.NDArray[np.complex128],
        Gamma: npt.NDArray[np.float64],
        ordering: _Ordering = "COLAMD",
        return_wavefunction: bool = False,
    ) -> _Sigma | tuple[_Sigma, _PsiOut]:
        """LCP dissociative-attachment cross section on this problem's NUCLEAR
        grid; see `qscat.core.lcp_da_cross_section`. The curve `(Vd, Gamma)`
        is per-call (compute it with `resonance_levels(..., return_curve=True)`
        -- see that docstring for why not `local_complex_potential` directly);
        `mu`/`eps`/`chi`/`v_init` come from the bundle, which is what pays down
        the functional signature's documented argument-order exception
        (docs/adr/0007). The LCP magnitude needs the FINE per-molecule nuclear
        deck -- construct the problem on it for physical numbers."""
        g_R = self.grid.grids[1]
        if return_wavefunction:
            return lcp_da_cross_section(
                g_R, self.model.mu, Vd, Gamma, self.eps, self.chi, self.v_init, E,
                ordering=ordering, return_wavefunction=True,
            )
        return lcp_da_cross_section(
            g_R, self.model.mu, Vd, Gamma, self.eps, self.chi, self.v_init, E,
            ordering=ordering,
        )

    # EXTENSION POINT: when `qscat.core.lcp.lcp_ve_cross_section` (the LCP
    # VIBRATIONAL-EXCITATION route, currently a validation-layer driver)
    # graduates into qscat.core.lcp, add the matching typed method here,
    # mirroring `lcp_da_cross_section` above: the bundle supplies
    # grid/mu/eps/chi/v_init, the curve arrives as `Vd`/`Gamma` keywords.

    @overload
    def resonance_levels(
        self,
        nuclear_grid_b: FemDvrEcsGrid,
        elec_grid_b: FemDvrEcsGrid,
        *,
        re_half_width: float = ...,
        im_half_width: float = ...,
        resid_tol: float = ...,
        window: _Window | None = ...,
        n_levels: int | None = ...,
        rel_tol: float = ...,
        atol: float = ...,
        golden_rule: bool = ...,
        return_curve: Literal[False] = ...,
    ) -> ResonanceLevels: ...

    @overload
    def resonance_levels(
        self,
        nuclear_grid_b: FemDvrEcsGrid,
        elec_grid_b: FemDvrEcsGrid,
        *,
        re_half_width: float = ...,
        im_half_width: float = ...,
        resid_tol: float = ...,
        window: _Window | None = ...,
        n_levels: int | None = ...,
        rel_tol: float = ...,
        atol: float = ...,
        golden_rule: bool = ...,
        return_curve: Literal[True],
    ) -> tuple[ResonanceLevels, npt.NDArray[np.complex128], npt.NDArray[np.float64]]: ...

    def resonance_levels(
        self,
        nuclear_grid_b: FemDvrEcsGrid,
        elec_grid_b: FemDvrEcsGrid,
        *,
        re_half_width: float = 0.05,
        im_half_width: float = 0.05,
        resid_tol: float = 1e-3,
        window: _Window | None = None,
        n_levels: int | None = None,
        rel_tol: float = 1e-4,
        atol: float = 1e-8,
        golden_rule: bool = True,
        return_curve: bool = False,
    ) -> ResonanceLevels | tuple[ResonanceLevels, npt.NDArray[np.complex128], npt.NDArray[np.float64]]:
        """Born-Oppenheimer quasi-bound levels of this problem's anion; see
        `qscat.core.resonance_levels`. This problem's own electronic/nuclear
        grids are the `_a` partners; `nuclear_grid_b`/`elec_grid_b` are the
        angle-moved partners (share every real node, differ only in the ECS
        tail angle -- `qscat.core.grids.ecs_angle_family` builds a valid
        family). `return_curve=True` also returns the `(Vd, Gamma)` curve the
        levels were computed in -- the input `lcp_da_cross_section` needs."""
        if return_curve:
            return resonance_levels(
                self.model, self.grid.grids[1], nuclear_grid_b,
                self.grid.grids[0], elec_grid_b,
                re_half_width=re_half_width, im_half_width=im_half_width,
                resid_tol=resid_tol, window=window, n_levels=n_levels,
                rel_tol=rel_tol, atol=atol, golden_rule=golden_rule,
                return_curve=True,
            )
        return resonance_levels(
            self.model, self.grid.grids[1], nuclear_grid_b,
            self.grid.grids[0], elec_grid_b,
            re_half_width=re_half_width, im_half_width=im_half_width,
            resid_tol=resid_tol, window=window, n_levels=n_levels,
            rel_tol=rel_tol, atol=atol, golden_rule=golden_rule,
        )

    def exact_resonance_states(
        self,
        grid_electronic: TensorGrid,
        grid_nuclear: TensorGrid,
        *,
        shifts: npt.ArrayLike,
        window: _Window,
        k: int = 8,
        rel_tol: float = 1e-4,
        atol: float = 1e-8,
    ) -> ExactResonanceStates:
        """Exact 2-D resonance states by two-angle ECS stability; see
        `qscat.core.exact_resonance_states`. This problem's grid is the base;
        `grid_electronic`/`grid_nuclear` are the one-angle-moved partner
        TensorGrids (`ecs_angle_family` builds all three consistently).
        Seeds (`shifts`) are passed in -- typically `resonance_levels`'s
        output -- so the exact solver never depends on the approximation it
        measures."""
        return exact_resonance_states(
            self.model, self.grid, grid_electronic, grid_nuclear,
            shifts=shifts, window=window, k=k, rel_tol=rel_tol, atol=atol,
        )

    # --- nonlocal resonance model (NRM) observables --------------------------

    def nrm_ve_cross_section(
        self,
        phi_d: DiscreteState,
        vprimes: list[int],
        E: float | npt.ArrayLike,
        *,
        ingredients: NrmIngredients | None = None,
        n_states: int | None = None,
        include_background: bool = True,
    ) -> npt.NDArray[np.float64]:
        """NRM vibrational-excitation cross section; see
        `qscat.core.nrm.nrm_ve_cross_section` (this problem's nuclear and
        electronic grids fill its leading NUCLEAR-grid-first pair)."""
        # Deferred import: `import qscat.core` must never pull `nrm` in
        # (the hard boundary documented in qscat.core.__init__).
        from .nrm import nrm_ve_cross_section

        return nrm_ve_cross_section(
            self.grid.grids[1], self.grid.grids[0], self.model, phi_d,
            self.eps, self.chi, self.v_init, vprimes, E,
            ingredients=ingredients, n_states=n_states,
            include_background=include_background,
        )

    def nrm_da_cross_section(
        self,
        phi_d: DiscreteState,
        E: float | npt.ArrayLike,
        *,
        ingredients: NrmIngredients | None = None,
        n_states: int | None = None,
    ) -> npt.NDArray[np.float64]:
        """NRM dissociative-attachment cross section; see
        `qscat.core.nrm.nrm_da_cross_section`."""
        from .nrm import nrm_da_cross_section  # deferred: see nrm_ve_cross_section

        return nrm_da_cross_section(
            self.grid.grids[1], self.grid.grids[0], self.model, phi_d,
            self.eps, self.chi, self.v_init, E,
            ingredients=ingredients, n_states=n_states,
        )
```

  3. Rewrite the module docstring's closing claim (currently ``This is the recommended entry point. The functional solvers remain public ...``) to state the coverage honestly:

```
This is the recommended entry point for every observable it carries: the TI
and TD cross sections (VE/DA/DR), the LCP DA cross section, the BO resonance
levels, the exact 2-D resonance states, and the NRM VE/DA cross sections. The
functional solvers remain public (they are the low-level layer this delegates
to, and each carries ADR 0004's *provisional* marker pending the pre-1.0
signature freeze); `ScatteringProblem` is the stable API. Deliberately NOT on
the facade: `lcp_resonance_levels` (its inputs are hand-built curves on
angle-paired grids -- `resonance_levels` here is the model-first route that
computes them internally), the `td_nrm_*` solvers (knob surface still
settling), and the curve/state builders (`local_complex_potential`,
`resonance_pole_walk`, `qscat.core.bo.*` -- ingredients, not observables).
```

  4. `docs/api/core.md`: update the intro paragraph's second sentence from "exposes every observable as a method" to "exposes the supported observables as methods (the class docstring lists what is deliberately functional-only)". No new autodoc entries needed — `ScatteringProblem` is documented via `autoclass ... :members:`, which picks up the new methods.
- [ ] Run the test file, expect pass; then the standard gate.
- [ ] Commit:
  ```bash
  git add libs/qscat/qscat/core/problem.py libs/qscat/tests/test_scattering_problem.py docs/api/core.md
  git commit -m "feat(core): ScatteringProblem gains lcp/resonance/exact/nrm observables"
  ```

---

## Task 9 — docs-N15: the provisional marker on every wide functional solver

**Files:** `libs/qscat/qscat/core/driven.py`, `libs/qscat/qscat/core/dissociation.py`, `libs/qscat/qscat/core/time_dependent.py`, `libs/qscat/qscat/core/lcp.py`, `libs/qscat/qscat/core/nrm/vibrational_excitation.py`, `libs/qscat/qscat/core/nrm/dissociation.py`

**Interfaces:** docstring-only. The marker goes on the IMPLEMENTATION def's docstring (overloads carry none), as the line immediately after the summary line, in each of these ten solvers: `ve_cross_section`, `da_cross_section`, `dr_cross_section`, `td_ve_cross_section`, `td_ve_cross_sections_all`, `td_da_cross_section`, `td_da_cross_sections_all`, `lcp_da_cross_section`, `nrm_ve_cross_section`, `nrm_da_cross_section`.

**Steps:**

- [ ] Re-verify the current marker count (this is the finding's own check):
  ```bash
  grep -rn "provisional" libs/qscat/qscat --include="*.py"
  ```
  Expect exactly one hit (`problem.py`'s module docstring, updated in Task 8). If the library-structure plan's lcp split landed first, `lcp_da_cross_section` lives under the `lcp` package — apply its edit there.
- [ ] Add the one-line marker to each of the ten docstrings, adapted per solver. Template (the `ScatteringProblem` method name matches each solver's facade counterpart, which exists for all ten after Task 8):

```
    *Provisional API* (docs/adr/0004-public-api-stability-policy.md): this wide
    functional signature is the layer the context-object refactor targets and
    may change in a minor release; `ScatteringProblem.<name>` is the stable route.
```

  Concretely, `<name>` per site: `ve_cross_section` (driven.py), `da_cross_section` / `dr_cross_section` (dissociation.py), `td_ve_cross_section` / `td_ve_cross_sections_all` / `td_da_cross_section` / `td_da_cross_sections_all` (time_dependent.py), `lcp_da_cross_section` (lcp.py), `nrm_ve_cross_section` / `nrm_da_cross_section` (nrm/).
- [ ] Verification (no failing-test step for docstring lines; the greps are the gate):
  ```bash
  grep -rn "Provisional API" libs/qscat/qscat --include="*.py" | wc -l   # expect 10
  ```
  Then the standard gate (the rendered API docs pick the lines up through autodoc; the coverage test is unaffected since no `__all__` changed).
- [ ] Commit:
  ```bash
  git add libs/qscat/qscat/core/driven.py libs/qscat/qscat/core/dissociation.py \
      libs/qscat/qscat/core/time_dependent.py libs/qscat/qscat/core/lcp.py \
      libs/qscat/qscat/core/nrm/vibrational_excitation.py libs/qscat/qscat/core/nrm/dissociation.py
  git commit -m "docs(core): apply ADR 0004's provisional marker to the ten wide solver signatures"
  ```

---

## Task 10 — Closeout: leftovers audit + full gates

**Files:** none expected to change (this task verifies; it edits only if an audit step finds something).

**Steps:**

- [ ] Leftover audit — each grep's hits must all fall in the allowed categories noted:
  ```bash
  # Old nrm name: only the shim, its test, and the two rewritten disambiguation notes may hit.
  grep -rn "free_hamiltonian" libs/qscat/qscat/core/nrm libs/qscat/tests/test_nrm_scattering.py
  # Old tuning name: only the two shims and their test may hit.
  grep -rn "\bresonance_curve\b" libs/qscat/qscat/tuning docs/api/tuning.md docs/physics/discretisation-tuning.md
  # No facade kwargs and no duck-typed incident remain.
  grep -n "kwargs" libs/qscat/qscat/core/problem.py            # expect empty
  grep -n "getattr(incident" libs/qscat/qscat/tuning/propose.py # expect empty
  ```
- [ ] Deprecation-cycle bookkeeping: confirm both shims cite ADR 0004 and that `CHANGELOG.md`, if the repo has one at execution time, records the two renames and the `Verdict` export (`ls CHANGELOG.md` — if absent, skip; ADR 0004's changelog requirement binds releases, and the release is deferred).
- [ ] The kernel-consolidation oracle still holds (cheap insurance that no task here accidentally touched numerics — every change was typing/wiring/docs):
  ```bash
  uv run --no-sync pytest libs/qscat/tests/test_kernel_consolidation_oracle.py -q -m "not slow"
  ```
- [ ] The standard gate, one last time, plus the two root-level docs gates:
  ```bash
  uv run --no-sync pytest tests/ -q
  ```
- [ ] No slow-tier run is required by this plan's own changes (no solver numerics touched; the delegation tests prove the facade routes bit-identically), but if the branch also carries the kernel-consolidation commits un-merged, follow THAT plan's Task 8 slow-tier requirements before merge.
- [ ] Commit only if an audit step changed a file:
  ```bash
  git add <exact paths changed by the audit, if any>
  git commit -m "chore(core): api-surface-pass closeout fixes from the leftovers audit"
  ```

---

## Self-review notes (kept for the executor)

- **Placeholder scan:** every code block above is complete and importable as written (module names, parameter lists, and defaults were read from the tree on 2026-08-25); the only deliberate `<...>` is Task 9's `<name>` template variable, resolved in the same bullet, and Task 10's `<exact paths>` commit placeholder, standard for audit tasks.
- **Name consistency with the kernel-consolidation plan:** that plan adds `s_vector_transform`, `sigma_from_s`, `correlation_channel_s`, `flux_channel_s`, `ChannelS` to `time_dependent.py` (not in `__all__`) and keeps every public solver name/signature. Nothing in this plan collides: `Method` is a new module-level name there; Task 2 edits the `build_extractor` closure that plan's Task 6 created (re-verified by grep before editing); Task 7/8's facade mirrors the unchanged public signatures.
- **Sequencing hooks honored:** `Ordering` export left to the library-structure plan (facade carries a pointered private copy; both plans carry adapt-if-landed-first notes); `lcp.py` split re-verify steps in Tasks 6, 8, 9; `lcp_ve_cross_section` left as a marked extension point for the experiment-lifecycle plan.
- **Coverage check against the findings:** lib-M14 → Task 1 (all four Literal sites + the runner consumer); lib-M18 → Task 2 (signature, guard, docstring, four dummy-passer tests); lib-M15 → Task 3 (annotation, direct calls, both docstrings, straggler grep); lib-C3 → Tasks 4-5 (both renames, both shim levels for the tuning name, all callers, docs/api + docs/physics, deprecation tests); lib-M3 → Task 6 (ADR 0007, skill, lcp docstring, decision recorded: document the exception, no model-accepting signature); lib-C2 → Tasks 7-8 (seven typed + five added methods, exclusions justified, honest docstring, signature-coverage + delegation tests); docs-N15 → Task 9 (ten markers, making `problem.py`'s pre-existing claim true).
- **Why argument-capture instead of real solves for dr/resonance_levels/exact_resonance_states/nrm:** dr is slow-tier (mpmath Coulomb), the pole walk and 2-D pole search are minutes-scale, and an NRM closed-channel "cheap" solve would compare zeros against zeros (vacuous). Argument capture proves the exact forwarding — which is the delegation contract — at zero cost; the solvers' own physics gates live in their dedicated test files, named in each test's docstring.
- **DeprecationWarning visibility:** pytest surfaces DeprecationWarning by default in this repo (no filterwarnings suppression in `pyproject.toml`'s pytest options was found); `pytest.warns` in Tasks 4-5 is therefore reliable. The shims use `stacklevel=2` so callers see their own import site.
