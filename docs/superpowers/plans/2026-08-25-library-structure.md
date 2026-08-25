# Library Structure Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Pay down the structural debt the 2026-08-25 release review found in `libs/qscat/qscat`: one shared `Ordering` type instead of three private copies, an `Extractor` protocol that matches its implementations, honest defaults (Padé order, backend context), removal of verified-dead surface, a settled numpydoc/ruff-D/ANN convention, a dataclass result for the DR solver, one row-per-state orientation (and one clamp/save-load story) across the two resonance-state containers, and `core/lcp.py` (999 lines, three capabilities) split into a three-module package with every import path preserved.

**Architecture:** All changes are internal restructurings of `libs/qscat/qscat` plus the thin layers that consume it (`apps/qscat-run`, `projects/n2_2d_td_cross_section`, `validation/`). No physics changes: every numeric output is either bit-identical or changed only at documented round-off level (the widths clamp, bounded by `2*atol = 2e-8 Ha`, far below the documented 1e-6 Ha width noise floor). The lcp split turns a module into a package so `from qscat.core.lcp import X` keeps working for every current name without shims.

**Tech Stack:** Python 3.12 / uv workspace, numpy/scipy, ruff (adding `D` + `ANN`), mypy --strict, pytest (fast tier as the gate).

**Spec:** the "Findings addressed" section below (self-contained; from the 2026-08-25 release review)

## Global Constraints

PyPI release DEFERRED until the peer-reviewed article publishes — repo-only distribution, no publishing tasks. This plan runs AFTER the kernel-consolidation and api-surface-pass plans (docs/superpowers/plans/2026-08-25-kernel-consolidation.md, ...-api-surface-pass.md); each task therefore begins by RE-VERIFYING its target still looks as described (line numbers will have moved; the structures will not). After every task: `uv run --no-sync pytest -m "not slow" -n auto --dist loadfile` green; `uv run --no-sync mypy libs/qscat/qscat apps/qscat-run/qscat_run` clean; `uv run --no-sync ruff check .` + `ruff format --check .` clean. Public-name moves keep a deprecated re-export for one cycle (module-level `__getattr__` emitting DeprecationWarning — show the code once, fully, in the first task that needs it). Never `git commit -a`.

## Findings addressed

- **lib-M11**: split core/lcp.py (999 lines as measured 2026-08-25 on this branch; the review counted 951, three capabilities) into three modules — (a) curve construction (`resonance_pole_walk`, `resonance_eigenstate`, `resonance_eigenstate_at_peak_width`, `_assemble_lcp`, `_walk_from_anion_seed`, `local_complex_potential`), (b) the cross-section solver (`lcp_da_cross_section`), (c) the BO nuclear eigenproblem (`ResonanceLevels`, `lcp_resonance_levels`, `resonance_levels` + helpers). Keep `qscat.core.lcp.X` working for all consumers (`core/nrm/discrete_state.py`, `tuning/resonance.py`, `validation/`, `apps/qscat-run`, docs).
- **lib-M9**: `lcp.ResonanceLevels.states` is `(n_levels, grid.n)` row-per-state while `resonance.ExactResonanceStates.states` is `(n, m)` column-per-state; widths clamped in one (`np.maximum(0.0, -2*imag)`) and not the other; save/load exists on one only. Reconcile orientation (justified choice), align the clamp behaviour explicitly, give `ResonanceLevels` save/load mirroring `ExactResonanceStates`'.
- **lib-m14**: replace `dr_cross_section`'s four-overload boolean-flag returns with a small frozen dataclass result, keeping the old signature working via the same flags for one deprecation cycle; same treatment for the two-overload ve/da/lcp_da only IF cheap, else scope to dr and note.
- **lib-m9**: the `Extractor` protocol's `sigma()` lacks the `n_steps: int | None = None` parameter all three implementations added — extend the protocol.
- **lib-M12**: export the ordering Literal from `qscat.linalg` (public name `Ordering` in `linalg/__init__.__all__`), delete the three private copies in driven.py/dissociation.py/lcp.py.
- **lib-M16**: dead surface — `tuning/propose.py` `del rtol` param; `core/correlation.py` `hankel_point_value`/`outgoing_surface_wave` discarding their `grid` param; `_CoordinateSpec.charge` written never read; `dvr/grid.py element_ranges` with zero consumers; `core/time_dependent.py _propagate`'s unused `eps` param.
- **lib-m10**: `make_pade_stepper`'s `order: int = 1` default is the setting its own module docstring warns under-converges; decide.
- **lib-m18**: `_DEFAULT_BACKEND` process-global — contextvars-based scoped override as the recommended API, `set_default_backend` keeps a docstring hazard note.
- **lib-m7**: private cross-module imports — `projects/n2_2d_td_cross_section/td_cross_section.py` imports `core.time_dependent._propagate` (and `_s_vector_one_energy`, `_sigma_one_energy`). (`td_extractors` importing `_quadrature_weights` is ALREADY resolved on this branch: it imports the public `quadrature_weights` today — re-verify and skip that half.)
- **lib-M13 + style-N2/N3**: settle numpydoc — `docs/conf.py` drops `napoleon_google_docstring = True`; ruff gains `D` with `convention = "numpy"` (per-file-ignores for tests) and the missing docstrings fixed; add `ANN` with the ANN401 decision made.

**Decided conventions (implemented by the tasks below):**

1. **States orientation: ROW-per-state everywhere.** `ExactResonanceStates.states` flips from `(n, m)` columns to `(m, n)` rows. Justification: rows are the repo-wide convention for "a set of states" — `chi[v]` (vibrational), `phi[c]` (`anion_electronic_states`), `ResonanceLevels.states`, and qscat-run's `EigenStates` are all row-per-state; only `ExactResonanceStates` follows the raw-eigensolver column layout. Flipping `ResonanceLevels` instead would churn `bo.bo_basis_from_levels` (documented `(n_levels, n_R)` contract, bo.py:568), its callers, and the qscat-run `resonance_levels` npz artifact layout — more consumers than the four `states[:, i]` readers of `ExactResonanceStates`. The one hazard — pre-flip `.npz` caches loading silently transposed — is closed by a shape guard in `ExactResonanceStates.load` (Task 12).
2. **Widths: UNCLAMPED in both.** `ResonanceLevels.widths` becomes raw `-2*Im E`, matching `ExactResonanceStates`. A negative width is a diagnostic and clamping hides it; after the existing `Im E <= atol` physicality filter the most negative representable width is `-2*atol = -2e-8 Ha`, i.e. round-off, far below the documented ~1e-6 Ha width noise floor. No `widths_clamped` property — it would differ from `widths` by at most 2e-8 and carry no information.
3. **Padé default: `order=3`.** Every direct caller in the repo passes `order` explicitly (tests, `validation/diatomic/td_nrm_figures.py`), every `qscat.core.time_dependent` entry point already defaults to `order=3`, and the module's own docstring records that order 1 under-converges the TD cross section ~100%. Crank-Nicolson stays reachable explicitly (`order=1`, or `make_sparse_cn_stepper`).
4. **Ruff `D` scope:** `libs/qscat/qscat` only (negated per-file-ignore), `convention = "numpy"`, additionally ignoring the five prose-style rules D205/D209/D400/D401/D403 (155 findings that would force rewriting deliberate long-summary house-style docstrings). Enforced set = docstring PRESENCE (25 sites to fix: 17 D102 + 2 D101 + 6 D105, listed in Tasks 8–9).
5. **Ruff `ANN`:** enabled with `ANN401` in the ignore list. Measured 2026-08-25: `--select ANN` yields exactly 31 findings, all ANN401, all at two deliberate `Any` boundaries — `core/problem.py`'s `**kwargs` pass-through facade (14) and the lazily-imported-matplotlib `viz` modules (17). mypy --strict polices the interior; a per-rule ignore is honest, fixing them would mean duplicating overload stacks or importing matplotlib types eagerly.
6. **DR result: scoped to `dr_cross_section`.** New `DrResult` dataclass + `dr_solve`; `ve_cross_section`/`da_cross_section`/`lcp_da_cross_section` keep their single-flag two-overload form — each has exactly one optional extra (`psi`), `da_cross_section` internally reuses `ve_cross_section(return_wavefunction=True)`, and converting all three would touch every solver consumer for no ambiguity gain (only DR had four overloads). Recorded here as the deliberate scope cut.
7. **Backend context: keep the existing name `default_backend`.** Re-verified: `linalg/sparse_lu.py` ALREADY has a `default_backend(name)` context manager (the finding's proposed `default_backend_ctx` predates it). The fix is to reimplement the trio over `contextvars.ContextVar` so concurrent threads can't race the global, keep the CM as the recommended API, and put the hazard note on `set_default_backend`.
8. **m7 promotion names:** `_propagate` → `propagate_wavepacket`, `_s_vector_one_energy` → `s_vector_one_energy`, `_sigma_one_energy` → `sigma_one_energy`, following the module's own `free_hamiltonian` promotion precedent; old private names live one cycle behind a module-level `__getattr__` DeprecationWarning shim.

---

## Task 1: Public `Ordering` in `qscat.linalg`; delete the three private copies (lib-M12)

**Files:**
- `libs/qscat/qscat/linalg/sparse_lu.py` (rename `_Ordering` → `Ordering`, ~line 61)
- `libs/qscat/qscat/linalg/__init__.py` (import + `__all__`)
- `libs/qscat/qscat/core/driven.py` (delete copy ~lines 59–62, import instead)
- `libs/qscat/qscat/core/dissociation.py` (delete copy ~lines 60–62, import instead)
- `libs/qscat/qscat/core/lcp.py` (delete copy ~lines 437–440, import instead)
- `libs/qscat/tests/test_sparse_lu.py` (new test)

**Interfaces:**
- Produces: `qscat.linalg.Ordering = Literal["NATURAL", "MMD_ATA", "MMD_AT_PLUS_A", "COLAMD"]` (public, in `__all__`).
- Consumes (unchanged behaviour): `SparseLU.__init__(..., ordering: Ordering = "COLAMD")`, and the `ordering:` kwargs of `ve_cross_section`, `da_cross_section`, `dr_cross_section`, `lcp_da_cross_section`.

**Steps:**

- [ ] Re-verify: `grep -n "_Ordering" libs/qscat/qscat/linalg/sparse_lu.py libs/qscat/qscat/core/driven.py libs/qscat/qscat/core/dissociation.py libs/qscat/qscat/core/lcp.py` still shows one definition per file (four total). If the api-surface-pass plan already exported an ordering type, adapt: keep ONE definition, follow its chosen name.
- [ ] Add to `libs/qscat/tests/test_sparse_lu.py`:

  ```python
  def test_ordering_literal_is_public_and_single_sourced() -> None:
      """lib-M12: one public Ordering type; the solver modules import it."""
      from typing import get_args

      from qscat import linalg
      from qscat.core import dissociation, driven, lcp

      assert "Ordering" in linalg.__all__
      assert get_args(linalg.Ordering) == ("NATURAL", "MMD_ATA", "MMD_AT_PLUS_A", "COLAMD")
      # the three former private copies are gone
      for mod in (driven, dissociation, lcp):
          assert not hasattr(mod, "_Ordering"), mod.__name__
  ```

- [ ] Run `uv run --no-sync pytest libs/qscat/tests/test_sparse_lu.py -k ordering_literal -q` — expect FAIL (`Ordering` not exported).
- [ ] In `sparse_lu.py`: rename `_Ordering` to `Ordering` (definition + the in-file uses at `__init__`/the `ordering: _Ordering = "COLAMD"` sites), give it a short trailing comment (`# scipy splu's permc_spec — the public name solver modules re-use`), add `"Ordering"` to the module `__all__`.
- [ ] In `linalg/__init__.py`: add `Ordering` to the `from .sparse_lu import (...)` block and to `__all__` (keep alphabetical order: it sorts after `""` — place as `"Ordering"` before `"ShiftInvertEigs"`), and mention it in the module-docstring Public API list.
- [ ] In `driven.py`, `dissociation.py`, `lcp.py`: delete the `_Ordering = Literal[...]` line and its "Mirrors ..." comment block; add `Ordering` to each file's existing `from qscat.linalg import ...` line; replace every `_Ordering` annotation with `Ordering`. Drop `Literal` from each file's `typing` import if now unused (lcp.py still uses `Literal` in its overloads — driven.py/dissociation.py also do; verify with ruff F401 rather than guessing).
- [ ] Run the new test — expect PASS. Run full gates (pytest fast tier, mypy, ruff check + format).
- [ ] Commit: `git add libs/qscat/qscat/linalg/sparse_lu.py libs/qscat/qscat/linalg/__init__.py libs/qscat/qscat/core/driven.py libs/qscat/qscat/core/dissociation.py libs/qscat/qscat/core/lcp.py libs/qscat/tests/test_sparse_lu.py` then commit `refactor(linalg): export Ordering, drop the three private copies`.

## Task 2: `Extractor.sigma` protocol gains `n_steps` (lib-m9)

**Files:**
- `libs/qscat/qscat/core/time_dependent.py` (protocol, ~lines 96–113)
- `libs/qscat/tests/test_td_extractors.py` (new test)

**Interfaces:**
- Produces (protocol):

  ```python
  class Extractor(Protocol):
      def record(self, psi: npt.NDArray[np.complex128]) -> None: ...
      def sigma(
          self,
          E: float | npt.ArrayLike,
          *,
          free: Extractor | None = None,
          n_steps: int | None = None,
      ) -> npt.NDArray[np.float64]: ...
  ```

- Consumes: the three implementations already carry exactly this signature (`TannorWeeks.sigma` td_extractors.py:371, `Dirac.sigma` :735, `Flux.sigma` :1123 — all `(self, E, *, free=None, n_steps=None)`), so widening the protocol breaks nothing typed against it (`propagate(extractors=...)` only calls `record`; `td_da_cross_section`'s `ext: Extractor` calls `sigma(E, free=...)` without `n_steps`, still valid).

**Steps:**

- [ ] Re-verify the three implementation signatures still read `def sigma(self, E, *, free=None, n_steps=None)` and the protocol still lacks `n_steps`.
- [ ] Add to `libs/qscat/tests/test_td_extractors.py`:

  ```python
  def test_extractor_protocol_declares_n_steps() -> None:
      """lib-m9: the protocol must carry the n_steps parameter all three
      implementations added, so a caller typed against Extractor can use it."""
      import inspect

      from qscat.core.time_dependent import Extractor

      params = inspect.signature(Extractor.sigma).parameters
      assert "n_steps" in params
      assert params["n_steps"].default is None
      assert params["n_steps"].kind is inspect.Parameter.KEYWORD_ONLY
  ```

- [ ] Run it — expect FAIL (`"n_steps" not in params`).
- [ ] Extend the protocol's `sigma` stub with `n_steps: int | None = None` (as above), and — since these two stubs are also two of the D102 sites Task 10 will enforce — give both protocol methods docstrings now:
  - `record`: `"""Accumulate this step's datum from the current `psi(t_n)` (called at every propagation step)."""`
  - `sigma`: `"""Transform the recorded series into a cross section, shape `(len(E), n_channels)`. `free` supplies the free-reference extractor for the elastic subtraction; `n_steps` truncates the transform to the first `n_steps` recorded samples (a convergence probe — both runs must share the step schedule)."""`
- [ ] Run the new test — expect PASS. Full gates. (mypy is the real check here: the three classes must still satisfy the protocol.)
- [ ] Commit: `git add libs/qscat/qscat/core/time_dependent.py libs/qscat/tests/test_td_extractors.py`; message `fix(core): declare n_steps on the Extractor protocol`.

## Task 3: `make_pade_stepper` defaults to `order=3` (lib-m10)

**Files:**
- `libs/qscat/qscat/evolution/pade.py` (signature ~line 84, docstring)
- `libs/qscat/tests/test_pade.py` (new test)

**Interfaces:**
- Produces: `def make_pade_stepper(H: sp.spmatrix, dt: float, order: int = 3) -> Callable[[npt.NDArray[np.complexfloating[Any, Any]]], npt.NDArray[np.complex128]]`

**Decision (do not punt):** default changes to 3. Rationale recorded in the docstring: the module's own docstring documents that order 1 (CN) under-converged the TD cross section (~100% accumulated error at dt=0.5–1.0), every `qscat.core.time_dependent` entry point already defaults to `order=3`, and — verified 2026-08-25 — every direct caller in the repo passes `order` explicitly (`libs/qscat/tests/test_pade.py:50,69,88,104`, `validation/diatomic/td_nrm_figures.py:1309`), so no call site changes behaviour implicitly. `order=1` remains available and remains bit-identical to `make_sparse_cn_stepper`.

**Steps:**

- [ ] Re-verify no repo caller relies on the implicit default: `grep -rn "make_pade_stepper(" libs projects validation apps benchmarks --include="*.py"` — every call passes `order` (or is the new test below).
- [ ] Add to `libs/qscat/tests/test_pade.py`:

  ```python
  def test_default_order_is_three() -> None:
      """lib-m10: the default must be the order the TD validation needed,
      not the CN order the module docstring warns under-converges."""
      import inspect

      assert inspect.signature(make_pade_stepper).parameters["order"].default == 3
  ```

- [ ] Run — expect FAIL (default is 1).
- [ ] Change the default to `order: int = 3`; append to the `make_pade_stepper` docstring: `Default order 3: order 1 (Crank-Nicolson) is documented above to under-converge long propagations; pass order=1 (or use make_sparse_cn_stepper) to get CN explicitly.`
- [ ] Run — expect PASS. Full gates.
- [ ] Commit: `git add libs/qscat/qscat/evolution/pade.py libs/qscat/tests/test_pade.py`; message `fix(evolution): default make_pade_stepper to the converged order 3`.

## Task 4: contextvars-based default backend (lib-m18)

**Files:**
- `libs/qscat/qscat/linalg/sparse_lu.py` (the `_DEFAULT_BACKEND` global ~line 99, `set_default_backend`/`get_default_backend`/`default_backend` ~lines 102–141, and the single `"auto"`-resolution read site inside `SparseLU`)
- `libs/qscat/tests/test_sparse_lu.py` (new test)

**Interfaces (names unchanged — see decided convention 7):**

```python
_DEFAULT_BACKEND: ContextVar[_Backend]  # replaces the module-global str
def set_default_backend(name: _Backend) -> None
def get_default_backend() -> _Backend
@contextmanager
def default_backend(name: _Backend) -> Iterator[None]  # the recommended API
```

**Steps:**

- [ ] Re-verify the current shape: a plain module-global `_DEFAULT_BACKEND: _Backend = "auto"`, the three functions, and exactly which line inside `SparseLU` reads it (`grep -n "_DEFAULT_BACKEND\|get_default_backend" libs/qscat/qscat/linalg/sparse_lu.py`).
- [ ] Add to `libs/qscat/tests/test_sparse_lu.py`:

  ```python
  def test_default_backend_is_context_local() -> None:
      """lib-m18: a default_backend(...) block in one thread must not leak
      into a concurrently running thread's "auto" resolution."""
      import threading

      from qscat.linalg import default_backend, get_default_backend

      seen: list[str] = []
      inside = threading.Event()
      release = threading.Event()

      def forcer() -> None:
          with default_backend("scipy"):
              inside.set()
              release.wait(timeout=10.0)

      t = threading.Thread(target=forcer)
      t.start()
      assert inside.wait(timeout=10.0)
      seen.append(get_default_backend())  # main thread, while forcer holds "scipy"
      release.set()
      t.join()
      assert seen == ["auto"]  # a process-global would have leaked "scipy"
      assert get_default_backend() == "auto"
  ```

- [ ] Run — expect FAIL (`seen == ["scipy"]` with the process-global). Note: a fresh `threading.Thread` copies the context at start; the forcer sets its var AFTER starting, inside its own context, so the main thread must not see it once the var is a ContextVar.
- [ ] Implement in `sparse_lu.py`:

  ```python
  from contextvars import ContextVar

  # Context-local default that backend="auto" resolves against. A ContextVar
  # rather than a module global so a scoped default_backend(...) block in one
  # thread (or async task) cannot leak into another's "auto" resolution.
  _DEFAULT_BACKEND: ContextVar[_Backend] = ContextVar(
      "qscat_sparse_lu_default_backend", default="auto"
  )


  def _validate_backend(name: _Backend) -> None:
      if name not in ("auto", "scipy", "mumps"):
          raise ValueError(f"unknown backend {name!r}; expected auto/scipy/mumps")


  def set_default_backend(name: _Backend) -> None:
      """Set the default backend `SparseLU(backend="auto")` resolves to.

      HAZARD: this mutates the CURRENT context for the rest of the process
      (or thread/task) lifetime and is easy to leave flipped — prefer the
      `default_backend` context manager, which restores the previous value
      on exit (including on exception). Only `"auto"` call sites consult
      this; an explicit `backend="scipy"`/`"mumps"` argument always wins.
      """
      _validate_backend(name)
      _DEFAULT_BACKEND.set(name)


  def get_default_backend() -> _Backend:
      """The current default backend (see `set_default_backend`)."""
      return _DEFAULT_BACKEND.get()


  @contextmanager
  def default_backend(name: _Backend) -> Iterator[None]:
      """Temporarily force the `"auto"` backend to `name` within a `with` block.

      The recommended way to steer internal `SparseLU(...)` construction
      (e.g. a backend-equivalence check over a whole computation): scoped,
      exception-safe, and context-local, so concurrent threads cannot race
      each other's defaults.
      """
      _validate_backend(name)
      token = _DEFAULT_BACKEND.set(name)
      try:
          yield
      finally:
          _DEFAULT_BACKEND.reset(token)
  ```

  Update the `"auto"` read site inside `SparseLU` to `get_default_backend()` (it may already call it; if it read the global directly, switch it). Update the old "Caveat: this is a plain process-global, not thread-local" comment block — that caveat is now false; replace it with the ContextVar rationale above.
- [ ] Run the new test — expect PASS. Also re-run `libs/qscat/tests/test_mumps_backend.py` (its skip logic touches backends). Full gates.
- [ ] Commit: `git add libs/qscat/qscat/linalg/sparse_lu.py libs/qscat/tests/test_sparse_lu.py`; message `fix(linalg): make the default-backend override context-local`.

## Task 5: dead surface, part 1 — `propose_grid` `rtol`, `_CoordinateSpec.charge`, `element_ranges` (lib-M16)

**Files:**
- `libs/qscat/qscat/tuning/propose.py` (`propose_grid` signature ~line 352, `del rtol` ~line 444; `_CoordinateSpec` ~lines 141–160 and its two constructor sites ~lines 168, 208)
- `libs/qscat/qscat/dvr/grid.py` (`element_ranges` attribute ~line 85, its build loop ~lines 134–145, assignment ~line 172)
- `.claude/skills/discretisation-tuner/SKILL.md` (lines 47 and 66 document/pass `rtol`)
- `libs/qscat/tests/test_tuning_propose.py`, `libs/qscat/tests/test_femdvr_ecs.py` (tests)

**Interfaces:**
- `propose_grid(model, coordinate, energy_range, *, rtol: float | None = None, incident=None, phase_coeff=None, channel="ve", elec_grids=None, resonance_n_dense=25) -> FemDvrEcsGrid` — `rtol` kept ONE cycle as an accepted-and-warned no-op (it is public API documented by the discretisation-tuner skill), removed next cycle.
- `_CoordinateSpec` loses its `charge: int` field (private dataclass; zero readers — verified: `grep -n "charge" libs/qscat/qscat/tuning/propose.py` shows only the field, two writes, and comments; the k-computation consumes `model.charge` via the adapters' other outputs, and `probes.py`'s `charge=` kwargs are its own separate API).
- `FemDvrEcsGrid` loses the `element_ranges` attribute (zero consumers — verified: `grep -rn "element_ranges" libs projects validation apps benchmarks --include="*.py"` matches only `dvr/grid.py` itself; all real consumers use `element_maps`).

**Steps:**

- [ ] Re-verify the three greps above still hold (kernel-consolidation may have touched `dvr/grid.py`).
- [ ] Add tests:

  In `test_tuning_propose.py`:

  ```python
  def test_rtol_is_deprecated_noop() -> None:
      """lib-M16: propose_grid never consumed rtol (`del rtol`); passing it
      now warns, omitting it is silent."""
      import warnings

      with pytest.warns(DeprecationWarning, match="rtol"):
          g1 = propose_grid(F2, "nuclear", (0.01, 0.05), rtol=1e-3)
      with warnings.catch_warnings():
          warnings.simplefilter("error")
          g2 = propose_grid(F2, "nuclear", (0.01, 0.05))
      assert np.array_equal(g1.points, g2.points)  # a no-op either way
  ```

  In `test_femdvr_ecs.py`:

  ```python
  def test_element_ranges_removed() -> None:
      """lib-M16: element_ranges had zero consumers; element_maps is the API."""
      g = FemDvrEcsGrid(GridSpec(elements=[ElementSpec(0.0, 4.0, 0.0)], quadrature=5))
      assert not hasattr(g, "element_ranges")
      assert len(g.element_maps) == 1
  ```

  (Match the surrounding tests' actual grid-construction helper — re-verify `GridSpec`/`ElementSpec` call shapes in that file and reuse whatever one-element fixture already exists there.)
- [ ] Run both — expect FAIL (no warning; attribute exists).
- [ ] Implement:
  - `propose_grid`: change `rtol: float = 1e-3` to `rtol: float | None = None`; delete `del rtol  # interface parity...` (~line 444); at the top of the body add:

    ```python
    if rtol is not None:
        warnings.warn(
            "propose_grid(rtol=...) was never consumed (the a-priori assembler "
            "has no eigensolve to converge) and will be removed; pass rtol to "
            "the probe/refine loop instead.",
            DeprecationWarning,
            stacklevel=2,
        )
    ```

    (`import warnings` if the module lacks it.) Delete the docstring paragraph "`rtol` is accepted for interface parity ..." and replace with one line noting the deprecation.
  - `_CoordinateSpec`: delete the `charge: int` field, the `charge=0,`/`charge=model.charge,` constructor lines, and the docstring sentence "`charge` is not consumed by this pipeline (...)".
  - `dvr/grid.py`: delete the `element_ranges` class-level annotation, the `element_ranges: list[tuple[int, int]] = []` build list, its `append`, its assignment `self.element_ranges = element_ranges`, and the `# element_ranges: half-open [start, stop) ...` comment.
  - `.claude/skills/discretisation-tuner/SKILL.md`: drop `rtol=1e-3` from the line-47 signature quote and `rtol=rtol,` from the line-66 call example (the skill's probe steps keep their own rtol usage — only `propose_grid` stops taking it).
- [ ] Run the new tests — expect PASS. Full gates.
- [ ] Commit: `git add libs/qscat/qscat/tuning/propose.py libs/qscat/qscat/dvr/grid.py .claude/skills/discretisation-tuner/SKILL.md libs/qscat/tests/test_tuning_propose.py libs/qscat/tests/test_femdvr_ecs.py`; message `refactor: remove dead tuning/dvr surface (rtol no-op, unused charge, element_ranges)`.

## Task 6: dead surface, part 2 — `correlation.py` unused `grid` parameters (lib-M16)

**Files:**
- `libs/qscat/qscat/core/correlation.py` (`hankel_point_value` ~line 189, `outgoing_surface_wave` ~line 225)
- `libs/qscat/qscat/core/td_extractors.py` (the four call sites: ~lines 483, 567, 865, 950)
- `libs/qscat/tests/test_correlation.py` (tests)

**Interfaces:**
- New: `hankel_point_value(z_position: float, k: float, l: int, charge: int = 0, *, mass: float = 1.0) -> complex`
- New: `outgoing_surface_wave(z_surface: float, k: float, l: int, charge: float = 0.0, *, mass: float = 1.0) -> tuple[complex, complex]`
- Both names are public in `qscat.core.__all__`, so the old grid-first call form is accepted for one cycle: a leading `FemDvrEcsGrid` first argument warns and shifts.

Both functions document their `grid` as "accepted (unused) to keep this call-compatible with `_regular_coeffs`/`_outgoing_coeffs`" — but no call site dispatches them interchangeably through a common variable (verified: all four consumers call them by name), so the parity argument is stylistic and the parameter is dead weight in a public signature.

**Steps:**

- [ ] Re-verify both signatures and that the only in-repo callers are the four td_extractors sites (`grep -rn "hankel_point_value(\|outgoing_surface_wave(" libs projects validation apps --include="*.py"`).
- [ ] Add to `test_correlation.py` (adapt the existing grid fixture in that file):

  ```python
  def test_point_value_functions_drop_grid_param() -> None:
      """lib-M16: the grid argument was documented-unused; the new signature
      drops it, the old grid-first form warns for one cycle."""
      k, l = 0.7, 0
      new = hankel_point_value(3.0, k, l)
      with pytest.warns(DeprecationWarning, match="grid"):
          old = hankel_point_value(_grid(), 3.0, k, l)  # legacy call form
      assert new == old

      pv, dv = outgoing_surface_wave(3.0, k, l)
      with pytest.warns(DeprecationWarning, match="grid"):
          pv2, dv2 = outgoing_surface_wave(_grid(), 3.0, k, l)
      assert (pv, dv) == (pv2, dv2)
  ```

- [ ] Run — expect FAIL (TypeError/no warning).
- [ ] Implement in `correlation.py` — the one-cycle positional-shift shim. Concretely for `hankel_point_value` (current body moves verbatim into the private worker, minus the grid parameter):

  ```python
  def _hankel_point_value(z_position: float, k: float, l: int, charge: int, mass: float) -> complex:
      # today's hankel_point_value body, unchanged, with the grid parameter gone
      if charge == 0:
          return riccati_hankel_en_mass(z_position, k, l, mass) / 2.0
      return coulomb_h1_en(z_position, k, charge, mass, l) / 2.0


  def hankel_point_value(
      z_position: float | FemDvrEcsGrid,
      k: float,
      l: float | int = 0,
      charge: float | int = 0,
      *,
      mass: float = 1.0,
      _legacy_charge: int = 0,
  ) -> complex:
      """<current docstring with the grid sentence replaced by the new
      signature and a one-line deprecation note for the old grid-first form>"""
      if isinstance(z_position, FemDvrEcsGrid):
          warnings.warn(
              "hankel_point_value no longer takes a grid argument (it was "
              "documented-unused); drop the first argument",
              DeprecationWarning,
              stacklevel=2,
          )
          # legacy form: (grid, z_position, k, l, charge=0, *, mass=1.0) --
          # every argument sits one slot to the right of its new home.
          return _hankel_point_value(float(k), float(l), int(charge), _legacy_charge, mass)
      return _hankel_point_value(float(z_position), float(k), int(l), int(charge), mass)
  ```

  The `_legacy_charge` keyword exists only so the legacy 5-positional call `(grid, z, k, l, charge)` still binds (its `charge` lands there); it is undocumented and dies with the shim next cycle. (Re-verify `_hankel_point_value`'s body against the actual current implementation before moving it — the two-branch return above is the 2026-08-25 shape.) Apply the same pattern to `outgoing_surface_wave` (`z_surface: float | FemDvrEcsGrid` first, `_legacy_charge: float = 0.0`, returning the `(phi_out, dphi_out)` tuple from its own private worker). The public signature documented in each docstring is the NEW one; both legacy tests above must pass. Next cycle both shims collapse to the plain new signatures.
- [ ] Update the four `td_extractors.py` call sites to the new form (delete the grid first argument: e.g. line 483 `f_i = hankel_point_value(z_position, kp, model.ell, model.charge)`, line 567 `f_c = hankel_point_value(R_position, k_r, 0, model.charge, mass=mu_r)`, lines 865/950 likewise for `outgoing_surface_wave`). Update the docstring references in `td_extractors.py` (lines ~28, 44, 65, 95, 632, 847) and `core/__init__.py` (lines ~97–99) that quote the old call shape.
- [ ] Run the new test + `test_td_extractors.py` — expect PASS, no DeprecationWarning from the library's own call sites (`-W error::DeprecationWarning` on the td_extractors run proves it). Full gates.
- [ ] Commit: `git add libs/qscat/qscat/core/correlation.py libs/qscat/qscat/core/td_extractors.py libs/qscat/qscat/core/__init__.py libs/qscat/tests/test_correlation.py`; message `refactor(core): drop the documented-unused grid parameter from the point-value helpers`.

## Task 7: promote the TD building blocks; drop `_propagate`'s dead `eps` (lib-m7 + lib-M16)

**Files:**
- `libs/qscat/qscat/core/time_dependent.py` (`_propagate` ~line 276, `_s_vector_one_energy` ~line 313, `_sigma_one_energy` ~line 350, `__all__`, new module `__getattr__`)
- `projects/n2_2d_td_cross_section/td_cross_section.py` (imports ~lines 18–23, the `_propagate` shim ~lines 40–66)
- `projects/n2_2d_td_cross_section/td_propagation.py` (~line 49), `projects/n2_2d_td_cross_section/test_td_cross_section.py` (~lines 92, 235, 247) — the shim's callers
- `libs/qscat/tests/test_core_td.py` (new test)

**Interfaces:**
- Produces (public, added to `__all__`):

  ```python
  def propagate_wavepacket(
      tgrid: TensorGrid,
      model: ResonanceModel,
      chi: npt.NDArray[np.complex128],
      v_init: int,
      vprimes: list[int],
      *,
      dt: float,
      n_steps: int,
      wp_in: _WpIn,
      wp_out: _WpOut,
      free: bool = False,
      order: int = 3,
  ) -> PropagationResult
  ```

  (= today's `_propagate` with the dead `eps: npt.NDArray[np.float64]` positional REMOVED — verified 2026-08-25: the body never reads `eps`; `psi0`/`out_channels`/`hamiltonian` use only `chi`, `wp_*`, `model`.)
  `s_vector_one_energy` and `sigma_one_energy` are `_s_vector_one_energy`/`_sigma_one_energy` renamed, signatures otherwise unchanged (both genuinely consume their `eps`).
- Keeps for one cycle: module-level `__getattr__` serving `_propagate`, `_s_vector_one_energy`, `_sigma_one_energy` with a DeprecationWarning. NOTE `_propagate` served this way still expects the OLD arg order (with `eps`) at its remaining external callers — the shim therefore adapts:

**Steps:**

- [ ] Re-verify: `td_extractors.py` imports the PUBLIC `quadrature_weights` (the review's `_quadrature_weights` half of lib-m7 is already resolved — confirm with `grep -n "quadrature_weights" libs/qscat/qscat/core/td_extractors.py`, expect no leading underscore) and the kernel-consolidation plan did not already move `_propagate` (`grep -n "_propagate" libs/qscat/qscat/core/time_dependent.py`).
- [ ] Add to `libs/qscat/tests/test_core_td.py`:

  ```python
  def test_td_building_blocks_are_public() -> None:
      """lib-m7: the names the n2_2d_td project consumes must be public;
      the old private names warn for one cycle."""
      from qscat.core import time_dependent as td

      for name in ("propagate_wavepacket", "s_vector_one_energy", "sigma_one_energy"):
          assert name in td.__all__
          assert callable(getattr(td, name))
      with pytest.warns(DeprecationWarning, match="_propagate"):
          legacy = td._propagate
      assert callable(legacy)
  ```

- [ ] Run — expect FAIL.
- [ ] Implement in `time_dependent.py`:
  - Rename `def _propagate(...)` → `def propagate_wavepacket(...)` deleting the `eps` parameter; rename `def _s_vector_one_energy` → `def s_vector_one_energy`, `def _sigma_one_energy` → `def sigma_one_energy`. Update every internal call site (the `td_ve_cross_section`/`td_da_cross_section` bodies) to the new names and to the eps-free `propagate_wavepacket` call. Extend each docstring's first line to note it is a public building block consumed by `projects/n2_2d_td_cross_section` (the module's `free_hamiltonian` precedent).
  - Add the three names to `__all__`.
  - Add the one-cycle shim — this is the canonical `__getattr__` deprecation pattern for the whole plan (referenced by later tasks, shown fully once here):

    ```python
    # One deprecation cycle for the pre-promotion private names (lib-m7). The
    # public defs above are the real objects; this only serves old imports.
    _DEPRECATED_ALIASES = {
        "_s_vector_one_energy": "s_vector_one_energy",
        "_sigma_one_energy": "sigma_one_energy",
    }


    def _propagate_legacy(
        tgrid: TensorGrid,
        model: ResonanceModel,
        eps: npt.NDArray[np.float64],
        chi: npt.NDArray[np.complex128],
        v_init: int,
        vprimes: list[int],
        *,
        dt: float,
        n_steps: int,
        wp_in: _WpIn,
        wp_out: _WpOut,
        free: bool = False,
        order: int = 3,
    ) -> PropagationResult:
        """Old `_propagate` call shape; `eps` was never read."""
        del eps
        return propagate_wavepacket(
            tgrid,
            model,
            chi,
            v_init,
            vprimes,
            dt=dt,
            n_steps=n_steps,
            wp_in=wp_in,
            wp_out=wp_out,
            free=free,
            order=order,
        )


    def __getattr__(name: str) -> object:
        if name == "_propagate":
            warnings.warn(
                "qscat.core.time_dependent._propagate is deprecated; use the "
                "public propagate_wavepacket (note: the unused eps argument "
                "is gone from the new signature)",
                DeprecationWarning,
                stacklevel=2,
            )
            return _propagate_legacy
        if name in _DEPRECATED_ALIASES:
            new = _DEPRECATED_ALIASES[name]
            warnings.warn(
                f"qscat.core.time_dependent.{name} is deprecated; use the "
                f"public {new}",
                DeprecationWarning,
                stacklevel=2,
            )
            return globals()[new]
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    ```

  - Check the existing `_free_hamiltonian = free_hamiltonian` alias still stands (its consumers include `libs/qscat/tests/test_td_extractors.py`) — leave it; it predates this plan and is a plain alias, not a finding.
- [ ] Update `projects/n2_2d_td_cross_section/td_cross_section.py`: imports become `from qscat.core.time_dependent import propagate_wavepacket as _core_propagate`, `s_vector_one_energy as _core_s_vector_one_energy`, `sigma_one_energy as _core_sigma_one_energy` (PropagationResult/sigma_from_correlations/td_ve_cross_section imports unchanged). Drop `eps` from the local `_propagate` shim's signature, its `_core_propagate(...)` forward, and its docstring; update the shim's callers (`td_propagation.py:49`, `test_td_cross_section.py:92,235,247`) to stop passing `eps`.
- [ ] Run the new test, `libs/qscat/tests/test_core_td.py`, `libs/qscat/tests/test_td_extractors.py`, and `uv run --no-sync pytest projects/n2_2d_td_cross_section -m "not slow" -q` — expect PASS with no DeprecationWarning raised by in-repo callers. Full gates.
- [ ] Commit: `git add libs/qscat/qscat/core/time_dependent.py projects/n2_2d_td_cross_section/td_cross_section.py projects/n2_2d_td_cross_section/td_propagation.py projects/n2_2d_td_cross_section/test_td_cross_section.py libs/qscat/tests/test_core_td.py`; message `refactor(core): promote the TD building blocks to public names, drop dead eps`.

## Task 8: missing docstrings, part 1 — `core/` (lib-M13 groundwork)

**Files:** `libs/qscat/qscat/core/assignment.py`, `core/bo.py`, `core/nrm/extended.py`, `core/problem.py`, `core/time_dependent.py`.

**Interfaces:** none change — docstrings only. The authoritative site list (measured 2026-08-25, `ruff check --select D101,D102,D105 libs/qscat/qscat` = 25 errors; line numbers will have drifted — re-run the command to refresh them):

| Site | Kind | Docstring to add |
|---|---|---|
| `assignment.py:157` `shift_mev` | D102 | `"""`shift_ev` in milli-eV -- the unit the reference tables quote."""` |
| `assignment.py:468` `__str__` | D105 | `"""One-line summary: marks vs peaks and the median distance in widths."""` |
| `bo.py:130` `n_curves` | D102 | `"""Number of electronic curves (rows of `energies`)."""` |
| `bo.py:134` `has_states` | D102 | `"""True when the curves were built `with_states=True`."""` |
| `bo.py:395` `__contains__` | D105 | `"""`(curve, v) in basis` membership test."""` |
| `bo.py:398` `__getitem__` | D105 | `"""The `BoState` stored under `(curve, v)`."""` |
| `bo.py:401` `__len__` | D105 | `"""Number of BO product states in the basis."""` |
| `bo.py:404` `items` | D102 | `"""All `((curve, v), BoState)` pairs, dict-style."""` |
| `bo.py:408` `has_states` | D102 | `"""True when the basis holds at least one state."""` |
| `nrm/extended.py:252` `rank` | D102 | `"""Number of separable-expansion vectors (columns of `vectors`)."""` |
| `problem.py:72` `__post_init__` | D105 | `"""Solve the vibrational basis once at construction (frozen dataclass)."""` |
| `time_dependent.py:124` `Snapshot` | D101 | `"""One observation of the propagated state at `time`: the two marginal densities, plus the full `psi` when requested (see `propagate`'s `snapshot_times`)."""` |
| `time_dependent.py:132` `PropagationResult` | D101 | `"""What `propagate` records per step: sample times, channel projections `c_{v'}(t_n)`, Hermitian norms, and any snapshots."""` |

(`time_dependent.py:109/111` — the `Extractor` protocol methods — were already covered in Task 2.)

**Steps:**

- [ ] Re-verify the list: `uv run --no-sync ruff check --select D101,D102,D105 libs/qscat/qscat --output-format concise` (expect the sites above minus Task 2's two, plus the Task 9 sites). Adjust wording where surrounding code moved.
- [ ] Failing check = that ruff command (non-empty output for these files).
- [ ] Add the docstrings from the table (adapt phrasing where a docstring reads wrong against the actual body after re-verification — the table text was written from the 2026-08-25 bodies).
- [ ] Re-run the ruff command — the `core/` sites are gone. Full gates (pytest fast tier especially: docstring insertion into a `Protocol`/`@dataclass` must not change behaviour).
- [ ] Commit: `git add libs/qscat/qscat/core/assignment.py libs/qscat/qscat/core/bo.py libs/qscat/qscat/core/nrm/extended.py libs/qscat/qscat/core/problem.py libs/qscat/qscat/core/time_dependent.py`; message `docs(core): add the missing public docstrings (ruff D groundwork)`.

## Task 9: missing docstrings, part 2 — `dvr/`, `linalg/`, `model/` (lib-M13 groundwork)

**Files:** `libs/qscat/qscat/dvr/spec.py`, `dvr/tensor.py`, `linalg/eigs.py`, `linalg/sparse_lu.py`, `model/diatomic.py`.

**Interfaces:** none change — docstrings only. Sites (same measurement; re-verify as in Task 8):

| Site | Kind | Docstring to add |
|---|---|---|
| `dvr/spec.py:52` `__post_init__` | D105 | `"""Validate the spec (quadrature >= 2, ordered elements) and derive `R0`."""` (read the body first; describe what it actually checks) |
| `dvr/tensor.py:47` `grids` | D102 | `"""The per-dimension `FemDvrEcsGrid`s, in tensor order."""` |
| `dvr/tensor.py:51` `ndim` | D102 | `"""Number of tensor dimensions `D`."""` |
| `dvr/tensor.py:55` `shape` | D102 | `"""Per-dimension point counts `(n_1, ..., n_D)`."""` |
| `dvr/tensor.py:59` `size` | D102 | `"""Total number of tensor-product points, `prod(shape)`."""` |
| `linalg/eigs.py:188` `shape` | D102 | `"""Operator shape `(n, n)`."""` |
| `linalg/sparse_lu.py:314` `shape` | D102 | `"""Shape of the factored matrix."""` |
| `linalg/sparse_lu.py:318` `ordering` | D102 | `"""The `permc_spec` column ordering this factorization was built with."""` |
| `model/diatomic.py:67` `mu` | D102 | `"""Nuclear reduced mass (a.u.)."""` |
| `model/diatomic.py:69` `ell` | D102 | `"""Resonance partial-wave angular momentum `l`."""` |

**Steps:**

- [ ] Re-verify sites (`ruff check --select D101,D102,D105` again).
- [ ] Add the docstrings (for the `model/diatomic.py` protocol properties the `...` body is replaced by the docstring — a docstring IS a valid protocol-stub body).
- [ ] `uv run --no-sync ruff check --select D101,D102,D105 libs/qscat/qscat` — expect ZERO findings now (Tasks 2+8+9 cover all 25).
- [ ] Full gates.
- [ ] Commit: `git add libs/qscat/qscat/dvr/spec.py libs/qscat/qscat/dvr/tensor.py libs/qscat/qscat/linalg/eigs.py libs/qscat/qscat/linalg/sparse_lu.py libs/qscat/qscat/model/diatomic.py`; message `docs: complete the missing dvr/linalg/model docstrings`.

## Task 10: settle numpydoc — ruff `D` + `ANN`, drop napoleon-google (lib-M13, style-N2/N3)

**Files:**
- `pyproject.toml` (`[tool.ruff.lint]` ~line 91, `[tool.ruff.lint.per-file-ignores]` ~line 98, new `[tool.ruff.lint.pydocstyle]`)
- `docs/conf.py` (line ~41 `napoleon_google_docstring = True`)

**Interfaces:** tool configuration only. Decided scope (conventions 4–5 above): `D` enforced for `libs/qscat/qscat` only; `ANN` enforced repo-wide with ANN401 ignored.

**Steps:**

- [ ] Re-verify current counts: `uv run --no-sync ruff check --select D --config 'lint.pydocstyle.convention="numpy"' --statistics libs/qscat/qscat` must show ONLY style rules (D205/D209/D400/D401/D403) after Tasks 8–9 (no D1xx), and `--select ANN --statistics libs/qscat/qscat` must show only ANN401.
- [ ] Edit `pyproject.toml`:

  ```toml
  [tool.ruff.lint]
  select = ["E", "F", "I", "UP", "B", "NPY", "RUF", "D", "ANN"]
  # RUF001-003 (ambiguous unicode) are off: the docstrings and comments here
  # deliberately use unicode maths characters (sigma, multiplication sign, en
  # dash) that these rules would flag as confusable with ASCII.
  # D205/D209/D400/D401/D403 are off: the library's house style opens with
  # multi-sentence summary paragraphs (measured 155 hits, all deliberate
  # prose); D here enforces docstring PRESENCE, not summary-line shape.
  # ANN401 is off: the two Any boundaries are deliberate -- core/problem.py's
  # **kwargs pass-through facade and the lazily-imported-matplotlib viz
  # modules -- and mypy --strict polices annotations everywhere else.
  ignore = ["RUF001", "RUF002", "RUF003", "D205", "D209", "D400", "D401", "D403", "ANN401"]

  [tool.ruff.lint.pydocstyle]
  convention = "numpy"
  ```

  and append to the existing `[tool.ruff.lint.per-file-ignores]` table:

  ```toml
  # Docstring presence (D) is a shipped-library contract; tests, research
  # projects, validation harnesses and app internals are exempt. ANN beyond
  # the annotations mypy --strict already enforces is likewise lib-only.
  "!libs/qscat/qscat/**" = ["D", "ANN"]
  ```

  (Negated per-file-ignore patterns are supported by the pinned ruff 0.15.21. If `ruff check .` still reports D/ANN outside the lib after this, fall back to enumerating: `"libs/qscat/tests/**" = ["D", "ANN"]`, `"projects/**" = [...]`, `"validation/**"`, `"apps/**"`, `"benchmarks/**"`, `"docs/**"`, `"native/**"`.)
- [ ] Edit `docs/conf.py`: change `napoleon_google_docstring = True` to `napoleon_google_docstring = False` with the comment `# numpydoc is the ONE docstring dialect (ruff D convention = "numpy"); google-style parsing off so a stray google-section renders wrong loudly instead of silently working.` Keep `napoleon_numpy_docstring = True`.
- [ ] Run `uv run --no-sync ruff check .` — expect clean. If new D/ANN findings surface inside `libs/qscat/qscat` that Tasks 8–9 didn't cover (rules beyond the measured set), fix them here if they are presence-rules, or add the specific rule to `ignore` with a one-line comment if they are style-rules — do not add per-file escapes for lib files.
- [ ] Build docs if the docs workflow builds locally (`uv run --no-sync sphinx-build -W docs docs/_build` or the repo's documented equivalent — re-verify the command in docs/README or the workflow file; skip if docs need extras absent locally, and note it in the commit message).
- [ ] Full gates.
- [ ] Commit: `git add pyproject.toml docs/conf.py`; message `chore: enforce numpydoc presence (ruff D) and ANN, retire google-docstring parsing`.

## Task 11: `DrResult` for the DR solver (lib-m14)

**Files:**
- `libs/qscat/qscat/core/dissociation.py` (`dr_cross_section` ~lines 294–469)
- `libs/qscat/qscat/core/__init__.py` (export)
- `libs/qscat/tests/test_dissociation.py` (tests; flag-call sites at lines ~189, 266, 285 migrate)

**Interfaces:**

```python
@dataclass(frozen=True)
class DrResult:
    """Result of the exact 2-D dissociative-recombination solve (`dr_solve`).

    - `sigma`: sigma_DR(E) in bohr^2 -- `(n_channels,)` for scalar `E`,
      `(len(E), n_channels)` for array `E`.
    - `psi`: the driven `Psi+` per energy (`None` when not stored, and `None`
      per energy below threshold) -- one array for scalar `E`, one list entry
      per energy for array `E`.
    - `amplitude`: the complex T-matrix amplitude, shaped like `sigma`
      (`None` when not stored). See `dr_solve` for the S-vs-T normalization
      note.
    """

    sigma: _Sigma
    psi: _PsiOut | None
    amplitude: _Amp | None


def dr_solve(
    tgrid: TensorGrid,
    model: ResonanceModel,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    E: float | npt.ArrayLike,
    *,
    n_channels: int = 3,
    ordering: Ordering = "COLAMD",
    store_wavefunction: bool = False,
    store_amplitude: bool = False,
) -> DrResult
```

`dr_cross_section` keeps its exact current signature AND its four overloads, becomes a thin delegate to `dr_solve`, and emits a DeprecationWarning ONLY when `return_wavefunction` or `return_amplitude` is passed True (the plain sigma-only call is the undisputed base case and stays silent and un-deprecated). `ve_cross_section`/`da_cross_section`/`lcp_da_cross_section` are deliberately NOT converted (decided convention 6).

**Steps:**

- [ ] Re-verify the four overloads and the tail-of-body flag dispatch (`if return_wavefunction and return_amplitude: ...`) still match the shape read on 2026-08-25; re-verify the flag consumers are only `libs/qscat/tests/test_dissociation.py` (lines ~189/266/285) plus the `**kwargs` facade `core/problem.py:100` and the flag-less `apps/qscat-run/qscat_run/runner.py:392`.
- [ ] Add to `test_dissociation.py` (reuse the small H2P fixture already used at line ~177):

  ```python
  def test_dr_solve_returns_dataclass_result() -> None:
      """lib-m14: one result object instead of four flag-shaped tuples."""
      from qscat.core.dissociation import DrResult, dr_solve

      res = dr_solve(tg, H2P, eps, chi, 0, E, n_channels=2)
      assert isinstance(res, DrResult)
      assert res.psi is None and res.amplitude is None

      full = dr_solve(
          tg, H2P, eps, chi, 0, E, n_channels=2,
          store_wavefunction=True, store_amplitude=True,
      )
      np.testing.assert_array_equal(full.sigma, res.sigma)
      assert full.psi is not None and full.amplitude is not None
      assert full.amplitude.shape == full.sigma.shape


  def test_dr_cross_section_flags_deprecated_but_working() -> None:
      s = dr_cross_section(tg, H2P, eps, chi, 0, E, n_channels=2)  # silent
      with pytest.warns(DeprecationWarning, match="dr_solve"):
          s2, psi = dr_cross_section(
              tg, H2P, eps, chi, 0, E, n_channels=2, return_wavefunction=True
          )
      np.testing.assert_array_equal(s, s2)
  ```

  (Build `tg/eps/chi/E` exactly as the existing `dr_cross_section` tests in that file do — copy their fixture lines, do not invent a new deck.)
- [ ] Run — expect FAIL (no `dr_solve`).
- [ ] Implement: move the current `dr_cross_section` body into `dr_solve` (loop unchanged; it already fills `out`, `amp`, `psi_list` unconditionally except that `psi_list` costs memory — gate `psi_list[ie] = psi_plus` on `store_wavefunction` and `amp[ie, n] = t` stays cheap and unconditional), ending with:

  ```python
  scalar = np.isscalar(E) or (isinstance(E, np.ndarray) and np.ndim(E) == 0)
  sigma = np.asarray(out[0] if scalar else out, dtype=np.float64)
  return DrResult(
      sigma=sigma,
      psi=(psi_list[0] if scalar else psi_list) if store_wavefunction else None,
      amplitude=(np.asarray(amp[0] if scalar else amp, dtype=np.complex128)
                 if store_amplitude else None),
  )
  ```

  Move the physics docstring to `dr_solve`; `dr_cross_section` keeps a short docstring pointing there plus the deprecation note, and its implementation becomes:

  ```python
  def dr_cross_section(..., return_wavefunction: bool = False, return_amplitude: bool = False):
      if return_wavefunction or return_amplitude:
          warnings.warn(
              "dr_cross_section's flag-shaped tuple returns are deprecated; "
              "call dr_solve(..., store_wavefunction=..., store_amplitude=...) "
              "and read the DrResult fields",
              DeprecationWarning,
              stacklevel=2,
          )
      res = dr_solve(tgrid, model, eps, chi, v_init, E, n_channels=n_channels,
                     ordering=ordering, store_wavefunction=return_wavefunction,
                     store_amplitude=return_amplitude)
      if return_wavefunction and return_amplitude:
          return res.sigma, res.psi, res.amplitude
      if return_amplitude:
          return res.sigma, res.amplitude
      if return_wavefunction:
          return res.sigma, res.psi
      return res.sigma
  ```

  Export `DrResult` and `dr_solve` from `qscat.core.dissociation.__all__` and `qscat/core/__init__.py` (imports + `__all__`).
- [ ] Migrate the in-repo flag callers: `test_dissociation.py` lines ~189/266/285 switch to `dr_solve(...)` field access (KEEP the two legacy-shape tests added above as the deprecation-cycle coverage). `runner.py:392` and `problem.py` facade are flag-less/pass-through — leave them.
- [ ] Run `test_dissociation.py` — expect PASS. Full gates.
- [ ] Commit: `git add libs/qscat/qscat/core/dissociation.py libs/qscat/qscat/core/__init__.py libs/qscat/tests/test_dissociation.py`; message `refactor(core): DrResult dataclass for the DR solve; deprecate the flag tuples`.

## Task 12: row-per-state orientation for `ExactResonanceStates` (lib-M9, part 1)

**Files:**
- `libs/qscat/qscat/core/resonance.py` (dataclass docstring ~lines 63–68, `save`/`load` ~lines 103–125, assembly ~lines 292–319)
- `validation/h2plus/bo_overlap.py` (`states[:, i]`-free already? no: cache + return at lines ~152–173; figure use)
- `validation/h2plus/resonance_state_figures.py` (~lines 173, 186)
- `validation/n2/exact_resonance_figures.py` (~line 140)
- `validation/n2/pole_verification.py` (~lines 117, 120, 148)
- `libs/qscat/tests/test_exact_resonance_states.py` (~line 182), `libs/qscat/tests/test_assignment.py` (~lines 144, 166, 187)

**Interfaces:**
- `ExactResonanceStates.states` becomes shape `(m, n)`: `states[i]` is the flattened, c-product-normalized 2-D eigenvector of `energies[i]` (decided convention 1 — row-per-state is the repo-wide container convention: `chi[v]`, `phi[c]`, `ResonanceLevels.states`, qscat-run `EigenStates`).
- `ExactResonanceStates.load` gains a stale-cache guard: an archive whose `states.shape[0] != energies.size` is a pre-flip column-layout cache and must FAIL LOUDLY (a silent transpose would feed garbage overlaps to hours-long H2+ analyses).

**Steps:**

- [ ] Re-verify the assembly block still ends `states = vecs_a[:, idx]` + per-column normalization + the `(vecs_a.shape[0], 0)` empty case, and that the consumer list above is complete: `grep -rn "\.states" libs projects validation apps --include="*.py" | grep -v test | grep -iv "bo\.\|curves\.\|EigenStates\|lv\."` and the test greps.
- [ ] Modify `libs/qscat/tests/test_exact_resonance_states.py`: change line ~182 `psi = res.states[:, i]` to `psi = res.states[i]`, and ADD:

  ```python
  def test_states_are_row_per_state_and_load_rejects_column_caches(tmp_path) -> None:
      """lib-M9: one orientation convention (rows, like chi/phi); pre-flip
      caches must fail loudly, not load transposed."""
      # `res` from the existing small fixture in this file
      assert res.states.shape[0] == res.energies.size
      p = tmp_path / "states.npz"
      res.save(p)
      loaded = type(res).load(p)
      np.testing.assert_array_equal(loaded.states, res.states)
      # forge a pre-flip archive: transpose states
      import numpy as np
      from dataclasses import fields
      legacy = {f.name: getattr(res, f.name) for f in fields(res)}
      legacy["states"] = legacy["states"].T
      np.savez(tmp_path / "legacy.npz", **legacy)
      with pytest.raises(ValueError, match="column-per-state"):
          type(res).load(tmp_path / "legacy.npz")
  ```

  (Requires the fixture's `m != n`, true for any real grid; assert `res.energies.size != res.states.shape[1]` in the test if the fixture is tiny.)
- [ ] Run the file — expect FAIL (shape assertions).
- [ ] Implement in `resonance.py`:
  - Assembly: `states = vecs_a[:, idx].T.copy()` then

    ```python
    for i in range(states.shape[0]):
        states[i] /= np.sqrt(c_product(states[i], states[i]))
    ```

    Empty case: `states=np.empty((0, vecs_a.shape[0]), dtype=np.complex128)`.
  - Dataclass docstring: `states : ndarray of complex128, shape (m, n)` / "`states[i]` is the 2-D eigenvector of `energies[i]` ... row-per-state, the same orientation as `chi`, `anion_electronic_states`, and `ResonanceLevels.states`."
  - `load`: after the missing-field check, add

    ```python
    if z["states"].shape[:1] != z["energies"].shape:
        raise ValueError(
            f"{path} stores states column-per-state (an archive from before "
            "the row-per-state layout); delete and regenerate it"
        )
    ```

- [ ] Update consumers, `states[:, i]` → `states[i]`:
  - `validation/h2plus/resonance_state_figures.py` ~173, 186; `validation/n2/exact_resonance_figures.py` ~140; `validation/n2/pole_verification.py` ~117, 120, 148; `libs/qscat/tests/test_assignment.py` ~144, 166, 187 (`pair_by_overlap(res.energies[i], res.states[i], basis)` — `pair_by_overlap` itself takes a 1-D vector and is unchanged).
  - `validation/h2plus/bo_overlap.py`: replace the hand-rolled cache (`np.load`/`np.savez` of `energies`+`states` only, lines ~152–155, 171) with the round-trip API this flip is exactly the hazard case for: write `res.save(cache)`; read via

    ```python
    if cache.exists():
        try:
            res = ExactResonanceStates.load(cache)
        except ValueError:
            cache.unlink()  # pre-flip or pre-save/load archive: regenerate
        else:
            return res.energies, res.states, base
    ```

    (old two-field archives fail the missing-field check, old full archives fail the orientation guard — both regenerate instead of silently transposing). Downstream in that file (and `main`) the returned `states` is now indexed `states[i]` — update its overlap loop accordingly (re-verify how `main` consumes the tuple before editing).
- [ ] Run: the new/updated tests, `uv run --no-sync pytest libs/qscat/tests/test_exact_resonance_states.py libs/qscat/tests/test_assignment.py -q`, and `uv run --no-sync python -m validation.n2.pole_verification` if it runs in fast-tier time (it is a driver, not a gate — if it needs the big deck, review the diff by eye instead and say so in the commit message). Full gates.
- [ ] Commit: `git add libs/qscat/qscat/core/resonance.py validation/h2plus/bo_overlap.py validation/h2plus/resonance_state_figures.py validation/n2/exact_resonance_figures.py validation/n2/pole_verification.py libs/qscat/tests/test_exact_resonance_states.py libs/qscat/tests/test_assignment.py`; message `refactor(core): row-per-state ExactResonanceStates with a loud stale-cache guard`.

## Task 13: unclamped widths + `ResonanceLevels.save`/`load` (lib-M9, part 2)

**Files:**
- `libs/qscat/qscat/core/lcp.py` (`ResonanceLevels` ~line 550, the `widths = np.maximum(0.0, -2.0 * energies.imag)` line ~line 780)
- `libs/qscat/tests/test_lcp_resonance_levels.py` (line ~245 `assert np.all(out.widths >= 0.0)  # clamped`, plus new tests)

**Interfaces:**
- `ResonanceLevels.widths` becomes raw `-2.0 * energies.imag` (decided convention 2). Docstring updated: "`widths`: `Gamma_v = -2 Im E_v` (Hartree), UNCLAMPED — after the `Im E <= atol` physicality filter the most negative representable value is `-2*atol`; a small negative width is a round-off diagnostic, and hiding it behind a clamp is how it goes unnoticed. Same convention as `ExactResonanceStates.widths`."
- New methods, mirroring `ExactResonanceStates`' exactly (fields-driven `.npz` round trip):

  ```python
  def save(self, path: str | os.PathLike[str]) -> None
  @classmethod
  def load(cls, path: str | os.PathLike[str]) -> ResonanceLevels
  ```

**Steps:**

- [ ] Re-verify the clamp line and that no in-repo consumer requires `widths >= 0` (`grep -rn "\.widths" libs apps validation projects --include="*.py"` — the known consumers use widths as magnitudes for line-widths/labels/complex reconstruction, all tolerant of a −2e-8 floor; `exact_resonance_figures.py:121` reconstructs `E - i*widths/2`, where the unclamped value is the MORE faithful one).
- [ ] Update/extend `test_lcp_resonance_levels.py`:
  - line ~245: `assert np.all(out.widths >= 0.0)  # clamped` → `assert np.all(out.widths >= -2.0 * 1e-8)  # unclamped: bounded below by the Im-filter's -2*atol` (use the test's actual `atol` variable if it passes one).
  - add:

    ```python
    def test_widths_are_raw_minus_two_im(tmp_path) -> None:
        """lib-M9: widths unclamped (match ExactResonanceStates), and the
        result round-trips through save/load."""
        # `out` from the existing analytic fixture in this file
        np.testing.assert_array_equal(out.widths, -2.0 * out.energies.imag)
        p = tmp_path / "levels.npz"
        out.save(p)
        back = ResonanceLevels.load(p)
        for name in ("energies", "widths", "states", "residuals", "real_weight", "golden_rule"):
            np.testing.assert_array_equal(getattr(back, name), getattr(out, name))
        with pytest.raises(ValueError, match="missing"):
            np.savez(tmp_path / "bad.npz", energies=out.energies)
            ResonanceLevels.load(tmp_path / "bad.npz")
    ```

- [ ] Run — expect FAIL (clamped equality breaks / no `save`).
- [ ] Implement in `lcp.py`: change the widths line to `widths = -2.0 * energies.imag`; update the `ResonanceLevels` docstring bullet as above; add to the dataclass (imports: `os`, and `fields` joins the existing `dataclass` import):

  ```python
  def save(self, path: str | os.PathLike[str]) -> None:
      """Write to a compressed `.npz` under the dataclass's own field names.

      Mirrors `ExactResonanceStates.save` -- the field names stay the
      dataclass's business, so a rename cannot silently desynchronize a
      hand-rolled cache.
      """
      np.savez(path, **{f.name: getattr(self, f.name) for f in fields(self)})

  @classmethod
  def load(cls, path: str | os.PathLike[str]) -> ResonanceLevels:
      """Read back a `save()` file, checking every field is present."""
      with np.load(path) as z:
          missing = [f.name for f in fields(cls) if f.name not in z]
          if missing:
              raise ValueError(f"{path} is not a ResonanceLevels archive: missing {missing}")
          return cls(**{f.name: z[f.name] for f in fields(cls)})
  ```

  (No orientation guard needed here: `ResonanceLevels.states` was ALWAYS row-per-state, and no pre-existing `ResonanceLevels.save` archives exist.)
- [ ] Run the file — expect PASS. Also run `libs/qscat/tests/test_plot_resonance_levels.py` and `test_bo.py` (widths flow into both). Full gates.
- [ ] Commit: `git add libs/qscat/qscat/core/lcp.py libs/qscat/tests/test_lcp_resonance_levels.py`; message `fix(core): unclamp ResonanceLevels.widths; add save/load round trip`.

## Task 14: split `core/lcp.py` into a three-module package (lib-M11) — LAST, biggest churn

**Files:**
- DELETE `libs/qscat/qscat/core/lcp.py` (999 lines pre-plan; Tasks 1, 11, 13 will have edited it — split whatever it then contains)
- NEW `libs/qscat/qscat/core/lcp/__init__.py`
- NEW `libs/qscat/qscat/core/lcp/curve.py`
- NEW `libs/qscat/qscat/core/lcp/cross_section.py`
- NEW `libs/qscat/qscat/core/lcp/levels.py`
- `libs/qscat/tests/test_lcp.py` (one import-surface test added; nothing else changes)

**Interfaces (all preserved — the package `__init__` re-exports the exact current `__all__`):**

```python
__all__ = [
    "ResonanceLevels",
    "lcp_da_cross_section",
    "lcp_resonance_levels",
    "local_complex_potential",
    "resonance_eigenstate",
    "resonance_eigenstate_at_peak_width",
    "resonance_levels",
    "resonance_pole_walk",
]
```

Verified consumer import surface (2026-08-25) — every one keeps working unchanged because a package named `qscat.core.lcp` answers the same `from qscat.core.lcp import X` statements:
`core/__init__.py:199` (5 names), `core/nrm/discrete_state.py:31` (`resonance_pole_walk`), `tuning/resonance.py:40` (`resonance_pole_walk`), `apps/qscat-run/qscat_run/runner.py:86-92` (4 names + `resonance_levels` alias), `validation/diatomic/{da_figure,nrm,td_nrm_figures,ve_nrm}.py`, and the seven test modules. No consumer imports a private `lcp._name` (verified: `grep -rn "lcp\._\|_walk_from_anion_seed\|_assemble_lcp\|_levels_from" libs/qscat/tests apps validation projects` is empty), so no `__getattr__` shim is needed — the split is invisible.

**Module split (each item lands with its constants, helpers, and docstrings):**

- `curve.py` — capability (a), the fixed-R electronic pole machinery: module docstring = the bulk of today's (the continuation-walk narrative, the freeze semantics, the Naming caution paragraph); `_FROZEN_TOL`, `_JUNCTION_GAMMA_TOL`, `_MIN_RESOLVABLE_GAMMA`; `_h_el`, `resonance_pole_walk`, `resonance_eigenstate`, `resonance_eigenstate_at_peak_width`, `_assemble_lcp`, `_walk_from_anion_seed`, `local_complex_potential`. Imports: `FemDvrEcsGrid, eigen, kinetic` from `qscat.dvr`; `find_resonance_pole` from `qscat.ecs`; `ConvergenceError`; `c_product`; `anion_electronic_states` from `..dissociation`.
- `cross_section.py` — capability (b), the 1-D TI resolvent solver: `_Sigma`/`_Psi`/`_PsiOut` aliases (with their comment block) and `lcp_da_cross_section` (post-Task-1 it imports `Ordering` from `qscat.linalg`). Imports: `kinetic_sparse`, `SparseLU`, scipy.sparse.
- `levels.py` — capability (c), the BO nuclear eigenproblem: `_C_NORM_TOL`, `_CLOSED_REGION_GAMMA_TOL`; `ResonanceLevels` (with Task 13's save/load), `_check_shared_real_nodes`, `_default_window`, `_levels_from`, `lcp_resonance_levels`, `_check_angle_bound`, `resonance_levels` (+ its two overloads). Imports: `kinetic`, `eigen`, `match_angle_stable` from `qscat.ecs`, `c_product`, `assert_shared_real_nodes` from `..grids`, and `from .curve import _assemble_lcp, _walk_from_anion_seed`.
- `__init__.py` — a SHORT module docstring (what the LCP approximation is, the research-program "approximation under test" framing, one pointer per submodule and to `docs/physics/lcp-resonance-levels.md` / `diatomic-ve-cross-sections.md`) followed by:

  ```python
  from .cross_section import lcp_da_cross_section
  from .curve import (
      local_complex_potential,
      resonance_eigenstate,
      resonance_eigenstate_at_peak_width,
      resonance_pole_walk,
  )
  from .levels import ResonanceLevels, lcp_resonance_levels, resonance_levels

  __all__ = [ ... exactly the list above ... ]
  ```

**Steps:**

- [ ] Re-verify: re-read the CURRENT `lcp.py` top to bottom (Tasks 1/11/13 edited it; the api-surface/kernel plans may have too) and re-run the consumer greps above. Adjust the split lists to what is actually in the file — the three-capability boundary (curve / cross-section / levels) is the invariant, the exact helper roster is not.
- [ ] Add to `test_lcp.py` FIRST (it passes before AND after the split — it is the no-regression contract, written to fail if the split loses a name):

  ```python
  def test_lcp_public_surface_survives_the_split() -> None:
      """lib-M11: every documented qscat.core.lcp name keeps its import path."""
      import qscat.core.lcp as lcp

      expected = {
          "ResonanceLevels", "lcp_da_cross_section", "lcp_resonance_levels",
          "local_complex_potential", "resonance_eigenstate",
          "resonance_eigenstate_at_peak_width", "resonance_levels",
          "resonance_pole_walk",
      }
      assert set(lcp.__all__) == expected
      for name in expected:
          assert callable(getattr(lcp, name)) or isinstance(getattr(lcp, name), type)
  ```

- [ ] Run it against the un-split module — expect PASS (baseline). Now perform the split: `git rm` is wrong here (never stage broadly) — do `git mv libs/qscat/qscat/core/lcp.py libs/qscat/qscat/core/lcp_tmp.py` is unnecessary ceremony; simply create the `lcp/` package files with content cut from `lcp.py`, then delete `lcp.py` in the same change. Move code verbatim — this task moves lines, it does not edit bodies (any body edit belongs to an earlier task or a follow-up).
- [ ] Run the surface test + the seven lcp-consuming test modules (`test_lcp.py test_lcp_resonance_levels.py test_assignment.py test_bo.py test_nrm_coupling.py test_nrm_dissociation.py test_nrm_extended.py test_nrm_td_cross_section.py`) — expect PASS. Run `uv run --no-sync pytest apps -m "not slow" -q` if qscat-run has fast tests (runner imports the 5 names).
- [ ] Check docs references: `grep -rn "core/lcp.py\|core\.lcp\b" docs/ --include="*.md" | head` — update any doc line that names the FILE `core/lcp.py` to name the package (`qscat.core.lcp` module references stay correct). CLAUDE.md's `qscat.core` entry describes `lcp` by capability, not by file — re-read the paragraph and leave it unless it names the single-file layout.
- [ ] Full gates (mypy will catch any import the move missed; ruff `I` will re-sort the new files' imports — run `ruff check --fix` on the four new files if needed, then `ruff format`).
- [ ] Commit: `git add libs/qscat/qscat/core/lcp.py libs/qscat/qscat/core/lcp/ libs/qscat/tests/test_lcp.py` (the deleted path stages the deletion; verify with `git status --short` that ONLY these paths are staged) plus any docs file updated; message `refactor(core): split lcp into curve / cross_section / levels package, import paths unchanged`.

---

## Self-review checklist (run before reporting done)

- [ ] Placeholder scan: `grep -n "TBD\|TODO\|similar to Task\|add error handling" docs/superpowers/plans/2026-08-25-library-structure.md` — only this line may match.
- [ ] Findings coverage: lib-M11 (Task 14), lib-M9 (Tasks 12–13), lib-m14 (Task 11), lib-m9 (Task 2), lib-M12 (Task 1), lib-M16 (Tasks 5–7), lib-m10 (Task 3), lib-m18 (Task 4), lib-m7 (Task 7), lib-M13/style-N2/N3 (Tasks 2, 8–10).
- [ ] Name consistency: `Ordering`, `DrResult`, `dr_solve`, `propagate_wavepacket`, `s_vector_one_energy`, `sigma_one_energy`, `default_backend`, `qscat.core.lcp.{curve,cross_section,levels}` — each introduced once and referenced with the same spelling throughout.
- [ ] Every task ends with the four gates and an explicit-path commit.
