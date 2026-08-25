# Experiment Lifecycle (Phase 2.B) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Finish the lifecycle for the N₂ experiment code that stalled at the toy-model stage: graduate the LCP vibrational-excitation solver into `qscat.core.lcp` (and wire it into qscat-run), single-source the N₂ potential on `qscat.model.N2`, delete the superseded projects-side copies (`cross_section.py`'s solver, `observation.py`, `nuclear_density.py`'s inlined solve), lock the qscat-run N₂ presets to the projects convergence studies, and make qscat-run config errors actionable (`ConfigError`, never `KeyError`).

**Architecture:** One new public solver, `qscat.core.lcp.lcp_ve_cross_section`, the VE sibling of the existing `lcp_da_cross_section` — same sparse `kinetic_sparse` + `SparseLU`/`refactor` energy-sweep structure, same `return_wavefunction` convention, model-independent (takes `Vd`/`Gamma` arrays, never a model). Consumers (`validation/n2/cross_section.py`, `validation/diatomic/ve_nrm.py`, the TD test's TI oracle, `nuclear_density.py`) rewire to it; the projects copy is deleted only after a differential parity test pins old-vs-new. qscat-run's `lcp` method gains the `ve` observable (F2/NO via their existing `lcp_grids`; N₂ gains an `lcp_grids` entry). The N₂ potential's three implementations collapse to delegation on `qscat.model.N2`, with identity checks in the `validation/n2/test_resonance.py` style replacing the numeric lockstep tests.

**Tech Stack:** Python 3.12 / uv, pytest (fast + `@slow` tiers), scipy.sparse + `qscat.linalg.SparseLU`, PyYAML/click (qscat-run), numpy.

**Spec:** the "Findings addressed" section below (self-contained; from the 2026-08-25 release review)

**Sequencing:** Runs AFTER the 2026-08-25 layering-close plan (Phase 1.4), in the same tree. This plan deliberately builds on its end state: layering-close first edits `projects/n2_ti_cross_section/test_cross_section.py` (removing its validation imports and the Houfek anchor duplicate) and creates `validation/n2/test_anchor_gate.py`; this plan then deletes `projects/n2_ti_cross_section/{cross_section.py,test_cross_section.py}` outright (Task 4). Do not reorder the two plans.

## Global Constraints

PyPI release DEFERRED until the peer-reviewed article publishes — repo-only distribution, no publishing tasks. After every task: `uv run --no-sync pytest -m "not slow" -n auto --dist loadfile` green; `uv run --no-sync mypy libs/qscat/qscat apps/qscat-run/qscat_run` clean; `uv run --no-sync ruff check .` + `ruff format --check .` clean. Physics-bearing moves are IDENTITY-PRESERVING (differential test pinning old vs new output before deleting the old). Tasks touching validation/n2 or validation/diatomic solver paths need a `validate:n2` / `validate:diatomic` labelled run before merge. Never `git commit -a`. Layering rule: **validation may import projects and qscat; projects may import qscat only; qscat imports neither.**

## Findings addressed

> **exp-M1.** Graduate `projects/n2_ti_cross_section/cross_section.py`'s LCP VE solver into `qscat.core.lcp` as `lcp_ve_cross_section`; consumers: `validation/n2/cross_section.py`, `validation/diatomic/ve_nrm.py:72`; then wire qscat-run `methods:[lcp]` + `ve` observable, currently rejected (apps/qscat-run README lines 55-57, `config.py` `validate_config`'s lcp block). Differential-test old vs new to 1e-12 before deleting the project copy. Give the graduated solver the per-energy reuse its docstring once falsely advertised OR document the dense choice.
>
> **exp-M10.** Delete `projects/n2_2d_cross_section/nuclear_density.py`'s inlined copy of the driven solve (`lcp_driven_solution`, lines 96-121), pointing it at the graduated function.
>
> **exp-M2.** Single-source the N₂ potential: make `validation/n2/model.py` and `projects/n2_resonance/potential.py` thin consumers of `qscat.model.N2`; note the parameter-provenance question (`config.json` vs `library.py` literals); turn the existing lockstep tests into identity checks (the `validation/n2/test_resonance.py` pattern).
>
> **exp-M11.** Lock test: qscat_run presets' `_n2_ti_grid`/`_n2_td_grid` byte-identical to `projects` convergence `WORKING_GRID`/`TD_WORKING_GRID` — mirror `validation/diatomic/test_da_grid.py::test_diatomic_decks_match_presets`.
>
> **exp-M4.** qscat-run config: missing sub-keys raise `ConfigError` not `KeyError` — wrap the raw[...] subscripts in `load_config`/`_load_td`/`_load_energies`/`_load_grid` with a helper `_require(raw, key, where)`; add the missing-subkey test class (three reproduced cases: `{kine: ve}` typo, `td` missing `n_steps`, `energies` missing `step`).
>
> **exp-M9.** Delete `projects/n2_2d_td_cross_section/observation.py` + its test IF nothing else imports them (verified: only `test_observation.py` imports it; `qscat/core/plot.py:12` and `qscat/core/time_dependent.py:187` mention it in docstrings only) — say what replaces each capability in `qscat_run/artifacts.py`; same delete-vs-keep decision for `nuclear_density.py` after M10.
>
> **lib-m15 / exp-m5.** `libs/qscat/tests/test_lcp.py` importorskips `projects` (`test_matches_n2_vres_oracle`, lines 72-89) — rework so the libs suite runs from the sdist.

**Design decisions taken here:**

- **Reuse, not dense.** The graduated `lcp_ve_cross_section` gets the per-energy symbolic-analysis reuse (`SparseLU` + `refactor`), exactly mirroring `lcp_da_cross_section` (`libs/qscat/qscat/core/lcp.py:475-546`) — the dense `np.linalg.solve` in the projects copy was a documented toy-model choice, and the graduated solver drives real energy sweeps (`validation/diatomic/ve_nrm.py`, qscat-run). Consequence: old-vs-new is dense-LAPACK-vs-SuperLU, so exact bit-identity is not on the table; the parity gate targets 1e-12 relative but measures first and, if the cross-solver floor is above that, gates at 10× the measured floor with the number recorded (the `validation/n2/test_anchors.py` tolerance-derivation pattern; see also the ci-test-portability rule: never pin sparse-solve outputs at rtol=1e-12 across architectures).
- **N₂ gets an `lcp_grids` preset entry** (its TI deck's own factors + a 44° electronic partner, the eMoScat 35/44 LCP pairing already used for F2/NO), because N₂ is the flagship VE molecule and `methods:[lcp]`+`ve` must work for it. The partner angle only gates two-angle pole *stability*; the pole values come from grid a's spectrum, so 44° vs the 40° partner `validation/diatomic/ve_nrm.py` measured with changes nothing about accepted pole values (both partners accept the same poles at `resid_tol=1e-3`).
- **`nuclear_density.py` is KEPT** (slimmed by M10): grep shows only its own test imports it, but its exact-2D-vs-LCP nuclear-density *shape* comparison (`compare_to_lcp`: unit-area normalization, centroid, RMS width) is a physics capability `qscat_run/artifacts.py` does not reproduce (artifacts emit TI/TD density marginals and LCP scattering states separately, never the normalized shape comparison), and `docs/physics/n2-2d-cross-section.md` (line 316) documents its findings. `observation.py` IS deleted — every capability has an artifacts.py replacement (mapping in Task 9).
- **The vres oracle test moves to validation.** `test_matches_n2_vres_oracle` compares `qscat.core.lcp.local_complex_potential` against the *independent* projects pole walk (`vres_on_grid`) — that independence is the point, so the test is not rewritten against qscat itself (which would be tautological); it moves to `validation/n2/` where importing projects is legal, and the libs suite becomes projects-free (sdist-runnable).

---

## Task 1 — `qscat.core.lcp.lcp_ve_cross_section` (sparse, refactor-swept, wavefunction-exposing)

**Files:**
- Modify: `libs/qscat/qscat/core/lcp.py` (new function + `__all__` entry), `libs/qscat/tests/test_lcp.py` (new tests), `libs/qscat/tests/test_nrm_td_cross_section.py` (docstring at lines 685-700 claims "There is no `lcp_ve_cross_section` in `qscat`" — now false)
- Test: `libs/qscat/tests/test_lcp.py`

**Interfaces:**
- Consumes: `qscat.dvr.kinetic_sparse(grid: FemDvrEcsGrid, mu: float) -> sp.spmatrix`, `qscat.linalg.SparseLU(A, ordering=...)` / `.refactor(A_new)` / `.solve(b)`, the module's existing `_Sigma` / `_Psi` / `_PsiOut` type aliases and `_Ordering`.
- Produces:

  ```python
  def lcp_ve_cross_section(
      nuclear_grid: FemDvrEcsGrid,
      mu: float,
      Vd: npt.NDArray[np.complex128],
      Gamma: npt.NDArray[np.float64],
      eps: npt.NDArray[np.float64],
      chi: npt.NDArray[np.complex128],
      v_init: int,
      vprimes: list[int],
      E: float | npt.ArrayLike,
      *,
      ordering: _Ordering = "COLAMD",
      return_wavefunction: bool = False,
  ) -> _Sigma | tuple[_Sigma, _PsiOut]
  ```

  Scalar `E` → shape `(len(vprimes),)`; array `E` → `(len(E), len(vprimes))` (the projects function's exact convention, so consumers rewire without call-site changes). `return_wavefunction` → also the 1-D driven solution `xi(R)` per energy (`None` for `E <= 0`), same convention as `lcp_da_cross_section`.

**Steps:**

- [ ] Write the failing tests in `libs/qscat/tests/test_lcp.py`, reusing the existing module-scoped `coarse_nuc`/`coarse_lcp_inputs` F2 fixtures (the shared ~3.5 s pole walk):

  ```python
  def test_lcp_ve_matches_dense_reference(coarse_nuc, coarse_lcp_inputs):
      # Differential oracle: the same driven equation assembled densely and
      # solved with np.linalg.solve -- the projects toy model's formulation.
      from qscat.core.lcp import lcp_ve_cross_section
      from qscat.dvr import kinetic

      g_R = coarse_nuc
      Vd, Gamma, eps, chi = coarse_lcp_inputs
      E, vprimes = 0.05, [0, 1, 2]
      sigma = lcp_ve_cross_section(g_R, F2.mu, Vd, Gamma, eps, chi, 0, vprimes, E)

      doorway = np.sqrt(Gamma / (2.0 * np.pi))[None, :] * chi
      H = kinetic(g_R, F2.mu) + np.diag(Vd - 0.5j * Gamma)
      M = (E + eps[0]) * np.eye(g_R.n, dtype=np.complex128) - H
      xi = np.linalg.solve(M, doorway[0])
      e_tot = E + eps[0]
      expected = np.zeros(len(vprimes))
      for k, vp in enumerate(vprimes):
          if e_tot - eps[vp] > 0.0:
              S = np.dot(doorway[vp], xi)  # c-product: no conjugate
              expected[k] = 4.0 * np.pi**3 * np.abs(S) ** 2 / (2.0 * E)
      np.testing.assert_allclose(sigma, expected, rtol=1e-9)


  def test_lcp_ve_shapes_closed_channels_and_sweep_reuse(coarse_nuc, coarse_lcp_inputs):
      from qscat.core.lcp import lcp_ve_cross_section

      g_R = coarse_nuc
      Vd, Gamma, eps, chi = coarse_lcp_inputs
      # scalar E -> (len(vprimes),); array E -> (len(E), len(vprimes))
      s1 = lcp_ve_cross_section(g_R, F2.mu, Vd, Gamma, eps, chi, 0, [0, 1], 0.05)
      assert s1.shape == (2,) and np.all(np.isfinite(s1)) and np.all(s1 >= 0.0)
      E = np.array([0.02, 0.05])
      s2 = lcp_ve_cross_section(g_R, F2.mu, Vd, Gamma, eps, chi, 0, [0, 1], E)
      assert s2.shape == (2, 2)
      # the refactor sweep must equal fresh per-energy solves
      for i, e in enumerate(E):
          fresh = lcp_ve_cross_section(g_R, F2.mu, Vd, Gamma, eps, chi, 0, [0, 1], float(e))
          np.testing.assert_allclose(s2[i], fresh, rtol=1e-12)
      # E <= 0 is closed -> exactly zero; a closed v' channel -> exactly zero
      assert np.all(
          lcp_ve_cross_section(g_R, F2.mu, Vd, Gamma, eps, chi, 0, [0, 1], -0.01) == 0.0
      )
      s3 = lcp_ve_cross_section(g_R, F2.mu, Vd, Gamma, eps, chi, 0, [2], 1e-4)
      assert s3[0] == 0.0  # eps[2]-eps[0] >> 1e-4: channel closed


  def test_lcp_ve_return_wavefunction_parity(coarse_nuc, coarse_lcp_inputs):
      from qscat.core.lcp import lcp_ve_cross_section

      g_R = coarse_nuc
      Vd, Gamma, eps, chi = coarse_lcp_inputs
      E = np.array([0.02, 0.05])
      s_plain = lcp_ve_cross_section(g_R, F2.mu, Vd, Gamma, eps, chi, 0, [0, 1], E)
      s2, xis = lcp_ve_cross_section(
          g_R, F2.mu, Vd, Gamma, eps, chi, 0, [0, 1], E, return_wavefunction=True
      )
      assert np.array_equal(s_plain, s2)
      assert isinstance(xis, list) and len(xis) == 2
      for xi in xis:
          assert xi is not None and xi.shape == (g_R.n,) and xi.dtype == np.complex128
      _s, xi0 = lcp_ve_cross_section(
          g_R, F2.mu, Vd, Gamma, eps, chi, 0, [0], -0.01, return_wavefunction=True
      )
      assert xi0 is None
  ```

- [ ] Run: `uv run --no-sync pytest libs/qscat/tests/test_lcp.py -q -m "not slow"` — expect the three new tests to FAIL with `ImportError: cannot import name 'lcp_ve_cross_section' from 'qscat.core.lcp'`.
- [ ] Implement in `libs/qscat/qscat/core/lcp.py` (place after `lcp_da_cross_section`; add `"lcp_ve_cross_section"` to `__all__`; add the two `@overload`s on `return_wavefunction`, mirroring `lcp_da_cross_section`'s at lines 443-472):

  ```python
  def lcp_ve_cross_section(
      nuclear_grid: FemDvrEcsGrid,
      mu: float,
      Vd: npt.NDArray[np.complex128],
      Gamma: npt.NDArray[np.float64],
      eps: npt.NDArray[np.float64],
      chi: npt.NDArray[np.complex128],
      v_init: int,
      vprimes: list[int],
      E: float | npt.ArrayLike,
      *,
      ordering: _Ordering = "COLAMD",
      return_wavefunction: bool = False,
  ) -> _Sigma | tuple[_Sigma, _PsiOut]:
      """LCP vibrational-excitation sigma_{v_init->v'}(E) (bohr^2), TI resolvent form.

      Solve `(E_tot I - H_res) xi = d_{v_init}`, `H_res = T_nuc(mu) + diag(V_d
      - i Gamma/2)`, doorway `d_v = sqrt(Gamma/2pi) chi_v`; S-matrix element
      `S_{v'<-v_init} = <d_{v'}|xi>` by the DVR c-product (no conjugate);
      `sigma = 4 pi^3 |S|^2 / 2E`, exactly zero for `E <= 0` and for a closed
      final channel (`E_tot - eps[v'] <= 0`).

      Graduated from `projects/n2_ti_cross_section/cross_section.py`'s
      `ve_cross_section` (the deliberately dense 1-D toy model). This version
      is SPARSE and sweep-reusing: `A(E) = E_tot I - H_res` has an
      E-independent sparsity pattern, so the symbolic analysis is done once
      and `SparseLU.refactor` re-runs only the numeric factor per energy --
      the same structure as `lcp_da_cross_section` and `driven.ve_cross_section`.
      `xi` depends only on `(E, v_init)`, so one solve per energy serves every
      channel in `vprimes`.

      If `return_wavefunction`, also returns `xi(R)` per energy (`None` when
      `E <= 0`): one array for scalar `E`, one list entry per energy for array
      `E` -- the driven solution `nuclear_density.lcp_driven_solution` consumes.
      """
      doorway = np.sqrt(Gamma / (2.0 * np.pi)).astype(np.complex128)[None, :] * chi
      H_res = (kinetic_sparse(nuclear_grid, mu) + sp.diags(Vd - 0.5j * Gamma)).tocsc()
      ident = sp.identity(nuclear_grid.n, format="csc", dtype=np.complex128)

      e_arr = np.atleast_1d(np.asarray(E, dtype=np.float64))
      out = np.zeros((e_arr.size, len(vprimes)), dtype=np.float64)
      psi_list: list[_Psi] = [None] * e_arr.size
      lu: SparseLU | None = None
      for ie, e in enumerate(e_arr):
          if float(e) <= 0.0:
              continue
          e_tot = float(e) + eps[v_init]
          a = (e_tot * ident - H_res).tocsc()
          if lu is None:
              lu = SparseLU(a, ordering=ordering)
          else:
              lu.refactor(a)
          xi = lu.solve(doorway[v_init])
          psi_list[ie] = np.asarray(xi, dtype=np.complex128)
          for k, vp in enumerate(vprimes):
              if e_tot - eps[vp] <= 0.0:
                  continue  # closed channel
              s_el = np.dot(doorway[vp], xi)  # c-product: no conjugate
              out[ie, k] = 4.0 * np.pi**3 * np.abs(s_el) ** 2 / (2.0 * float(e))

      scalar = np.isscalar(E) or (isinstance(E, np.ndarray) and np.ndim(E) == 0)
      sigma = np.asarray(out[0] if scalar else out, dtype=np.float64)
      if return_wavefunction:
          return sigma, (psi_list[0] if scalar else psi_list)
      return sigma
  ```

- [ ] Update the docstring of `libs/qscat/tests/test_nrm_td_cross_section.py::test_markovian_ve_reproduces_the_local_cross_section` (lines ~685-700): replace "There is no `lcp_ve_cross_section` in `qscat` -- the repository's LCP VE route lives in `projects/...`, which `libs/qscat` must not import" with: "`qscat.core.lcp.lcp_ve_cross_section` exists, but the reference here stays hand-assembled from `solve_nuclear` + the doorway on purpose: an INDEPENDENT assembly of the same formula is an oracle; comparing the shipped function to itself would not be." No code change to that test.
- [ ] Run: `uv run --no-sync pytest libs/qscat/tests/test_lcp.py libs/qscat/tests/test_nrm_td_cross_section.py -q -m "not slow"` — expect PASS.
- [ ] Run the full gate (fast tier + mypy + ruff, per Global Constraints).
- [ ] Commit: `git add libs/qscat/qscat/core/lcp.py libs/qscat/tests/test_lcp.py libs/qscat/tests/test_nrm_td_cross_section.py && git commit -m "feat(qscat.core.lcp): graduate the LCP VE cross section as lcp_ve_cross_section (sparse, refactor-swept)"`

---

## Task 2 — Differential parity gate: projects copy vs graduated solver (pin before delete)

The identity-preservation gate the Global Constraints require. It lives in `validation/n2/` (validation may import projects) and is deleted in Task 4 together with the old copy, its purpose fulfilled.

**Files:**
- Create: `validation/n2/test_lcp_ve_parity.py`
- Test: itself

**Interfaces:**
- Consumes: `projects.n2_ti_cross_section.cross_section.ve_cross_section` (old), `qscat.core.lcp.lcp_ve_cross_section` (new), `validation.n2.cross_section.build_system() -> (grid, eps, chi, Vd, Gamma)` (lru_cached, ~7 s once per process), `qscat.model.N2.mu`.
- Produces: the pinned old-vs-new differential record.

**Steps:**

- [ ] Create `validation/n2/test_lcp_ve_parity.py`:

  ```python
  """TEMPORARY parity gate for the exp-M1 graduation (deleted with the old copy).

  Pins `projects/n2_ti_cross_section/cross_section.py::ve_cross_section` (dense
  np.linalg.solve) against `qscat.core.lcp.lcp_ve_cross_section` (sparse
  SparseLU + refactor) on the real N2 system, before the projects copy is
  deleted. Dense LAPACK vs SuperLU on the same ~300x300 complex-symmetric
  matrix: the target is 1e-12 relative; if the measured cross-solver floor is
  above that, gate at 10x the measured maximum and record it below (the
  test_anchors.py tolerance-derivation pattern; cf. the ci-test-portability
  rule against pinning sparse-solve outputs at 1e-12 cross-arch).

  MEASURED (fill in at implementation time): max |sigma_new/sigma_old - 1| =
  <value> over the grid below.
  """

  from __future__ import annotations

  import numpy as np
  from qscat.core.lcp import lcp_ve_cross_section
  from qscat.model import N2

  from projects.n2_ti_cross_section.cross_section import ve_cross_section
  from validation.n2.cross_section import build_system

  RTOL = 1e-12  # raise to 10x the measured floor if dense-vs-sparse exceeds it


  def test_graduated_solver_reproduces_the_projects_copy():
      grid, eps, chi, Vd, Gamma = build_system()
      E = np.array([0.02, 0.05, 0.1, 0.15, 0.2])
      vprimes = [0, 1, 2, 3]
      old = ve_cross_section(grid, N2.mu, Vd, Gamma, eps, chi, 0, vprimes, E)
      new = lcp_ve_cross_section(grid, N2.mu, Vd, Gamma, eps, chi, 0, vprimes, E)
      assert old.shape == new.shape == (5, 4)
      dev = np.abs(new - old) / np.maximum(np.abs(old), 1e-300)
      print(f"max relative deviation old-vs-new: {dev.max():.3e}")
      np.testing.assert_allclose(new, old, rtol=RTOL)
  ```

- [ ] Run: `uv run --no-sync pytest validation/n2/test_lcp_ve_parity.py -q -s` — read the printed max deviation. If it exceeds 1e-12, set `RTOL = 10 * <measured>` and record the measured value in the docstring; rerun. Expect PASS with the honest tolerance recorded.
- [ ] Commit: `git add validation/n2/test_lcp_ve_parity.py && git commit -m "test(validation/n2): pin the projects LCP VE solver against the graduated qscat.core.lcp copy"`

---

## Task 3 — Rewire every consumer to the graduated solver (incl. exp-M10)

**Files:**
- Modify: `validation/n2/cross_section.py` (import at line 60 + docstring), `validation/diatomic/ve_nrm.py` (import at lines 72-74 + docstring's route table at lines 29-31), `projects/n2_td_cross_section/test_td_cross_section.py` (import at line 41 + docstring mentions), `projects/n2_2d_cross_section/nuclear_density.py` (M10: `lcp_driven_solution` delegates; module docstring lines 39-47 rewritten)
- Test: existing suites of all four modules (`validation/n2/test_anchor_gate.py`, `validation/diatomic/test_ve_nrm.py` (slow), `projects/n2_td_cross_section/test_td_cross_section.py`, `projects/n2_2d_cross_section/test_nuclear_density.py`), plus the Task 2 parity test still green.

**Interfaces:**
- Consumes: `qscat.core.lcp.lcp_ve_cross_section` (Task 1 signature).
- Produces: `nuclear_density.lcp_driven_solution(grid: FemDvrEcsGrid, mu: float, eps, chi, v_init: int, E: float) -> npt.NDArray[np.complex128]` — signature unchanged, body delegated.

**Steps:**

- [ ] `validation/n2/cross_section.py`: replace `from projects.n2_ti_cross_section.cross_section import ve_cross_section` with `from qscat.core.lcp import lcp_ve_cross_section as ve_cross_section` (call sites at lines 148-150 are signature-compatible; keep the local alias so the rest of the module is untouched). Update the docstring's "the object under test for C5 *is* that project's resolvent/driven-equation solver" sentence: the solver now ships as `qscat.core.lcp.lcp_ve_cross_section`; the `V_d`/`Gamma` inputs still come from the projects pole walk (`vres_on_grid`).
- [ ] `validation/diatomic/ve_nrm.py`: replace lines 72-74 with `from qscat.core.lcp import lcp_ve_cross_section` (the call at line 331 already uses that name; delete the alias import). Update the docstring route table (line 30): `lcp -- qscat.core.lcp.local_complex_potential + qscat.core.lcp.lcp_ve_cross_section`.
- [ ] `projects/n2_td_cross_section/test_td_cross_section.py`: replace line 41 with `from qscat.core.lcp import lcp_ve_cross_section as ve_cross_section`; update the module docstring's oracle reference (lines 8-11) to name `qscat.core.lcp.lcp_ve_cross_section`.
- [ ] `projects/n2_2d_cross_section/nuclear_density.py` (M10): delete the body of `lcp_driven_solution` (the inlined `kinetic` + `np.diag` + `np.linalg.solve` block, lines 112-121) and delegate:

  ```python
  def lcp_driven_solution(
      grid: FemDvrEcsGrid,
      mu: float,
      eps: npt.NDArray[np.float64],
      chi: npt.NDArray[np.complex128],
      v_init: int,
      E: float,
  ) -> npt.NDArray[np.complex128]:
      """The 1-D LCP driven solution `xi(R)` at collision energy `E` -- now the
      graduated solver's own `return_wavefunction` output (`qscat.core.lcp.
      lcp_ve_cross_section`), no longer an inlined re-solve. `vres_on_grid`
      still supplies this project's `(V_d, Gamma)`."""
      Vd, Gamma = vres_on_grid(grid)
      _sigma, xi = lcp_ve_cross_section(
          grid, mu, Vd, Gamma, eps, chi, v_init, [v_init], float(E),
          return_wavefunction=True,
      )
      if xi is None:
          raise ValueError(f"lcp_driven_solution: E={E} Ha <= 0, no driven solve")
      return np.asarray(xi, dtype=np.complex128)
  ```

  Remove the now-unused `kinetic` import if nothing else in the module uses it (check: `nuclear_density` doesn't); add `from qscat.core.lcp import lcp_ve_cross_section`. Rewrite the module docstring's "No helper anywhere already exposes `xi`, so `lcp_driven_solution` below inlines..." paragraph (lines 39-47): the graduated function exposes `xi` directly, so this module no longer inlines anything.
- [ ] Run: `uv run --no-sync pytest projects/n2_td_cross_section projects/n2_2d_cross_section/test_nuclear_density.py validation/n2/test_lcp_ve_parity.py validation/n2/test_anchor_gate.py -q -m "not slow"` — expect PASS.
- [ ] Run the full gate. This task touches `validation/n2` and `validation/diatomic` solver paths — the branch needs `validate:n2` AND `validate:diatomic` labelled runs before merge.
- [ ] Commit: `git add validation/n2/cross_section.py validation/diatomic/ve_nrm.py projects/n2_td_cross_section/test_td_cross_section.py projects/n2_2d_cross_section/nuclear_density.py && git commit -m "refactor: rewire every LCP VE consumer to qscat.core.lcp.lcp_ve_cross_section (exp-M1, exp-M10)"`

---

## Task 4 — Delete the projects copy; re-home its internal tests

`projects/n2_ti_cross_section/cross_section.py` now has zero importers outside its own test (Task 3 rewired all four; verify with grep). The three surviving internal tests in `test_cross_section.py` (post-layering-close: `test_sigma_real_and_nonnegative`, `test_closed_channel_is_exactly_zero`, `test_v0_to_v1_resonance_enhancement`) are genuine N₂ physics coverage — they move to validation, driven by the graduated solver on the cached `build_system()`.

**Files:**
- Create: `validation/n2/test_lcp_ve.py`
- Delete: `projects/n2_ti_cross_section/cross_section.py`, `projects/n2_ti_cross_section/test_cross_section.py`, `validation/n2/test_lcp_ve_parity.py` (purpose fulfilled — it pinned old-vs-new; the old side no longer exists)
- Modify: `projects/n2_td_cross_section/td_cross_section.py` (docstring references to the deleted module at lines 8, 28, 158), `CLAUDE.md` (the `n2_ti_cross_section` repo-map line names `cross_section.py`; the `qscat.core` blurb gains `lcp_ve_cross_section`), `docs/physics/n2-cross-section.md` + `docs/physics/diatomic-ve-cross-sections.md` (grep for `n2_ti_cross_section/cross_section.py` / `ve_cross_section` project references and update to the graduated name)
- Test: `validation/n2/test_lcp_ve.py`, full fast tier

**Interfaces:**
- Consumes: `validation.n2.cross_section.build_system()`, `qscat.core.lcp.lcp_ve_cross_section`, `qscat.model.N2.mu`.
- Produces: no new API — deletions plus the re-homed tests.

**Steps:**

- [ ] Create `validation/n2/test_lcp_ve.py` with the three internal tests ported verbatim in physics, driven by `build_system()` + `lcp_ve_cross_section` (module docstring: "the model-independent internal checks that shipped with the projects toy model, re-homed onto the graduated solver"):

  ```python
  from __future__ import annotations

  import numpy as np
  from qscat.core.lcp import lcp_ve_cross_section
  from qscat.model import N2

  from validation.n2.cross_section import build_system


  def test_sigma_real_and_nonnegative():
      grid, eps, chi, Vd, Gamma = build_system()
      for E in (0.02, 0.05, 0.1, 0.15, 0.2):
          sigma = lcp_ve_cross_section(grid, N2.mu, Vd, Gamma, eps, chi, 0, [0, 1, 2, 3], E)
          assert sigma.shape == (4,)
          assert np.all(np.isfinite(sigma)) and np.all(sigma >= 0.0)


  def test_closed_channel_is_exactly_zero():
      grid, eps, chi, Vd, Gamma = build_system()
      E = 0.001
      assert E + eps[0] - eps[3] < 0  # sanity: v'=3 is closed at this E
      assert lcp_ve_cross_section(grid, N2.mu, Vd, Gamma, eps, chi, 0, [3], E)[0] == 0.0


  def test_v0_to_v1_resonance_enhancement():
      grid, eps, chi, Vd, Gamma = build_system()
      near = lcp_ve_cross_section(grid, N2.mu, Vd, Gamma, eps, chi, 0, [1], 0.02)[0]
      res = lcp_ve_cross_section(grid, N2.mu, Vd, Gamma, eps, chi, 0, [1], 0.1)[0]
      assert res > near  # the ~2-3 eV Pi_g resonance enhances 0->1
  ```

- [ ] Run: `uv run --no-sync pytest validation/n2/test_lcp_ve.py -q` — expect 3 PASS.
- [ ] Verify no importer remains: `grep -rn "n2_ti_cross_section.cross_section\|n2_ti_cross_section import cross_section" --include="*.py" projects validation libs apps benchmarks` — expect hits only in the two files being deleted plus `projects/n2_td_cross_section/td_cross_section.py`'s docstring prose (updated below); no `import` statements elsewhere.
- [ ] Delete: `git rm projects/n2_ti_cross_section/cross_section.py projects/n2_ti_cross_section/test_cross_section.py validation/n2/test_lcp_ve_parity.py`
- [ ] Update the docstring references: `projects/n2_td_cross_section/td_cross_section.py` lines 8/28/158 → name `qscat.core.lcp.lcp_ve_cross_section`; CLAUDE.md's repo-map line for `n2_ti_cross_section` (drop `cross_section.py` from the file list, note the solver graduated to `qscat.core.lcp`) and its `qscat.core` lcp blurb (add `lcp_ve_cross_section`); `docs/physics/n2-cross-section.md` and `docs/physics/diatomic-ve-cross-sections.md` wherever `grep -rn "n2_ti_cross_section/cross_section" docs/physics` hits.
- [ ] Run the full gate (fast tier catches any missed importer).
- [ ] Commit: `git add -u projects/n2_ti_cross_section validation/n2/test_lcp_ve_parity.py && git add validation/n2/test_lcp_ve.py projects/n2_td_cross_section/td_cross_section.py CLAUDE.md docs/physics/n2-cross-section.md docs/physics/diatomic-ve-cross-sections.md && git commit -m "refactor(lifecycle): delete the graduated LCP VE toy model; re-home its internal tests onto validation"`

---

## Task 5 — qscat-run: `methods: [lcp]` + `ve` observable

**Files:**
- Modify: `apps/qscat-run/qscat_run/presets.py` (N₂ `lcp_grids`), `apps/qscat-run/qscat_run/config.py` (`validate_config` lcp block), `apps/qscat-run/qscat_run/runner.py` (`_run_lcp` ve branch + n_vib), `apps/qscat-run/README.md` (lcp method bullet, lines 51-57, History's `ve_nrm_figure.py` bullet), `apps/qscat-run/tests/test_config.py` (three tests change meaning), `apps/qscat-run/tests/test_runner_lcp.py` (new `@slow` end-to-end)
- Test: `apps/qscat-run/tests/test_config.py` (fast gate), `apps/qscat-run/tests/test_runner_lcp.py` (slow)

**Interfaces:**
- Consumes: `qscat.core.lcp.lcp_ve_cross_section`; `presets.MoleculePreset.lcp_grids: Callable[[], tuple[FemDvrEcsGrid, FemDvrEcsGrid, FemDvrEcsGrid]] | None`.
- Produces: `_n2_lcp_grids() -> tuple[FemDvrEcsGrid, FemDvrEcsGrid, FemDvrEcsGrid]` in presets; cross-section result keys `"lcp:ve:v{v_init}->{vp}"`.

**Steps:**

- [ ] Failing tests first, in `apps/qscat-run/tests/test_config.py`:
  - Replace `test_lcp_without_da_observable_rejected` (lines 71-84) with `test_lcp_with_ve_observable_accepted_for_f2`: same F2 config (`methods: [lcp]`, `observables: [{kind: ve, channels: 2}]`), now expecting `validate_config(cfg)` to NOT raise.
  - Replace `test_lcp_without_da_is_still_rejected_when_no_levels_requested` (lines ~458-468) with `test_lcp_with_ve_only_and_energies_accepted` (same shape, no raise).
  - Replace `test_lcp_rejected_for_n2` (lines 55-68) with two tests: `test_lcp_ve_accepted_for_n2` (`molecule: N2`, `methods: [lcp]`, `observables: [{kind: ve, channels: 2}]` → no raise) and `test_lcp_rejected_for_h2p` (`molecule: H2P`, `methods: [lcp]`, `observables: [{kind: dr, channels: 1}]` → `pytest.raises(ConfigError, match=r"lcp.*not available")`).
- [ ] Run: `uv run --no-sync pytest apps/qscat-run/tests/test_config.py -q` — expect the new/changed tests to FAIL: F2+ve raises `ConfigError: methods includes 'lcp' but no 'da' or 'resonance_levels' observable is requested...`; N2 raises `ConfigError: the 'lcp' method is not available for N2...`.
- [ ] Implement, `apps/qscat-run/qscat_run/presets.py`:

  ```python
  def _n2_lcp_grids() -> tuple[FemDvrEcsGrid, FemDvrEcsGrid, FemDvrEcsGrid]:
      """(nuclear, elec_a, elec_b) for N2 LCP -- the TI deck's own factors plus
      the electronic deck rebuilt at the eMoScat LCP partner angle. 35/44 is
      the pairing F2/NO use here AND the pairing the projects N2 pole walk
      itself uses (`projects/n2_ti_cross_section/vres.py`: _ANGLE_A_DEG=35,
      _ANGLE_B_DEG=44). N2's LCP observable is VE (its DA is
      closed-in-range); the partner angle only gates two-angle pole
      STABILITY -- accepted pole values come from grid a's spectrum."""
      return (
          nuclear_grid(angle_deg=35.0, r_max=20.0, n_complex=5, quadrature=10),
          electronic_grid(**_N2_ELEC),  # type: ignore[arg-type]
          electronic_grid(**{**_N2_ELEC, "angle_deg": _LCP_ANGLE_B}),  # type: ignore[arg-type]
      )
  ```

  and add `lcp_grids=_n2_lcp_grids,` to the `"N2:emoscat"` preset entry. Update `resolve_lcp_grids`'s error message ("LCP is defined only for the DA molecules (F2, NO)" → "LCP is wired for the neutral diatomics (N2, F2, NO); H2P's observable is DR").
- [ ] Implement, `apps/qscat-run/qscat_run/config.py` `validate_config` lcp block (lines 578-603): DELETE the `if "da" not in kinds and "resonance_levels" not in kinds:` check with a replacing comment — after this task every lcp-capable molecule's validity matrix (`N2: {ve, da}`, `NO`/`F2`: `{ve, da, resonance_levels}`) guarantees any observable that survived the per-observable validity check is one LCP can serve, so the branch is unreachable; keep the explicit-grid rejection and the preset-availability check, updating the availability message to match `resolve_lcp_grids`'s. Update the docstring's check list accordingly.
- [ ] Implement, `apps/qscat-run/qscat_run/runner.py` `_run_lcp`:
  - imports: add `lcp_ve_cross_section` to the `qscat.core.lcp` import block (line 86-91).
  - `wants_sigma = "da" in kinds` (line 789) → `wants_sigma = bool(kinds & {"da", "ve"})`.
  - n_vib (line 799): mirror `_run_ti`'s widening —

    ```python
    required = cfg.v_init + 1
    for obs in cfg.observables:
        if obs.kind == "ve":
            required = max(required, max(_vprimes(obs), default=-1) + 1)
    n_vib = _n_vib(cfg, required)
    ```

  - in the observables loop (line 863-869), change the guard to `if obs.kind not in ("da", "ve"): continue` and add the ve branch:

    ```python
    if obs.kind == "ve":
        vprimes = _vprimes(obs)
        sigma_ve = lcp_ve_cross_section(
            g_R, model.mu, v_d, gamma, eps, chi, cfg.v_init, vprimes, energies
        )
        for j, vp in enumerate(vprimes):
            cross_sections[f"lcp:ve:v{cfg.v_init}->{vp}"] = np.asarray(
                sigma_ve, dtype=np.float64
            )[:, j]
        timings["lcp:ve"] = timings.get("lcp:ve", 0.0) + (time.time() - t0)
        continue
    ```

    (keep the existing da branch as-is). Update `_run_lcp`'s docstring ("LCP only produces the DA cross section" is no longer true).
- [ ] Add the end-to-end test in `apps/qscat-run/tests/test_runner_lcp.py`:

  ```python
  @pytest.mark.slow
  def test_lcp_run_produces_ve_cross_section_for_n2(tmp_path: Path) -> None:
      out_dir = tmp_path / "out"
      cfg = load_config(
          _write(
              tmp_path,
              f"""
          molecule: N2
          methods: [lcp]
          observables: [{{kind: ve, channels: [0, 1]}}]
          energies: {{values: [0.05, 0.1]}}
          output_dir: {out_dir}
      """,
          )
      )
      validate_config(cfg)
      result = run_experiment(cfg)
      for key in ("lcp:ve:v0->0", "lcp:ve:v0->1"):
          assert key in result.cross_sections
          sigma = result.cross_sections[key]
          assert sigma.shape == (2,)
          assert np.all(np.isfinite(sigma)) and np.all(sigma >= 0.0)
      # the Pi_g resonance region beats near-threshold on 0->1
      assert result.cross_sections["lcp:ve:v0->1"][1] > result.cross_sections["lcp:ve:v0->1"][0]
  ```

- [ ] Run: `uv run --no-sync pytest apps/qscat-run/tests/test_config.py -q` — expect PASS; run `uv run --no-sync pytest apps/qscat-run/tests/test_runner_lcp.py -q -m slow` once locally — expect PASS (the N₂ pole walk is ~7 s + two sparse solves).
- [ ] `apps/qscat-run/README.md`: (1) lcp method bullet (line 33): "the local-complex-potential *approximation* of DA **and VE** (N2/F2/NO)"; (2) delete the sentence "The LCP's own VE route is *not* exposed here (it lives in `projects/n2_ti_cross_section`), so `methods: [lcp]` with a `ve` observable is rejected." and update "`lcp` needs a `da` or `resonance_levels` observable" → "`lcp` needs a `ve`, `da` or `resonance_levels` observable"; (3) History's `ve_nrm_figure.py` bullet: its justification is now ONLY that the figure needs both `include_background` settings (two runs from one flag) — the "LCP's VE route is not a qscat-run method" half is obsolete; (4) Registry table N₂ row: valid observables note unchanged, but the Methods section now lists N₂ under lcp.
- [ ] Run the full gate.
- [ ] Commit: `git add apps/qscat-run/qscat_run/presets.py apps/qscat-run/qscat_run/config.py apps/qscat-run/qscat_run/runner.py apps/qscat-run/README.md apps/qscat-run/tests/test_config.py apps/qscat-run/tests/test_runner_lcp.py && git commit -m "feat(qscat-run): the lcp method serves the ve observable (N2/F2/NO) via lcp_ve_cross_section"`

---

## Task 6 — Single-source the N₂ potential on `qscat.model.N2` (exp-M2)

After layering-close, `projects/n2_resonance/potential.py` already takes its PARAMETERS from `N2` but still re-implements the formulas; `validation/n2/model.py` still parses `config.json` and re-implements them too. This task delegates the function BODIES and converts the numeric lockstep tests into identity checks (the `validation/n2/test_resonance.py::test_resonance_reexports_the_toy_model_grid_factory` pattern, adapted for bound methods).

**Parameter provenance, resolved:** `validation/n2/config.json` stays as the validation-side provenance record (its `provenance` field names the eMoScat deck; its `note` field carries the model-vs-reality caveat) and keeps feeding `model.PARAMS`/`model_checks()`; `qscat.model.library.N2`'s literals are the single RUNTIME source (cited there to PRA 73 Table I). A new exact-equality test (`test_params_match_the_library_fields_exactly`) makes the two un-driftable, which is what "single source" has to mean while both files exist.

**Numeric caution:** `N2`'s methods compute in `complex128` (ECS-safe). Two consequences the edits below must handle: (1) `float()` on a 0-d complex array raises `TypeError` — every `float(v0(...))`-style coercion in `model_checks()` becomes `float(np.real(...))`; (2) complex division can differ from real division in the last ulp, so `test_potential.py`'s two EXACT `==` assertions on `lam` become `abs(...) < 1e-12`.

**Files:**
- Modify: `validation/n2/model.py`, `validation/n2/test_model.py`, `projects/n2_resonance/potential.py`, `projects/n2_resonance/test_potential.py`
- Test: `validation/n2/test_model.py`, `projects/n2_resonance/test_potential.py`, plus consumers `validation/n2/resonance.py` (via `validation/n2` suite) and `projects/n2_ti_cross_section/test_vres.py`

**Interfaces:**
- Consumes: `qscat.model.N2` (bound methods `v0`, `lam`, `v_int`; fields as in the layering-close plan).
- Produces: `validation.n2.model.v0/lam/v_int` and `projects.n2_resonance.potential.v0/lam/v_int` are the SAME bound-method objects as `N2.v0/lam/v_int` (assignable, `== N2.v0` holds); each module keeps a local `v_eff_el(r, R)` (`N2.v_int(r, R) + N2.ell*(N2.ell+1)/(2*r**2)`, complex-safe — `N2` has no v_eff_el because `surface` includes `v0`).

**Steps:**

- [ ] Write the failing identity tests. In `validation/n2/test_model.py` add:

  ```python
  from qscat.model import N2


  def test_model_is_the_library_model():
      """Identity, not lockstep: re-introducing a local copy makes these stop
      BEING the library methods and fails immediately (the test_resonance.py
      pattern; bound-method == compares function and instance)."""
      assert model.v0 == N2.v0
      assert model.lam == N2.lam
      assert model.v_int == N2.v_int


  def test_params_match_the_library_fields_exactly():
      p = model.PARAMS
      assert p["reduced_mass"] == N2.mu
      assert p["impulsemomentum"] == N2.ell
      assert p["potential"] == {
          "D_0": N2.D0, "alpha_0": N2.alpha0, "R_0": N2.R0,
          "lambda_inf": N2.lambda_inf, "lambda_1": N2.lambda_1,
          "R_lambda": N2.R_lambda, "lambda_c": N2.lambda_c,
          "R_c": N2.R_c, "alpha_c": N2.alpha_c,
      }
  ```

  In `projects/n2_resonance/test_potential.py`, replace `test_matches_library_model_to_1e_12` (the layering-close version) with the same-shaped identity test over `potential.v0/lam/v_int`, keep a numeric `v_eff_el` check at 1e-12 against `N2.v_int(r, R) + N2.ell*(N2.ell+1)/(2*r**2)`, and relax the two exact assertions: `test_lambda_at_Rc_matches_config` → `abs(potential.lam(p["R_c"]) - p["lambda_c"]) < 1e-12`; `test_v0_at_R0_is_minus_D0` → `abs(potential.v0(p["R_0"]) + p["D_0"]) < 1e-12`.
- [ ] Run: `uv run --no-sync pytest validation/n2/test_model.py projects/n2_resonance/test_potential.py -q` — expect the identity tests to FAIL (`model.v0 == N2.v0` is False while `model.v0` is still a local function).
- [ ] Implement `validation/n2/model.py`: keep `PARAMS = json.loads((Path(__file__).parent / "config.json").read_text())` (provenance record) and `model_checks()`; delete the local `v0`/`lam`/`v_int` bodies and bind `v0 = N2.v0`, `lam = N2.lam`, `v_int = N2.v_int` (with `from qscat.model import N2`); reimplement `v_eff_el` as

  ```python
  def v_eff_el(r, R):
      """Fixed-R electronic effective potential incl. the centrifugal term.
      N2.surface includes v0(R); this deliberately does NOT, so it is
      N2.v_int + centrifugal (complex-safe -- see N2.v_int's docstring)."""
      rc = np.asarray(r, dtype=np.complex128)
      return N2.v_int(rc, R) + N2.ell * (N2.ell + 1) / (2 * rc**2)
  ```

  and in `model_checks()` replace every `float(v0(...))`/`float(lam(...))`/`float(v_int(...))`/`float(v_eff_el(...))` with `float(np.real(...))` (values are real to round-off for real inputs; `float()` on complex raises). Update the module docstring: this module is now a thin consumer of `qscat.model.N2` plus the provenance `PARAMS` and the `model_checks()` harness hook.
- [ ] Implement `projects/n2_resonance/potential.py` the same way: `v0 = N2.v0`, `lam = N2.lam`, `v_int = N2.v_int`, local `v_eff_el` as above (drop the `PARAMS["impulsemomentum"]` lookup for `ell` in favor of `N2.ell`); keep the Task-2-of-layering-close `PARAMS` dict; update the docstring ("re-implements the same formulas" is no longer true — it now IS the library model).
- [ ] Run: `uv run --no-sync pytest validation/n2 projects/n2_resonance projects/n2_ti_cross_section/test_vres.py -q -m "not slow"` — expect PASS (delegation is value-preserving; `resonance.py`'s `model.v_eff_el` consumer and `vres.py`'s `v0` consumer both handle complex output already).
- [ ] Run the full gate; the branch needs a `validate:n2` labelled run before merge (this touches the potential every `validation/n2` solver consumes).
- [ ] Commit: `git add validation/n2/model.py validation/n2/test_model.py projects/n2_resonance/potential.py projects/n2_resonance/test_potential.py && git commit -m "refactor(n2): single-source the N2 potential on qscat.model.N2; lockstep tests become identity checks"`

---

## Task 7 — Lock test: N₂ presets ↔ projects convergence grids (exp-M11)

`qscat_run.presets._n2_ti_grid`/`_n2_td_grid` transcribe `projects.n2_2d_cross_section.convergence.WORKING_GRID` and `projects.n2_2d_td_cross_section.convergence.TD_WORKING_GRID` by value (presets.py docstring, lines 6-9) because qscat_run must not import projects. Mirror `validation/diatomic/test_da_grid.py::test_diatomic_decks_match_presets`: the guarded-duplication lock lives in validation, which may import both sides.

**Files:**
- Create: `validation/n2/test_presets_lock.py`
- Test: itself

**Interfaces:**
- Consumes: `qscat_run.presets._n2_ti_grid() -> TensorGrid`, `qscat_run.presets._n2_td_grid() -> TensorGrid`, `projects.n2_2d_cross_section.convergence.working_tgrid() -> TensorGrid`, `projects.n2_2d_td_cross_section.convergence.td_working_tgrid() -> TensorGrid`.
- Produces: the drift lock.

**Steps:**

- [ ] Create `validation/n2/test_presets_lock.py`:

  ```python
  """Guard: qscat_run's N2 preset decks and the projects convergence grids must
  stay identical. The two exist separately because layering forbids qscat_run
  importing projects (see qscat_run/presets.py's docstring); this test is what
  keeps the transcription from drifting -- the N2 sibling of
  validation/diatomic/test_da_grid.py::test_diatomic_decks_match_presets.
  Byte-identical (np.array_equal), not allclose: both sides build through the
  same qscat.core.grids factories from the same literals."""

  from __future__ import annotations

  import numpy as np
  from qscat.dvr import TensorGrid

  from projects.n2_2d_cross_section.convergence import working_tgrid
  from projects.n2_2d_td_cross_section.convergence import td_working_tgrid


  def _assert_identical(a: TensorGrid, b: TensorGrid) -> None:
      for ga, gb in zip(a.grids, b.grids, strict=True):
          assert ga.n == gb.n
          assert ga.R0 == gb.R0
          assert np.array_equal(ga.points, gb.points)
          assert np.array_equal(ga.weights, gb.weights)
          assert np.array_equal(ga.real_points, gb.real_points)


  def test_n2_ti_preset_matches_projects_working_grid():
      import qscat_run.presets as presets

      _assert_identical(presets._n2_ti_grid(), working_tgrid())


  def test_n2_td_preset_matches_projects_td_working_grid():
      import qscat_run.presets as presets

      _assert_identical(presets._n2_td_grid(), td_working_tgrid())
  ```

- [ ] Run: `uv run --no-sync pytest validation/n2/test_presets_lock.py -q` — expect 2 PASS (the transcription is currently faithful: TI = r_max 16/order 7/n_complex 5 electronic × r_max 20/quad 10/n_complex 5 nuclear at 35°; TD = r_max 50/order 8/n_complex 6 × r_max 22/quad 10/n_complex 5). If either FAILS, that is a real pre-existing drift: STOP and report it — do not "fix" the test.
- [ ] Detection check (mirrors the layering-close canary discipline): temporarily change `_n2_ti_grid`'s `quadrature=10` to `11`, run — expect FAIL on `ga.n == gb.n`; revert, rerun — PASS.
- [ ] Run the full gate.
- [ ] Commit: `git add validation/n2/test_presets_lock.py && git commit -m "test(validation/n2): lock the qscat-run N2 preset decks to the projects convergence grids"`

---

## Task 8 — qscat-run config: missing sub-keys raise `ConfigError`, never `KeyError` (exp-M4)

**Files:**
- Modify: `apps/qscat-run/qscat_run/config.py` (`_require` helper; `_load_observables`, `_load_energies`, `_load_td`, `_load_segment`/`_load_grid` wrap their raw subscripts)
- Test: `apps/qscat-run/tests/test_config.py` (new `TestMissingSubkeys` class)

**Interfaces:**
- Consumes: nothing new.
- Produces: `_require(raw: dict[str, Any], key: str, where: str) -> Any` (module-private); `_load_segment(raw: dict[str, Any], where: str) -> SegmentSpec` (gains the `where` parameter; `_load_grid` passes `"grid.electronic"` / `"grid.nuclear"`).

**Steps:**

- [ ] Write the failing tests in `apps/qscat-run/tests/test_config.py` — the three cases the release review reproduced, plus the incident one for coverage of every wrapped loader:

  ```python
  class TestMissingSubkeys:
      """A typo'd or missing sub-key must produce an actionable ConfigError
      naming the block and the key -- never a bare KeyError (exp-M4)."""

      def _load(self, tmp_path: Path, body: str):
          return load_config(_write(tmp_path, body))

      def test_observable_missing_kind_is_actionable(self, tmp_path: Path) -> None:
          # the reproduced `{kine: ve}` typo
          with pytest.raises(ConfigError, match=r"observables\[0\].*'kind'"):
              self._load(
                  tmp_path,
                  """
              molecule: F2
              methods: [ti]
              observables: [{kine: ve}]
              output_dir: out
          """,
              )

      def test_td_missing_n_steps_is_actionable(self, tmp_path: Path) -> None:
          with pytest.raises(ConfigError, match=r"'td'.*'n_steps'"):
              self._load(
                  tmp_path,
                  """
              molecule: F2
              methods: [ti, td]
              observables: [{kind: ve, channels: 2}]
              td: {dt: 0.5}
              output_dir: out
          """,
              )

      def test_energies_missing_step_is_actionable(self, tmp_path: Path) -> None:
          with pytest.raises(ConfigError, match=r"'energies'.*'step'"):
              self._load(
                  tmp_path,
                  """
              molecule: F2
              methods: [ti]
              observables: [{kind: ve, channels: 2}]
              energies: {min: 0.01, max: 0.05}
              output_dir: out
          """,
              )

      def test_td_incident_missing_sigma_is_actionable(self, tmp_path: Path) -> None:
          with pytest.raises(ConfigError, match=r"'td.incident'.*'sigma'"):
              self._load(
                  tmp_path,
                  """
              molecule: F2
              methods: [ti, td]
              observables: [{kind: ve, channels: 2}]
              td: {dt: 0.5, n_steps: 10, incident: {r0: 20.0, p0: -0.4}}
              output_dir: out
          """,
              )
  ```

- [ ] Run: `uv run --no-sync pytest apps/qscat-run/tests/test_config.py::TestMissingSubkeys -q` — expect 4 FAIL, each with a raw `KeyError` (`'kind'`, `'n_steps'`, `'step'`, `'sigma'`) instead of a `ConfigError`.
- [ ] Implement in `apps/qscat-run/qscat_run/config.py`, above `_load_observables`:

  ```python
  def _require(raw: dict[str, Any], key: str, where: str) -> Any:
      """`raw[key]`, or a ConfigError naming the block and the missing key --
      a hand-written config's KeyError names the key but not where it was
      expected, which is the difference between actionable and not."""
      if key not in raw:
          raise ConfigError(f"'{where}' is missing required key '{key}'")
      return raw[key]
  ```

  Then wrap every raw subscript in the loaders:
  - `_load_observables`: `enumerate(raw)`; reject a non-mapping item (`raise ConfigError(f"observables[{i}] must be a mapping with a 'kind' key")`); `kind=str(_require(item, "kind", f"observables[{i}]"))`.
  - `_load_energies` (non-`values` branch): `min`/`max`/`step` via `_require(raw, k, "energies")`.
  - `_load_td`: `dt`/`n_steps` via `_require(raw, k, "td")`; incident `r0`/`p0`/`sigma` via `_require(ir, k, "td.incident")`; flat test-function `p0_out`/`sigma_out` via `_require(tr, k, "td.test_function")` (`"r0_out" in tr` stays the branch discriminator); per-kind blocks via `_require(block, k, f"td.test_function.{kind}")` for all three keys.
  - `_load_segment(raw, where)`: add the `where` parameter; `real`/`ecs` via `_require(raw, k, where)`; `angle`/`elements`/`quadrature` via `_require(ecs_raw, k, f"{where}.ecs")`; `_load_grid` calls `_load_segment(raw["electronic"], "grid.electronic")` / `..., "grid.nuclear")`.
- [ ] Run: `uv run --no-sync pytest apps/qscat-run/tests/test_config.py -q` — expect PASS (the whole file: the new class green, no existing test broken — the happy paths never hit `_require`'s raise).
- [ ] Run the full gate.
- [ ] Commit: `git add apps/qscat-run/qscat_run/config.py apps/qscat-run/tests/test_config.py && git commit -m "fix(qscat-run): missing config sub-keys raise actionable ConfigError, never KeyError"`

---

## Task 9 — Delete `observation.py` + its test; keep `nuclear_density.py` (exp-M9)

Verified importer census (2026-08-25): `projects/n2_2d_td_cross_section/observation.py` is imported ONLY by `projects/n2_2d_td_cross_section/test_observation.py`; `libs/qscat/qscat/core/plot.py:12` and `libs/qscat/qscat/core/time_dependent.py:187` mention it in docstrings only. `nuclear_density.py` is imported only by its own test, but is KEPT — see "Design decisions" (its exact-vs-LCP density-shape comparison has no artifacts.py replacement, and `docs/physics/n2-2d-cross-section.md` documents its findings); after Task 3 it no longer inlines any solve.

Capability → replacement mapping (goes into `docs/physics/n2-2d-td-cross-section.md`):

| `observation.py` capability | Replacement in `qscat_run` |
|---|---|
| `save_numeric_outputs` (`t`/`c`/`sigma_E` npz + dt/wp metadata) | `artifacts.py`: `cross_section.{csv,npz}` (sigma per key), `correlations.npz` (raw per-step `t`/`c` — `artifacts.correlations: true`), `config.resolved.yaml` + `manifest.json` (dt/wp/grid provenance) |
| `plot_snapshots` (rho(R,t)/rho(r,t) panels, real-region-masked, R0 marked) | `artifacts.py`: `wavefunction/psi_*.{npz,png}` via `artifacts.wavefunction_snapshots.td_times` |
| `plot_correlation` (per-channel c(t) figure) | `correlations.npz` raw series (npz only — no committed figure; accepted, the arrays are the deliverable) |
| `plot_sigma_vs_ti` (TD-vs-TI overlay + usable-window shading) | `cross_section.png` overlays `ti:`/`td:` keys from a `methods: [ti, td]` run + `reference:` overlays. NOT replaced: the `usable_window` shading (the function itself survives in `projects/n2_2d_td_cross_section/convergence.py`) and the `norm(t)` decay panel — both recorded as accepted losses |

**Files:**
- Delete: `projects/n2_2d_td_cross_section/observation.py`, `projects/n2_2d_td_cross_section/test_observation.py`
- Modify: `libs/qscat/qscat/core/plot.py` (docstring line 12), `libs/qscat/qscat/core/time_dependent.py` (docstring line 187), `docs/physics/n2-2d-td-cross-section.md` (lines 5, 229, 407 reference `observation.py`/`test_observation.py`; add the mapping table above), `CLAUDE.md` if its `n2_2d_td_cross_section` repo-map line names `observation.py`
- Test: full fast tier (catches any missed importer)

**Interfaces:**
- Consumes/Produces: deletions only; no API change.

**Steps:**

- [ ] Re-verify the census at execution time: `grep -rn "observation" --include="*.py" projects validation libs apps benchmarks | grep "n2_2d_td_cross_section"` — expect hits only in the two files being deleted (plus the two library docstrings).
- [ ] `git rm projects/n2_2d_td_cross_section/observation.py projects/n2_2d_td_cross_section/test_observation.py`
- [ ] Update the two library docstrings: `qscat/core/plot.py:12` (drop the "unlike `projects/n2_2d_td_cross_section/observation.py`" comparison — name `qscat_run.artifacts` as the serialization layer instead) and `qscat/core/time_dependent.py:187` (same substitution).
- [ ] Update `docs/physics/n2-2d-td-cross-section.md`: remove `observation.py` from the module list (line 5) and the test inventory (line 407); at the line-229 discussion note the module was retired into `qscat_run.artifacts` and add the capability→replacement table above, including the two named accepted losses (usable-window shading, norm(t) panel). Check `grep -n "observation" CLAUDE.md` and update the repo-map line if it names the file.
- [ ] Run the full gate (fast tier proves nothing imported it).
- [ ] Commit: `git add -u projects/n2_2d_td_cross_section && git add libs/qscat/qscat/core/plot.py libs/qscat/qscat/core/time_dependent.py docs/physics/n2-2d-td-cross-section.md CLAUDE.md && git commit -m "refactor(lifecycle): retire observation.py -- every capability now lives in qscat_run.artifacts (losses documented)"`

---

## Task 10 — The libs test suite runs from the sdist (lib-m15 / exp-m5)

Move `libs/qscat/tests/test_lcp.py::test_matches_n2_vres_oracle` (lines 72-89, the only `projects` import left in `libs/qscat/tests` after Tasks 1-9) to `validation/n2/`, keeping it verbatim: it is a genuine cross-implementation oracle (`qscat.core.lcp.local_complex_potential` vs the independent projects pole walk `vres_on_grid`) and must not be rewritten against qscat itself.

**Files:**
- Create: `validation/n2/test_lcp_vres_parity.py`
- Modify: `libs/qscat/tests/test_lcp.py` (delete the moved test + its `importorskip` preamble)
- Test: `validation/n2/test_lcp_vres_parity.py` (`@slow`), `libs/qscat/tests/test_lcp.py`

**Interfaces:**
- Consumes: `qscat.core.lcp.local_complex_potential(model, nuclear_grid, elec_a, elec_b) -> (Vd, Gamma)`, `projects.n2_ti_cross_section.vres.vres_on_grid(grid) -> (Vd, Gamma)`, `projects.n2_ti_cross_section.nuclear_grid.n2_nuclear_grid()`, `qscat.core.grids.electronic_grid`, `qscat.model.N2`.
- Produces: a projects-free `libs/qscat/tests` (sdist-runnable).

**Steps:**

- [ ] Create `validation/n2/test_lcp_vres_parity.py` — the moved test verbatim, with a local `_elec_grids()` (copied from `test_lcp.py:12-16`), plain top-level imports (no `importorskip` — projects always exists in the monorepo, and this file ships only in the repo), `pytestmark = pytest.mark.slow`, and a module docstring: "moved from libs/qscat/tests/test_lcp.py so the library suite runs from the sdist; the projects pole walk is the independent oracle for `local_complex_potential`, and validation is the layer allowed to import both."
- [ ] In `libs/qscat/tests/test_lcp.py`: delete `test_matches_n2_vres_oracle` (lines 72-89). Update `qscat/core/lcp.py`'s module docstring line 22 (it names `test_lcp.py::test_matches_n2_vres_oracle`) to the new location.
- [ ] Run: `grep -rn "projects" libs/qscat/tests --include="*.py"` — expect docstring/comment prose at most, no imports.
- [ ] Run: `uv run --no-sync pytest libs/qscat/tests/test_lcp.py -q -m "not slow"` — expect PASS; run `uv run --no-sync pytest validation/n2/test_lcp_vres_parity.py -q -m slow` once locally — expect 1 PASS.
- [ ] Run the full gate.
- [ ] Commit: `git add validation/n2/test_lcp_vres_parity.py libs/qscat/tests/test_lcp.py libs/qscat/qscat/core/lcp.py && git commit -m "test: move the vres cross-implementation oracle to validation; libs suite is projects-free"`

---

## Completion checklist

- [ ] `grep -rn "n2_ti_cross_section.cross_section" --include="*.py" .` → no matches (module deleted, all consumers on `qscat.core.lcp.lcp_ve_cross_section`).
- [ ] `grep -rn "import projects\|from projects" libs/qscat/tests apps/qscat-run` → no matches.
- [ ] `tests/test_layering.py` (from the layering-close plan) still green — none of this plan's edits re-introduced a projects→validation edge.
- [ ] Fast tier, mypy, ruff, ruff-format clean after EVERY task (Global Constraints), plus one local run of the new `@slow` tests (Tasks 5, 10).
- [ ] `validate:n2` and `validate:diatomic` labelled slow-tier runs green before merge (Tasks 3 and 6 touch those solver paths).
- [ ] Findings ledger: exp-M1 (Tasks 1-5), exp-M10 (Task 3), exp-M2 (Task 6), exp-M11 (Task 7), exp-M4 (Task 8), exp-M9 (Task 9), lib-m15/exp-m5 (Task 10) — all closed.
