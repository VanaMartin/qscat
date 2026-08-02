# qscat Hardening & Optimization Roadmap

**Date:** 2026-08-02
**Author:** Martin (martin@qscat.com) with Claude
**Status:** Review complete — roadmap for approval
**Scope:** A harsh structural/logic audit of `libs/qscat` and the monorepo, turned
into a prioritized plan to (1) clean up code + docs to top-tier-package standard,
(2) publish `qscat` on PyPI, (3) establish CPU-optimized and GPU-ready base Docker
images, and (4) narrow down where to invest optimization effort.

## Method

A panel review: the **packaging/PyPI dimension** was audited by a separate Fable
model reviewer that *empirically* built the wheel (`uv build`), installed it into a
pristine venv, and inspected METADATA/RECORD — findings are verified, not
speculative. The **structure/logic, documentation, build/Docker, and optimization**
dimensions were audited directly against the source. (A five-agent Fable panel was
intended; the session's subagent cap was hit after the first, so four dimensions
were reviewed inline.)

## Verdict

**Structure: B. Distribution: F.** The monorepo layout, the `qscat.core`/`qscat.model`
split, the clean internal dependency DAG (`units/linalg/special` → `dvr/ecs/evolution`
→ `core` → `tuning`, with `core → model` TYPE_CHECKING-only and test-enforced), and
the differential-oracle discipline are genuinely better than most research code. But
the distribution layer barely exists, and the gap is stark: the wheel that builds
today is **legally unusable (no license), has a blank PyPI page (no metadata/README),
a meaningless version (0.0.0 ×4), ships no types, has no CI — and two of its eight
submodules crash on import in a clean environment.** Nobody has installed this package
outside its own workspace venv, and it shows. Every blocker is cheap; that they are
all still open means no release-readiness pass has ever been run.

---

## Part 0 — Blockers (must fix before ANY release)

Ordered by how much they disqualify a release. All are small.

| # | Blocker | Evidence | Fix |
|---|---------|----------|-----|
| B1 | **`import qscat.core` / `import qscat.tuning` crash on clean install** — `core/plot.py` imports matplotlib at module scope, `core/__init__.py` eagerly re-exports it, `tuning/*` imports `core`; matplotlib is not a declared dependency (works only via dev-group leakage). | Verified in a clean venv: `ModuleNotFoundError: No module named 'matplotlib'`. | Make plotting lazy (import matplotlib *inside* `plot_cross_sections`; drop/guard the eager re-export in `core/__init__.py`), add a `plot = ["matplotlib>=3.9"]` extra, and add a **clean-venv import smoke test** to CI so it can't recur. |
| B2 | **No LICENSE anywhere** → all-rights-reserved → downstream cannot lawfully use it. | No `LICENSE`/`license` key at root or `libs/qscat`; wheel METADATA has no `License-Expression`. | Add root `LICENSE` (BSD-3-Clause, the scientific-Python norm) + `license = "BSD-3-Clause"` + `license-files = ["LICENSE"]` (PEP 639) in `libs/qscat/pyproject.toml`. |
| B3 | **Metadata is a 300-byte stub; no `libs/qscat/README.md`.** PyPI page would be blank. | No readme/authors/classifiers/keywords/urls. | Write `libs/qscat/README.md` (what it is, install, 10-line example, MUMPS caveat, citation) + fill `[project]`: readme, authors, classifiers (`Topic :: Scientific/Engineering :: Physics`, `Intended Audience :: Science/Research`, `License ::`, `Programming Language :: Python :: 3.12`), keywords, `[project.urls]`. |
| B4 | **Version is four hardcoded `0.0.0`s and zero git tags.** No single source of truth; release = editing files in lockstep. | `libs/qscat/pyproject.toml:3`, `qscat/__init__.py:3`, native + apps + root. `git tag` empty. | Single-source via hatch (`dynamic = ["version"]`, `[tool.hatch.version] path = "qscat/__init__.py"`) or hatch-vcs from tags. Adopt `qscat-v0.1.0` tags. Add `CHANGELOG.md` (Keep-a-Changelog). First upload is a deliberate `0.1.0`. |
| B5 | **No `.github/` at all — no CI, no build check, no publish.** No automated process has ever run the suite or built artifacts. | Directory absent. | Stand up workflows: test matrix (3.12/3.13); `uv build` + `uvx twine check dist/*` on every PR; tag-triggered PyPI **Trusted Publishing** (OIDC, no tokens). |
| B6 | **No `py.typed`** → every `mypy --strict` pass is discarded at the wheel boundary; downstream gets `import-untyped`. | Marker absent; confirmed by the recurring CLI mypy notices. | Add empty `qscat/py.typed`; verify it lands in RECORD. Cheapest high-value fix in the audit. |

**Part 0 deliverable:** a green `0.1.0` on **TestPyPI** installable into a clean venv,
importable in full, typed, licensed, with a real PyPI page.

---

## Part 1 — Library structure & logic (the primary focus)

The bones are right; these raise it to top-tier. Ranked by leverage.

1. **Context-object the wide, array-threading signatures.** `ve_cross_section`,
   `da_cross_section`, `propagate`, `td_ve_cross_section` take 7+ loose positional
   args (`tgrid, model, eps, chi, v_init, vprimes, E, …`). `eps`/`chi`/`v_init` are
   one thing — a vibrational basis — recomputed and re-passed everywhere. Introduce a
   small frozen `VibrationalBasis` (or `ScatteringProblem`) dataclass bundling
   `tgrid, model, eps, chi, v_init`; solvers take `(problem, vprimes, E, *, ...)`.
   Biggest single readability/ergonomics win; also shrinks the CLI runner.
2. **Add a typed exception hierarchy.** 72 bare `ValueError`/`KeyError`/`RuntimeError`
   in `qscat`. Top packages raise domain errors (`scipy.linalg.LinAlgError`). Add
   `qscat/exceptions.py`: `QscatError(Exception)` base + `GridError`, `ConvergenceError`,
   `BackendError`, `ModelError`. Makes failures catchable and self-documenting; the
   CLI's `ConfigError` can subclass it.
3. **Split the two oversized modules.** `core/td_extractors.py` (1233 LOC — three
   extractor classes × electronic/nuclear axes) and `core/time_dependent.py` (879 LOC
   — propagate engine + Tannor-Weeks transform + td_ve/td_da orchestration) each do too
   much. Split `time_dependent.py` into `propagation.py` (the engine) + `td_cross_section.py`
   (the observable orchestration); consider one file per extractor family in
   `td_extractors/`. Improves your own editability and review reliability.
4. **Lazy top-level namespace (PEP 562).** `qscat/__init__.py` is 3 lines — `import qscat;
   qscat.dvr` raises AttributeError. Add `__getattr__` lazily exposing the eight
   submodules + `__all__` (scipy's pattern) so the namespace is discoverable without
   forcing heavy imports. (Ties to B1/B6.)
5. **Move `ecs_map` out of `ecs/__init__.py`.** Implementation (38 LOC + numpy) lives
   in the package `__init__`, unlike every sibling. Move to `ecs/map.py`, re-export;
   keep `__init__` files as pure API manifests.
6. **Write the public-API stability policy before v0.1.0.** `core/__init__.py`
   re-exports ~30 names; each freezes as public API at first publish. Decide
   public-vs-provisional now; mark provisional APIs; document deprecation policy.
7. **Resolve the `hypothesis` gap.** Declared as a dev dep, **zero usage** in lib tests
   though the lifecycle claims property-based testing. Either add property tests where
   they pay (grid invariants, `kron_sum`, ECS map round-trips, quadrature-weight sums)
   or drop the dependency. Prefer adding them — cheap coverage on the D-general layer.

---

## Part 2 — Documentation & developer experience

Currently: prose docstrings (only ~8 structured numpydoc/google section markers in all
of `qscat` → autodoc renders no param/return tables), no rendered docs site, no
contributor/citation surface. The physics notes in `docs/physics/` are excellent
*domain* content but are not user docs.

1. **Adopt numpydoc docstrings** (Parameters/Returns/Raises/Examples/Notes) across the
   public API. Enforce with the `numpydoc` validation hook + `ruff`'s `D` (pydocstyle)
   rules. This is the precondition for autodoc.
2. **Stand up a docs site: Sphinx + MyST + `numpydoc`/napoleon + `sphinx.ext.autodoc` +
   `intersphinx`** (cross-linking numpy/scipy), hosted on Read the Docs or GH Pages.
   (Sphinx is the scientific-Python standard — numpy/scipy/astropy; mkdocs-material +
   mkdocstrings is the lighter alternative if preferred.) Structure:
   - **Getting started** (install incl. the MUMPS caveat, first cross section in 10 lines)
   - **Tutorials** (VE cross section; TI vs TD; DA/DR; using the tuner; the CLI)
   - **API reference** (autodoc per submodule)
   - **Theory** (link/adapt `docs/physics/` — FEM-DVR-ECS, ECS, the extractors)
   - **Design/ADRs** (the core/model split, backend dispatch)
3. **Add the OSS front matter:** `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, issue/PR
   templates, README badges (CI/PyPI/docs/license), and a **`CITATION.cff`** — this is
   scientific software; citability is a feature.
4. **Separate agent-ops from user docs.** `CLAUDE.md` is the agent operating manual and
   should stay that (and never be published); user-facing narrative belongs in the docs
   site and `libs/qscat/README.md`.
5. **Define units/notation once.** Atomic units and FEM-DVR-ECS notation should be a
   single "conventions" page the API docs link to, not restated per module.

---

## Part 3 — Packaging & PyPI release engineering

Beyond the Part-0 blockers:

- **`mumps` extra is a laptop landmine** — `pip install "qscat[mumps]"` fails without a
  system MUMPS. Document the prerequisite in the README; point at conda-forge's prebuilt
  `mumps` as the supported route.
- **De-duplicate dependency pins.** Runtime deps are re-pinned in the root dev group —
  the exact duplication that masked B1. Dev group carries tooling only.
- **Lowest-resolution CI leg** (`uv run --resolution lowest-direct pytest`) so
  `numpy>=2`/`scipy>=1.14` floors are real, not just whatever `uv.lock` resolves.
- **Name risk:** `qscat` is free on PyPI *today* but collides in mindshare (a QGIS
  shoreline tool; NASA QuikSCAT). Register early; make metadata clearly distinguish it.
- **Decide the Rust coupling policy now (see Part 5).** Today `qscat` never imports
  `qscat_kernels`, so the wheel is accidentally clean pure-Python (`py3-none-any`).
  Lock in: **pure-Python core forever + optional binary via a `qscat[kernels]` extra +
  try-import-else-Python-oracle** in every kernel. Ship `py.typed` in the same PR that
  first wires a kernel; add `.pyi` stubs when `qscat_kernels` publishes.
- **Housekeeping:** commit `native/qscat-kernels/Cargo.lock` (binary-wheel repro);
  explicit `[tool.hatch.build.targets.sdist]` include/exclude (drop the stray
  `.gitignore` from the sdist); comment the root `[project]` as "virtual root, never
  published"; pin `qscat` with a compatible range in `apps/qscat-run` at release.

---

## Part 4 — Docker: CPU-optimized + GPU-ready base images

Today: one layered CPU base (`base.Dockerfile`, OpenBLAS/LAPACKE/FFTW3 + system MUMPS,
ARM-friendly) → `Dockerfile` (build/test/runtime). Good foundation. Target:

1. **Two CPU base flavors, one swappable ABI.** Keep **`base-cpu-portable`**
   (OpenBLAS/FFTW3, ARM/Graviton-friendly, current default) and add
   **`base-cpu-mkl`** (Intel MKL + MKL PARDISO on x86-64) — the "CPU-optimized" image.
   Code already targets the standard CBLAS/LAPACKE/FFTW3 ABIs, so this is a
   `BASE_IMAGE`/`ARG VENDOR` swap. This also seeds the PARDISO optimization path (Part 5).
2. **GPU base scaffold now, fill later.** Add `base-gpu.Dockerfile` `FROM nvidia/cuda:*-runtime`
   with cuBLAS/cuSPARSE/cuSOLVER (and cuDSS for sparse), wired through the same
   `BASE_IMAGE` ARG so app layers are unchanged. Do not implement GPU kernels yet — just
   fix the layer/ARG seam so the future GPU sparse backend (Part 5) slots in like MUMPS did.
3. **Clean up:** remove the obsolete `docker/run-n2.sh` (superseded by the general
   `docker/run.sh`); document the three base variants in `docker/README.md`.
4. **Docker is for run/test, not wheels.** Build PyPI artifacts with
   `cibuildwheel`/`maturin`-generated CI (manylinux/macos/windows), not these images.

---

## Part 5 — Optimization: narrowed directions

**Current state:** the "optimize-in-Rust" lifecycle stage is *unstarted* — `qscat_kernels`
is a 13-line `l2_norm` stub, not imported anywhere. So this is greenfield; invest where
the 2-D/ion decks actually spend time. Reference bar (from prior work): eMoScat's MKL
PARDISO did all of N2/NO/F2 in <1 hr.

**Step 0 — profile before cutting.** There is no profiling harness. Extend `benchmarks/`
(already has `mumps_vs_superlu`, `sweep_reuse`) with a `py-spy`/`scalene` pass on the
**H2P full deck (~1.15M unknowns)** and a converged **F2 2-D** run to *confirm* the split
below with numbers, not intuition. Gate every optimization against its Python
differential oracle (the lifecycle already mandates keeping it alive).

**Priority-ranked targets:**

1. **Sparse factorization (`qscat.linalg.SparseLU`) — the #1 hot path.** Dominates the
   2-D driven solves and is the one-time cost of every propagation (factor-once). MUMPS
   already beats SuperLU 72×/9× and `refactor` reuses the symbolic analysis across energy
   sweeps. Next: (a) benchmark **MKL PARDISO** on x86 (the eMoScat bar) as a third backend
   behind the existing dispatch; (b) profile whether ordering (SCOTCH/METIS) or numeric
   factor dominates at the largest decks; (c) design a **GPU sparse backend (cuDSS/cuSOLVER)**
   as a fourth dispatch option — the backend protocol already exists, so this is additive.
2. **The propagation inner loop (Padé/CN steppers + extractor `record`).** Thousands of
   steps, each a cached-factor triangular solve plus per-step `c_product` projections in
   the extractors. Candidates for the **first real Rust kernel**: the `record`
   accumulation (tight complex BLAS-1/2 over the recorded rows) and batched RHS solves.
   High call count, small hot code, clean differential oracle — ideal kernel shape.
3. **`mpmath` Coulomb special functions (`qscat.special.coulomb`).** The H2P DR path loops
   over Rydberg channels calling `coulomb_f/g/h1_en`; mpmath is slow. Vectorize/cache, or
   move to a fast implementation (GSL/Boost via a kernel). Real sink specifically for ions.
4. **Hamiltonian assembly (`kinetic_nd`/`potential_nd`/`kron_sum`).** One-time but scales
   with grid size; sparse assembly in Rust could help only the very largest decks — verify
   with the profiler before investing.
5. **Keep the D-general layer backend-agnostic.** `linalg`/`dvr`/`tensor` must not bake in
   CPU-only assumptions, so the GPU sparse backend and array-API (cupy) paths stay a
   dispatch choice, not a rewrite. This is the structural precondition for the GPU phase.

**Not worth optimizing now:** the tuner (`qscat.tuning`) is setup-time; the CLI is I/O.

---

## Suggested sequencing

- **Phase A (days): Part 0 blockers** → `0.1.0` on TestPyPI. Highest value, lowest effort.
- **Phase B: Part 1 structure/logic** (exceptions, context object, module splits, lazy
  namespace) — do before the API freezes at first public release.
- **Phase C: Part 2 docs site + OSS front matter + numpydoc.**
- **Phase D: Part 3 release engineering** → first real **PyPI** `0.1.0` (or `0.2.0`).
- **Phase E: Part 4 Docker** CPU-MKL + GPU scaffold.
- **Phase F: Part 5 optimization** — profile → PARDISO backend + first propagation-loop
  Rust kernel + mpmath Coulomb → GPU sparse backend design.

Each phase is its own spec/plan; Phase A should start immediately. Phases B and C can run
in parallel with A once the blockers are green.

## Detailed sub-project specs

Parts 0, and the safe subsets of 1–5, are done (see the CHANGELOG and the git
history). The remaining implementation work is captured as focused specs:

- `2026-08-02-sparse-lu-pardiso-backend.md` — MKL PARDISO backend (Part 5 #1).
- `2026-08-02-gpu-sparse-backend.md` — GPU cuDSS backend (design-only; deferred).
- `2026-08-02-first-rust-kernel.md` — the optimize-in-Rust starter (Part 5 #2).
- `2026-08-02-scattering-context-refactor.md` — full functional-signature
  refactor onto the context object (the deferred Part 1 piece; SDD-suited).
- `2026-08-02-docs-numpydoc-and-hosting.md` — numpydoc conversion + Pages (Part 2).

Also outstanding (no separate spec needed — mechanical): the `td_extractors.py`/
`time_dependent.py` module splits; extending the exception migration to deeper
internal raises; dev-group dependency de-dup + a lowest-resolution CI leg.
