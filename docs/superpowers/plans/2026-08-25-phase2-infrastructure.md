# Phase 2 Infrastructure Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close the infrastructure findings from the 2026-08-25 release review:
coverage reporting in the fast gate, an OS/version test matrix, honest Rust-crate
metadata with the crate's CI cost made conditional, honest qscat-run metadata,
community/security files, a CITATION.cff version gate, a mypy config that stops
overstating its scope, a correct `validate:n2` path map, and digest-pinned Docker
base images.

**Architecture:** All changes are CI/workflow/metadata-level; no solver code moves.
`ci.yml` gains a `changes` job (dorny/paths-filter) whose output gates the Rust
toolchain + kernel build, a wider test matrix, coverage on one matrix cell, and a
CITATION version assertion. `validation.yml`'s SUITES map is narrowed and gains a
`factory` suite (the cover test forces one — `projects/potential_factory` holds
slow tests that today ride on the over-broad `projects` prefix). New files:
`.github/dependabot.yml`, `SECURITY.md`, `.github/ISSUE_TEMPLATE/bug.yml`,
`.github/PULL_REQUEST_TEMPLATE.md`, `docs/adr/0006-*.md`.

**Tech Stack:** GitHub Actions, uv, pytest-cov, dorny/paths-filter@v3, maturin/PyO3, Docker.

**Spec:** the "Findings addressed" section below (self-contained; from the 2026-08-25 release review)

## Global Constraints

- **PyPI release DEFERRED until the peer-reviewed citation article publishes — distribution is repo-only; NO task may register publishers, tag releases, or claim pip-installability.** When the article is out, CITATION.cff gains it as preferred-citation (out of scope here). `publish.yml` stays as-is (tag-triggered, and nobody may push a `qscat-v*` tag under this policy); no task edits it.
- After every task: `uv run --no-sync pytest -m "not slow" -n auto --dist loadfile` green; `uv run --no-sync mypy libs/qscat/qscat apps/qscat-run/qscat_run` clean; `uv run --no-sync ruff check .` + `ruff format --check .` clean; `uv run --no-sync pytest tests/test_docs_portability.py -q` green.
- Workflow YAML edits validated with `yaml.safe_load` + a pushed draft run where noted.
- Never `git commit -a`.
- Tasks 1–4 edit `.github/workflows/ci.yml`; run them in order, not in parallel, or the YAML merges will conflict.

## Findings addressed

- **style-M1 (coverage):** no coverage measurement exists anywhere. Add pytest-cov to the dev group; the fast-tier CI job gains `--cov=qscat --cov-report=term` and a summary in `$GITHUB_STEP_SUMMARY`. Decision required on fail-under.
- **style-M5 (matrix):** CI tests only ubuntu × {3.12, 3.13}. Add a macos-latest row (respecting ci.yml:58-71's own finding that macOS scipy links Accelerate, so the Linux OpenBLAS pins must not be copied blindly), a Python 3.14 row, and a minimum-versions job pinning the floors from `libs/qscat/pyproject.toml`.
- **style-M3/N6 (Rust crate):** `native/qscat-kernels` is a stub (`l2_norm`) whose toolchain install + maturin build is paid on every CI job, and its metadata (no license, no repository, edition 2021, pyo3 0.22) misrepresents it. DEFAULT option: keep the crate in the workspace, make its CI cost conditional on `native/**` changes, fix the metadata. The alternative — implement a real kernel now — is explicitly deferred (see Task 3's ADR).
- **style-M6 (qscat-run):** `apps/qscat-run/pyproject.toml` has placeholder metadata (version 0.0.0, no license, no authors, no readme). Given repo-only distribution, give it honest metadata and a README line stating it is repo-only; no publish path.
- **style-M7 (community files):** no dependabot config, no SECURITY.md, no issue/PR templates.
- **style-M4 remainder:** nothing asserts `CITATION.cff`'s `version:` equals `qscat.__version__` (the file's own comment says "keep in step" by hand).
- **style-N11 (mypy scope):** root `[tool.mypy] strict = true` claims the whole repo while only `libs/qscat/qscat` + `apps/qscat-run/qscat_run` are actually clean; `projects/`, `validation/`, `benchmarks/` would fail. Encode the real scope as per-module overrides.
- **recent-m5 (validate:n2 paths):** `validation.yml`'s `"n2": "validation/n2 projects"` claims ALL of `projects/` for the n2 suite — including `projects/potential_factory`, which is not N₂. Narrow it to the actual n2 project dirs without orphaning potential_factory's slow tests (the cover test `tests/test_validation_suites.py` forbids orphans).
- **style-N7 (promoted rider):** the two Docker base images (`docker/base.Dockerfile` line 5 and `docker/Dockerfile`'s runtime stage) are pinned by mutable tag, not digest.

---

## Task 1: Coverage in the fast gate (style-M1)

**Files:**
- Modify: `pyproject.toml` (dev group)
- Modify: `.github/workflows/ci.yml` (test job)

**Interfaces:** none (CI + dev-deps only).

**Decision (state it, then implement it): report-only first, threshold later.**
No baseline exists, so any `--fail-under` today would be a number invented on the
spot — either vacuous or an instant red gate. Ship report-only now; revisit around
2026-09-25 with a month of numbers and set `fail_under` in a `[tool.coverage.report]`
block then. Record this decision in the ci.yml comment added below so the next
reader knows the threshold is pending by choice, not oversight.

**Steps:**

- [ ] Add pytest-cov to the dev group in `pyproject.toml`, directly after the
  `pytest-xdist` entry (pytest-cov's subprocess support covers the xdist workers):

  ```toml
      # Coverage on the fast tier (CI reports it in the step summary).
      # Report-only for now -- no fail-under until a baseline month of numbers
      # exists (decision 2026-08-25, revisit ~2026-09-25).
      "pytest-cov>=6",
  ```

- [ ] Run `uv sync --all-packages` so the lockfile picks it up (commit `uv.lock` with the task).

- [ ] In `.github/workflows/ci.yml`, give the test job's matrix a coverage
  variable on exactly one cell (paying the coverage overhead once, not per row)
  and append it to the pytest invocation. The test job's `strategy`/`steps`
  become (context lines shown; this is the exact end state of the touched
  region — the matrix `include` is also where Task 4 adds rows):

  ```yaml
      strategy:
        fail-fast: false
        matrix:
          python-version: ["3.12", "3.13"]
          include:
            # Coverage is collected on one cell only: it costs measurable
            # overhead and the number is the same on every row.
            - python-version: "3.12"
              coverage: "--cov=qscat --cov-report=term"
  ```

  and the run step + a new summary step:

  ```yaml
        - name: Test (fast suite; MUMPS tests skip without a system MUMPS)
          env:
            OMP_NUM_THREADS: "1"
            OPENBLAS_NUM_THREADS: "1"
            MKL_NUM_THREADS: "1"
          run: uv run --no-sync pytest -q -m "not slow" -n auto --dist loadfile ${{ matrix.coverage }}
        - name: Coverage summary
          if: matrix.coverage != ''
          run: |
            {
              echo "### Coverage — qscat, fast tier (report-only; no fail-under until a baseline exists)"
              echo
              uv run --no-sync python -m coverage report --format=markdown
            } >> "$GITHUB_STEP_SUMMARY"
  ```

- [ ] Validate: `python3 -c "import yaml; yaml.safe_load(open('.github/workflows/ci.yml'))"`.

- [ ] Verify locally that the flag combination works with xdist before pushing:
  `uv run --no-sync pytest -q -m "not slow" -n auto --dist loadfile --cov=qscat --cov-report=term tests/` (a subset is fine locally; the point is that pytest-cov + xdist + loadfile coexist and a non-zero total prints).

- [ ] Add `.coverage` and `htmlcov/` to `.gitignore` if not already present.

- [ ] Run the four Global Constraints checks. Push the branch and confirm the CI
  draft run's step summary shows the markdown coverage table (noted verification).

## Task 2: Honest Rust crate metadata + current pyo3 (style-M3/N6, metadata half)

**Files:**
- Modify: `native/qscat-kernels/Cargo.toml`
- Modify: `native/qscat-kernels/pyproject.toml`

**Interfaces:** the compiled module name `qscat_kernels` and `l2_norm` are unchanged.

This task comes before the matrix task because Python 3.14 (Task 4) needs a pyo3
with 3.14 support; pyo3 0.22 predates it.

**Steps:**

- [ ] Replace `native/qscat-kernels/Cargo.toml` `[package]` and `[dependencies]` with:

  ```toml
  [package]
  name = "qscat-kernels"
  version = "0.0.0"
  edition = "2024"
  rust-version = "1.85"
  description = "QSCAT compiled kernels (Rust/PyO3). Currently a stub: l2_norm exists to validate the maturin/PyO3 toolchain end-to-end; real kernels arrive via the qm-method-lifecycle optimize stage."
  license = "BSD-3-Clause"
  repository = "https://github.com/VanaMartin/qscat"

  [lib]
  name = "qscat_kernels"
  crate-type = ["cdylib"]

  [dependencies]
  pyo3 = { version = "0.26", features = ["extension-module"] }
  ```

  If `cargo build` under pyo3 0.26 reports API changes against `src/lib.rs`,
  fix them (the `#[pymodule] fn(m: &Bound<'_, PyModule>)` signature used there
  is already the modern form; expect zero or trivial changes). If 0.26 is no
  longer the latest at implementation time, use the latest release that
  supports Python 3.14 — the requirement is "current + 3.14-capable", not the
  literal number.

- [ ] In `native/qscat-kernels/pyproject.toml`, add under `[project]`:

  ```toml
  license = "BSD-3-Clause"
  ```

  (description already exists there; maturin takes the rest from Cargo.toml).

- [ ] Rebuild and differential-test:
  `uv run maturin develop --manifest-path native/qscat-kernels/Cargo.toml`
  then `uv run --no-sync pytest native/qscat-kernels/tests -q`.

- [ ] `cargo clippy --manifest-path native/qscat-kernels/Cargo.toml -- -D warnings`
  and `cargo fmt --check --manifest-path native/qscat-kernels/Cargo.toml`.

- [ ] Run the four Global Constraints checks.

## Task 3: Stop paying for the stub crate on every CI job (style-M3/N6, CI half) + ADR 0006

**Files:**
- Modify: `.github/workflows/ci.yml` (new `changes` job; conditional Rust in `lint` and `test`)
- Modify: `native/qscat-kernels/tests/test_l2_norm.py`
- Create: `docs/adr/0006-rust-kernels-stay-a-stub-until-a-proven-hot-path.md`

**Interfaces:** none.

**Mechanism.** `uv sync --all-packages` builds the crate because it is a
workspace member; `uv sync --all-packages --no-install-package qscat-kernels`
resolves the same lock but skips installing (and therefore building) the
extension, so no Rust toolchain is needed. The kernel's own test then has to
skip rather than die at collection when the module is absent.

**Steps:**

- [ ] In `native/qscat-kernels/tests/test_l2_norm.py`, replace the bare
  `import qscat_kernels` with:

  ```python
  import pytest

  # CI builds the kernel only when native/** changes (docs/adr/0006); on jobs
  # that skipped the build this module must skip, not fail collection.
  qscat_kernels = pytest.importorskip("qscat_kernels")
  ```

- [ ] Add a `changes` job to `ci.yml` (before `lint`):

  ```yaml
    changes:
      name: what changed
      runs-on: ubuntu-latest
      permissions:
        pull-requests: read
      outputs:
        native: ${{ steps.filter.outputs.native }}
      steps:
        - uses: actions/checkout@v5
        - uses: dorny/paths-filter@v3
          id: filter
          with:
            # On push events (main) the filter diffs against the previous push
            # to this base; on PRs, against the PR base.
            base: ${{ github.ref_name }}
            filters: |
              native:
                - 'native/**'
                - 'uv.lock'
                - '.github/workflows/ci.yml'
  ```

  (`uv.lock` and the workflow itself are in the filter so a dependency bump or
  a change to this very mechanism still exercises the build path.)

- [ ] In the `lint` and `test` jobs: add `needs: changes`, make the Rust
  toolchain step conditional, and pick the sync flags from the filter output:

  ```yaml
        - name: Install Rust toolchain (only when native/** changed — docs/adr/0006)
          if: needs.changes.outputs.native == 'true'
          uses: dtolnay/rust-toolchain@stable
  ```

  ```yaml
        - name: Sync workspace (kernel build skipped unless native/** changed — docs/adr/0006)
          run: uv sync --all-packages ${{ needs.changes.outputs.native != 'true' && '--no-install-package qscat-kernels' || '' }}
  ```

- [ ] Leave `validation.yml` untouched: slow-tier runs are minutes long and
  explicitly requested; a one-minute kernel build there is noise, and the slow
  decks may legitimately exercise the kernel.

- [ ] Validate the YAML with `yaml.safe_load`; verify locally that
  `uv sync --all-packages --no-install-package qscat-kernels` followed by
  `uv run --no-sync pytest -q -m "not slow" -n auto --dist loadfile` is green
  with `native/qscat-kernels/tests/test_l2_norm.py` reported as skipped, then
  `uv sync --all-packages` to restore the kernel.

- [ ] Write `docs/adr/0006-rust-kernels-stay-a-stub-until-a-proven-hot-path.md`:

  ```markdown
  # 6. The Rust kernel crate stays a stub until a proven hot path exists

  Date: 2026-08-25

  ## Status

  Accepted

  ## Context

  `native/qscat-kernels` exists to keep the PyO3/maturin toolchain proven
  end-to-end (build, import, differential test, Docker), but its only content
  is an `l2_norm` placeholder. Meanwhile every CI job paid a Rust toolchain
  install and a maturin build for it, and its metadata (no license, edition
  2021, pyo3 0.22) misrepresented the project. The measured hot path is the
  sparse LU factorization and its per-step triangular solves
  (docs/physics/optimization-targets.md), which are already served by the
  MUMPS backend — the profile explicitly found NO pure-Python hot loop that a
  first Rust kernel could win on the current direct-solver architecture.

  ## Decision

  Keep the crate in the workspace with honest metadata, and stop paying for
  it on unrelated CI runs: ci.yml builds the kernel only when `native/**`,
  `uv.lock`, or the workflow itself changes (dorny/paths-filter); otherwise
  `uv sync --all-packages --no-install-package qscat-kernels` skips the build
  and the kernel's tests skip via `importorskip`. Implementing a real kernel
  now — the qm-method-lifecycle stage-4 path, e.g. a sparse-LU-adjacent or
  propagation-inner-loop kernel — is explicitly deferred until a profile
  shows a hot path the MUMPS backend does not already cover.

  ## Consequences

  - CI jobs untouched by `native/**` skip the toolchain + build entirely.
  - The Docker images and local `uv sync --all-packages` still always build
    the kernel, so the toolchain never rots unexercised.
  - A future kernel starts from a crate whose metadata is already honest.
  ```

- [ ] Run the four Global Constraints checks. Push and confirm on the CI draft
  run that a docs-only commit skips the toolchain step while a `native/**`
  commit runs it (noted verification; two throwaway commits on the draft
  branch are the cheapest way to see both paths).

## Task 4: Widen the test matrix — macOS, Python 3.14, minimum floors (style-M5)

**Files:**
- Modify: `.github/workflows/ci.yml` (test job matrix + a new `min-versions` job)
- Modify: `libs/qscat/pyproject.toml` (3.14 classifier)

**Interfaces:** none.

**BLAS-pin rationale (required by the finding):** ci.yml:58-71 records that the
Linux pins exist because OpenBLAS starts one thread per core in each xdist
worker on the 4-vCPU runner (16 threads on 4 cores; the unpinned run was slower
than serial), and that this does NOT reproduce on macOS because scipy there
links Accelerate — pinning measured as a no-op locally (76.1 s pinned vs 74.6 s
unpinned at `-n 12`). Decision: keep the three env pins on ALL rows, macOS
included. On macOS they are a measured no-op, and if a future macOS scipy wheel
ever links OpenBLAS the pin is already in place. This keeps the workflow free
of per-OS env conditionals; extend the existing env comment with one line
saying the pins are a deliberate no-op on the macOS row.

**Steps:**

- [ ] Extend the test job (which after Tasks 1 and 3 already has the coverage
  include and `needs: changes`):

  ```yaml
    test:
      name: test (${{ matrix.os }}, py${{ matrix.python-version }})
      needs: changes
      runs-on: ${{ matrix.os }}
      strategy:
        fail-fast: false
        matrix:
          os: [ubuntu-latest]
          python-version: ["3.12", "3.13", "3.14"]
          include:
            - os: ubuntu-latest
              python-version: "3.12"
              coverage: "--cov=qscat --cov-report=term"
            # One macOS row: catches Accelerate-vs-OpenBLAS numeric divergence
            # and BSD-vs-GNU path assumptions. The BLAS pins below are a
            # measured no-op here (Accelerate; see the env comment) and are
            # kept for uniformity.
            - os: macos-latest
              python-version: "3.13"
  ```

  The checkout/uv/sync/test steps need no per-OS branches (uv and dtolnay's
  action are cross-platform).

- [ ] Add the minimum-versions job. With `--resolution lowest-direct` uv
  installs every *direct* dependency at its declared floor; from
  `libs/qscat/pyproject.toml` and `apps/qscat-run/pyproject.toml` those floors
  are (write them into the job comment so the run is self-explaining):
  numpy 2.0.0, scipy 1.14.0, mpmath 1.3.0, matplotlib 3.9.0 (plot extra) /
  3.0 (qscat-run — expect the resolver to reconcile upward), click 8.0,
  pyyaml 6.0, pytest 8.0, pytest-xdist 3.0, hypothesis 6.0.

  ```yaml
    min-versions:
      name: test (lowest direct deps)
      needs: changes
      runs-on: ubuntu-latest
      steps:
        - uses: actions/checkout@v5
        - name: Install Rust toolchain (only when native/** changed — docs/adr/0006)
          if: needs.changes.outputs.native == 'true'
          uses: dtolnay/rust-toolchain@stable
        - name: Install uv
          uses: astral-sh/setup-uv@v7
          with:
            python-version: "3.12"
            enable-cache: true
        # Floors under test (from the pyprojects): numpy 2.0.0, scipy 1.14.0,
        # mpmath 1.3.0, matplotlib 3.9.0, click 8.0, pyyaml 6.0. If this job
        # fails while the normal rows pass, either fix the code or RAISE the
        # floor -- never let the floor claim a version the suite cannot run on.
        - name: Sync at the declared floors
          run: uv sync --all-packages --resolution lowest-direct ${{ needs.changes.outputs.native != 'true' && '--no-install-package qscat-kernels' || '' }}
        - name: Test (fast suite at the floors)
          env:
            OMP_NUM_THREADS: "1"
            OPENBLAS_NUM_THREADS: "1"
            MKL_NUM_THREADS: "1"
          run: uv run --no-sync pytest -q -m "not slow" -n auto --dist loadfile
  ```

- [ ] Add `"Programming Language :: Python :: 3.14"` to the classifiers in
  `libs/qscat/pyproject.toml`, next to the existing 3.13 entry — but only in
  the same commit where the 3.14 CI row is green; the classifier is a claim
  the matrix must back.

- [ ] Validate with `yaml.safe_load`; push and watch the draft run (noted
  verification). Known risks to check in that run, with the planned response:
  Python 3.14 wheels for numpy/scipy (expected fine by 2026); the macOS row's
  MUMPS/ffmpeg-gated tests must SKIP (they are `@skipif`-guarded already —
  a fail there is a bug to fix, not to mask); `--resolution lowest-direct`
  may surface a floor that is simply wrong — raise that floor in the
  pyproject rather than pinning around it.

- [ ] Run the four Global Constraints checks.

## Task 5: CITATION.cff version gate (style-M4 remainder)

**Files:**
- Modify: `.github/workflows/ci.yml` (lint job)

**Interfaces:** none.

**Steps:**

- [ ] Append to the `lint` job (it already syncs the whole workspace, and
  pyyaml is present via qscat-run's dependency):

  ```yaml
        - name: CITATION.cff version matches qscat.__version__
          run: |
            uv run --no-sync python -c "import sys, yaml, qscat; v = yaml.safe_load(open('CITATION.cff'))['version']; sys.exit(0 if v == qscat.__version__ else f'CITATION.cff version {v!r} != qscat.__version__ {qscat.__version__!r} -- update CITATION.cff (both top-level version: and preferred-citation.version)')"
  ```

- [ ] Run the one-liner locally first and confirm it passes on the current tree
  (`CITATION.cff` says 0.1.0.dev0; `qscat.__version__` is 0.1.0.dev0), then
  temporarily edit one side to confirm it fails with the message, then revert.

- [ ] Note the deliberate scope: the gate checks the top-level `version:` only.
  `preferred-citation.version` is mentioned in the failure message so the
  fixer updates both, but asserting it too is deferred — when the article
  publishes, `preferred-citation` becomes the article (per the Global
  Constraints) and a hard assert on its `version` key would then be wrong.

- [ ] Validate YAML; run the four Global Constraints checks.

## Task 6: Make the mypy config say what is actually checked (style-N11)

**Files:**
- Modify: `pyproject.toml` (`[tool.mypy]`)

**Interfaces:** none.

**Steps:**

- [ ] Replace the bare `[tool.mypy]` block at the bottom of `pyproject.toml` with:

  ```toml
  [tool.mypy]
  python_version = "3.12"
  strict = true

  # Strict-clean scope is libs/qscat/qscat + apps/qscat-run/qscat_run -- the
  # two trees CI checks. projects/, validation/ and benchmarks/ are research
  # code with untyped helpers; running strict mypy there reports hundreds of
  # pre-existing errors that nothing gates on. These overrides make the config
  # say so, instead of `strict = true` overstating a repo-wide claim.
  [[tool.mypy.overrides]]
  module = ["projects.*", "validation.*", "benchmarks.*"]
  ignore_errors = true
  ```

- [ ] Confirm the CI-checked scope is unchanged:
  `uv run --no-sync mypy libs/qscat/qscat apps/qscat-run/qscat_run` still clean.

- [ ] Confirm the overrides do what they claim:
  `uv run --no-sync mypy projects/n2_resonance` now exits 0 (errors ignored)
  where it previously reported errors.

- [ ] Run the four Global Constraints checks.

## Task 7: Narrow `validate:n2`, add `validate:factory` (recent-m5)

**Files:**
- Modify: `.github/workflows/validation.yml` (SUITES map, header comment, advise-job label line)
- Modify: `CONTRIBUTING.md` (label list)

**Interfaces:** the SUITES dict is machine-read by `tests/test_validation_suites.py`
(regex + `ast.literal_eval`), which also enforces: every slow test covered, no
empty suite, every path exists, header labels == SUITES keys. All four gates
constrain this edit.

**Why a new suite:** `projects/potential_factory` holds slow tests
(`test_fit.py`, `test_roundtrip.py`, `test_target.py`, `test_tracker.py`) that
today are covered only by the over-broad `"n2": "validation/n2 projects"`
prefix. Narrowing n2 without a new home orphans them and
`test_every_slow_test_is_selected_by_some_suite` fails. Potential-factory work
is not N₂, so it gets its own label rather than squatting in someone else's.

**Steps:**

- [ ] In the selection script, replace the SUITES map with (one line per suite —
  the cover test's regex/literal_eval parse handles long lines, and one-per-line
  keeps the diff honest):

  ```python
          # suite -> the paths that suite owns
          SUITES = {
              "core":     "libs/qscat/tests",
              "n2":       "validation/n2 projects/n2_resonance projects/n2_ti_cross_section projects/n2_td_cross_section projects/n2_2d_cross_section projects/n2_2d_td_cross_section",
              "diatomic": "validation/diatomic",
              "factory":  "projects/potential_factory",
              "h2plus":   "validation/h2plus",
              "tuning":   "validation/tuning",
              "run":      "apps/qscat-run/tests",
          }
  ```

  (All five n2 project dirs are listed even though only the two 2-D ones hold
  slow tests today — the paths exist, and a slow test added to the 1-D
  projects tomorrow is then covered without another map edit.)

- [ ] Update the workflow header's label list — the cover test asserts it
  matches SUITES exactly (regex `^#\s+validate:(\w+)\s{2,}`), so add:

  ```
  #   validate:factory   the potential-factory project (fits, round trips, targets)
  ```

  and reword the `validate:n2` line to "the N2 sub-projects (validation/n2 +
  projects/n2_*)".

- [ ] Update the advise job's suggested-labels line to include `validate:factory`:

  ```
            echo '`validate:core` `validate:n2` `validate:diatomic` `validate:factory` `validate:h2plus` `validate:tuning` `validate:run` `validate:all`'
  ```

- [ ] Update CONTRIBUTING.md's label list ("Asking CI to run the production
  tier") to the same eight labels.

- [ ] Create the GitHub label so applying it is possible:
  `gh label create "validate:factory" --color 0E8A16 --description "Run the potential-factory slow suite"`.

- [ ] Verify: `uv run --no-sync pytest tests/test_validation_suites.py -q` green
  (this single command exercises the cover, non-empty, paths-exist,
  header-match, and selection-script gates).

- [ ] Validate YAML; on the draft PR, apply `validate:factory` and confirm the
  matrix spawns exactly the `validate:factory` job (noted verification).

- [ ] Run the four Global Constraints checks.

## Task 8: Community + security files (style-M7)

**Files:**
- Create: `.github/dependabot.yml`
- Create: `SECURITY.md`
- Create: `.github/ISSUE_TEMPLATE/bug.yml`
- Create: `.github/PULL_REQUEST_TEMPLATE.md`

**Interfaces:** none.

**Steps:**

- [ ] Write `.github/dependabot.yml`. The Python ecosystem entry uses `"uv"`
  (the repo's lockfile is `uv.lock`, which Dependabot's uv support reads;
  the `pip` ecosystem would see only the pyprojects and miss the lock). If the
  repository's Dependabot rejects the `uv` ecosystem at implementation time,
  fall back to `pip` with the same directory and note it in the file:

  ```yaml
  version: 2
  updates:
    - package-ecosystem: "uv"
      directory: "/"
      schedule:
        interval: "weekly"
      groups:
        python-deps:
          patterns: ["*"]
    - package-ecosystem: "github-actions"
      directory: "/"
      schedule:
        interval: "weekly"
    - package-ecosystem: "cargo"
      directory: "/native/qscat-kernels"
      schedule:
        interval: "weekly"
    # Keeps the digest-pinned base images (docker/base.Dockerfile,
    # docker/Dockerfile -- see the style-N7 task) from rotting.
    - package-ecosystem: "docker"
      directory: "/docker"
      schedule:
        interval: "weekly"
  ```

- [ ] Write `SECURITY.md` at the repo root:

  ```markdown
  # Security Policy

  ## Supported versions

  | Version | Supported |
  |---|---|
  | `main` (this repository) | yes |

  qscat is not on PyPI; the only distribution is this repository, and only
  the tip of `main` is supported. There are no maintained release branches.

  ## Reporting a vulnerability

  Report privately through GitHub's advisory form:
  <https://github.com/VanaMartin/qscat/security/advisories/new>
  (repository **Security** tab → "Report a vulnerability"). Please do not
  open a public issue for a suspected vulnerability. You can expect an
  acknowledgement within one week.

  Scope worth knowing: `qscat-run` executes YAML configuration files that
  select solvers and write artifacts to paths named in the config. Treat
  configs from untrusted sources as untrusted input.
  ```

- [ ] Write `.github/ISSUE_TEMPLATE/bug.yml`:

  ```yaml
  name: Bug report
  description: A reproducible defect in qscat, qscat-run, or the build.
  labels: ["bug"]
  body:
    - type: textarea
      id: repro
      attributes:
        label: Minimal reproducer
        description: >-
          The smallest input that shows the problem — a small grid config is
          ideal (see CONTRIBUTING.md). For numeric issues, include the number
          you got and the number you expected, with the tolerance you judge
          them by.
        placeholder: |
          uv run qscat-run run bug.yaml --output /tmp/bug
          # bug.yaml: ...
      validations:
        required: true
    - type: textarea
      id: expected
      attributes:
        label: Expected vs actual
      validations:
        required: true
    - type: input
      id: version
      attributes:
        label: Commit SHA
        description: "`git rev-parse HEAD` — qscat is distributed repo-only, so the SHA is the version."
      validations:
        required: true
    - type: input
      id: platform
      attributes:
        label: Platform
        placeholder: "macOS 15 arm64 / Ubuntu 24.04 x86-64 / Docker qmodeling-base"
      validations:
        required: true
  ```

- [ ] Write `.github/PULL_REQUEST_TEMPLATE.md`:

  ```markdown
  ## What

  ## Why

  ## Checks

  - [ ] `uv run pytest -m "not slow" -n auto --dist loadfile` green
  - [ ] `uv run ruff check .` and `uv run mypy libs/qscat/qscat apps/qscat-run/qscat_run` clean

  ## Validation tier

  Does this change touch calculation-bearing source (`libs/qscat/qscat/`,
  `native/`, `projects/`, `apps/qscat-run/qscat_run/`)?

  - **No** — say so; no label needed.
  - **Yes** — name the `validate:*` label that covers what changed and apply
    it, or state why the change cannot move a number (a pure refactor, a
    docstring, plumbing asserted on toy decks). The Validation workflow's
    advisory note will call out the paths either way; the judgement is yours,
    not the path filter's.
  ```

- [ ] Validate both YAML files with `yaml.safe_load`. Push; confirm on GitHub
  that the issue form renders and Dependabot accepts the config (Insights →
  Dependency graph → Dependabot shows the four ecosystems, or an error to fix)
  (noted verification).

- [ ] Run the four Global Constraints checks (`test_docs_portability.py` does
  not scan these files, but run the block as always).

## Task 9: Honest qscat-run metadata, repo-only stated (style-M6)

**Files:**
- Modify: `apps/qscat-run/pyproject.toml`
- Modify: `apps/qscat-run/README.md`

**Interfaces:** none (the `qscat-run` entry point and package name are unchanged).

**Steps:**

- [ ] In `apps/qscat-run/pyproject.toml`, replace the `[project]` header fields with:

  ```toml
  [project]
  name = "qscat-run"
  # Kept in step with qscat.__version__ by hand (the CLI is a thin surface
  # over the library and is versioned with it). Repo-only: this package is
  # deliberately NOT published to PyPI -- and will not be until the qscat
  # citation article is out; install it via `uv sync --all-packages` from a
  # clone.
  version = "0.1.0.dev0"
  description = "qscat-run — a config-driven CLI for 2-D electron-diatomic model experiments"
  readme = "README.md"
  requires-python = ">=3.12"
  license = "BSD-3-Clause"
  authors = [{ name = "Martin Vana", email = "martin@qscat.com" }]
  ```

  (dependencies, scripts, sources, build-system unchanged. No `license-files`:
  the LICENSE lives at `libs/qscat/LICENSE`; since this package is never built
  for distribution, referencing it from here is not needed and hatchling would
  reject a path outside the project root.)

- [ ] Add to `apps/qscat-run/README.md`, directly under the opening bold line:

  ```markdown
  qscat-run is **repo-only**: it is not published to PyPI (and will not be
  until the qscat citation article is out). Install it from a clone with
  `uv sync --all-packages`.
  ```

- [ ] Confirm `uv sync --all-packages` still installs it and
  `uv run qscat-run list` still runs.

- [ ] Run the four Global Constraints checks.

## Task 10: Digest-pin the Docker base images (style-N7)

**Files:**
- Modify: `docker/base.Dockerfile` (line 5's `FROM`)
- Modify: `docker/Dockerfile` (the runtime stage's `FROM`)

**Interfaces:** none. The two other base variants (`base-cpu-mkl.Dockerfile`,
`base-gpu.Dockerfile`) inherit different images and are optional/experimental;
pin them the same way only if they build cleanly today — do not block this task
on them.

**Steps:**

- [ ] Resolve the current digest of the shared base image (the multi-arch
  manifest-list digest, so the same line works on the arm Mac and the x86
  `sadaharu` host):

  ```bash
  docker buildx imagetools inspect ghcr.io/astral-sh/uv:python3.12-bookworm-slim \
    | head -5   # the "Digest: sha256:..." line of the manifest list
  ```

- [ ] In `docker/base.Dockerfile` replace line 5 with (using the digest actually
  resolved above — the `sha256:...` below is a placeholder for the resolved
  value, the ONE spot in this plan where a value cannot be known until run time):

  ```dockerfile
  # Tag kept for humans, digest is what builds: `docker buildx imagetools
  # inspect ghcr.io/astral-sh/uv:python3.12-bookworm-slim` to refresh, and
  # dependabot's docker ecosystem watches it.
  FROM ghcr.io/astral-sh/uv:python3.12-bookworm-slim@sha256:<resolved-digest>
  ```

- [ ] Same edit for `docker/Dockerfile`'s runtime stage
  (`FROM ghcr.io/astral-sh/uv:python3.12-bookworm-slim AS runtime` gains the
  same `@sha256:` suffix, same comment). Both files must carry the SAME digest —
  they intend the same image.

- [ ] Verify the pinned build still works: `docker/build.sh test-deps` (the
  cheapest target that exercises the base image; a full `test` run is not
  required for a FROM-line change, and the slow tier inside it costs ~13 min).

- [ ] Run the four Global Constraints checks.

---

## Explicitly-deferred alternative (style-M3/N6): implement the real kernel now

The alternative to Task 2/3 is to give `native/qscat-kernels` a real kernel
immediately — the qm-method-lifecycle stage-4 path: profile, pick the hot loop,
mirror the validated Python API, differential-test against it. It is deferred
because the profile already answered the question: ~98% of TI cost is the sparse
factorization and ~82% of TD cost is the per-step solve, both inside
SuperLU/MUMPS, and `docs/physics/optimization-targets.md` records that the
pure-Python loops a first kernel would target are ~0.1% of runtime. The next
real kernel opportunity appears only with an iterative-solver or
very-large-scale-assembly path, neither of which exists yet. No tasks; ADR 0006
(Task 3) records this so the decision is in the clone, not only in this plan.

## Final verification

- [ ] All four Global Constraints checks green on the finished branch.
- [ ] One pushed draft-PR run showing: coverage table in the step summary; the
  matrix rows ubuntu {3.12, 3.13, 3.14} + macos 3.13 + min-versions all green;
  the Rust toolchain step skipped on a docs-only commit; the CITATION gate
  green; `validate:factory` selectable.
- [ ] `grep -rn "TBD\|TODO\|placeholder" docs/superpowers/plans/2026-08-25-phase2-infrastructure.md`
  finds only the one documented `<resolved-digest>` runtime value (Task 10).
