# Docs Structure Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the physics notes readable results-first and structurally
consistent: split the ~387-line vibrational-excitation section out of
`nonlocal-resonance-model.md` into its own note, retrofit results-first openings
onto four buried-headline notes, record the heading scheme and the
backticked-maths policy in the conventions skill, sweep the offending notes, and
extend the portability test with the one zero-false-positive detector that
exists.

**Architecture:** Documentation + one test file + one skill file. No solver
code, no numbers change anywhere — every content move is checked by a
numeric-token multiset assertion (every numeric token of the original text must
survive, verbatim, into the successor text). Policy/skill edits land first, the
mechanical sweeps second, the enforcement detector third, and the two big
content moves (M10 rewrites, M11 split) last, so the tree passes the four
checks after every task.

**Tech Stack:** Markdown/MyST, `tests/test_docs_portability.py` (pytest),
Sphinx (`sphinx-build -W`), python one-liners for the numeric checks.

**Spec:** the "Findings addressed" section below (self-contained; from the 2026-08-25 release review)

## Global Constraints

- **PyPI release DEFERRED until the peer-reviewed citation article publishes — distribution is repo-only; NO task may register publishers, tag releases, or claim pip-installability.** When the article is out, CITATION.cff gains it as preferred-citation (out of scope here).
- After every task: `uv run --no-sync pytest -m "not slow" -n auto --dist loadfile` green; `uv run --no-sync mypy libs/qscat/qscat apps/qscat-run/qscat_run` clean; `uv run --no-sync ruff check .` + `ruff format --check .` clean; `uv run --no-sync pytest tests/test_docs_portability.py -q` green.
- Workflow YAML edits validated with `yaml.safe_load` + a pushed draft run where noted (this plan touches no workflow YAML; the clause is carried for uniformity).
- Never `git commit -a`.
- **Numbers never change.** Every task that moves or rewrites note prose runs the numeric-token multiset check given in that task and pastes its output into the task's completion note. The check asserts the BEFORE token multiset is a sub-multiset of the AFTER one (openings and stubs may *repeat* a number; nothing may lose or alter one).
- After any task that adds/removes a docs page or edits a toctree:
  `uv sync --package qscat --group docs && uv run sphinx-build -b html -W --keep-going docs docs/_build/html` clean, then `uv sync --all-packages` to restore the full workspace before the pytest run.

## Findings addressed

- **docs-M11:** `nonlocal-resonance-model.md` §8 (Vibrational excitation, lines 635–1021, subsections 8.1–8.10) is a second note trapped inside a first: DA readers scroll past ~387 lines of VE, and the VE result (the note's strongest) is unfindable from the section indexes. Split it into `docs/physics/nrm-vibrational-excitation.md`, filed under the ENGINE section index (`docs/physics/engine.md`), leaving a pointer paragraph; update `dissociation.md`'s blurb, the toctree, and every inbound link; byte-identical numbers.
- **docs-M10:** four notes bury their headline: `discretisation-tuning.md` ("Key result" at line 261 of 283), `mumps-sparse-backend.md` (the production speedup mid-note, behind two mechanism sections), `h2plus-dr.md` (the delivered σ_DR curve at line 334, behind a 201-line discrepancy narrative), `optimization-targets.md` (the "MUMPS is the correct default" conclusion mid-note). Retrofit a results-first opening on each; prose bodies move, numbers never change.
- **docs-M9:** no shared heading scheme across notes. Pick Key result → Physical picture → Method → Validation, record it in the `qscat-conventions` skill as the rule for NEW notes with opportunistic migration; no big-bang re-heading. Already conforming (per the review): n2-resonance, n2-cross-section, n2-td-cross-section, nd-tensor-hamiltonian, femdvr-ecs, exact-2d-resonances.
- **recent-M6:** `validation-harnesses.md`'s "Where the @slow boundary falls" predates ADR 0005 and the `validate:*` workflow; it explains one harness group's budget but not the actual boundary rule. Rewrite it to cite ADR 0005's cost rule and the label workflow.
- **docs-M8 remainder:** backticked unicode maths (`` `Γ(R)` ``, `` `μ = 918.25` ``) blurs the code-vs-maths line the conventions skill draws. Decide the policy (decided below: outlaw it — a symbol is either code or maths), sweep the offenders, and extend `test_docs_portability.py` only if a zero-false-positive detector exists (evaluated below: one does, for the Greek-letter subset).
- **docs-N16/N19 (riders):** the conventions skill's subpackage map is missing `core`/`model`/`tuning`/`viz`; its tolerance guidance is missing the cross-arch BLAS floor (sparse-solve goldens at 1e-9; 1e-12 failed in CI).

### Decided policies (referenced by the tasks)

**Heading scheme (M9):** a note opens with the header block (Location / Origin
/ Units), then **Key result** (5–10 lines, the measured numbers first), then
**Physical picture**, **Method**, **Validation** as the top-level spine. New
notes must follow it; existing notes migrate opportunistically.

**Backticked-maths policy (M8):** outlawed. A symbol is either a code
identifier (`R0` the model attribute, backticks) or mathematics ($\Gamma(R)$,
dollars) — never Greek/maths dressed in backticks. The existing skill
exception for short level labels (`v = 0`, `Ry₄`) stands: those contain no
Greek and survive table cells without the `\vert` hazard.

**Detector (M8c):** enforce the **Greek-letter subset only**. A backtick span
containing a Greek letter (U+0370–U+03FF, U+1F00–U+1FFF) is always
maths-in-backticks in this repository: Python/Rust identifiers here are ASCII,
so no genuine code span can contain one. The broader review-time heuristic
(any superscript/subscript unicode in a span) is NOT enforceable: it flags
molecule names in real paths (`reference/eMoScat/input/{NO,F₂}/grids.txt`),
reaction equations, and the skill's own level-label exception (`Ry₄`). Those
remain convention-by-review. Re-verified counts with the Greek detector
(2026-08-25, this worktree — they differ from the review's counts, which used
the broader heuristic, and three notes have been added since the review):

| note | Greek-in-backtick spans |
|---|---|
| potential-factory.md | 104 |
| nrm-time-dependent.md | 87 |
| potential-factory-options.md | 67 |
| h2plus-resonance-states.md | 40 |
| exact-2d-resonances.md | 14 |
| nonlocal-resonance-model.md | 14 |
| diatomic-ve-cross-sections.md | 10 |
| h2plus-dr.md | 7 |
| n2-resonance.md | 4 |
| **total** | **347** |

---

## Task 1: Record the heading scheme in the conventions skill (docs-M9)

**Files:**
- Modify: `.claude/skills/qscat-conventions/SKILL.md`

**Steps:**

- [ ] Add a new section after "Mathematics in Documentation":

  ```markdown
  ## Physics-note structure

  A note under `docs/physics/` opens with the standard header block
  (**Location** / **Origin** or **Source** / **Units**), then follows this
  top-level spine:

  1. **Key result** — 5–10 lines, the measured numbers and the one-sentence
     claim first. A reader who stops here leaves with the result and where
     it was measured.
  2. **Physical picture** — what the method is and why it exists here.
  3. **Method** — the equations, conventions, and implementation choices.
  4. **Validation** — what was measured against what, with the evidence
     (and the negative results; the notes record limitations too).

  This is the rule for NEW notes. Existing notes migrate opportunistically —
  when a note is being edited anyway — never as a big-bang re-heading.
  Already conforming: `n2-resonance`, `n2-cross-section`,
  `n2-td-cross-section`, `nd-tensor-hamiltonian`, `femdvr-ecs`,
  `exact-2d-resonances`; results-first retrofits (Key result prepended,
  bodies untouched): `discretisation-tuning`, `mumps-sparse-backend`,
  `h2plus-dr`, `optimization-targets`, and the split-out
  `nrm-vibrational-excitation`.
  ```

  (The retrofit list names notes this plan's Tasks 7–11 produce; land this
  task first anyway — the list documents intent and is corrected in Task 12's
  final pass if any of those tasks changes shape.)

- [ ] Run the Global Constraints checks.

## Task 2: Record the backticked-maths policy + the two rider fixes (docs-M8a, N16, N19)

**Files:**
- Modify: `.claude/skills/qscat-conventions/SKILL.md`

**Steps:**

- [ ] In "Mathematics in Documentation", extend the first bullet ("Backticks
  are for code identifiers only…") to name the unicode form of the defect
  explicitly:

  ```markdown
  - **Backticks are for code identifiers only** — `ve_cross_section`,
    `SparseLU`, `backend="mumps"`. A backticked `sigma` that means σ is a
    defect; so is a backticked Greek letter — `` `Γ(R)` `` and
    `` `μ = 918.25` `` are maths and are written `$\Gamma(R)$`,
    `$\mu = 918.25$`. A symbol is either code or maths, never Greek dressed
    in backticks. (Enforced: `tests/test_docs_portability.py` flags any
    backtick span containing a Greek letter — zero false positives, because
    this repository's identifiers are ASCII.)
  ```

- [ ] Rewrite the closing sentence of "What is enforced, and what is
  convention-by-review" (currently: "Backticked unicode that is really maths …
  is deliberately **not** detected — the false-positive rate … is too high to
  gate on"), which the new detector makes wrong. Replacement:

  ```markdown
  Backticked maths is detected for the **Greek-letter subset** only
  (`find_greek_in_backticks`): a Greek letter inside a backtick span is
  always maths here, since the repository's identifiers are ASCII. The
  broader family — superscript/subscript unicode in a span — stays
  convention-by-review, because it collides with legitimate spans: molecule
  names in real paths (`reference/eMoScat/input/{NO,F₂}/grids.txt`),
  reaction equations, and the level-label exception above (`Ry₄`).
  ```

  and add the row to the enforcement table:

  ```markdown
  | no Greek letters inside backtick spans | `find_greek_in_backticks` |
  ```

  (The detector itself lands in Task 6, AFTER the sweeps — the skill may
  briefly describe enforcement that is one task away; the alternative,
  landing a red test first, violates the Global Constraints.)

- [ ] N16 — extend the subpackage map table with the four missing rows:

  ```markdown
  | `core`      | Model-independent electron–diatomic scattering engine (driven/TI, TD, DA/DR, LCP, NRM, resonance verification) |
  | `model`     | The `ResonanceModel` protocol + the per-molecule registry (`N2`, `NO`, `F2`, `H2P`) |
  | `tuning`    | Automatic FEM-DVR-ECS discretisation tuner                |
  | `viz`       | Wavefunction rendering and animation (plot extra)         |
  ```

- [ ] N19 — add to "Tolerance Defaults":

  ```markdown
  - **Cross-arch BLAS floor:** a golden value produced *through a sparse
    solve* is compared at `rtol=1e-9`, never `1e-12`. The identical solve on
    a different BLAS (CI's Linux OpenBLAS vs a Mac's Accelerate) differs at
    ~1e-10–1e-11, and a `1e-12` gate on such a golden has already failed in
    CI for exactly this reason. The `1e-12` band above is for two
    implementations of the *same deterministic arithmetic* on the same
    machine, not for cross-architecture goldens.
  ```

- [ ] Run the Global Constraints checks.

## Task 3: Sweep the six review-flagged notes (docs-M8b, part 1)

**Files:**
- Modify: `docs/physics/n2-resonance.md` (4 spans), `docs/physics/h2plus-dr.md` (7),
  `docs/physics/diatomic-ve-cross-sections.md` (10), `docs/physics/exact-2d-resonances.md` (14),
  `docs/physics/nonlocal-resonance-model.md` (14), `docs/physics/h2plus-resonance-states.md` (40)

**Steps:**

- [ ] Generate the exact worklist per note with the detector prototype (also
  the acceptance check — rerun after the sweep and expect zero):

  ```bash
  python3 - <<'PY'
  import re, pathlib
  fence = re.compile(r'^```.*?^```', re.DOTALL | re.MULTILINE)
  span, greek = re.compile(r'`([^`\n]+)`'), re.compile(r'[Ͱ-Ͽἀ-῿]')
  for name in ["n2-resonance", "h2plus-dr", "diatomic-ve-cross-sections",
               "exact-2d-resonances", "nonlocal-resonance-model",
               "h2plus-resonance-states"]:
      t = fence.sub("", pathlib.Path(f"docs/physics/{name}.md").read_text())
      hits = [m.group(1) for m in span.finditer(t) if greek.search(m.group(1))]
      print(f"{name}: {len(hits)}")
      for h in hits: print("   ", h)
  PY
  ```

- [ ] Convert each span, by category (worked examples — these exact spans are
  in the worklist):
  - pure maths → dollars: `` `Γ(R)` `` → `$\Gamma(R)$`; `` `μ = 918.25` `` →
    `$\mu = 918.25$`; `` `ω_i^j` `` → `$\omega_i^j$`; `` `φ⊗χ` `` →
    `$\phi \otimes \chi$`.
  - mixed code + maths → split at the boundary:
    `` `c_product(φ_d(R_j), φ_d(R_{j+1}))` `` → "`c_product` of
    $\phi_d(R_j)$ with $\phi_d(R_{j+1})$";
    `` `ε_e(R_∞) − eps[0]` `` → "$\varepsilon_e(R_\infty)$ minus `eps[0]`".
  - formula transcriptions →  dollars with the unicode operators promoted to
    LaTeX: `` `G_l(0,ρ)=−ρ y_l(ρ)` `` → `$G_l(0,\rho) = -\rho\, y_l(\rho)$`.
  - Rules that bind during conversion: inside a Markdown table cell any
    `\lvert`/`\rvert`/`\vert`, never a bare `|`; NEVER move a `$...$` into a
    heading (the portability test fails the build); a span inside a fenced
    code block is out of scope (the detector strips fences).

- [ ] Rerun the detector script: all six notes report 0.

- [ ] Numbers-preserved check per note (paste the output into the completion note):

  ```bash
  for n in n2-resonance h2plus-dr diatomic-ve-cross-sections exact-2d-resonances \
           nonlocal-resonance-model h2plus-resonance-states; do
    git show HEAD:docs/physics/$n.md > /tmp/before-$n.md
    python3 - "$n" <<'PY'
  import re, sys, collections, pathlib
  n = sys.argv[1]
  tok = lambda p: collections.Counter(re.findall(r'\d+(?:\.\d+)?(?:[eE][+-]?\d+)?', pathlib.Path(p).read_text()))
  lost = tok(f"/tmp/before-{n}.md") - tok(f"docs/physics/{n}.md")
  assert not lost, f"{n}: numeric tokens lost: {dict(lost)}"
  print(f"{n}: no numeric token lost")
  PY
  done
  ```

- [ ] Eyeball the built pages for the two heaviest notes
  (`h2plus-resonance-states`, `exact-2d-resonances`) in
  `docs/_build/html/physics/` after the sphinx build — MathJax rendering of
  ~50 new inline spans is where a stray unescaped `_` or `^` would show.

- [ ] Run the Global Constraints checks (including the sphinx build, per the
  toctree clause — no toctree changed, but the maths must build warning-free).

## Task 4: Sweep the three post-review notes (docs-M8b, part 2)

**Files:**
- Modify: `docs/physics/nrm-time-dependent.md` (87 spans),
  `docs/physics/potential-factory.md` (104),
  `docs/physics/potential-factory-options.md` (67)

These three notes postdate the review (potential-factory landed 2026-08-24)
and carry more offenders than the six reviewed notes combined. Without this
task the Task 6 detector cannot be enforced.

**Steps:**

- [ ] Same procedure as Task 3 (detector worklist → categorized conversion →
  detector reports 0 → per-note numeric multiset check → sphinx build). The
  dominant patterns here are wavefunction/curve symbols (`Ψ_d(R;E)` →
  `$\Psi_d(R;E)$`, `λ(R)` → `$\lambda(R)$`) and Morse-form transcriptions
  (`` `v0(R) = D_e [(1 − e^{−β(R)(R − R_e)})² − 1]` `` →
  `$v_0(R) = D_e[(1 - e^{-\beta(R)(R - R_e)})^2 - 1]$` — note `v0` the code
  attribute stays backticked where the prose means the callable, and becomes
  $v_0$ only inside a formula).
- [ ] These notes are long; convert in 2–3 commits per note if that keeps
  review sane, but the task is done only when the detector reports 0 for all
  three and the numeric checks pass.
- [ ] Run the Global Constraints checks.

## Task 5: Rewrite "Where the @slow boundary falls" (recent-M6)

**Files:**
- Modify: `docs/physics/validation-harnesses.md` (the `## Where the @slow boundary falls` section only)

**Steps:**

- [ ] Replace the section body with the following (the F1/T-scan example and
  every number are retained — the numeric multiset check applies):

  ```markdown
  ## Where the `@slow` boundary falls

  The boundary is ADR 0005's cost rule, not a judgement of importance: a test
  needing more than a few seconds or ~0.5 GB belongs in the `slow` tier — or
  wants a smaller deck, which is often the better answer
  (`docs/adr/0005-test-tiers-fast-and-slow.md`, points 2 and 7). The default
  tier is toy-scale and is the CI gate on every push; the `slow` tier is
  production-scale physics and runs in three places: locally
  (`uv run pytest -m slow`, serial — the decks are sized in gigabytes), in
  the Docker `test` image, and on demand in CI, where a reviewer applies the
  `validate:*` label that covers the change and
  `.github/workflows/validation.yml` runs that suite (or a manual
  `workflow_dispatch` names it). The label is a human judgement — a path
  filter cannot tell whether a change can move a number — and the workflow's
  advisory note exists only to point out when that judgement has not been
  taken.

  The harness groups in `validation/n2/experiment.py` sit under the same cost
  rule with a per-group budget of roughly 60 seconds. A full time-dependent
  2-D propagation at `TD_WORKING_GRID` costs `~210-250s`
  (`validation/n2/experiment.py:189`) — even the *shortest* configuration on
  the sub-project's own T-scan (T=600) costs `~85s`, over budget, and is also
  the least-converged point that scan measured (sigma_TD/sigma_TI = 0.760
  there vs 0.931 at the converged T=1500). So Group F1
  (`validation/n2/td_exact2d.py`) does not run a TD propagation at all: it
  reports the already-validated $\sigma_\mathrm{TD}$ as a cited, literal
  constant, recomputing only the cheap $\sigma_\mathrm{TI}$ side live every
  harness run. The genuine PASS/FAIL gate on the TD-vs-TI agreement lives in
  the `slow` tier where ADR 0005 puts costs of this size:
  `projects/n2_2d_td_cross_section/test_td_cross_section.py`'s
  `@pytest.mark.slow` tests, covered by the `validate:n2` label. The same
  pattern — a cheap, always-live check in the harness paired with an
  expensive, opt-in `@slow` pytest gate elsewhere — recurs in
  `validation/tuning`, whose 2-D spot-checks and convergence tests
  (`test_emoscat_decks.py`, `test_resonance_aware.py`) are `@pytest.mark.slow`
  and run under the `validate:tuning` label.
  ```

- [ ] Numeric multiset check on `validation-harnesses.md` (same one-liner
  pattern as Task 3, single file).

- [ ] Cross-check the label names against `.github/workflows/validation.yml`'s
  SUITES map AS OF THE IMPLEMENTATION DATE — the parallel
  2026-08-25-phase2-infrastructure plan renames nothing but ADDS
  `validate:factory`; this section names only `validate:n2`/`validate:tuning`,
  which exist under both states, so no coupling — verify and move on.

- [ ] Run the Global Constraints checks.

## Task 6: The Greek-in-backticks detector (docs-M8c)

**Files:**
- Modify: `tests/test_docs_portability.py`

**Prototype and its false-positive test, evaluated honestly:** the regex pair

```python
_BACKTICK_SPAN_RE = re.compile(r"`([^`\n]+)`")
_GREEK_RE = re.compile(r"[Ͱ-Ͽἀ-῿]")
```

applied to fence-stripped text flags exactly 347 spans across 9 notes on the
pre-sweep tree (table above) and, hand-audited on this worktree, every one of
them is maths — no identifier, path, or literal in this repository contains a
Greek letter (Python/Rust source is ASCII; verify once during implementation
with `grep -rlP '[\x{0370}-\x{03ff}]' libs/qscat/qscat apps/qscat-run/qscat_run --include='*.py' | xargs -I{} sh -c 'echo {}'` and
confirm any hits are docstrings/comments, not identifiers a note would
backtick). What the detector deliberately does NOT flag — and why gating on it
would false-positive: superscript/subscript-only spans (`Ry₄` is the skill's
level-label exception; `F₂` appears inside the real path
`reference/eMoScat/input/{NO,F₂}/grids.txt`; `e⁻ + H₂⁺(v) → H + H` is a
reaction equation). Those stay convention-by-review.

**Steps:**

- [ ] Add to `tests/test_docs_portability.py`, following the existing
  detector-function pattern:

  ```python
  _BACKTICK_SPAN_RE = re.compile(r"`([^`\n]+)`")
  _GREEK_RE = re.compile(r"[Ͱ-Ͽἀ-῿]")


  def find_greek_in_backticks(text: str) -> list[str]:
      """Return every backtick span containing a Greek letter.

      A Greek letter inside backticks is always maths dressed as code in this
      repository -- identifiers are ASCII -- so this subset gates at zero
      false positives. Superscript/subscript unicode is deliberately NOT
      flagged: it appears in legitimate spans (level labels like ``Ry₄``,
      molecule names inside real paths) and stays convention-by-review.
      """
      return [
          m.group(1)
          for m in _BACKTICK_SPAN_RE.finditer(_strip_fences(text))
          if _GREEK_RE.search(m.group(1))
      ]
  ```

- [ ] Detector unit tests, same style as the existing ones:

  ```python
  def test_find_greek_in_backticks_flags_a_backticked_gamma():
      assert find_greek_in_backticks("the width `Γ(R)` is frozen\n") == ["Γ(R)"]


  def test_find_greek_in_backticks_passes_code_and_labels():
      text = "the attribute `R0`, the level `Ry₄`, and `$\\Gamma(R)$` prose\n"
      assert find_greek_in_backticks(text) == []


  def test_find_greek_in_backticks_ignores_fenced_blocks():
      text = "```python\n# σ_DA printed here\n```\n"
      assert find_greek_in_backticks(text) == []
  ```

- [ ] The tree-scan test, parametrized over `_notes()` like the others:

  ```python
  @pytest.mark.parametrize("note", _notes(), ids=lambda p: p.name)
  def test_note_has_no_greek_in_backticks(note: Path):
      found = find_greek_in_backticks(note.read_text())
      assert not found, (
          f"{note.name} backticks maths {found}. A symbol is either code "
          f"(`R0`) or maths ($\\Gamma(R)$) -- never Greek in backticks. "
          f"See the qscat-conventions skill, Mathematics in Documentation."
      )
  ```

- [ ] `uv run --no-sync pytest tests/test_docs_portability.py -q` — green only
  because Tasks 3–4 already swept every note; if it is red here, a sweep
  missed a span: fix the note, not the detector.

- [ ] Run the Global Constraints checks.

## Task 7: Results-first — `discretisation-tuning.md` (docs-M10)

**Files:**
- Modify: `docs/physics/discretisation-tuning.md`

**Reorder outline:** the existing `## Key result` (line 261) MOVES to directly
after the header block, lightly reworded to stand as an opening (below); the
trailing "See also" paragraph stays at the end of the file. Resulting H2 order:
Key result → Why → The hybrid approach → Task 8: calibrating C → the two
"Genuine finding" sections → The gate → See also. Body prose is otherwise
untouched.

**Steps:**

- [ ] Insert as the first section (this replaces the current line-261 section,
  which is deleted where it stands — it is a move, not a duplication):

  ```markdown
  ## Key result

  Calibrated once against F₂'s genuinely-open dissociative-attachment channel
  ($C = 0.10$), the tuner reproduces-and-beats the hand-tuned eMoScat F₂ deck
  on the 1-D probes (37% fewer points, clean absolute convergence) and
  correctly flags the coarse shared N₂-style grid's historical
  under-resolution of the same K≈58–78 wave — the two things this
  sub-project set out to prove — and the same clean result holds for H₂⁺'s
  proxy deck. The resonance-aware `channel="dissociation"` nuclear path
  converges F₂'s 2-D $\sigma_\mathrm{DA}$ on the FIRST a-priori pass
  (1.6562 bohr², matching the eMoScat deck and finding #3's own refine²
  value) at deck-parity size (1.027×), and sizes H₂⁺'s resonant grid ~4%
  under its proxy deck. Two honest caveats carry the same weight as the
  result: reaching convergence costs approximately deck-sized resolution —
  the deliverable is convergence + automation at deck-competitive size, not
  the "10–20% smaller" grid originally hoped for — and N₂/NO's proposed
  nuclear grids cost more points than their decks (root-caused to a fixed
  real-region extent default, not to $C$), a documented limitation for a
  follow-on.
  ```

- [ ] Numeric multiset check on the file (Task 3's one-liner pattern).
- [ ] Run the Global Constraints checks (sphinx build included — heading moves
  change anchors; `-W` catches a broken internal reference if one exists).

## Task 8: Results-first — `mumps-sparse-backend.md` (docs-M10)

**Files:**
- Modify: `docs/physics/mumps-sparse-backend.md`

**Reorder outline:** a NEW `## Key result` section is inserted directly after
the header block; every existing section keeps its current order (Why the
sparse LU is the hot path → Why every matrix is complex-symmetric → The MUMPS
SYM=2 backend → The benchmark → The mechanism → Honest caveats → Provisioning
→ Lifecycle → Reproduce → See also). Nothing else moves — the note's body
reads bottom-up fine once the headline exists.

**Steps:**

- [ ] Insert:

  ```markdown
  ## Key result

  On the real complex-symmetric N₂ 2-D matrices, MUMPS `SYM=2` beats SuperLU
  at every size, and by more the larger the deck: at the 143k-unknown
  production deck the factorization is **81.3× faster (3.2 s vs 258 s)** and
  peak RSS is **11.9× smaller (0.6 GB vs 7.4 GB)**; the 47k TD deck measures
  23.4× / 6.6×, the 27k working deck 12.0× / 3.3×. Residuals agree across
  backends to ~1e-12 — the two engines compute the same solution, only the
  cost differs — which puts a multi-hundred-point energy sweep comfortably
  inside the "under an hour" bar. `SparseLU(backend="auto")` picks MUMPS
  whenever it is importable; SuperLU stays the fallback and the differential
  oracle. (Earlier citations of 72.6× / 9.2× are the pre-fix `SYM=0` numbers
  — see the historical note under The benchmark.)
  ```

  The parenthetical matters: `CLAUDE.md` and `optimization-targets.md` still
  quote 72×/9×, and a reader reconciling the two must land on the historical
  note, not conclude the numbers disagree. (Updating those two citations to
  81.3×/11.9× is a worthwhile rider IF touched anyway; not required here.)

- [ ] Numeric multiset check; Global Constraints checks.

## Task 9: Results-first — `h2plus-dr.md` (docs-M10)

**Files:**
- Modify: `docs/physics/h2plus-dr.md`

**Reorder outline:** two changes. (1) A NEW `## Key result` after the header
block. (2) The section `## The converged full-size σ_DR(E) curve (delivered)`
(currently line 334, AFTER the 201-line "Against a previously obtained
reference sweep" narrative) MOVES up to directly after
"## Discretization (…)" + its cost subsection — so the resulting order is:
Key result → What this is → Coulomb special functions → The ionic model → The
DR cross section → Discretization (+ cost) → **The converged full-size σ_DR(E)
curve (delivered)** → Against a previously obtained reference sweep (the whole
discrepancy narrative, all four subsections, order unchanged) → The resonance
positions behind these peaks → Follow-ons. The discrepancy narrative loses no
prose — it just stops standing between the reader and the delivered curve.

**Steps:**

- [ ] Insert:

  ```markdown
  ## Key result

  The exact-2D TI $\sigma_\mathrm{DR}(E)$ for e⁻ + H₂⁺ is delivered at full
  deck size — ~1.15 M unknowns on the 1300-bohr Coulomb electronic grid —
  through `apps/qscat-run` under Docker/MUMPS at ~8 s/energy: the DR1 (n=0)
  channel peaks at E ≈ 6.31×10⁻³ Ha, σ ≈ 1.54×10⁻³ bohr² above a ~10⁻¹⁰
  background; DR2 (n=1) is ~10⁻⁶; DR3 (n=2) is closed in the window
  (threshold ≈ 0.0426 Ha). Against a previously obtained reference sweep,
  DR₀ agrees once a measured 2π normalization convention is accounted for
  (geometric mean 1.001), while DR₁ carries a systematic ~1.3× deficit that
  remains an OPEN, documented discrepancy — the narrative below records what
  was ruled out, and why the Born-Oppenheimer levels are where the agreement
  is genuinely readable.
  ```

- [ ] Move the delivered-curve section as outlined; verify the committed
  figure reference
  (`figures/h2plus-dr-cross-section-shortrange.png`) still resolves from its
  new position (relative path unchanged — it will).
- [ ] Numeric multiset check; Global Constraints checks.

## Task 10: Results-first — `optimization-targets.md` (docs-M10)

**Files:**
- Modify: `docs/physics/optimization-targets.md`

**Reorder outline:** a NEW `## Key result` after the two-line preamble;
existing sections keep their order (Measured hot path TI → Measured hot path
TD → Measured: does MUMPS fix it → Ranked plan → Thread oversubscription →
its two subsections).

**Steps:**

- [ ] Insert:

  ```markdown
  ## Key result

  Three measured findings set the agenda. **(1) The sparse LU is the whole
  story:** ~98% of the TI cost is the SuperLU numeric factorization, and the
  TD propagation is ~82% per-step triangular solves. **(2) MUMPS fixes
  both:** 8.2× on the factor and 4.6× on the per-solve at N=20328, so MUMPS
  is the correct default backend for TI and TD alike — and
  `SparseLU(backend="auto")` already selects it where provisioned. On the
  current direct-solver architecture that is near the practical optimum;
  PARDISO and GPU cuDSS are incremental beyond it, not step changes.
  **(3) There is NO pure-Python hot loop worth a first Rust kernel:** the
  `c_product`/extractor loops the first-kernel spec targeted are ~0.1% of
  runtime, so that spec's premise is invalidated. Separately: thread
  oversubscription cost ~300× on a concurrent sweep — pin BLAS threads per
  worker whenever processes are multiplied.
  ```

- [ ] Numeric multiset check; Global Constraints checks.

## Task 11: Split §8 out of the NRM note (docs-M11)

**Files:**
- Create: `docs/physics/nrm-vibrational-excitation.md`
- Modify: `docs/physics/nonlocal-resonance-model.md` (§8 → pointer stub; §9–§13 keep their numbers)
- Modify: `docs/physics/engine.md` (toctree + bullet)
- Modify: `docs/physics/dissociation.md` (blurb cross-link)
- Modify: `docs/physics/nrm-time-dependent.md`, `docs/physics/diatomic-ve-cross-sections.md`,
  `docs/molecules/n2.md`, `validation/diatomic/test_ve_nrm.py` (inbound §8.x references)

**Numbering decision (load-bearing):** the remaining sections §9–§13 KEEP their
numbers. `libs/qscat/qscat/core/nrm/ingredients.py` cites "Sec. 5 and Sec. 11"
in a raised error message and a comment; renumbering would silently break
those citations. The stub keeps the `## 8. Vibrational excitation` heading and
says explicitly that §8.1–§8.10 map to §1–§10 of the new note.

**Steps:**

- [ ] Snapshot for the numeric check:
  `git show HEAD:docs/physics/nonlocal-resonance-model.md > /tmp/nrm-before.md`.

- [ ] Create `docs/physics/nrm-vibrational-excitation.md`:
  - Title + header block (Key result opening per the Task 1 scheme):

    ```markdown
    # NRM vibrational excitation

    **Location:** `qscat.core.nrm.vibrational_excitation` (`j_dk`, `t_resonant`,
    `t_background`, `nrm_ve_cross_section`) and `qscat.core.nrm.scattering`
    (`scattering_state_minus`); `libs/qscat/tests/test_nrm_ve.py`,
    `test_nrm_scattering.py`; `validation/diatomic/ve_nrm.py` + `test_ve_nrm.py` +
    `ve_nrm_figure.py`; `apps/qscat-run`'s `nrm` method on a `ve` observable.
    **Source:** K. Houfek, T. N. Rescigno, C. W. McCurdy, *Phys. Rev. A* **77**,
    012710 (2008) — `reference/literature/houfek-2008-pra77-012710.md`. Split out
    of [`nonlocal-resonance-model.md`](nonlocal-resonance-model.md), which holds
    the method core (the kernel $F(E)$, the discrete-state choices, the
    ingredients) and the dissociative-attachment side; this note's §1–§10 were
    that note's §8.1–§8.10.
    **Units:** atomic units throughout (hartree, bohr).

    ## Key result

    Vibrational excitation is the channel PRA 77 plots for every molecule in
    its study, so it is where the nonlocal model can be checked against the
    exact solver most broadly — and choice B plus the background term
    reproduces the exact `driven.ve_cross_section` oracle to better than 0.7%
    on both N₂ (11 energies, 0.06–0.16 Ha) and F₂ (6 energies, 0.02–0.09 Ha),
    elastic and first-inelastic alike (worst ratio 0.99623–1.00692 over all
    four molecule/transition pairs), while choice A degrades to 0.565–1.140.
    The reason B is that good is physics, not luck: an R-independent
    $\phi_d$ carries no derivative couplings, so the model is formally exact
    and the residual is discretization error. The comparison is differential
    (both routes on the same grids): it validates the model reduction, not
    the grid.
    ```

    (The Key result paragraph's claims and numbers are §8's own — restated
    from 8.4 for the opening; the multiset check permits repetition.)
  - Body: current lines 647–1021 of the old note (subsections 8.1–8.10),
    verbatim, with ONLY these mechanical edits: headings `### 8.N Title` →
    `## N. Title` (and `#### 8.N.M` → `### N.M` if any exist — there are
    none today); internal references "§8.N"/"Sec. 8.N" → "§N"; references to
    the old note's other sections ("§2.3", "§6.2", "Eq. (60) … §6.2") →
    "[`nonlocal-resonance-model.md`](nonlocal-resonance-model.md) §2.3" form.
    Find them all with
    `grep -n '§[0-9]\|Sec\. [0-9]' docs/physics/nrm-vibrational-excitation.md`
    after the paste and fix each hit deliberately.
  - The intro lines of old §8 (635–645: "DA is the observable PRA 77 plots
    for **one** molecule…" + the Location block) fold into the header block
    above rather than repeating as a section.

- [ ] Reduce old §8 (lines 635–1021) to the stub, keeping the heading and
  number:

  ```markdown
  ## 8. Vibrational excitation

  Moved to its own note:
  [`nrm-vibrational-excitation.md`](nrm-vibrational-excitation.md) — the
  Eq. (28)/(31)/(37)/(38) two-potential decomposition, the φ⁻ gate, the VE
  state sum, the figure, and the measured result: choice B + background
  reproduces the exact `driven.ve_cross_section` oracle to better than 0.7%
  on both N₂ and F₂ (worst ratio 0.99623–1.00692 over all four
  molecule/transition pairs), while choice A degrades to 0.565–1.140. That
  note's §1–§10 are this section's former §8.1–§8.10; the sections below
  keep their original numbers (§9–§13) so existing references to them —
  including `qscat.core.nrm.ingredients`'s "Sec. 5 and Sec. 11" — stay
  valid.
  ```

- [ ] Register the new note in `docs/physics/engine.md` — toctree entry after
  `n2-2d-cross-section`, plus the bullet:

  ```markdown
  - {doc}`nrm-vibrational-excitation` — the nonlocal resonance model's
    vibrational-excitation route: two-potential background + resonant
    T-matrices on the shared nonlocal kernel, reproducing the exact solver
    to better than 0.7 % on N₂ and F₂.
  ```

- [ ] Update `docs/physics/dissociation.md`'s NRM bullet: after the
  "…reproduces the exact oracle to better than 0.7 % on both N₂ and F₂…"
  sentence, insert "(now its own note, {doc}`nrm-vibrational-excitation`,
  filed under the engine section)". The old note stays in this section's
  toctree — it still holds the method core and the DA side.

- [ ] Update every inbound §8.x reference (the complete list, from
  `grep -rn "§8\.\|Sec\. 8\." --include='*.py' --include='*.md' libs apps validation projects docs/physics docs/molecules`):
  - `validation/diatomic/test_ve_nrm.py:62` "Sec. 8.6" →
    "docs/physics/nrm-vibrational-excitation.md Sec. 6"; `:201` "Sec. 8.9" →
    "…nrm-vibrational-excitation.md Sec. 9".
  - `docs/physics/nrm-time-dependent.md` "§8.5" / "§8.4" (lines ~657, ~663) →
    "`nrm-vibrational-excitation.md` §5" / "§4".
  - `docs/physics/diatomic-ve-cross-sections.md:245`
    "`docs/physics/nonlocal-resonance-model.md` §8.4 publishes…" →
    "`docs/physics/nrm-vibrational-excitation.md` §4 publishes…".
  - `docs/molecules/n2.md:~138` "(§8.9 of the note quantifies it)" → "(§9 of
    {doc}`../physics/nrm-vibrational-excitation` quantifies it)"; and the
    figure paragraph's "Full derivation and validation ladder in
    {doc}`../physics/nonlocal-resonance-model`" →
    "{doc}`../physics/nrm-vibrational-excitation`". The "Where to read more"
    line gains the new note alongside the old one.
  - `docs/molecules/no-f2.md` — verify and leave: its NRM references are the
    DA story, which stays in the old note.
  - Plain-word references without a § number (`qscat.core.nrm.__init__`,
    `runner.py`, example YAML comments) still point at a note that exists
    and still describes the model core — verify each grep hit is not
    §8-specific and leave them.

- [ ] Numbers-preserved check across the split (paste output into the
  completion note):

  ```bash
  python3 - <<'PY'
  import re, collections, pathlib
  tok = lambda p: collections.Counter(re.findall(r'\d+(?:\.\d+)?(?:[eE][+-]?\d+)?', pathlib.Path(p).read_text()))
  before = tok('/tmp/nrm-before.md')
  after = tok('docs/physics/nonlocal-resonance-model.md') + tok('docs/physics/nrm-vibrational-excitation.md')
  lost = before - after
  assert not lost, f"numeric tokens lost in the split: {dict(lost)}"
  print("split preserves all", sum(before.values()), "numeric tokens")
  PY
  ```

  (Sub-multiset, not equality: heading renumbering 8.N → N legitimately
  changes tokens like "8.4", and the stub repeats a few results. The check
  direction — nothing LOST — is the invariant that matters; the "8.N" tokens
  survive via the stub's own mapping sentence and the old note's remaining
  sections, but if the check flags one, account for it by hand before
  overriding anything.)

- [ ] `uv run --no-sync pytest tests/test_docs_portability.py -q` — the new
  note is a normal note (NOT added to `SITE_FIRST_PAGES`), so it must pass
  every note rule including Task 6's detector; the ≥20-notes floor gains one.

- [ ] Sphinx build per the Global Constraints toctree clause — `-W` proves the
  new page is in exactly one toctree and every `{doc}` reference resolves.

- [ ] Run the Global Constraints checks.

## Task 12: Self-review pass

**Files:** none new.

- [ ] `grep -rn "8\.\(1\|2\|3\|4\|5\|6\|7\|8\|9\|10\)" docs/physics/nonlocal-resonance-model.md`
  — no live reference into the moved section remains inside the old note.
- [ ] Reconcile Task 1's conforming/retrofit lists in the skill against what
  Tasks 7–11 actually produced.
- [ ] Full sweep detector run over `docs/physics/*.md`: zero hits.
- [ ] All Global Constraints checks green; one full
  `uv run --no-sync pytest -m "not slow" -n auto --dist loadfile` on the
  finished branch.
