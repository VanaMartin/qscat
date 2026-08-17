# Documentation Accuracy Audit — Corrections Plan

**Goal:** Remove documentation content that the current code contradicts, and recompute the figures that are too coarse to support the claims made about them.

**Architecture:** The audit is already done (findings below). Phase 1 is text-only corrections, no compute. Phase 2 recomputes figures — three from a single cheap TD propagation, the rest as dense energy sweeps under Docker/MUMPS on `sadaharu`.

**Tech Stack:** the existing `validation.n2.experiment` harness, `qscat-run` configs, Docker `test` image (MUMPS), matplotlib.

---

## The audit result

Every quantitative claim in `docs/physics/` was cross-checked against current code. The
`validation.n2.experiment` harness was re-run on 2026-08-17 (**23 PASS / 0 FAIL / 6 NOTE**)
and reproduces, exactly, the numbers in:

| note | verified |
|---|---|
| `n2-2d-cross-section.md` | E1 anchors 0.125645, 0.055934, 0.012030, 0.009313, 0.0021926, 0.0018121, 0.4452, 0.7742, 0.8265, 1.000267 |
| `n2-td-cross-section.md` | D1 anchors 1.002, 0.988, 1.035, 1.017, 0.446, 0.765, 0.855 |
| `n2-cross-section.md` | C5 anchors 0.445, 0.774, 0.826, 1.010, 8133.082 |
| `n2-resonance.md` | B1 `E_res = 2.445 eV`, `Γ = 0.455 eV`, residual 3.30e-06 |
| `ti-energy-sweep-reuse.md` | 1.000267 |

`diatomic-ve-cross-sections.md`'s model quantities were recomputed directly and match to the
quoted precision (N₂ −0.2432/+0.5016, NO −0.0598/+0.1719, F₂ −0.1269/−0.0691, spacings
0.0124/0.0091/0.0039). `n2-resonance.md`'s pole is independently corroborated to 0.2% by
`shift-invert-eigensolver.md`. `td-extractors.md`, `lcp-resonance-levels.md`, `h2plus-dr.md`
post-date the fixes. The tooling notes (`mumps-sparse-backend`, `discretisation-tuning`,
`optimization-targets`) quote machine measurements and grid-point counts that no physics fix
invalidates.

**Exactly one note is substantively wrong: `n2-2d-td-cross-section.md`.** Plus a set of
under-resolved figures.

Two traps this audit had to avoid, both of which nearly produced wrong deletions:

- **"Absent from any `.py`" ≠ wrong.** ~150 numbers are absent from the code; almost all are
  correct-but-unpinned measurements. Deleting them mechanically would have destroyed correct
  physics — `diatomic-ve-cross-sections.md`'s thresholds are the proof.
- **Date-matching is not enough.** `n2-td-cross-section.md`'s numbers predate both fixes
  (2026-07-22) but describe the **1-D** solver, which legitimately still uses Crank–Nicolson.
  The Padé change was a 2-D change. The harness confirms those numbers are current.

## Global Constraints

- Do not delete a number without first checking whether current code reproduces it.
- Atomic units throughout.
- Never edit anything under `docs/superpowers/` except this plan file.
- Recomputed figures must be reproducible: each one records the command that made it.
- Preserve the repo's commit trailers on every commit.

---

## Phase 1 — corrections (no compute)

### Task 1: Fix `docs/physics/n2-2d-td-cross-section.md`

The note presents order-1 Crank–Nicolson measurements as current validation, while also
carrying the post-Padé capstone figure that contradicts them.

**Authoritative current values** (harness group F1, 2026-08-17):

| E (Ha) | σ_TD | σ_TI | ratio | the note currently claims |
|---|---|---|---|---|
| 0.10 | 5.9595 | 6.1230 | **0.9733** | 0.9305 |
| 0.15 | 0.61850 | 0.62576 | **0.9884** | 1.1033 |

Gates, from `projects/n2_2d_td_cross_section/test_td_cross_section.py`: `rel=0.06` at both
anchors (lines 110, 124); elastic TD/TI ~1.01 (E=0.14) / ~0.99 (E=0.15) at `rel=0.08`
(line 243).

- [ ] **D1** — Replace the validation-ladder table (`5.6973 / 6.1228 / 0.9305`, `1.1033`)
      with the values above, and the "rtol 0.10 / 0.15 … window-edge loosening" prose with
      the actual `rel=0.06` gate.
- [ ] **D2** — Delete the T-scan table (0.760/1.204 … 0.945/1.145) **and** the claim that
      E=0.15 "oscillates in [1.10, 1.20] … a floor, not a still-converging transient,
      diagnostic of the usable-spectral-window effect." That is an invented physical limit
      explaining a solver bug; `test_td_cross_section.py:115` says the opposite.
- [ ] **D3** — Delete the entire "A SEPARATE, finite-T resolution limit" section (ratios
      0.229/0.575/0.348, "off by 2-4×", "not fixable by re-tuning the wavepacket"). Same
      invented-limit error, and directly contradicted by the capstone figure's "~1-2% median
      across 0.04-0.18 Ha, boomerang oscillations resolved point-by-point" in the same file.
- [ ] **D4** — Delete the `n2-2d-td-sigma.png` figure block and its "three honestly-marked
      regions" caption; the figure renders D3's fictitious regions. Replaced in Task 3.
- [ ] **D5** — Replace "elastic … within ~15% (elastic TD/TI = 1.03–1.14 for E=0.13–0.17,
      from ~500× wrong)" with the measured ~1.01/~0.99. Keep the "from ~500× wrong" history
      — that part is true and explains why the free reference exists.
- [ ] **D6** — In the `F_out` discriminator table, keep the structural argument (regular vs
      Hankel differs by five to six orders of magnitude — that survives any propagator) but
      delete the order-1 CN ratio column values (1.101/1.144/1.113/1.234), which no longer
      describe the code. State the qualitative discriminator without the stale numbers.
- [ ] Also re-check the norm-decay profile (`1.000 → 0.024`) quoted as "physics fact 2"; it
      is an order-1 CN measurement. Re-measure in Task 3 or delete.
- [ ] Verify no remaining number in the file is contradicted by the harness.

### Task 2: Fix the same stale numbers where they propagated

The 0.93/1.10 pair leaked into three more files.

- [ ] `CLAUDE.md:285` — "σ_TD/σ_TI = 0.93 at E=0.10, 1.10 at E=0.15" → 0.9733 / 0.9884.
      Note `CLAUDE.md` currently **contradicts itself**: lines 86 and 294 already say the TD
      route matches to "~1-2%". Make the three consistent.
- [ ] `projects/n2_2d_td_cross_section/observation.py:346` — docstring "ratio 0.9305".
- [ ] `projects/n2_2d_td_cross_section/test_td_convergence.py:70-71` — the comment carrying
      the whole stale set (0.10→0.9305, 0.14→0.8648, 0.15→1.1033, 0.16→1.2249, 0.17→1.2777).
      If the surrounding test still depends on those values being true, fix the test too;
      if it is only a comment, correct the comment.
- [ ] `uv run --no-sync pytest projects/n2_2d_td_cross_section -q -m "not slow"` stays green.

---

## Phase 2 — recomputation

### Task 3: Recompute the three pre-fix TD figures

`n2-2d-td-sigma.png`, `n2-2d-td-correlation.png`, `n2-2d-td-snapshots.png` are all 2026-07-25,
i.e. order-1 CN and no free-reference. **All three derive from one propagation** (~210–250 s at
`TD_WORKING_GRID`), so this is cheap and can run locally.

- [ ] Re-run the propagation at the current defaults (order-3 Padé, `dt=1.0`, free reference).
- [ ] Regenerate all three figures from that single run.
- [ ] Re-measure the norm-decay profile for D5/"physics fact 2".
- [ ] If the σ-vs-E figure no longer shows distinct "regions" (expected — D3's regions were an
      artefact), plot it plainly against the TI oracle instead of inventing new regions.

### Task 4: Recompute the under-resolved figures densely (~0.001 Ha) on `sadaharu`

Current sampling vs what the physics needs. `diatomic-ve-cross-sections.md` states F₂'s
boomerang features are **~0.004 Ha wide**; its own F₂ figures sample at 0.004 and 0.0033 Ha —
about one point per feature. The note refutes its own figures.

| figure | now | target |
|---|---|---|
| `no-2d-ti-da-cross-section` | 12 pts @ 0.0136 Ha | 0.15–0.30 @ 0.001 |
| `no-2d-da-lcp-vs-exact` | 12 pts @ 0.0136 Ha | 0.15–0.30 @ 0.001 |
| `f2-2d-ti-da-cross-section` | 13 pts @ 0.0033 Ha | 0.01–0.05 @ 0.001 |
| `f2-2d-da-lcp-vs-exact` | 13 pts @ 0.0033 Ha | 0.01–0.05 @ 0.001 |
| `f2-2d-ti-cross-section` | 25 pts @ 0.004 Ha | 0.004–0.10 @ 0.001 |
| `no-2d-ti-cross-section` | 30 pts @ 0.004 Ha | 0.004–0.12 @ 0.001 |
| `n2-2d-ti-cross-section` | 60 pts @ 0.0033 Ha | 0.005–0.20 @ 0.001 |
| `diatomic-ve-comparison` | derived | rebuild from the above |

- [ ] Build the Docker `test` image on `sadaharu` (repo is at `/home/kooza/qModeling`; it has
      `qmodeling-base:cpu` and `qmodeling:runtime-cpu` but **not** the `test` image, which is
      the one carrying MUMPS).
- [ ] Write a `qscat-run` config per figure at the dense grid. The N₂ one uses the new
      `reference:` block to overlay Houfek directly.
- [ ] Run them under Docker/MUMPS. `SparseLU.refactor` reuse makes the per-energy cost small
      after the first factorization, so a ~200-energy sweep is far cheaper than 200 solves.
- [ ] Regenerate each figure and its `.npz` from the run output.
- [ ] Record, per figure, the exact command that produced it.

### Task 5: Reconcile the notes with the new figures

- [ ] Update any claim whose supporting figure changed — in particular the LCP-vs-exact
      "~11% away from threshold" figure, which rests on the 12–13 point curves.
- [ ] If a denser curve changes a stated conclusion, say so plainly rather than restating the
      old conclusion over a new figure. A conclusion that moves is a finding, not an
      embarrassment.
- [ ] Re-run `python -m validation.n2.experiment`; confirm still 23 PASS / 0 FAIL.
- [ ] `uv run sphinx-build -b html -W --keep-going docs docs/_build/html` clean.
- [ ] Full non-slow suite green.

---

## Definition of done

- No number in `docs/physics/` is contradicted by current code.
- No invented physical limitation remains in the documentation.
- Every recomputed figure resolves the structure its own note claims exists.
- The stale 0.93/1.10 pair is gone from all four files, and `CLAUDE.md` no longer contradicts
  itself.
- Harness 23 PASS / 0 FAIL; docs build `-W` clean; suite green.

## Explicitly not in this plan

The showcase gallery (its own plan). Adding new physics. Re-deriving the tooling notes'
machine-measured timings, which are true as recorded.
