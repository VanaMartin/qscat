# TI Energy-Sweep Symbolic Reuse + Cross-Section Display — Design Spec (sub-project #9)

**Date:** 2026-07-26
**Author:** Martin (martin@qscat.com) with Claude
**Status:** Approved design — implementation pending
**Lifecycle:** `qm-method-lifecycle` stage 4 (optimize) for the reuse, plus a validation/display
deliverable. Builds on #6 (`ve_cross_section_2d`, the exact 2-D TI solver) and #8
(`qscat.linalg.SparseLU` + the MUMPS backend).

## Context

The exact 2-D time-independent solver (#6) computes σ(E) by factoring `A(E) = E_tot·I − H`
**once per energy**: `ve_cross_section_2d` loops over energies, building a fresh `SparseLU` each
time — a full symbolic analysis (ordering) **and** numeric factorization at every energy. But
`A(E)` has an **identical sparsity pattern at every energy** — only the diagonal shift `E_tot·I`
changes (`H` is fixed; the identity keeps the diagonal populated). MUMPS separates the
**symbolic analysis** (ordering + symbolic factorization, pattern-only) from the **numeric
factorization** (values), so the analysis can be done **once** and reused across all energies —
skipping the SCOTCH ordering `N−1` times for an `N`-energy sweep.

**Verified (Docker spike):** MUMPS `analyze()` once on the pattern, then
`factor(reuse_analysis=True)` per shifted matrix (SYM=2), matches a fresh SuperLU solve to
~1e-15 across different diagonal shifts. The reuse is correct.

**The paired deliverable:** with the dense sweep now cheap (and clean), compute the **full exact
2-D σ(E) curve** for N₂ — which we do not yet have (#6 gave anchors only; #7's TD curve was
finite-T-limited) — and **display it** against Houfek's `CSVE.V00.J00` golden data. The reuse
makes the dense curve affordable; the display is what one does with it, and validating the curve
against Houfek exercises the reuse end-to-end.

## Deliverable 1 — symbolic/numeric reuse

**API: `SparseLU.refactor(A_new)`** — re-run only the **numeric** factorization, reusing the
symbolic analysis from construction.
- **MUMPS backend:** `factor(reuse_analysis=True)` on the same Context — skips re-ordering.
- **scipy backend:** no clean symbolic reuse in SuperLU, so `refactor` re-runs `splu(A_new)` —
  correct, but no speedup. This is the honest fallback; the reuse win materializes only on MUMPS.
- **Pattern guard:** `A_new` must share the analyzed matrix's structure (shape + nnz, and for
  the symmetric path the `triu` pattern). A cheap guard raises on a mismatch — `reuse_analysis`
  is only valid for an identical pattern. (`E_tot·I − H` always satisfies this.)
- `backend_used`/`symmetric`/diagnostics behave as in #8; `refactor` keeps the same backend and
  symmetry decision as the original factorization.

**Wire `ve_cross_section_2d` to sweep with reuse:** build `H` once, construct the solver once at
the first energy (analyze + factor), then `refactor(E_tot·I − H)` + `solve` per subsequent
energy — no change to the returned σ or its shape contract. Physics is untouched; only the
factorization path changes.

## Deliverable 2 — cross-section display

**A reusable plotting utility** `plot_cross_sections(E_grid, sigma, *, reference=None, thresholds=None, ...) -> saves a PNG` that takes computed cross sections `sigma[E, channel]` (bohr²) on an
energy grid and renders the σ_{0→v'}(E) curves per VE channel, optionally overlaid on reference
golden data and with channel thresholds marked. Experiment-agnostic (N₂ now; F₂/NO reuse it).

**Applied to N₂:** compute the dense exact-2D σ_{0→v'}(E) curve over the resonance region
(e.g. [0, 0.2] Ha) for the VE channels, via the reuse-enabled sweep, and plot it against
Houfek's `CSVE.V00.J00`. Commit the figure and the numeric `.npz` (the σ(E) data), matching the
#6/#7 deliverable pattern.

## Interface

```
libs/qscat/qscat/linalg/sparse_lu.py     SparseLU.refactor(A_new) + (MUMPS reuse / scipy fallback)
projects/n2_2d_cross_section/
  cross_section_2d.py                     ve_cross_section_2d sweeps with analyze-once/refactor
  sigma_curve.py (or extend)              dense exact-2D sigma_{0->v'}(E) over an energy grid
  cross_section_plot.py                   plot_cross_sections(...) reusable utility (+ N2 figure driver)
docs/physics/figures/n2-2d-ti-cross-section.png   the committed curve vs Houfek
benchmarks/                               (extend) the energy-sweep reuse speedup benchmark
```

Reuses `qscat.linalg.SparseLU` (the reuse lives here), `projects.n2_2d_cross_section`
(`build_h2d`, `ve_cross_section_2d`, `vibrational_states`), `validation.n2.loader` (Houfek data
— note the import direction: the *display driver* may live in `validation/` or read Houfek via a
path, since `projects/` must not import `validation/`).

## Validation

**V1 — reuse is bit-identical.** `SparseLU.refactor(A_new)` then `solve(b)` gives the same result
(to round-off) as a fresh `SparseLU(A_new).solve(b)` — on both backends, on complex-symmetric
matrices, over several distinct shifts. The pattern guard rejects a mismatched matrix. (MUMPS
path `@skipif` in Docker; scipy path on the Mac.)

**V2 — physics unchanged.** `ve_cross_section_2d` with the reuse sweep returns σ identical (to
round-off) to the pre-reuse per-energy result, at the benchmark anchors and across a small grid.
The N₂ harness (23/0/6/0) is unchanged.

**V3 — the sweep speedup (measured).** Benchmark a multi-energy sweep (e.g. 50–100 energies on
the working grid) with reuse vs without (per-energy fresh factorization), on the MUMPS backend,
reporting the wall-time saving. **Measure it** — the analysis is ~10–30% of a MUMPS
factorization, so the expected saving is bounded and grows with the energy count; if it is
smaller than expected, report honestly. The scipy fallback shows no speedup (expected).

**V4 — the σ(E) curve matches Houfek.** The dense exact-2D σ(E) curve reproduces Houfek's
`CSVE.V00.J00` across the resonance region within the documented cross-model bound (the same
tolerance #6's anchors met — this is the same solver, now on a dense grid). The figure honestly
shows the agreement (and any near-threshold / elastic-background limitations already documented
for the LCP, which the exact solver largely closes — as #6 showed at the anchors).

## Out of scope

- `complex64` + iterative refinement; MKL PARDISO; the publishing pipeline; a Rust non-LU kernel
  (each a later sub-project, per #8's deferred list).
- Interactive/HTML display (a possible later bonus; this ships the static figure + utility).
- The DA channel; higher partial waves; the TD route (the display is the TI curve).
- Reusing the factorization across *right-hand sides* (already done — `solve` takes many RHS);
  this is reuse across *energies* (the matrix's diagonal shift).

## Verification

- `uv run pytest libs/qscat projects/n2_2d_cross_section validation/n2 -q` → all pass (MUMPS
  reuse tests skip on the Mac, run in Docker).
- `uv run mypy libs/qscat` → 0; `uv run ruff check .` → clean.
- `uv run python -m validation.n2.experiment` → 23 PASS / 0 PENDING / 6 NOTE / 0 FAIL, exit 0,
  unchanged (physics is backend/path-independent).
- `docker/build.sh test` → passes (the reuse benchmark + differential run in-container).
- The committed σ(E)-vs-Houfek figure + the `.npz`; the measured sweep-reuse speedup and the
  curve-vs-Houfek agreement recorded in `docs/physics/`; `CLAUDE.md` updated (`SparseLU.refactor`,
  the reuse sweep, the display utility).
