# 5. Two test tiers: a toy-scale gate and a production-scale suite

Date: 2026-08-24

## Status

Accepted

## Context

A `slow` marker has existed since early in the project, and CI has always run
`pytest -m "not slow"`. But nothing recorded what qualified a test for the
marker, so it was applied case by case and drifted. An audit measured the
resulting state.

Profiling the `-m "not slow"` suite serially, with per-test wall time and peak
RSS, gave 879 tests in 19m56s locally and 24.7 min on a GitHub runner — for the
tier that is supposed to be the quick one. The distribution was the finding:

| threshold | tests | share of test time |
|---|---|---|
| >= 30 s | 13 | 51% |
| >= 10 s | 30 | 73% |
| < 1 s | 784 | 4% |

Thirty tests carried 73% of the gate. Some were doing real production-scale
physics and simply had not been marked. Others were the opposite problem —
tests whose own docstrings describe them as cheap plumbing checks, running on
production decks:

- `test_public_api_shape_contract` asserts only a return *shape*, is documented
  "Cheap (not `slow`)", and cost 49.9 s and 3.83 GB — the largest single
  allocation in the gate.
- `test_td_propagation.py` spent 206 s across six snapshot-bookkeeping tests,
  each paying its own sparse factorization of a 43674-dimension system.

Memory turned out to matter more than time. Peak RSS in the "fast" gate reached
3.83 GB, and the marked tier is far heavier — the H2+ DR example was measured at
~19 GB peak on its own, on a deck its own config calls a "laptop-sized
reduction". Meanwhile parallelism was configured nowhere: `pytest-xdist` was
installed but no workflow, Dockerfile, or `addopts` ever passed `-n`. The one
place `-n 8` was recommended (CLAUDE.md, as "the fast way to iterate") omitted
`-m "not slow"`, so following that advice ran the 19 GB test alongside seven
other workers. That combination exhausted memory on a 64 GB development machine
during ordinary work.

## Decision

1. **Two tiers, split by cost.** The default tier is toy-scale: every test
   proves a piece of library behaviour on the smallest deck that can still fail
   for the right reason. The `slow` tier is production-scale: real molecule
   decks, converged grids, multi-thousand-step propagations, and comparisons
   against published values.

2. **The threshold is a few seconds, or ~0.5 GB.** A test exceeding either
   belongs in `slow` — *or* it wants a smaller deck, which is often the better
   answer. Ask which of the two applies before reaching for the marker.

3. **The boundary is cost, not importance.** A `slow` test is frequently the
   more valuable one. Nothing about this tiering says the fast tier is what
   matters.

4. **CI runs the fast tier only.** GitHub runs `-m "not slow"`. The slow tier
   runs locally (`pytest -m slow`) and in the Docker `test` image. This is
   deliberate: these are research-grade physics runs whose cost is measured in
   minutes and gigabytes, and a hosted runner is the wrong place for them.

5. **Parallel runs use `--dist loadfile` and always pair with `-m "not slow"`.**
   Memory is the binding constraint, not CPU. Several modules build their grid
   at module scope; the default `--dist load` distributes per test, so every
   worker touching the module builds its own copy. `loadfile` bounds that to one
   copy per file. The Docker image runs the fast tier parallel and the slow tier
   serially, because the slow decks would OOM a container long before
   concurrency saved wall-clock.

6. **Pin BLAS to one thread per worker wherever `-n` is used on Linux.** The
   first parallel CI run was *slower* than the serial one it replaced — the step
   ran past 22 min, against 24.7 min serially and 123 s for the identical
   command locally at `-n 4`. The runner has 4 vCPU, `-n auto` starts 4 workers,
   and OpenBLAS defaults to one thread per core in each: 16 threads for 4 cores.
   With `OMP_NUM_THREADS=1` / `OPENBLAS_NUM_THREADS=1` / `MKL_NUM_THREADS=1` the
   step takes 186 s.

   This does not reproduce on macOS, where scipy links Accelerate: pinning there
   changed nothing (74.6 s vs 76.1 s at a saturating `-n 12`). A local
   measurement of threading behaviour does not transfer to the Linux runners,
   and CLAUDE.md's older note that pinning "is not worth bothering with" was one
   of those.

7. **When a test crosses the boundary, prefer shrinking the deck.** If what is
   under test is plumbing — shapes, cadences, dispatch, bookkeeping — it should
   assert that on a toy deck and stay in the gate. Reserve the marker for tests
   whose *assertion* genuinely depends on the production deck. A tolerance
   calibrated to one grid's round-off floor is the clearest example: it cannot
   move to a smaller grid without silently becoming vacuous.

## Consequences

- The CI gate drops from 24.7 min to a measured 186 s (3.12) / 239 s (3.13) for
  the pytest step, and no longer depends on a developer noticing that a new test
  is expensive.
- Both Python versions keep running the *full* fast tier. The earlier argument
  for trimming the matrix was the 50 runner-minutes it cost; at ~4 min a job
  that argument no longer holds.
- The slow tier is not run by any automation on GitHub. It is therefore possible
  to merge a change that breaks it. That is an accepted trade: the tier exists
  because these runs are expensive, and the mitigation is to run
  `pytest -m slow`, or `docker/build.sh test`, before merging anything touching
  the solvers.
- Point 6 creates ongoing work. Every test that gets marked `slow` should first
  be examined for whether a smaller deck would do; the audit found roughly half
  of the offenders were in that category.
