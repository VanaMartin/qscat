---
name: qm-method-lifecycle
description: Use when adding or porting any quantum-mechanics method or numerical capability to qModeling — enforces the design → Python toy model → validate → optimize-in-Rust → promote-to-qscat lifecycle.
---

# qm-method-lifecycle

## Overview

Every method or numerical capability added to qModeling — new DVR grid, ECS
contour, propagator, special function, etc. — moves through the same five
stages. Skipping a stage (e.g. hand-writing a Rust kernel before a Python toy
model exists, or promoting unvalidated code into `libs/qscat`) is the most
common way this repo accumulates unverified physics. Turn the checklist below
into todos at the start of the work and check them off in order.

## When to Use

- Starting work on a new quantum-mechanics method, algorithm, or numerical
  routine anywhere in the repo.
- Porting an algorithm out of a `reference/` oracle (`reference/eMoScat`,
  `reference/libXcuda`) into qModeling.
- Deciding whether code is "done" enough to live in `libs/qscat/`.

## The Lifecycle Checklist

1. **Design** — Write or point to a spec under `docs/superpowers/specs/`.
   State: the physics being modeled, the unit system (atomic units — Hartree
   for energy, Bohr for length, unless the spec says otherwise), the public
   interface (function signatures, inputs/outputs), and concrete success
   criteria (what numerical result proves correctness). Use
   `superpowers:brainstorming` first if the approach itself is still open.

2. **Toy model** — Implement a pure-Python version under `projects/<name>/`.
   Optimize for correctness and readability, not speed. This is throwaway
   scaffolding; it does not need to be production quality, but it must be
   numerically right — it becomes the differential-test oracle in step 4.

3. **Validate** — Invoke the `numerical-validation` skill. Do not proceed to
   optimization until validation passes: analytic benchmarks, convergence
   studies, conservation checks, and/or differential tests vs. `reference/`
   or `mpmath`, as appropriate to the method.

4. **Optimize (only if warranted)** — Only move to Rust if profiling
   (`pytest-benchmark` or `cProfile`) shows an actual hot path. If so, invoke
   the `python-to-rust-kernel` skill. If nothing is measurably slow, skip this
   stage — pure Python that meets its accuracy bar is a valid end state.

5. **Promote** — Move the validated, reusable implementation into
   `libs/qscat/qscat/<subpackage>/` (`special`, `dvr`, `ecs`, `evolution`,
   `linalg`, or `units`, per `qscat-conventions`). Bring its tests with it,
   and update `CLAUDE.md`/docs to describe the new capability and where it
   lives.

## Rule: the Python oracle is permanent

Once a Rust kernel exists (step 4), the pure-Python implementation from step
2 is never deleted. It stays in the codebase as the differential-test oracle
that the Rust kernel is continuously checked against — see
`native/qscat-kernels/tests/test_l2_norm.py` for the pattern: the Python/NumPy
computation and the Rust-backed call are asserted equal (within tolerance) on
the same inputs. Removing the Python version removes the ability to catch a
regression in the compiled kernel.

## Common Mistakes

- Writing the Rust kernel first "because it'll be faster anyway" — skips
  steps 2–3, so there is no oracle and no evidence optimization was needed.
- Promoting straight from the toy model in `projects/<name>/` into
  `libs/qscat/` without running `numerical-validation`.
- Deleting the Python implementation after a Rust kernel ships.
- Treating this skill as a substitute for `superpowers:test-driven-development`
  — TDD governs how each stage's code gets written; this skill governs which
  stage you're in and when to advance.
