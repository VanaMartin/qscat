---
name: python-to-rust-kernel
description: Use when a validated Python method has a proven hot path worth moving to a compiled kernel — profile, scaffold a PyO3/maturin crate in native/, mirror the Python API, benchmark, and keep Python as the differential oracle.
---

# python-to-rust-kernel

## Overview

qModeling only moves code to Rust when profiling proves it's worth it.
`native/qscat-kernels` is the reference pattern for every new kernel crate:
PyO3 bindings built with maturin, a Rust function that mirrors its Python
counterpart's signature, and a differential test that checks the Rust output
against the retained Python implementation on the same inputs.

**REQUIRED BACKGROUND:** the Python implementation must already exist and
have passed `numerical-validation`. This skill does not create new physics —
it recompiles proven physics for speed.

## When to Use

- Step 4 ("Optimize") of `qm-method-lifecycle`, after profiling shows a real
  hot path.
- Adding a new kernel alongside `qscat_kernels`, or extending it.
- Someone proposes "let's just write this in Rust" without profiling data —
  this skill's step 1 blocks that.

## Steps

1. **Profile first — no kernel without evidence.** Run `pytest-benchmark`
   (for function-level timing in tests) or `cProfile` on the Python
   implementation under realistic input sizes. Only proceed if a specific
   function is the measured bottleneck; "it's Python, it must be slow" is not
   evidence.

2. **Scaffold a crate in `native/<name>/`, mirroring `native/qscat-kernels`:**
   - `Cargo.toml` — `crate-type = ["cdylib"]`, `pyo3` dependency with the
     `extension-module` feature (match the pinned version in
     `native/qscat-kernels/Cargo.toml`).
   - `pyproject.toml` — `build-backend = "maturin"`, `[tool.maturin]
     module-name = "<crate>_kernels"` (or the appropriate module name), same
     shape as `native/qscat-kernels/pyproject.toml`.
   - `src/lib.rs` — a `#[pymodule]` function that registers one
     `#[pyfunction]` per exposed operation, following
     `native/qscat-kernels/src/lib.rs`.
   - Add the new crate to the workspace: it's picked up automatically by
     `[tool.uv.workspace] members = ["libs/*", "native/*"]` in the root
     `pyproject.toml`.

3. **API parity.** The Rust function's signature (argument order, types,
   units) mirrors the Python function it replaces exactly — same inputs
   in atomic units, same return shape. Callers should not need to know
   whether they're calling the Python or Rust version.

4. **Build it.** `uv run maturin develop` (add `--release` for a
   performance-representative build; see `containerize-and-run` for the
   `--release` build used in Docker/CI). Rebuilding after any `src/lib.rs`
   change is required before Python can see the new symbols.

5. **Differential test vs. the retained Python implementation.** Follow the
   pattern in `native/qscat-kernels/tests/test_l2_norm.py`: generate inputs
   (ideally with a seeded RNG across many samples), call both the Python
   oracle and the Rust-backed function, and assert agreement within an
   explicit tolerance (see `numerical-validation` for tolerance
   conventions). The Python implementation is never deleted — it's the
   permanent oracle.

6. **Add a `criterion` benchmark in Rust** for the hot function, so
   regressions in the compiled kernel's performance are caught independent
   of the Python-level `pytest-benchmark` numbers from step 1.

7. **Record the measured speedup** (e.g. in the PR description or the spec
   under `docs/superpowers/specs/`) — the before/after numbers from step 1
   vs. step 6 are the evidence that justified steps 2–6 in the first place.

## Common Mistakes

- Skipping step 1 and writing Rust for a function that was never actually
  the bottleneck.
- Changing the function's argument order or units relative to the Python
  version "since we're rewriting it anyway" — breaks API parity and silently
  invalidates the differential test's assumption that inputs are equivalent.
- Forgetting to rebuild (`uv run maturin develop`) after editing `src/lib.rs`
  and testing against a stale compiled extension.
- Deleting or skipping the Python oracle test once the Rust kernel is in
  place — it must keep running as the differential-test baseline.
