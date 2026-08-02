# First Rust kernel (optimize-in-Rust starter) — Design Spec

**Date:** 2026-08-02
**Status:** Approved design (roadmap Part 5 #2)
**Lifecycle:** optimize stage — the first real `native/qscat-kernels` kernel
(the crate is currently a `l2_norm` stub, unused). Python stays the differential
oracle, per the lifecycle.

## Goal

Port the single hottest **pure-Python** inner loop on the TD propagation path to
a PyO3/Rust kernel in `native/qscat-kernels`, wired into `qscat` behind a
try-import-else-Python fallback, with a differential test and a benchmark.

## Why the TD path (not the factorization)

Profiling shows the TI cost is ~98% SuperLU `gstrf` — that is C already, owned by
the sparse backend (see the PARDISO/GPU specs), not a Rust target. The Rust
opportunity is the **per-step Python work** in the propagation loop: the extractor
`record` projections (`c_product(Φ_v', ψ)` per channel, every step, thousands of
steps) and the wavepacket/`eta` construction. High call count + small hot code +
a clean oracle = the ideal first-kernel shape.

## Step 0 — profile to SELECT the target (do not pre-commit)

Run `uv run python -m benchmarks.profile_hotpaths --td --top 30` and pick the top
pure-Python cumulative-time function that is NOT a scipy/numpy C call. Expected
candidates, in likelihood order: `td_extractors` `record` / `c_product`
accumulation; `correlation.eta_*`; `wavepacket.initial_state`. The spec's kernel
is "the measured winner"; the rest of this doc assumes it is the `c_product`
projection loop and adapts if the profile says otherwise.

## Non-goals

- No API change to the Python callers — the kernel is an internal accelerator.
- Not the factorization (backends), not GPU.
- No behavior change: the kernel must reproduce the Python result to round-off.

## Design

- Add the kernel to `native/qscat-kernels/src/lib.rs` (e.g. `c_product_rows`: the
  bilinear non-conjugated ECS inner product over a stack of complex rows). Accept
  numpy arrays via `numpy`/`ndarray` PyO3 bindings; return complex128.
- Wire into the Python hot path with a guarded import:
  `try: from qscat_kernels import <fn>; except ImportError: <python impl>`. The
  Python implementation stays as both fallback and oracle (the wheel is
  pure-Python by default; the kernel ships via the `qscat[kernels]` extra — see
  the PyPI coupling decision in the roadmap).
- `abi3-py312` on the crate so one wheel covers 3.12+ (a packaging fix noted in
  the audit).

## Validation

- **Differential test:** `native/`'s Python test asserts the kernel matches the
  Python oracle to `rtol=1e-12` on random complex inputs + the real N2 vectors
  (hypothesis-generated shapes). This is the lifecycle-required oracle test.
- **Benchmark:** pytest-benchmark comparing kernel vs Python on realistic sizes;
  and a before/after `profile_hotpaths --td` showing the loop dropped out of the
  top functions.
- Full suite green with the kernel built (`maturin develop`) AND with it absent
  (fallback path) — both must pass.

## Deliverables

- The kernel in `native/qscat-kernels` (+ `abi3-py312`, committed `Cargo.lock`).
- Guarded wiring in the qscat hot path + the Python oracle retained.
- Differential test + benchmark + a note in `docs/physics/optimization-targets.md`.

## Verification

`uv run maturin develop --manifest-path native/qscat-kernels/Cargo.toml` then
`uv run pytest native/ libs/qscat/tests -k "kernel or extractor"` green; the
fallback path green with the extension uninstalled; benchmark recorded.
