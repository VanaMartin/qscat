# GPU sparse backend for `qscat.linalg.SparseLU` (cuDSS) — Design Spec

**Date:** 2026-08-02
**Status:** Approved design — DESIGN ONLY (implementation deferred; GPU is a
"later" roadmap item). Written now so the interface is fixed before it lands.
**Lifecycle:** optimize stage — a fourth backend behind the existing dispatch.

## Goal

Define the GPU sparse-direct backend `backend="cuda"` for `SparseLU`, using
**NVIDIA cuDSS** (the GPU sparse direct solver; cuSOLVER as fallback) for the
complex-symmetric ECS matrices, so the largest decks (H₂⁺ ~1.15M unknowns, and
future higher-dimensional models) factor on-GPU. The `base-gpu.Dockerfile` scaffold
and `docker/build.sh gpu` already stage the environment.

## Why now (design only)

The point of writing this before implementing: the `SparseLU` backend dispatch,
`refactor` symbolic-reuse contract, and `BackendError` behavior must not have to
change when GPU arrives. Fixing the interface now means the GPU backend "slots in
like MUMPS did," and the D-general layer (`linalg`/`dvr`/`tensor`) stays
backend-agnostic (no CPU-only assumptions baked in).

## Non-goals

- No GPU kernels for the rest of qscat (propagation, assembly) — this is the
  solver backend only. Those are separate optimize tasks.
- No multi-GPU / distributed. Single-GPU first.
- Not implemented in this spec — this is the contract + a phased plan.

## Design

- New `libs/qscat/qscat/linalg/_cuda_backend.py`, same factor/`solve` shape as
  `_mumps_backend.py`. Matrix + RHS transferred to device via **cupy**; cuDSS
  performs analyze/factor/solve. `cuda_available()` guards import (cupy + a
  visible device).
- `refactor(reuse_analysis=True)` maps to cuDSS's analyze-once/refactor-many
  (the exact analog of MUMPS reuse); a pattern guard raises `BackendError` on a
  structure mismatch (same contract as today).
- Data-movement policy: keep the factor resident on-device across a `refactor`
  sweep; only RHS/solution cross the PCIe boundary per energy. Document the
  host↔device transfer as the cost model input.
- `auto` remains CPU-only (`mumps → scipy`); GPU is always explicit
  (`backend="cuda"`) — a GPU may be absent or smaller than the problem.
- Precision: complex128 on-device; fall back to `BackendError` (not a silent
  precision drop) if a device lacks fp64 capacity for the deck.

## Validation (when implemented)

- **Differential oracle:** cuDSS `solve` matches SuperLU to `rtol=1e-8`
  (fp64 GPU vs CPU) on the N2 2-D matrix and a complex-symmetric RHS.
- **refactor equivalence** across a diagonal-shift sweep, matching CPU.
- **Scale test:** the H₂⁺ full deck factors on-GPU; benchmark factor+solve+
  transfer vs MUMPS/PARDISO (extend the three-way benchmark to four-way).
- `@skipif` no CUDA device; runs only in the `gpu` Docker image with `--gpus all`.

## Deliverables (when implemented)

- `_cuda_backend.py` + dispatch wiring + `cuda_available()`; a `cuda` extra
  (cupy + cuDSS bindings); `base-gpu.Dockerfile` completed (cuDSS + MUMPS `.pc`
  + sanity gate — the TODOs in that file); a `docs/physics/gpu-sparse-backend.md`.

## Phased plan

1. **Now:** this interface contract (no code).
2. Complete `base-gpu.Dockerfile` (cuDSS libs, sanity gate).
3. Implement `_cuda_backend.py` against the contract; differential + refactor
   tests in the `gpu` image.
4. Four-way benchmark; decide whether `auto` should prefer GPU when present.
